"""
exonerate_utils.py
------------------
Exonerate-based, frame-safe locus / gene-of-interest extraction.

Why this exists
===============
The BLAST HSP-as-exon path (``blast_loci_utils.extract_from_hsps``) stitches
tblastn / blastn HSPs together as exons with no reading-frame check. A 1-2 bp
HSP-boundary error frameshifts the CDS silently (see PLANNING.md → RM-002 risk
register). Exonerate's *spliced-alignment* models solve this directly: they align
a reference protein (``protein2genome``) or coding sequence (``coding2genome``)
against genomic DNA using an explicit intron + splice-site model, yielding
accurate exon/intron boundaries and a translatable reading frame.

Two jobs (D-008):
  1. Frame-safe CDS for the curated protein-coding loci (TEF1, RPB2, TUB2, …).
  2. "Locate any gene of interest" — give it an orthologous protein/CDS and a
     genome, and get back just the CDS + exon model for that gene.

Hybrid pipeline
===============
Whole-genome Exonerate is slow, so by default we *narrow first*: tblastn/blastn
finds the single best contig (reusing ``blast_loci_utils.select_best_locus_group``),
then Exonerate runs against that one contig. Coordinates therefore stay in contig
space — identical to the rest of phylofetch — with no offset arithmetic. When
BLAST finds nothing (common for distant orthologs where seeding fails) we fall
back to running Exonerate against the whole assembly.

Output (per locus per strain) mirrors the BLAST path plus a translated protein:
  LOCUS_CDS.fasta       — intron-free CDS, rich provenance header
  LOCUS_protein.fasta   — translated CDS (Exonerate is protein-aware)
  LOCUS_genomic.fasta   — gene span (exons + introns)
  LOCUS_introns.fasta   — individual introns with splice sites
  LOCUS.gff3            — exon/CDS annotation with proper phase (reused writer)
  LOCUS_partition.nex   — codon-position partitions (reused writer)
  LOCUS_extraction.log  — self-contained log incl. frame/stop/splice QC

All external commands route through RunManager when one is supplied, so every
Exonerate (and narrowing-BLAST) call is logged with its command + tool version.
"""

from __future__ import annotations

import os
import re
import shutil
import subprocess
import tempfile
from datetime import datetime, timezone
from pathlib import Path
from typing import Optional

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# Reuse the BLAST-path primitives — no duplication.
from .blast_loci_utils import (
    _load_contig,
    detect_fasta_type,
    parse_blast_hsps,
    run_blast,
    select_best_locus_group,
    write_codon_partition,
    write_gff3,
    write_region_gff3,
    _build_description,
)

# Query alphabet → Exonerate model. protein2genome is the workhorse (protein
# orthologs map across wider divergence and are frame-safe); coding2genome takes
# a nucleotide CDS / cDNA query.
MODEL_FOR_QUERYTYPE = {
    "protein": "protein2genome",
    "nucleotide": "coding2genome",
}

# RYO ("roll your own") format. %tcs is Exonerate's authoritative spliced target
# coding sequence — we rebuild the CDS from GFF coords and cross-check it against
# %tcs. Fields after the marker: query id, target id, score, %id, query length,
# query aln begin/end (0-based), target coding begin/end (0-based; unreliable for
# extraction — we use the GFF cds features instead, but keep them for the log).
_RYO_FIELDS = ["qi", "ti", "s", "pi", "ql", "qab", "qae", "tcb", "tce"]
_RYO = ">RYO\\t" + "\\t".join(f"%{f}" for f in _RYO_FIELDS) + "\\n%tcs\\n"


def _now_utc() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat()


def _probe_exonerate_version(exonerate_bin: str = "exonerate") -> str:
    try:
        r = subprocess.run([exonerate_bin, "--version"],
                           capture_output=True, text=True, timeout=10)
        for line in (r.stdout + r.stderr).splitlines():
            if "exonerate" in line.lower() and any(c.isdigit() for c in line):
                return line.strip()[:80]
    except Exception:
        pass
    return "unknown"


# ── Exonerate runner ───────────────────────────────────────────────────────

def run_exonerate(
    query_fasta: str,
    target_fasta: str,
    model: str,
    output_dir: str,
    *,
    exonerate_bin: str = "exonerate",
    minintron: int = 20,
    maxintron: int = 2000,
    bestn: int = 1,
    score: Optional[float] = None,
    refine: str = "none",
    geneticcode: int = 1,
    manager=None,
    module: str = "loci_extraction",
    action: str = "exonerate",
    timeout: int = 900,
) -> tuple[int, str, str]:
    """
    Run Exonerate and return ``(returncode, raw_stdout_text, gff_path)``.

    ``maxintron`` defaults to 2000 bp: Exonerate's own default (200 kb) invites
    spurious giant introns, whereas fungal introns are typically 50-300 bp.
    When ``manager`` (a RunManager) is given the command is executed and logged
    through it (command + exonerate version captured in the run folder).
    ``timeout`` (seconds) bounds the call so a wedged Exonerate can't hang the app;
    a timeout or launch failure is surfaced as a non-zero rc, not an exception (D-033).
    """
    os.makedirs(output_dir, exist_ok=True)
    cmd = [
        exonerate_bin,
        "--model", model,
        "--query", query_fasta,
        "--target", target_fasta,
        "--showtargetgff", "yes",
        "--showalignment", "no",
        "--showvulgar", "no",
        "--bestn", str(max(1, bestn)),
        "--minintron", str(minintron),
        "--maxintron", str(maxintron),
        "--geneticcode", str(geneticcode),
        "--ryo", _RYO,
    ]
    if refine and refine != "none":
        cmd += ["--refine", refine]
    if score is not None:
        cmd += ["--score", str(score)]

    if manager is not None:
        rr = manager.run(
            cmd, module=module, action=action,
            tool_version_keys=["exonerate"],
            inputs={"query": query_fasta, "target": target_fasta},
            params={"model": model, "minintron": minintron, "maxintron": maxintron,
                    "bestn": bestn, "refine": refine, "geneticcode": geneticcode},
            timeout=timeout,
        )
        try:
            raw = Path(rr.stdout_path).read_text(encoding="utf-8")
        except OSError:
            raw = ""
        rc = rr.returncode
    else:
        try:
            proc = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)
            raw = proc.stdout or ""
            rc = proc.returncode
        except subprocess.TimeoutExpired:
            raw, rc = "", 124
        except (FileNotFoundError, OSError):
            raw, rc = "", 127

    # Save Exonerate's verbatim stdout for provenance/debugging — but as `.txt`, NOT `.gff`:
    # the raw output begins with `Command line: [...]` / `Hostname:` lines and is NOT valid GFF,
    # so a `.gff` name invited importing it into Geneious, which then choked on `--model` where it
    # expected a coordinate. The clean, importable annotation is the per-locus `LOCUS.gff3`
    # (written by `write_gff3`); this file is only the raw tool log. (D-029)
    raw_path = os.path.join(output_dir, "exonerate_raw.txt")
    Path(raw_path).write_text(raw, encoding="utf-8")
    return rc, raw, raw_path


# ── Output parsing ───────────────────────────────────────────────────────────

def parse_exonerate_gff(text: str) -> list[dict]:
    """
    Parse Exonerate ``--showtargetgff yes`` + ``--ryo`` output into a list of
    gene models (one dict per alignment result), highest-scoring first as emitted.

    Exonerate prints, *per result*, a ``# --- START OF GFF DUMP ---`` … END block
    immediately followed by that result's ``>RYO`` line + wrapped %tcs sequence.
    Splitting on the START marker therefore isolates each result with its own RYO.

    Each model dict:
      contig, strand ('+'/'-'), model, query_id, score, pident,
      query_len, q_aln_begin, q_aln_end (0-based), tcs (spliced CDS, upper),
      gene_start, gene_end, cds (list of (start,end), 1-based, start<=end),
      exons, introns  (same coord convention)
    """
    models: list[dict] = []
    # Drop the leading banner; each subsequent chunk is one result.
    chunks = text.split("# --- START OF GFF DUMP ---")
    for chunk in chunks[1:]:
        gff_part, _, ryo_part = chunk.partition("# --- END OF GFF DUMP ---")

        model: dict = {
            "contig": "", "strand": "+", "model": "", "query_id": "",
            "score": 0.0, "pident": 0.0, "query_len": 0,
            "q_aln_begin": 0, "q_aln_end": 0, "tcs": "",
            "gene_start": 0, "gene_end": 0,
            "cds": [], "exons": [], "introns": [],
        }

        for line in gff_part.splitlines():
            if not line or line.startswith("#"):
                continue
            cols = line.split("\t")
            if len(cols) < 8:
                continue
            seqname, source, feature = cols[0], cols[1], cols[2]
            try:
                start, end = int(cols[3]), int(cols[4])
            except ValueError:
                continue
            if start > end:
                start, end = end, start
            strand = cols[6]
            if feature == "gene":
                model["contig"] = seqname
                model["model"] = source
                model["strand"] = strand if strand in ("+", "-") else "+"
                model["gene_start"], model["gene_end"] = start, end
                try:
                    model["score"] = float(cols[5])
                except ValueError:
                    pass
            elif feature == "cds":
                model["cds"].append((start, end))
            elif feature == "exon":
                model["exons"].append((start, end))
            elif feature == "intron":
                model["introns"].append((start, end))

        # RYO line + wrapped %tcs (sequence lines until blank / next marker).
        tcs_lines: list[str] = []
        in_tcs = False
        for line in ryo_part.splitlines():
            if line.startswith(">RYO"):
                parts = line.split("\t")[1:]
                for key, val in zip(_RYO_FIELDS, parts):
                    if key == "qi":
                        model["query_id"] = val
                    elif key == "s":
                        try:
                            model["score"] = float(val)
                        except ValueError:
                            pass
                    elif key == "pi":
                        try:
                            model["pident"] = float(val)
                        except ValueError:
                            pass
                    elif key == "ql":
                        model["query_len"] = int(val) if val.isdigit() else 0
                    elif key == "qab":
                        model["q_aln_begin"] = int(val) if val.lstrip("-").isdigit() else 0
                    elif key == "qae":
                        model["q_aln_end"] = int(val) if val.lstrip("-").isdigit() else 0
                in_tcs = True
                continue
            if in_tcs:
                s = line.strip()
                if s and all(c in "ACGTNacgtn" for c in s):
                    tcs_lines.append(s)
                elif s.startswith(">"):
                    break
                elif s:
                    break
        model["tcs"] = "".join(tcs_lines).upper()

        if model["contig"] and (model["cds"] or model["exons"]):
            models.append(model)

    return models


def select_best_model(models: list[dict],
                      min_score: float = 0.0) -> tuple[Optional[dict], list[dict]]:
    """
    Pick the highest-scoring model; return ``(best, others)``.

    ``others`` are the remaining candidates (e.g. paralogs / tandem duplicates),
    surfaced so the UI can let the user disambiguate (RM-002 paralog awareness).
    """
    eligible = [m for m in models if m.get("score", 0.0) >= min_score]
    if not eligible:
        return None, []
    ranked = sorted(eligible, key=lambda m: m.get("score", 0.0), reverse=True)
    return ranked[0], ranked[1:]


# ── CDS quality control (RM-002) ───────────────────────────────────────────

def validate_cds(cds_seq: str, geneticcode: int = 1) -> dict:
    """
    Translate the CDS and report frame / internal-stop QC.

    Returns: protein, length, len_mod3, n_internal_stops, starts_with_atg,
    ends_with_stop, clean (len%3==0 and no internal stops).
    """
    cds = (cds_seq or "").upper().replace("-", "")
    length = len(cds)
    # Translate complete codons only (a trailing partial codon can't be translated);
    # length / len_mod3 below are still measured on the full CDS.
    try:
        protein = str(Seq(cds[: length // 3 * 3]).translate(table=geneticcode))
    except Exception:
        protein = ""
    ends_with_stop = protein.endswith("*")
    n_total_stops = protein.count("*")
    n_internal_stops = n_total_stops - (1 if ends_with_stop else 0)
    return {
        "protein": protein,
        "length": length,
        "len_mod3": length % 3,
        "n_internal_stops": n_internal_stops,
        "starts_with_atg": cds.startswith("ATG"),
        "ends_with_stop": ends_with_stop,
        "clean": (length % 3 == 0) and (n_internal_stops == 0),
    }


# ── Genomic soft-masking (exon/intron boundary annotation) ──────────────────

def soft_mask_genomic(contig_seq: str, g_min: int, g_max: int,
                      intron_coords: list[tuple[int, int]], is_minus: bool) -> str:
    """
    The genomic gene span ``contig_seq[g_min..g_max]`` (1-based, inclusive), oriented to the
    coding strand and written **exons UPPERCASE / introns lowercase**.

    Soft-masking the introns by case carries Exonerate's exon/intron boundaries into the
    sequence itself, so they stay visible — and pinned to the right bases as gaps are inserted —
    when the gene is aligned/curated by hand, while remaining inert to aligners and tree tools
    (nucleotide case is ignored). Used for every Exonerate genomic product so the user's own
    extracted loci and the comparison tips are annotated identically (D-022). Returns "" for an
    empty/invalid span; Biopython's reverse_complement preserves case, so the masking survives
    the minus-strand flip.
    """
    if g_min <= 0 or g_max <= 0 or g_max < g_min:
        return ""
    chars = [c.upper() for c in contig_seq[g_min - 1:g_max]]
    for s, e in intron_coords:
        lo, hi = (s, e) if s <= e else (e, s)
        for p in range(max(lo, g_min), min(hi, g_max) + 1):
            chars[p - g_min] = chars[p - g_min].lower()
    seq = "".join(chars)
    return str(Seq(seq).reverse_complement()) if is_minus else seq


# ── Build the shared result dict from an Exonerate model ────────────────────

def build_result_from_model(
    model: dict,
    contig_seq: str,
    strain_id: str,
    locus_name: str,
    *,
    ref_accession: str = "",
    exonerate_version: str = "",
    geneticcode: int = 1,
    n_other_models: int = 0,
) -> dict:
    """
    Convert a parsed Exonerate model + the contig sequence into the same result
    dict shape ``blast_loci_utils.extract_from_hsps`` produces (so write_gff3 /
    write_codon_partition are reused unchanged), augmented with QC fields.

    CDS/exon segments are rebuilt from genomic coordinates so that the CDS,
    genomic span and introns are mutually consistent. On the minus strand,
    segments are taken in transcript order (high genomic coord first) and each is
    reverse-complemented — matching how the BLAST path orders minus-strand exons.
    """
    is_minus = model["strand"] == "-"

    def _segments(feats: list[tuple[int, int]]) -> list[tuple[int, int]]:
        ordered = sorted(feats, key=lambda t: t[0])      # genomic ascending
        return list(reversed(ordered)) if is_minus else ordered

    def _oriented(seq: str) -> str:
        return str(Seq(seq).reverse_complement()) if is_minus else seq

    # Prefer cds features; exonerate emits cds==exon for these coding models.
    cds_feats = model["cds"] or model["exons"]

    exons: list[dict] = []
    cds_parts: list[str] = []
    for (g_start, g_end) in _segments(cds_feats):
        seg = contig_seq[g_start - 1:g_end]
        cds_parts.append(_oriented(seg))
        exons.append({
            "g_start": g_start, "g_end": g_end,
            "length": g_end - g_start + 1,
            "pident": model.get("pident", 0.0),
            "q_start": 0, "q_end": 0,
        })
    cds_seq = "".join(cds_parts).upper()

    # Genomic gene span (exons + introns), soft-masked: exons UPPERCASE / introns lowercase so
    # the boundaries are visible for by-hand alignment (case is inert to aligners/tree tools).
    g_min = model["gene_start"] or min((e["g_start"] for e in exons), default=0)
    g_max = model["gene_end"] or max((e["g_end"] for e in exons), default=0)
    amplicon = soft_mask_genomic(contig_seq, g_min, g_max, model["introns"], is_minus)

    # Introns: Exonerate gives them explicitly (no gap inference needed).
    introns: list[dict] = []
    for i, (g_start, g_end) in enumerate(_segments(model["introns"]), start=1):
        seq = _oriented(contig_seq[g_start - 1:g_end]).upper()
        introns.append({
            "intron_num": i,
            "g_start": g_start, "g_end": g_end,
            "length": g_end - g_start + 1,
            "sequence": seq,
            "splice_5": seq[:2] if len(seq) >= 2 else "",
            "splice_3": seq[-2:] if len(seq) >= 2 else "",
        })

    qc = validate_cds(cds_seq, geneticcode=geneticcode)
    q_cov = 0.0
    if model.get("query_len"):
        q_cov = (model["q_aln_end"] - model["q_aln_begin"]) / model["query_len"] * 100.0

    return {
        "strain_id": strain_id,
        "locus_name": locus_name,
        "contig_id": model["contig"],
        "strand": "-" if is_minus else "+",
        "blast_type": model.get("model", "exonerate"),
        "ref_accession": ref_accession or model.get("query_id", ""),
        "amplicon_start": g_min,
        "amplicon_end": g_max,
        "cds_seq": cds_seq,
        "amplicon_seq": amplicon,
        "exons": exons,
        "introns": introns,
        "n_exons": len(exons),
        "n_introns": len(introns),
        "cds_length": len(cds_seq),
        "amplicon_length": len(amplicon),
        # ── Exonerate / QC extras ──
        "tool": "exonerate",
        "tool_version": exonerate_version,
        "exonerate_score": model.get("score", 0.0),
        "query_pident": model.get("pident", 0.0),
        "query_coverage": round(q_cov, 1),
        "protein": qc["protein"],
        "n_internal_stops": qc["n_internal_stops"],
        "len_mod3": qc["len_mod3"],
        "qc_clean": qc["clean"],
        "tcs_matches": (model.get("tcs", "") == cds_seq) if model.get("tcs") else None,
        "n_other_models": n_other_models,
    }


# ── Output writers ───────────────────────────────────────────────────────────

def write_exonerate_fastas(result: dict, output_dir: str,
                           extracted_at: Optional[str] = None) -> dict:
    """
    Write CDS, translated protein, genomic span and introns with rich
    Exonerate-appropriate headers (model / score / %id / coverage rather than the
    BLAST e-value). Returns {file_type: path}.
    """
    os.makedirs(output_dir, exist_ok=True)
    sid, locus = result["strain_id"], result["locus_name"]
    ts = extracted_at or _now_utc()
    strand, contig = result["strand"], result["contig_id"]
    coords = f"{result['amplicon_start']}..{result['amplicon_end']}"
    paths: dict[str, str] = {}

    cds_desc = _build_description({
        "sample": sid, "locus": locus, "type": "CDS", "contig": contig,
        "coords": coords, "strand": strand, "n_exons": result["n_exons"],
        "cds_len": f"{result['cds_length']}bp", "ref": result.get("ref_accession"),
        "tool": "exonerate", "model": result.get("blast_type"),
        "score": result.get("exonerate_score"),
        "pid": f"{result.get('query_pident', 0.0):.1f}",
        "qcov": f"{result.get('query_coverage', 0.0):.1f}",
        "internal_stops": result.get("n_internal_stops"),
        "extracted": ts,
    })
    paths["cds"] = os.path.join(output_dir, f"{locus}_CDS.fasta")
    SeqIO.write([SeqRecord(Seq(result["cds_seq"]), id=f"{sid}_{locus}_CDS",
                           description=cds_desc)], paths["cds"], "fasta")

    protein = result.get("protein", "")
    if protein:
        prot_desc = _build_description({
            "sample": sid, "locus": locus, "type": "protein", "contig": contig,
            "coords": coords, "strand": strand,
            "aa_len": f"{len(protein.rstrip('*'))}aa",
            "tool": "exonerate", "model": result.get("blast_type"),
            "internal_stops": result.get("n_internal_stops"), "extracted": ts,
        })
        paths["protein"] = os.path.join(output_dir, f"{locus}_protein.fasta")
        SeqIO.write([SeqRecord(Seq(protein), id=f"{sid}_{locus}_protein",
                               description=prot_desc)], paths["protein"], "fasta")

    gen_desc = _build_description({
        "sample": sid, "locus": locus, "type": "genomic", "contig": contig,
        "coords": coords, "strand": strand,
        "amplicon_len": f"{result['amplicon_length']}bp",
        "n_introns": result["n_introns"], "tool": "exonerate", "extracted": ts,
    })
    paths["genomic"] = os.path.join(output_dir, f"{locus}_genomic.fasta")
    SeqIO.write([SeqRecord(Seq(result["amplicon_seq"]), id=f"{sid}_{locus}_genomic",
                           description=gen_desc)], paths["genomic"], "fasta")

    if result["introns"]:
        recs = []
        for intr in result["introns"]:
            desc = _build_description({
                "sample": sid, "locus": locus, "type": "intron",
                "intron_num": intr["intron_num"],
                "coords": f"{intr['g_start']}..{intr['g_end']}", "strand": strand,
                "length": f"{intr['length']}bp",
                "splice5": intr["splice_5"], "splice3": intr["splice_3"],
                "tool": "exonerate", "extracted": ts,
            })
            recs.append(SeqRecord(Seq(intr["sequence"]),
                                  id=f"{sid}_{locus}_intron{intr['intron_num']}",
                                  description=desc))
        paths["introns"] = os.path.join(output_dir, f"{locus}_introns.fasta")
        SeqIO.write(recs, paths["introns"], "fasta")

    return paths


def write_exonerate_log(result: dict, output_dir: str,
                        *, exonerate_cmd: str = "", reference_fasta: str = "",
                        narrowed: bool = True, run_dir: str = "") -> str:
    """Self-contained per-locus log: model, command, QC verdict, run provenance."""
    sid, locus = result["strain_id"], result["locus_name"]
    n_introns = result["n_introns"]
    splices_canonical = sum(
        1 for intr in result["introns"]
        if intr["splice_5"] == "GT" and intr["splice_3"] == "AG"
    )
    splice_note = (f"{splices_canonical}/{n_introns} canonical GT-AG"
                   if n_introns > 0 else "n/a (no introns)")
    frame_ok = result.get("len_mod3", 0) == 0
    stops = result.get("n_internal_stops", 0)
    qc_verdict = "PASS" if (frame_ok and stops == 0) else "REVIEW"
    frame_note = ("len % 3 == 0 (OK)" if frame_ok
                  else f"len % 3 == {result.get('len_mod3', 0)} (FRAMESHIFT?)")

    lines = [
        "# phylofetch Exonerate extraction log",
        f"# Generated: {_now_utc()}",
        "",
        f"[sample]        {sid}",
        f"[locus]         {locus}",
        f"[type]          CDS (spliced alignment)",
        f"[model]         {result.get('blast_type', 'exonerate')}",
        f"[tool_version]  {result.get('tool_version', 'unknown')}",
        "",
        f"[reference]     {result.get('ref_accession', 'unknown')}",
        f"[ref_fasta]     {reference_fasta or '(not recorded)'}",
        f"[narrowed]      {'tblastn/blastn best contig' if narrowed else 'whole assembly'}",
        f"[command]       {exonerate_cmd or '(not recorded)'}",
        "",
        "[result]",
        f"  contig        {result['contig_id']}",
        f"  coords        {result['amplicon_start']}..{result['amplicon_end']}",
        f"  strand        {result['strand']}",
        f"  n_exons       {result['n_exons']}",
        f"  cds_len       {result['cds_length']} bp",
        f"  amplicon_len  {result['amplicon_length']} bp",
        f"  n_introns     {result['n_introns']}",
        f"  splice_sites  {splice_note}",
        f"  score         {result.get('exonerate_score', 0.0)}",
        f"  query_pident  {result.get('query_pident', 0.0)}",
        f"  query_cov     {result.get('query_coverage', 0.0)} %",
        f"  paralog_cands {result.get('n_other_models', 0)} other model(s)",
        "",
        "[QC]",
        f"  geneticcode    {result.get('geneticcode', 1)}",
        f"  reading_frame {frame_note}",
        f"  internal_stops {stops}",
        f"  tcs_crosscheck {result.get('tcs_matches')}",
        f"  refine_used    {result.get('refine_used', 'none')}",
        f"  verdict        {qc_verdict}",
        "",
    ]
    if run_dir:
        lines += [
            f"[run_folder]    {run_dir}",
            f"[terminal_log]  {os.path.join(run_dir, 'terminal.log')}",
        ]

    log_path = os.path.join(output_dir, f"{locus}_extraction.log")
    Path(log_path).write_text("\n".join(lines) + "\n", encoding="utf-8")
    return log_path


# ── Top-level pipeline ───────────────────────────────────────────────────────

def extract_locus_exonerate(
    assembly_fasta: str,
    reference_fasta: str,
    output_dir: str,
    strain_id: str,
    locus_name: str,
    *,
    exonerate_bin: str = "exonerate",
    blastn_bin: str = "blastn",
    tblastn_bin: str = "tblastn",
    narrow: bool = True,
    minintron: int = 20,
    maxintron: int = 2000,
    bestn: int = 1,
    score: Optional[float] = None,
    refine: str = "none",
    geneticcode: int = 1,
    min_pident: float = 70.0,
    evalue: float = 1e-20,
    threads: int = 4,
    strict_qc: bool = False,
    escalate_refine: bool = True,
    escalate_ceiling: str = "region",
    manager=None,
    run_dir: str = "",
) -> tuple[Optional[dict], str]:
    """
    Hybrid BLAST-narrow → Exonerate-refine locus extraction.

    Returns ``(result_or_None, status)``. On success writes CDS / protein /
    genomic / introns FASTAs, GFF3, codon partition and a self-contained log to
    ``output_dir``. ``status`` starts with "OK" on success (with a "(QC: …)"
    suffix when the CDS frame/stop QC warrants review).

    ``strict_qc=True`` rejects a model whose CDS has internal stop codons or is
    not a multiple of three. Default is False: imperfect models are still written
    but flagged, consistent with D-007 (partial/type sequences are not silently
    dropped).
    """
    os.makedirs(output_dir, exist_ok=True)

    if not shutil.which(exonerate_bin):
        return None, f"ERROR: exonerate not found on PATH ({exonerate_bin})"

    ref_type = detect_fasta_type(reference_fasta)
    model = MODEL_FOR_QUERYTYPE.get(ref_type, "protein2genome")
    narrow_bin = tblastn_bin if ref_type == "protein" else blastn_bin

    # ── Step 1: narrow to the best contig (fast) or use the whole assembly ──
    target_fasta = assembly_fasta
    ref_acc = ""
    contig_seq_cache: dict[str, str] = {}
    narrowed = False

    narrow_warning = ""
    if narrow and shutil.which(narrow_bin):
        rc, _stderr, tsv = run_blast(
            reference_fasta, assembly_fasta, output_dir,
            blast_bin=narrow_bin, task=None, evalue=evalue, threads=threads,
            manager=manager, module="blast", action=f"narrow_locus:{locus_name}",
        )
        if rc == 0:
            best, ref_acc = select_best_locus_group(
                parse_blast_hsps(tsv), min_pident=min_pident)
            if best:
                contig_id = best[0]["sseqid"]
                contig_seq = _load_contig(assembly_fasta, contig_id)
                if contig_seq is not None:
                    contig_seq_cache[contig_id] = contig_seq
                    target_fasta = os.path.join(output_dir, "_exonerate_target.fasta")
                    SeqIO.write([SeqRecord(Seq(contig_seq), id=contig_id,
                                           description="")], target_fasta, "fasta")
                    narrowed = True
        else:
            # Don't silently fall back: a BLAST *error* (vs. a genuine no-hit) means Exonerate
            # runs on the WHOLE assembly — slower, and it can pick a paralog. Flag it so the
            # riskier path is a visible decision, not a hidden one. (D-033)
            narrow_warning = (f" [WARN: narrowing {narrow_bin} failed (exit {rc}) — "
                              f"ran Exonerate on whole assembly]")
    # If narrowing produced nothing, target_fasta stays = whole assembly.

    # ── Steps 2-3: run Exonerate + build result, escalating boundary refinement when the CDS
    # comes out frameshifted / with internal stops (D-030). Exonerate's protein2genome can misplace
    # a splice boundary by a base or two on some sequences (a couple of substitutions near a splice
    # site are enough), frameshifting the reconstructed CDS even though the genomic DNA is clean.
    # `--refine region`/`full` re-optimizes the boundaries and recovers the frame in most such cases.
    # The first pass uses the requested `refine`; only if that isn't clean (and `escalate_refine`)
    # do we retry at higher refinement, keeping the cleanest result — so the common clean case stays
    # a single fast pass. Escalation stops at `escalate_ceiling` (default "region"): `--refine full`
    # re-runs *exhaustive* DP over the whole narrowed contig (minutes on a multi-Mb contig / long
    # multi-intron gene like RPB2) and — verified on the aff. eureka strains — rescued nothing that
    # `region` didn't, so it is opt-in (Thorough effort / the page's targeted deep-refine pass). (D-035)
    REFINE_ORDER = ["none", "region", "full"]
    start = REFINE_ORDER.index(refine) if refine in REFINE_ORDER else 0
    if escalate_refine:
        ceil = REFINE_ORDER.index(escalate_ceiling) if escalate_ceiling in REFINE_ORDER else 1
        levels = REFINE_ORDER[start:max(ceil, start) + 1]
    else:
        levels = [refine]

    def _run_build(refine_level: str):
        rc, raw, _gff = run_exonerate(
            reference_fasta, target_fasta, model, output_dir,
            exonerate_bin=exonerate_bin, minintron=minintron, maxintron=maxintron,
            bestn=bestn, score=score, refine=refine_level, geneticcode=geneticcode,
            manager=manager, action=f"exonerate_{locus_name}_{strain_id}",
        )
        if rc != 0:
            return None, f"Exonerate failed (exit {rc})"
        ms = parse_exonerate_gff(raw)
        if not ms:
            return None, "Exonerate found no gene model (try whole-assembly target or lower thresholds)"
        bm, others = select_best_model(ms, min_score=(score or 0.0))
        if bm is None:
            return None, f"No Exonerate model scored ≥ {score}"
        cseq = contig_seq_cache.get(bm["contig"]) or _load_contig(assembly_fasta, bm["contig"])
        if cseq is None:
            return None, f"Contig '{bm['contig']}' not found in assembly"
        res = build_result_from_model(
            bm, cseq, strain_id, locus_name,
            ref_accession=ref_acc, exonerate_version=_probe_exonerate_version(exonerate_bin),
            geneticcode=geneticcode, n_other_models=len(others),
        )
        return res, "OK"

    result = None
    best_q = None
    used_refine = refine
    last_err = "Exonerate found no gene model"
    for lvl in levels:
        res, st = _run_build(lvl)
        if res is None:
            last_err = st
            continue
        q = (res["len_mod3"] == 0, -res["n_internal_stops"])   # frame-OK first, then fewer stops
        if best_q is None or q > best_q:
            result, best_q, used_refine = res, q, lvl
        if res["n_internal_stops"] == 0 and res["len_mod3"] == 0:
            break                                              # clean — stop escalating
    if result is None:
        return None, last_err
    result["refine_used"] = used_refine
    result["refine_escalated"] = used_refine != refine
    result["geneticcode"] = geneticcode   # recorded so a later re-validation uses the SAME table

    if strict_qc and (result["n_internal_stops"] > 0 or result["len_mod3"] != 0):
        return None, (
            f"strict QC failed: {result['n_internal_stops']} internal stop(s), "
            f"len % 3 = {result['len_mod3']} (disable strict QC to keep flagged output)"
        )

    # ── Step 4: write outputs (reuse GFF3 + codon-partition writers) ──
    ts = _now_utc()
    exo_cmd = (
        f"{exonerate_bin} --model {model} --query {reference_fasta} "
        f"--target {os.path.basename(target_fasta)} --showtargetgff yes "
        f"--minintron {minintron} --maxintron {maxintron} --bestn {bestn}"
    )
    write_exonerate_fastas(result, output_dir, extracted_at=ts)
    write_gff3(result, output_dir)
    write_region_gff3(result, output_dir)
    write_codon_partition(result["cds_length"], output_dir, locus_name)
    write_exonerate_log(result, output_dir, exonerate_cmd=exo_cmd,
                        reference_fasta=reference_fasta, narrowed=narrowed,
                        run_dir=run_dir)

    status = "OK"
    if result["n_internal_stops"] > 0 or result["len_mod3"] != 0:
        status = (f"OK (QC: review — {result['n_internal_stops']} internal stop(s), "
                  f"len % 3 = {result['len_mod3']})")
    if result.get("refine_escalated"):
        status += f" [boundary-refined: refine={used_refine}]"
    status += narrow_warning
    return result, status


def scan_flagged_cds(per_strain_dir: str, geneticcode: int = 1) -> list[dict]:
    """
    Stateless scan of extracted coding loci for a CDS that needs attention (frameshift / internal
    stop), so a targeted ``--refine full`` rescue pass can be offered at **any** time — including a
    future session, e.g. when downstream analysis fails — without relying on the run's session
    state. Reads each ``<strain>/<locus>/<locus>_CDS.fasta`` and re-validates its reading frame
    (`validate_cds`), independent of how it was produced (D-035).

    Returns one row per flagged ``(strain, locus)``::

        {strain, locus, locus_dir, cds_length, len_mod3, n_internal_stops, refine_used, geneticcode}

    The CDS is re-validated under the genetic code recorded in its extraction log — **not** the
    `geneticcode` argument — whenever the log records one, so a later change to the page's genetic-code
    widget can't by itself flip a clean CDS to flagged or vice versa (D-036 review). The argument is
    only the fallback for loci extracted before the code was logged. ``refine_used`` is parsed from
    the log when present, so the user can see which refinement level already ran; "" when unknown.
    Clean CDS are omitted.
    """
    rows: list[dict] = []
    root = Path(per_strain_dir)
    if not root.is_dir():
        return rows
    for strain_dir in sorted(p for p in root.iterdir() if p.is_dir()):
        for locus_dir in sorted(p for p in strain_dir.iterdir() if p.is_dir()):
            locus = locus_dir.name
            cds_fa = locus_dir / f"{locus}_CDS.fasta"
            if not cds_fa.exists():
                continue
            try:
                recs = list(SeqIO.parse(str(cds_fa), "fasta"))
            except Exception:
                continue
            if not recs:
                continue
            log = locus_dir / f"{locus}_extraction.log"
            log_text = log.read_text(errors="replace") if log.exists() else ""
            # Re-validate under the code the locus was ORIGINALLY extracted with (logged), not the
            # caller's current widget — otherwise a code change alone could flip the verdict.
            mgc = re.search(r"geneticcode[:= ]+(\d+)", log_text)
            gc = int(mgc.group(1)) if mgc else geneticcode
            qc = validate_cds(str(recs[0].seq), geneticcode=gc)
            if qc["len_mod3"] == 0 and qc["n_internal_stops"] == 0:
                continue   # clean — nothing to rescue
            m = re.search(r"refine[_ ]?used[:= ]+(\w+)|boundary-refined: refine=(\w+)", log_text)
            refine_used = (m.group(1) or m.group(2) or "") if m else ""
            rows.append({
                "strain": strain_dir.name, "locus": locus, "locus_dir": str(locus_dir),
                "cds_length": qc["length"], "len_mod3": qc["len_mod3"],
                "n_internal_stops": qc["n_internal_stops"], "refine_used": refine_used,
                "geneticcode": gc,
            })
    return rows


def resolve_guide_path(locus: str, project_dir: str, goi_dir, ref_dir,
                       is_goi: bool = False) -> Optional[str]:
    """
    Find the guide FASTA that would extract ``locus``, for a deep-refine re-run (D-035): the saved
    bundled/both/protein-library query (``scratch/guides``), a gene-of-interest ref
    (``_goi_refs``), or the nucleotide library FASTA — whichever exists, with **no directory side
    effects**. For a gene of interest, the user's uploaded ortholog takes precedence over a
    same-named catalogue guide (D-036 review: a GOI named e.g. "RPB2" must not silently resolve to
    the bundled RPB2 guide). Returns ``None`` when no guide is on disk.
    """
    goi = Path(goi_dir) / f"{locus}.fasta"
    scratch = Path(project_dir) / "scratch" / "guides" / f"{locus}_guide.fasta"
    lib = Path(ref_dir) / locus / f"{locus}_refs.fasta"
    order = [goi, scratch, lib] if is_goi else [scratch, goi, lib]
    for c in order:
        if c.exists() and c.stat().st_size > 0:
            return str(c)
    return None
