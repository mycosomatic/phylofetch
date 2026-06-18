"""
blast_loci_utils.py
-------------------
BLAST-based locus extraction with HSP-as-exon parsing.

Core insight: when a CDS query is BLASTed against a genomic target, BLAST
splits the alignment at intron boundaries because the query has no intronic
sequence. Each HSP corresponds to one exon; gaps between consecutive HSPs
(in target coordinates, ordered by query position) are introns.

LXD-002 fix: HSPs are grouped first by qseqid (reference accession), then by
(sseqid, strand). A single reference accession is chosen as the best query
before selecting the best contig/strand group. This prevents HSPs from
different reference accessions being stitched into one false CDS.

Outputs per locus per strain:
  LOCUS_CDS.fasta        — intron-free CDS with rich provenance header
  LOCUS_genomic.fasta    — amplicon-style (exons + introns) with provenance
  LOCUS_introns.fasta    — individual intron sequences
  LOCUS.gff3             — exon/CDS annotation with proper phase
  LOCUS_partition.nex    — codon position partitions for IQ-TREE / RAxML
  LOCUS_extraction.log   — self-contained log: command, version, result
"""

import os
import shutil
import subprocess
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Optional

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


# ── Reference type detection ─────────────────────────────────────────────────

def detect_fasta_type(fasta_path: str) -> str:
    """Return 'nucleotide' or 'protein' by inspecting the first sequence."""
    nuc_chars = set("ACGTUNWSMKRYBDHVacgtunwsmkrybdhv-.")
    with open(fasta_path) as f:
        in_seq = False
        sample: list[str] = []
        for line in f:
            if line.startswith(">"):
                if in_seq and sample:
                    break
                in_seq = True
                continue
            if in_seq:
                sample.append(line.strip())
                if sum(len(s) for s in sample) >= 200:
                    break
    seq_text = "".join(sample)
    if not seq_text:
        return "nucleotide"
    nuc_count = sum(1 for c in seq_text if c in nuc_chars)
    return "nucleotide" if nuc_count / len(seq_text) > 0.9 else "protein"


# ── BLAST runner ─────────────────────────────────────────────────────────────

_OUTFMT_FIELDS = [
    "qseqid", "sseqid", "pident", "length", "mismatch",
    "gapopen", "qstart", "qend", "sstart", "send",
    "evalue", "bitscore", "sstrand", "qlen", "slen",
]
_OUTFMT_STR = "6 " + " ".join(_OUTFMT_FIELDS)


def run_blast(query_fasta: str, target_fasta: str, output_dir: str,
              blast_bin: str = "blastn",
              task: Optional[str] = None,
              evalue: float = 1e-20,
              threads: int = 4) -> tuple[int, str, str]:
    """Run blastn or tblastn. Returns (returncode, stderr, tsv_path)."""
    os.makedirs(output_dir, exist_ok=True)
    out_tsv = os.path.join(output_dir, "blast_hsps.tsv")

    cmd_parts = [
        blast_bin,
        "-query", query_fasta,
        "-subject", target_fasta,
        "-evalue", str(evalue),
        "-num_threads", str(threads),
        "-outfmt", _OUTFMT_STR,
        "-out", out_tsv,
    ]
    if blast_bin == "blastn" and task:
        cmd_parts += ["-task", task]

    result = subprocess.run(cmd_parts, capture_output=True, text=True)
    return result.returncode, result.stderr, out_tsv


def run_blast_alignment(query_fasta: str, target_fasta: str, output_dir: str,
                        blast_bin: str = "blastn",
                        task: Optional[str] = None,
                        evalue: float = 1e-20,
                        threads: int = 4) -> tuple[int, str, str]:
    """Run BLAST with standard pairwise alignment output (outfmt 0) for viewing."""
    os.makedirs(output_dir, exist_ok=True)
    out_txt = os.path.join(output_dir, "blast_alignment.txt")

    cmd_parts = [
        blast_bin,
        "-query", query_fasta,
        "-subject", target_fasta,
        "-evalue", str(evalue),
        "-num_threads", str(threads),
        "-outfmt", "0",
        "-out", out_txt,
    ]
    if blast_bin == "blastn" and task:
        cmd_parts += ["-task", task]

    result = subprocess.run(cmd_parts, capture_output=True, text=True)
    return result.returncode, result.stderr, out_txt


# ── HSP parsing ──────────────────────────────────────────────────────────────

def parse_blast_hsps(tsv_path: str) -> list[dict]:
    hsps = []
    if not os.path.exists(tsv_path):
        return hsps
    with open(tsv_path) as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 15:
                continue
            try:
                hsps.append({
                    "qseqid":   parts[0],
                    "sseqid":   parts[1],
                    "pident":   float(parts[2]),
                    "length":   int(parts[3]),
                    "qstart":   int(parts[6]),
                    "qend":     int(parts[7]),
                    "sstart":   int(parts[8]),
                    "send":     int(parts[9]),
                    "evalue":   float(parts[10]),
                    "bitscore": float(parts[11]),
                    "sstrand":  parts[12],
                    "qlen":     int(parts[13]),
                    "slen":     int(parts[14]),
                })
            except ValueError:
                continue
    return hsps


def select_best_locus_group(hsps: list[dict],
                             min_pident: float = 70.0) -> tuple[list[dict], str]:
    """
    Filter by % identity, then pick the best single reference accession
    (highest cumulative bitscore per qseqid), then pick the best
    (contig, strand) group from that reference's HSPs.

    Returns (sorted_hsps_in_gene_order, best_ref_accession).

    LXD-002: grouping by qseqid first prevents HSPs from different reference
    accessions from being stitched into one extracted CDS.
    """
    filtered = [h for h in hsps if h["pident"] >= min_pident]
    if not filtered:
        return [], ""

    # Step 1: best reference accession
    per_query: dict[str, list] = defaultdict(list)
    for h in filtered:
        per_query[h["qseqid"]].append(h)
    best_qseqid = max(per_query, key=lambda q: sum(x["bitscore"] for x in per_query[q]))
    query_hsps = per_query[best_qseqid]

    # Step 2: best (contig, strand) group within that reference
    per_target: dict[tuple, list] = defaultdict(list)
    for h in query_hsps:
        strand = "plus" if h["send"] >= h["sstart"] else "minus"
        per_target[(h["sseqid"], strand)].append(h)
    best_key = max(per_target, key=lambda k: sum(x["bitscore"] for x in per_target[k]))
    best_group = sorted(per_target[best_key], key=lambda h: h["qstart"])

    return best_group, best_qseqid


# ── Sequence extraction ───────────────────────────────────────────────────────

def _load_contig(fasta_path: str, contig_id: str) -> Optional[str]:
    for rec in SeqIO.parse(fasta_path, "fasta"):
        if rec.id == contig_id:
            return str(rec.seq)
    return None


def extract_from_hsps(hsps: list[dict], assembly_fasta: str,
                      strain_id: str, locus_name: str,
                      blast_type: str = "blastn",
                      ref_accession: str = "") -> Optional[dict]:
    """
    Extract CDS, genomic amplicon, and intron sequences from a set of exon HSPs.
    Returns a result dict, or None if the contig cannot be found.
    """
    if not hsps:
        return None

    contig_id = hsps[0]["sseqid"]
    is_minus  = hsps[0]["send"] < hsps[0]["sstart"]
    contig_seq = _load_contig(assembly_fasta, contig_id)
    if contig_seq is None:
        return None

    exons = []
    for h in hsps:
        g_start, g_end = (h["send"], h["sstart"]) if is_minus else (h["sstart"], h["send"])
        exons.append({
            "g_start": g_start, "g_end": g_end,
            "length":  g_end - g_start + 1,
            "pident":  h["pident"],
            "q_start": h["qstart"], "q_end": h["qend"],
        })

    # CDS: concatenated exons in gene order, reverse-comp if minus
    cds_segments = []
    for ex in exons:
        seg = contig_seq[ex["g_start"] - 1 : ex["g_end"]]
        if is_minus:
            seg = str(Seq(seg).reverse_complement())
        cds_segments.append(seg)
    cds_seq = "".join(cds_segments).upper()

    # Genomic amplicon: first exon start → last exon end
    g_min = min(ex["g_start"] for ex in exons)
    g_max = max(ex["g_end"]   for ex in exons)
    amplicon = contig_seq[g_min - 1 : g_max]
    if is_minus:
        amplicon = str(Seq(amplicon).reverse_complement())
    amplicon = amplicon.upper()

    # Introns: gaps between consecutive exons in gene order
    introns = []
    for i in range(len(exons) - 1):
        if is_minus:
            intron_g_start = exons[i + 1]["g_end"] + 1
            intron_g_end   = exons[i]["g_start"] - 1
        else:
            intron_g_start = exons[i]["g_end"] + 1
            intron_g_end   = exons[i + 1]["g_start"] - 1

        if intron_g_end < intron_g_start:
            continue

        intron_seq = contig_seq[intron_g_start - 1 : intron_g_end]
        if is_minus:
            intron_seq = str(Seq(intron_seq).reverse_complement())
        introns.append({
            "intron_num": i + 1,
            "g_start":    intron_g_start,
            "g_end":      intron_g_end,
            "length":     intron_g_end - intron_g_start + 1,
            "sequence":   intron_seq.upper(),
            "splice_5":   intron_seq[:2].upper() if len(intron_seq) >= 2 else "",
            "splice_3":   intron_seq[-2:].upper() if len(intron_seq) >= 2 else "",
        })

    return {
        "strain_id":       strain_id,
        "locus_name":      locus_name,
        "contig_id":       contig_id,
        "strand":          "-" if is_minus else "+",
        "blast_type":      blast_type,
        "ref_accession":   ref_accession,
        "amplicon_start":  g_min,
        "amplicon_end":    g_max,
        "cds_seq":         cds_seq,
        "amplicon_seq":    amplicon,
        "exons":           exons,
        "introns":         introns,
        "n_exons":         len(exons),
        "n_introns":       len(introns),
        "cds_length":      len(cds_seq),
        "amplicon_length": len(amplicon),
    }


# ── Rich FASTA header builder ─────────────────────────────────────────────────

def _now_utc() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat()


def _build_description(fields: dict) -> str:
    """Build NCBI-style bracket key=value description for FASTA headers."""
    return " ".join(f"[{k}={v}]" for k, v in fields.items() if v is not None)


# ── Output writers ────────────────────────────────────────────────────────────

def write_fastas(result: dict, output_dir: str,
                 evalue: float = 1e-20,
                 min_pident: float = 70.0,
                 extracted_at: Optional[str] = None) -> dict:
    """
    Write CDS, genomic amplicon, and introns to FASTA files with rich headers.
    Returns {file_type: path}.
    """
    os.makedirs(output_dir, exist_ok=True)
    sid    = result["strain_id"]
    locus  = result["locus_name"]
    ts     = extracted_at or _now_utc()
    ref    = result.get("ref_accession", "")
    strand = result["strand"]
    contig = result["contig_id"]
    pident = max((ex["pident"] for ex in result["exons"]), default=0.0)

    paths = {}

    # CDS
    cds_desc = _build_description({
        "sample": sid, "locus": locus, "type": "CDS",
        "contig": contig,
        "coords": f"{result['amplicon_start']}..{result['amplicon_end']}",
        "strand": strand,
        "n_exons": result["n_exons"], "cds_len": f"{result['cds_length']}bp",
        "ref_acc": ref, "pident": f"{pident:.1f}", "evalue": evalue,
        "extracted": ts,
    })
    paths["cds"] = os.path.join(output_dir, f"{locus}_CDS.fasta")
    SeqIO.write([SeqRecord(
        Seq(result["cds_seq"]),
        id=f"{sid}_{locus}_CDS",
        description=cds_desc,
    )], paths["cds"], "fasta")

    # Genomic amplicon
    gen_desc = _build_description({
        "sample": sid, "locus": locus, "type": "genomic",
        "contig": contig,
        "coords": f"{result['amplicon_start']}..{result['amplicon_end']}",
        "strand": strand,
        "amplicon_len": f"{result['amplicon_length']}bp",
        "n_introns": result["n_introns"],
        "extracted": ts,
    })
    paths["genomic"] = os.path.join(output_dir, f"{locus}_genomic.fasta")
    SeqIO.write([SeqRecord(
        Seq(result["amplicon_seq"]),
        id=f"{sid}_{locus}_genomic",
        description=gen_desc,
    )], paths["genomic"], "fasta")

    # Introns
    if result["introns"]:
        intron_recs = []
        for intr in result["introns"]:
            desc = _build_description({
                "sample": sid, "locus": locus, "type": "intron",
                "intron_num": intr["intron_num"],
                "coords": f"{intr['g_start']}..{intr['g_end']}",
                "strand": strand,
                "length": f"{intr['length']}bp",
                "splice5": intr["splice_5"], "splice3": intr["splice_3"],
                "extracted": ts,
            })
            intron_recs.append(SeqRecord(
                Seq(intr["sequence"]),
                id=f"{sid}_{locus}_intron{intr['intron_num']}",
                description=desc,
            ))
        paths["introns"] = os.path.join(output_dir, f"{locus}_introns.fasta")
        SeqIO.write(intron_recs, paths["introns"], "fasta")

    return paths


def write_gff3(result: dict, output_dir: str) -> str:
    sid    = result["strain_id"]
    locus  = result["locus_name"]
    contig = result["contig_id"]
    strand = result["strand"]
    exons  = result["exons"]
    gene_id = f"{sid}_{locus}"
    mrna_id = f"{gene_id}.t1"

    lines = ["##gff-version 3",
             f"##sequence-region {contig} {result['amplicon_start']} {result['amplicon_end']}",
             f"{contig}\tblast_loci\tgene\t{result['amplicon_start']}\t{result['amplicon_end']}"
             f"\t.\t{strand}\t.\tID={gene_id};Name={locus}",
             f"{contig}\tblast_loci\tmRNA\t{result['amplicon_start']}\t{result['amplicon_end']}"
             f"\t.\t{strand}\t.\tID={mrna_id};Parent={gene_id}"]

    cumulative = 0
    phase_by_idx: dict[int, int] = {}
    for i, ex in enumerate(exons):
        phase_by_idx[i] = (3 - (cumulative % 3)) % 3
        cumulative += ex["length"]

    for gff_order, (gene_idx, ex) in enumerate(
        sorted(enumerate(exons), key=lambda t: t[1]["g_start"]), start=1
    ):
        lines.append(
            f"{contig}\tblast_loci\texon\t{ex['g_start']}\t{ex['g_end']}"
            f"\t{ex['pident']:.1f}\t{strand}\t.\t"
            f"ID={mrna_id}.exon{gff_order};Parent={mrna_id}"
        )
        lines.append(
            f"{contig}\tblast_loci\tCDS\t{ex['g_start']}\t{ex['g_end']}"
            f"\t{ex['pident']:.1f}\t{strand}\t{phase_by_idx[gene_idx]}\t"
            f"ID={mrna_id}.cds;Parent={mrna_id}"
        )

    gff_path = os.path.join(output_dir, f"{locus}.gff3")
    Path(gff_path).write_text("\n".join(lines) + "\n")
    return gff_path


def write_codon_partition(cds_length: int, output_dir: str, locus_name: str) -> str:
    path = os.path.join(output_dir, f"{locus_name}_partition.nex")
    content = f"""#nexus
[ Codon position partition for {locus_name} ]
[ IQ-TREE: iqtree2 -s {locus_name}_CDS.fasta -p {locus_name}_partition.nex -m MFP+MERGE -B 1000 ]
begin sets;
    charset {locus_name}_pos1 = 1-{cds_length}\\3;
    charset {locus_name}_pos2 = 2-{cds_length}\\3;
    charset {locus_name}_pos3 = 3-{cds_length}\\3;
end;
"""
    Path(path).write_text(content)
    return path


def write_extraction_log(result: dict, output_dir: str,
                         blast_cmd: str = "",
                         blast_version: str = "",
                         evalue: float = 1e-20,
                         min_pident: float = 70.0,
                         run_dir: str = "") -> str:
    """
    Write a self-contained extraction log alongside the FASTA files.
    Every output directory can be understood independently of the project workspace.
    """
    ts = _now_utc()
    sid   = result["strain_id"]
    locus = result["locus_name"]

    splices_canonical = sum(
        1 for intr in result["introns"]
        if intr["splice_5"] == "GT" and intr["splice_3"] == "AG"
    )
    n_introns = result["n_introns"]
    splice_note = (
        f"{splices_canonical}/{n_introns} canonical GT-AG"
        if n_introns > 0 else "n/a (no introns)"
    )

    lines = [
        "# phylofetch extraction log",
        f"# Generated: {ts}",
        "",
        f"[sample]       {sid}",
        f"[locus]        {locus}",
        f"[type]         CDS ({result['blast_type']})",
        "",
        f"[command]      {blast_cmd or '(not recorded)'}",
        f"[tool_version] {blast_version or 'unknown'}",
        f"[reference]    {result.get('ref_accession', 'unknown')}",
        f"[e_value]      {evalue}",
        f"[min_pident]   {min_pident}",
        "",
        "[result]",
        f"  contig       {result['contig_id']}",
        f"  coords       {result['amplicon_start']}..{result['amplicon_end']}",
        f"  strand       {result['strand']}",
        f"  n_exons      {result['n_exons']}",
        f"  cds_len      {result['cds_length']} bp",
        f"  amplicon_len {result['amplicon_length']} bp",
        f"  n_introns    {result['n_introns']}",
        f"  splice_sites {splice_note}",
        "",
    ]
    if run_dir:
        lines += [
            f"[run_folder]   {run_dir}",
            f"[terminal_log] {os.path.join(run_dir, 'terminal.log')}",
            f"[command_json] {os.path.join(run_dir, 'command.json')}",
        ]

    log_path = os.path.join(output_dir, f"{locus}_extraction.log")
    Path(log_path).write_text("\n".join(lines) + "\n")
    return log_path


# ── Length acceptability ──────────────────────────────────────────────────────

def min_acceptable_cds_length(reference_fasta: str, blast_type: str,
                              min_cds_pct_of_ref: float) -> float:
    """
    Minimum CDS length (bp) for a hit to count as 'complete enough', based on the
    shortest reference. tblastn references are in amino acids → ×3 for bp.
    Returns 0.0 if no references can be read.
    """
    ref_lengths = [len(r.seq) for r in SeqIO.parse(reference_fasta, "fasta")]
    if not ref_lengths:
        return 0.0
    scale = 3 if blast_type == "tblastn" else 1
    return min(ref_lengths) * scale * (min_cds_pct_of_ref / 100)


# ── Top-level pipeline ────────────────────────────────────────────────────────

def extract_locus(
    assembly_fasta: str,
    reference_fasta: str,
    output_dir: str,
    strain_id: str,
    locus_name: str,
    min_pident: float = 70.0,
    min_cds_pct_of_ref: float = 50.0,
    evalue: float = 1e-20,
    blastn_task: str = "dc-megablast",
    threads: int = 4,
    blastn_bin: str = "blastn",
    tblastn_bin: str = "tblastn",
    run_dir: str = "",
    require_complete_cds: bool = True,
) -> tuple[Optional[dict], str]:
    """
    Top-level locus extraction. Returns (result_dict_or_None, status_message).
    Also writes FASTA files, GFF3, codon partition, and extraction log to output_dir.

    require_complete_cds: when True (coding/ortholog strategy), reject hits whose
    stitched CDS is shorter than ``min_cds_pct_of_ref`` of the reference. Set False
    for PCR-amplicon references (the usual NCBI nucleotide case), where partial /
    intron-containing matches are expected and the genomic amplicon output is used.
    """
    os.makedirs(output_dir, exist_ok=True)

    ref_type = detect_fasta_type(reference_fasta)
    if ref_type == "protein":
        blast_bin  = tblastn_bin
        blast_type = "tblastn"
        task       = None
    else:
        blast_bin  = blastn_bin
        blast_type = "blastn"
        task       = blastn_task

    if not shutil.which(blast_bin):
        return None, f"ERROR: {blast_bin} not found on PATH"

    rc, stderr, tsv = run_blast(
        reference_fasta, assembly_fasta, output_dir,
        blast_bin=blast_bin, task=task, evalue=evalue, threads=threads,
    )
    if rc != 0:
        return None, f"BLAST failed (exit {rc}): {stderr.strip()[:200]}"

    hsps = parse_blast_hsps(tsv)
    if not hsps:
        return None, "No BLAST hits at e-value threshold"

    best, ref_acc = select_best_locus_group(hsps, min_pident=min_pident)
    if not best:
        return None, f"No HSPs ≥ {min_pident}% identity"

    result = extract_from_hsps(
        best, assembly_fasta, strain_id, locus_name,
        blast_type=blast_type, ref_accession=ref_acc,
    )
    if result is None:
        return None, "Sequence extraction failed (contig not found)"

    # CDS-completeness gate — only enforced for the coding/ortholog strategy.
    # PCR-amplicon references (require_complete_cds=False) are expected to be
    # partial or intron-containing, so the gate is skipped and the genomic
    # amplicon output is the primary product.
    if require_complete_cds:
        min_acceptable = min_acceptable_cds_length(
            reference_fasta, blast_type, min_cds_pct_of_ref)
        if result["cds_length"] < min_acceptable:
            return None, (
                f"CDS {result['cds_length']} bp < {min_acceptable:.0f} bp threshold "
                f"(uncheck 'require complete CDS' for PCR-amplicon refs, or lower min identity)"
            )

    ts = _now_utc()
    blast_cmd = (
        f"{blast_bin} -query {reference_fasta} -subject {assembly_fasta} "
        f"-evalue {evalue} -outfmt '{_OUTFMT_STR}'"
        + (f" -task {task}" if task else "")
    )

    write_fastas(result, output_dir, evalue=evalue, min_pident=min_pident, extracted_at=ts)
    write_gff3(result, output_dir)
    write_codon_partition(result["cds_length"], output_dir, locus_name)
    write_extraction_log(result, output_dir,
                         blast_cmd=blast_cmd, evalue=evalue,
                         min_pident=min_pident, run_dir=run_dir)

    return result, "OK"


# ── Combined multi-FASTA merger ───────────────────────────────────────────────

def merge_per_strain_outputs(per_strain_dir: str, combined_dir: str,
                             locus_name: str) -> dict[str, str]:
    """Concatenate per-strain FASTA files into combined multi-FASTAs."""
    os.makedirs(combined_dir, exist_ok=True)
    outputs: dict[str, list] = {"CDS": [], "genomic": [], "introns": []}

    for strain_dir in sorted(Path(per_strain_dir).iterdir()):
        if not strain_dir.is_dir():
            continue
        locus_dir = strain_dir / locus_name
        if not locus_dir.is_dir():
            continue
        for key, suffix in [
            ("CDS",     f"{locus_name}_CDS.fasta"),
            ("genomic", f"{locus_name}_genomic.fasta"),
            ("introns", f"{locus_name}_introns.fasta"),
        ]:
            fp = locus_dir / suffix
            if fp.exists():
                outputs[key].extend(list(SeqIO.parse(str(fp), "fasta")))

    combined_paths: dict[str, str] = {}
    for key, recs in outputs.items():
        if recs:
            path = os.path.join(combined_dir, f"{locus_name}_{key}_combined.fasta")
            SeqIO.write(recs, path, "fasta")
            combined_paths[key] = path
    return combined_paths
