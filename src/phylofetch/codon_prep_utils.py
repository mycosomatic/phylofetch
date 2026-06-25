"""
codon_prep_utils.py
-------------------
RM-008 component 2 (D-022): turn coding-locus **amplicon** sequences — the
comparison tips pulled from NCBI (and, by the same logic, any genomic amplicon) —
into frame-consistent products ready for a *hand-checked* codon alignment. This
is the step a researcher used to do by hand in Mesquite.

Why this exists
===============
Your own isolates' coding loci come out of the Exonerate extraction path (D-008)
already clean: intron-free, frame-pinned, strand-correct CDS. But the GenBank
comparison tips (component 3, ``tips_utils``) are raw nucleotide barcode
amplicons — for fungal coding markers these are usually *partial-CDS genomic*
records that carry **introns**, sit in an **arbitrary reading frame**, and may be
on **either strand** (the D-017 finding). You cannot codon-align those against
intron-free isolate CDS as-is.

The fix is to run each tip amplicon through Exonerate against the **same bundled
protein guide** (component 1, D-020) with ``protein2genome``: that strips the
introns, pins the reading frame to the guide ORF, and orients to the coding
strand — exactly the manual-Mesquite work, automated and logged.

What it produces (per coding locus, isolates + tips together)
=============================================================
  <locus>_CDS_combined.fasta      intron-stripped, codon-phased CDS
                                  → MACSE / hand-alignment (codon-aware)
  <locus>_genomic_combined.fasta  the full gene (exons + introns, oriented),
                                  exons UPPERCASE / introns lowercase so the
                                  exon-intron boundaries stay visible — and move
                                  with the gaps — while you hand-align in
                                  AliView / Geneious / SeaView / a text editor
  <locus>_protein_combined.fasta  the translation (AA tree / hand-align guide)

Design constraints (from the user)
==================================
* **Runs no aligner** and adds **no new external dependency** beyond Exonerate
  (already required for extraction). Alignment stays a separate, hand-curated
  step; MACSE / AliView / Geneious are optional companions, never dependencies.
* **rDNA loci (ITS/LSU/SSU) are out of scope** here — they are non-coding, have
  no protein guide, and go straight to MAFFT.
* QC is **write-and-flag** by default (D-007/D-008); ``strict_qc`` excludes
  frameshift / internal-stop CDS.
* **Nucleotide fallback (D-027):** some standard coding barcodes — fungal TEF1
  (EF1-728F/986R) above all — amplify a largely **intronic** region: the tip
  matches the genomic gene continuously but has almost no exon, so
  ``protein2genome`` finds no model and it cannot be codon-framed. Rather than
  drop it, such a tip is oriented (blastn vs the isolates' genomic, which also
  confirms the locus) and written to the **genomic** matrix only, flagged
  ``nucleotide_only`` — so the intron-inclusive tree still carries the
  comparison taxon. The published phylogenies for these markers align the
  introns too; the CDS/protein matrices stay isolate-only for that locus. A tip
  that orients to nothing (wrong locus / contamination) is still reported, not
  included.
"""

from __future__ import annotations

import os
import tempfile
from pathlib import Path
from typing import Optional

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from .blast_loci_utils import (
    _build_description,
    detect_fasta_type,
    parse_blast_hsps,
    run_blast,
)
from .exonerate_utils import (
    MODEL_FOR_QUERYTYPE,
    build_result_from_model,
    parse_exonerate_gff,
    run_exonerate,
    select_best_model,
    _probe_exonerate_version,
)
from .ncbi_utils import count_refs, load_ref_meta, load_ref_records
from .protein_guide_utils import get_guides, guide_loci, write_guide_fasta

# Product keys, in the order they are written / merged. Mirrors the extraction
# path's per-strain product names so the merged matrices line up with the
# isolate combined files (``blast_loci_utils.merge_per_strain_outputs``).
PRODUCTS = ("CDS", "genomic", "protein")


# ── Boundary-aware sequence rendering ─────────────────────────────────────────

# The full gene (exons + introns) is already returned soft-masked — exons UPPERCASE, introns
# lowercase — by ``exonerate_utils.build_result_from_model`` (``result['amplicon_seq']``), the
# *same* code path the user's own extracted loci go through, so isolates and tips are annotated
# identically (D-022). Case is inert to aligners / tree tools; it just keeps the exon-intron
# boundaries visible (and pinned to the right bases as gaps move) for by-hand work.

def exon_marked_cds(result: dict) -> str:
    """
    The spliced CDS with **case alternating per exon** (exon 1 upper, exon 2
    lower, …), so exon *junctions* are visible in the codon alignment too. The
    bases are unchanged (case only); ``exon_marked_cds(...).upper() ==
    result['cds_seq']``.
    """
    cds = result["cds_seq"]
    out, pos = [], 0
    for i, ex in enumerate(result["exons"]):
        chunk = cds[pos:pos + ex["length"]]
        out.append(chunk.upper() if i % 2 == 0 else chunk.lower())
        pos += ex["length"]
    out.append(cds[pos:])          # defensive: any remainder kept verbatim
    return "".join(out)


# ── Single-amplicon framing ───────────────────────────────────────────────────

def frame_consistent_amplicon(
    seq: str,
    seq_id: str,
    guide_fasta: str,
    locus: str,
    *,
    organism: str = "",
    exonerate_bin: str = "exonerate",
    minintron: int = 20,
    maxintron: int = 2000,
    geneticcode: int = 1,
    mark_exons: bool = False,
    manager=None,
) -> tuple[Optional[dict], str]:
    """
    Align one amplicon to the protein guide and return its frame-consistent
    products.

    Returns ``(record_or_None, status)``. ``record`` is a dict with ``id``,
    ``organism`` and the three rendered sequences (``cds`` codon-phased,
    ``genomic`` exon-upper/intron-lower, ``protein``) plus QC fields
    (``n_exons``, ``n_introns``, ``cds_length``, ``len_mod3``,
    ``n_internal_stops``, ``qc_clean``, ``verdict``). ``None`` means the amplicon
    could not be aligned to the guide (reported, not framed).
    """
    model_kind = detect_fasta_type(guide_fasta)
    model = MODEL_FOR_QUERYTYPE.get(model_kind, "protein2genome")

    with tempfile.TemporaryDirectory() as scratch:
        target = os.path.join(scratch, "amplicon.fasta")
        SeqIO.write([SeqRecord(Seq(seq), id=seq_id, description="")], target, "fasta")
        rc, raw, _gff = run_exonerate(
            guide_fasta, target, model, scratch,
            exonerate_bin=exonerate_bin, minintron=minintron, maxintron=maxintron,
            bestn=1, geneticcode=geneticcode, manager=manager,
            module="codon_prep", action=f"codonprep_{locus}_{seq_id}",
        )
    if rc != 0:
        return None, f"exonerate failed (exit {rc})"

    best, _others = select_best_model(parse_exonerate_gff(raw))
    if best is None:
        return None, "no gene model (amplicon does not align to the guide for this locus)"

    result = build_result_from_model(
        best, seq, seq_id, locus,
        exonerate_version=_probe_exonerate_version(exonerate_bin),
        geneticcode=geneticcode,
    )
    verdict = "PASS" if result["qc_clean"] else "REVIEW"
    record = {
        "id": seq_id,
        "organism": organism,
        "cds": exon_marked_cds(result) if mark_exons else result["cds_seq"],
        "genomic": result["amplicon_seq"],          # already exon-upper / intron-lower
        "protein": result["protein"],
        "n_exons": result["n_exons"],
        "n_introns": result["n_introns"],
        "cds_length": result["cds_length"],
        "len_mod3": result["len_mod3"],
        "n_internal_stops": result["n_internal_stops"],
        "qc_clean": result["qc_clean"],
        "verdict": verdict,
    }
    status = "OK" if result["qc_clean"] else (
        f"OK (QC: review — {result['n_internal_stops']} internal stop(s), "
        f"len % 3 = {result['len_mod3']})"
    )
    return record, status


# ── Nucleotide fallback for un-framable (intronic / short) amplicons ──────────

def orient_amplicon(
    seq: str,
    seq_id: str,
    orient_ref_fasta: str,
    *,
    blastn_bin: str = "blastn",
    min_pident: float = 65.0,
    evalue: float = 1e-5,
    manager=None,
) -> Optional[tuple[str, str, float, int]]:
    """
    blastn ``seq`` against ``orient_ref_fasta`` (the isolates' genomic locus) to BOTH confirm the
    amplicon belongs to this locus and orient it to the reference strand (D-027).

    Returns ``(oriented_seq, strand, pident, aln_len)`` for the best HSP at ``>= min_pident``, or
    ``None`` when nothing matches — a genuine failure (wrong locus / contamination), which is kept
    out of the matrices rather than silently included. This is the path for the standard
    intron-rich coding barcodes (e.g. fungal TEF1 EF1-728F/986R), which are largely **intronic**:
    they match the genomic gene continuously but have almost no exon for ``protein2genome`` to
    model, so they cannot be codon-framed — but they are still valid comparison data for the
    **nucleotide** (intron-inclusive) tree.
    """
    with tempfile.TemporaryDirectory() as scratch:
        q = os.path.join(scratch, "amplicon.fasta")
        SeqIO.write([SeqRecord(Seq(seq), id=seq_id, description="")], q, "fasta")
        rc, _err, tsv = run_blast(q, orient_ref_fasta, scratch,
                                  blast_bin=blastn_bin, task="blastn", evalue=evalue,
                                  manager=manager, module="blast",
                                  action=f"orient_amplicon:{seq_id}")
        if rc != 0:
            return None
        hsps = [h for h in parse_blast_hsps(tsv) if h["pident"] >= min_pident]
    if not hsps:
        return None
    best = max(hsps, key=lambda h: h["bitscore"])
    oriented = seq if best["sstrand"] == "plus" else str(Seq(seq).reverse_complement())
    return oriented, best["sstrand"], best["pident"], best["length"]


def _nt_tip_seqrecord(nt: dict, locus: str) -> SeqRecord:
    """Wrap a nucleotide-only (un-framed, oriented) tip for the **genomic** matrix only."""
    desc = _build_description({
        "tip": "yes", "locus": locus, "type": "genomic", "framed": "no",
        "nucleotide_only": "yes", "organism": nt.get("organism", ""),
        "amplicon_len": f"{nt['tip_len']}bp",
        "orient": f"{nt['strand']} vs isolate genomic",
        "orient_pident": f"{nt['pident']:.1f}",
        "note": "intronic/short amplicon — not codon-framed",
    })
    return SeqRecord(Seq(nt["genomic"]), id=nt["id"], description=desc)


def _tip_seqrecord(rec: dict, locus: str, product: str) -> SeqRecord:
    """Wrap a framed-tip product as a SeqRecord with a rich, NCBI-bracket header."""
    seq = {"CDS": rec["cds"], "genomic": rec["genomic"], "protein": rec["protein"]}[product]
    desc = _build_description({
        "tip": "yes", "locus": locus, "type": product,
        "organism": rec.get("organism", ""),
        "exons": rec["n_exons"], "introns": rec["n_introns"],
        "cds_len": f"{rec['cds_length']}bp",
        "internal_stops": rec["n_internal_stops"], "qc": rec["verdict"],
    })
    return SeqRecord(Seq(seq), id=rec["id"], description=desc)


# ── Per-locus orchestration ───────────────────────────────────────────────────

def coding_loci_with_tips(tips_dir, include_user: bool = True) -> list[str]:
    """
    Loci that (a) have at least one imported tip and (b) have a protein guide —
    i.e. the coding loci this page can process. rDNA tips (no guide) are excluded.
    """
    guided = set(guide_loci(include_user=include_user))
    return sorted(loc for loc in guided if count_refs(loc, ref_dir=tips_dir) > 0)


def _isolate_records(combined_dir, locus: str, product: str) -> list[SeqRecord]:
    """Read the isolate combined file for one product (empty list if absent)."""
    fp = Path(combined_dir) / f"{locus}_{product}_combined.fasta"
    if not fp.exists():
        return []
    return list(SeqIO.parse(str(fp), "fasta"))


def prepare_codon_locus(
    locus: str,
    tips_dir,
    out_dir,
    *,
    include_isolates: bool = True,
    isolate_combined_dir=None,
    include_user_guides: bool = True,
    exonerate_bin: str = "exonerate",
    blastn_bin: str = "blastn",
    minintron: int = 20,
    maxintron: int = 2000,
    geneticcode: int = 1,
    mark_exons: bool = False,
    strict_qc: bool = False,
    nt_fallback: bool = True,
    manager=None,
) -> dict:
    """
    Frame every tip for ``locus`` and write the three merged matrices (isolates +
    tips) to ``out_dir``.

    Returns a summary dict::

      {locus, n_tips, n_framed, n_flagged, n_failed, n_nt_only, n_isolates,
       outputs: {CDS, genomic, protein},   # written file paths
       rows: [ {id, organism, n_exons, n_introns, cds_length, verdict, product,
                included, status}, ... ],  # one per tip, for the UI table
       guide: <n guides used> or 0,
       status}

    ``n_nt_only`` tips were not codon-framable (intronic/short amplicons) but oriented to the
    isolate genomic and added to the **genomic** matrix as nucleotide-only (D-027).

    A locus with no protein guide is skipped (``status='no guide'``); a locus with
    no tips is skipped (``status='no tips'``).
    """
    rows: list[dict] = []
    summary = {
        "locus": locus, "n_tips": 0, "n_framed": 0, "n_flagged": 0,
        "n_failed": 0, "n_nt_only": 0, "n_isolates": 0, "outputs": {}, "rows": rows,
        "guide": 0, "status": "ok",
    }

    guides = get_guides(locus, include_user=include_user_guides)
    if not guides:
        summary["status"] = "no guide"
        return summary
    summary["guide"] = len(guides)

    tip_records = load_ref_records(locus, ref_dir=tips_dir)
    summary["n_tips"] = len(tip_records)
    if not tip_records:
        summary["status"] = "no tips"
        return summary

    meta = load_ref_meta(locus, ref_dir=tips_dir)
    os.makedirs(out_dir, exist_ok=True)

    # Orientation reference for the nucleotide fallback: the isolates' genomic locus. A tip that
    # cannot be codon-framed is oriented (and locus-confirmed) by blastn against this; with no
    # isolate genomic there is nothing to orient against, so the fallback is disabled (D-027).
    iso_dir = isolate_combined_dir if (include_isolates and isolate_combined_dir) else None
    orient_ref = None
    if nt_fallback and iso_dir:
        cand = Path(iso_dir) / f"{locus}_genomic_combined.fasta"
        if cand.exists() and cand.stat().st_size > 0:
            orient_ref = str(cand)

    with tempfile.TemporaryDirectory() as gdir:
        guide_fasta = write_guide_fasta(locus, os.path.join(gdir, f"{locus}_guide.fasta"),
                                        include_user=include_user_guides)
        framed: list[dict] = []
        nt_only: list[dict] = []
        for rec in tip_records:
            m = meta.get(rec.id) or meta.get(rec.id.split(".")[0])
            organism = m.organism if m else ""
            fr, status = frame_consistent_amplicon(
                str(rec.seq), rec.id, guide_fasta, locus, organism=organism,
                exonerate_bin=exonerate_bin, minintron=minintron, maxintron=maxintron,
                geneticcode=geneticcode, mark_exons=mark_exons, manager=manager,
            )
            included = fr is not None
            nt = None
            if fr is not None:
                if strict_qc and not fr["qc_clean"]:
                    included = False
                    status = "excluded by strict QC (" + status + ")"
                    summary["n_flagged"] += 1
                elif not fr["qc_clean"]:
                    summary["n_flagged"] += 1
                if included:
                    framed.append(fr)
                    summary["n_framed"] += 1
            else:
                # Could not codon-frame: keep it for the nucleotide (intron-inclusive) matrix
                # if it orients to the isolate genomic — the intron-rich-barcode path (D-027).
                if orient_ref:
                    o = orient_amplicon(str(rec.seq), rec.id, orient_ref,
                                        blastn_bin=blastn_bin, manager=manager)
                    if o:
                        oriented, strand, pident, aln_len = o
                        nt = {"id": rec.id, "organism": organism, "genomic": oriented,
                              "strand": strand, "pident": pident, "aln_len": aln_len,
                              "tip_len": len(rec.seq)}
                        nt_only.append(nt)
                        summary["n_nt_only"] += 1
                        included = True
                        status = (f"nucleotide-only — intronic/short amplicon, not codon-framed; "
                                  f"oriented {strand} vs isolate genomic "
                                  f"({pident:.0f}% / {aln_len} bp)")
                if nt is None:
                    summary["n_failed"] += 1
            rows.append({
                "id": rec.id, "organism": organism,
                "n_exons": fr["n_exons"] if fr else 0,
                "n_introns": fr["n_introns"] if fr else 0,
                "cds_length": fr["cds_length"] if fr else (nt["tip_len"] if nt else 0),
                "verdict": fr["verdict"] if fr else ("NT-ONLY" if nt else "FAIL"),
                "product": ("CDS+genomic+protein" if fr else
                            ("genomic (nt)" if nt else "—")),
                "included": included, "status": status,
            })

    # ── merge isolates (as-is) + framed tips (+ nucleotide-only tips, genomic only) ──
    for product in PRODUCTS:
        iso = _isolate_records(iso_dir, locus, product) if iso_dir else []
        tips = [_tip_seqrecord(fr, locus, product) for fr in framed]
        nt_tips = [_nt_tip_seqrecord(nt, locus) for nt in nt_only] if product == "genomic" else []
        if not iso and not tips and not nt_tips:
            continue
        path = os.path.join(out_dir, f"{locus}_{product}_combined.fasta")
        SeqIO.write(iso + tips + nt_tips, path, "fasta")
        summary["outputs"][product] = path
        if product == "CDS":
            summary["n_isolates"] = len(iso)
    return summary
