"""
primer_utils.py
---------------
PCR primer-based locus extraction from genome assemblies.

Strategy: run blastn-short to locate forward and reverse primer binding sites
in assembly contigs with edit-distance tolerance (max_mismatches), then extract
the amplicon sequence between them on the same contig and strand.

Primer orientation conventions (consistent with wetlab usage):
  fwd — 5'→3' sequence that anneals to the minus (antisense) strand
         → blastn reports a PLUS-strand hit (sstart < send)
  rev — 5'→3' sequence that anneals to the plus (sense) strand
         → blastn reports a MINUS-strand hit (sstart > send)

When the amplicon itself lies on the minus strand the roles are reversed:
  fwd anneals to plus strand (sstart > send), rev anneals to minus strand (sstart < send)
"""

from __future__ import annotations

import os
import subprocess
import tempfile
from collections import defaultdict
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Optional

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


@dataclass
class PrimerPair:
    name: str
    locus: str
    fwd: str          # 5'→3' forward primer (including IUPAC degenerate codes)
    rev: str          # 5'→3' reverse primer (including IUPAC degenerate codes)
    min_amplicon: int = 100
    max_amplicon: int = 5000


# Common fungal phylogenetic primer pairs (curated catalogue)
PRIMER_CATALOGUE: dict[str, PrimerPair] = {
    "ITS1/ITS4": PrimerPair(
        name="ITS1/ITS4", locus="ITS",
        fwd="TCCGTAGGTGAACCTGCGG",
        rev="TCCTCCGCTTATTGATATGC",
        min_amplicon=400, max_amplicon=800,
    ),
    "ITS1F/ITS4": PrimerPair(
        name="ITS1F/ITS4", locus="ITS",
        fwd="CTTGGTCATTTAGAGGAAGTAA",
        rev="TCCTCCGCTTATTGATATGC",
        min_amplicon=400, max_amplicon=800,
    ),
    "NL1/NL4": PrimerPair(
        name="NL1/NL4", locus="LSU",
        fwd="GCATATCAATAAGCGGAGGAAAAG",
        rev="GGTCCGTGTTTCAAGACGG",
        min_amplicon=600, max_amplicon=1100,
    ),
    "NL1/NL4B": PrimerPair(
        name="NL1/NL4B", locus="LSU",
        fwd="GCATATCAATAAGCGGAGGAAAAG",
        rev="CAAGGCATGCAAGTTGACTCGAG",
        min_amplicon=900, max_amplicon=1400,
    ),
    "EF1-728F/EF1-986R": PrimerPair(
        name="EF1-728F/EF1-986R", locus="TEF1",
        fwd="CATCGAGAAGTTCGAGAAGG",
        rev="TACTTGAAGGAACCCTTACC",
        min_amplicon=250, max_amplicon=350,
    ),
    "EF1-526F/EF1-1567R": PrimerPair(
        name="EF1-526F/EF1-1567R", locus="TEF1",
        fwd="GTCGTYGTYATYGGHCAYGT",
        rev="ACHGTRCCRATACCACCRATCTT",
        min_amplicon=900, max_amplicon=1200,
    ),
    "RPB2-5F2/RPB2-7cR": PrimerPair(
        name="RPB2-5F2/RPB2-7cR", locus="RPB2",
        fwd="GAYGAYMGWGATCAYTTYGG",
        rev="CCCATRGCTTGTYYTTGCAT",
        min_amplicon=700, max_amplicon=1100,
    ),
    "BT2a/BT2b": PrimerPair(
        name="BT2a/BT2b", locus="TUB2",
        fwd="GGTAACCAAATCGGTGCTGCTTTC",
        rev="ACCCTCAGTGTAGTGACCCTTGGC",
        min_amplicon=400, max_amplicon=600,
    ),
    "ACT-512F/ACT-783R": PrimerPair(
        name="ACT-512F/ACT-783R", locus="ACT",
        fwd="ATGTGCAAAGCCGGCTTCGCGGGCGACGATGCCCC",
        rev="TACGAGTCCTTCTGGCCCAT",
        min_amplicon=250, max_amplicon=350,
    ),
    "CAL-228F/CAL-737R": PrimerPair(
        name="CAL-228F/CAL-737R", locus="CAL",
        fwd="GAGTTCAAGGAGGCCTTCTCCC",
        rev="CATCTTTCTGGCCATCATGG",
        min_amplicon=450, max_amplicon=600,
    ),
    "GS1/GS2r": PrimerPair(
        name="GS1/GS2r", locus="GS",
        fwd="ATGGCCACCGTCGAGTTCTA",
        rev="TTGTCCACCACCGACTTCTT",
        min_amplicon=400, max_amplicon=700,
    ),
    "HIS3-1F/HIS3-1R": PrimerPair(
        name="HIS3-1F/HIS3-1R", locus="HIS3",
        fwd="TGGCTTGTCGACTTCTGGTG",
        rev="TGGCGATGTTGACGTCTTCC",
        min_amplicon=400, max_amplicon=700,
    ),
}

# Group catalogue by locus for UI lookup
LOCUS_PRIMER_MAP: dict[str, list[str]] = {}
for _pp in PRIMER_CATALOGUE.values():
    LOCUS_PRIMER_MAP.setdefault(_pp.locus, []).append(_pp.name)


def _run_blastn_short(
    primer_fa: str,
    assembly_fa: str,
    perc_identity: float,
    blastn_bin: str = "blastn",
) -> list[dict]:
    """Run blastn-short; return list of parsed HSP dicts."""
    fields = [
        "qseqid", "sseqid", "pident", "length", "mismatch",
        "qstart", "qend", "sstart", "send", "evalue", "bitscore",
    ]
    cmd = [
        blastn_bin, "-task", "blastn-short",
        "-query", primer_fa,
        "-subject", assembly_fa,
        "-perc_identity", f"{perc_identity:.1f}",
        "-outfmt", "6 " + " ".join(fields),
        "-word_size", "7",
        "-num_alignments", "1000",
    ]
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode != 0:
        return []

    rows: list[dict] = []
    for line in proc.stdout.splitlines():
        parts = line.split("\t")
        if len(parts) < len(fields):
            continue
        try:
            rows.append({
                "qseqid":   parts[0],
                "sseqid":   parts[1],
                "pident":   float(parts[2]),
                "length":   int(parts[3]),
                "mismatch": int(parts[4]),
                "sstart":   int(parts[7]),
                "send":     int(parts[8]),
                "evalue":   float(parts[9]),
                "bitscore": float(parts[10]),
            })
        except (ValueError, IndexError):
            continue
    return rows


def find_primer_amplicons(
    assembly_fasta: str,
    primer_pair: PrimerPair,
    max_mismatches: int = 2,
    blastn_bin: str = "blastn",
) -> list[dict]:
    """
    Search an assembly for amplicons bounded by primer_pair.

    Returns candidates sorted by ascending total-mismatch count.
    Each candidate dict:
      contig, strand, amp_start, amp_end, amp_len,
      fwd_mismatch, rev_mismatch, fwd_hit, rev_hit
    """
    fwd_len = len(primer_pair.fwd.replace(" ", ""))
    rev_len = len(primer_pair.rev.replace(" ", ""))
    min_len = min(fwd_len, rev_len)
    perc_id = max(70.0, (min_len - max_mismatches) / min_len * 100)

    with tempfile.TemporaryDirectory() as tdir:
        primer_fa = os.path.join(tdir, "primers.fasta")
        with open(primer_fa, "w") as f:
            f.write(f">FWD\n{primer_pair.fwd}\n>REV\n{primer_pair.rev}\n")
        hits = _run_blastn_short(primer_fa, assembly_fasta, perc_id, blastn_bin)

    fwd_hits = [h for h in hits if h["qseqid"] == "FWD"]
    rev_hits  = [h for h in hits if h["qseqid"] == "REV"]

    fwd_by_contig: dict[str, list] = defaultdict(list)
    rev_by_contig: dict[str, list] = defaultdict(list)
    for h in fwd_hits:
        fwd_by_contig[h["sseqid"]].append(h)
    for h in rev_hits:
        rev_by_contig[h["sseqid"]].append(h)

    candidates: list[dict] = []
    for contig in set(fwd_by_contig) & set(rev_by_contig):
        for fh in fwd_by_contig[contig]:
            for rh in rev_by_contig[contig]:
                # Amplicon on + strand: F hit is plus (sstart < send), R hit is minus (sstart > send)
                if fh["sstart"] < fh["send"] and rh["sstart"] > rh["send"]:
                    amp_start, amp_end, strand = fh["sstart"], rh["sstart"], "+"
                # Amplicon on - strand: F hit is minus (sstart > send), R hit is plus (sstart < send)
                elif fh["sstart"] > fh["send"] and rh["sstart"] < rh["send"]:
                    amp_start, amp_end, strand = rh["sstart"], fh["sstart"], "-"
                else:
                    continue

                amp_len = amp_end - amp_start + 1
                if primer_pair.min_amplicon <= amp_len <= primer_pair.max_amplicon:
                    candidates.append({
                        "contig":       contig,
                        "strand":       strand,
                        "amp_start":    amp_start,
                        "amp_end":      amp_end,
                        "amp_len":      amp_len,
                        "fwd_mismatch": fh["mismatch"],
                        "rev_mismatch": rh["mismatch"],
                        "fwd_hit":      fh,
                        "rev_hit":      rh,
                    })

    candidates.sort(key=lambda c: c["fwd_mismatch"] + c["rev_mismatch"])
    return candidates


def extract_primer_amplicon(
    assembly_fasta: str,
    hit: dict,
    strain_id: str,
    locus_name: str,
    output_dir: str,
    primer_pair: PrimerPair,
) -> Optional[SeqRecord]:
    """Extract and write the amplicon FASTA; return SeqRecord or None."""
    contig_seq: Optional[Seq] = None
    for rec in SeqIO.parse(assembly_fasta, "fasta"):
        if rec.id == hit["contig"]:
            contig_seq = rec.seq
            break
    if contig_seq is None:
        return None

    # BLAST coords are 1-based, inclusive
    extracted = contig_seq[hit["amp_start"] - 1 : hit["amp_end"]]
    if hit["strand"] == "-":
        extracted = extracted.reverse_complement()

    ts = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
    description = (
        f"[sample={strain_id}] [locus={locus_name}] [type=amplicon] "
        f"[primer={primer_pair.name}] [contig={hit['contig']}] "
        f"[coords={hit['amp_start']}..{hit['amp_end']}] [strand={hit['strand']}] "
        f"[amplicon_len={hit['amp_len']}bp] "
        f"[fwd_mm={hit['fwd_mismatch']}] [rev_mm={hit['rev_mismatch']}] "
        f"[extracted={ts}]"
    )
    out_rec = SeqRecord(
        Seq(str(extracted)),
        id=f"{strain_id}_{locus_name}_amplicon",
        description=description,
    )
    out_path = Path(output_dir) / f"{locus_name}_amplicon.fasta"
    out_path.parent.mkdir(parents=True, exist_ok=True)
    SeqIO.write([out_rec], str(out_path), "fasta")
    return out_rec


def run_primer_extraction(
    assembly_fasta: str,
    primer_pair: PrimerPair,
    output_dir: str,
    strain_id: str,
    locus_name: str,
    max_mismatches: int = 2,
    blastn_bin: str = "blastn",
) -> tuple[Optional[dict], str]:
    """
    Top-level entry point for primer-based locus extraction.

    Returns (result_dict, status).
    result_dict keys: contig, strand, amp_start, amp_end, amp_len,
                      fwd_mismatch, rev_mismatch, primer_name,
                      output_fasta, n_candidates.
    status is "ok" on success or a human-readable failure message.
    """
    candidates = find_primer_amplicons(
        assembly_fasta, primer_pair,
        max_mismatches=max_mismatches, blastn_bin=blastn_bin,
    )
    if not candidates:
        return None, (
            f"No amplicon found with {primer_pair.name} "
            f"(max_mismatches={max_mismatches}, "
            f"size {primer_pair.min_amplicon}–{primer_pair.max_amplicon} bp)"
        )

    best = candidates[0]
    rec = extract_primer_amplicon(
        assembly_fasta, best, strain_id, locus_name, output_dir, primer_pair,
    )
    if rec is None:
        return None, f"Contig '{best['contig']}' not found in assembly FASTA"

    out_fasta = str(Path(output_dir) / f"{locus_name}_amplicon.fasta")
    return {
        "contig":       best["contig"],
        "strand":       best["strand"],
        "amp_start":    best["amp_start"],
        "amp_end":      best["amp_end"],
        "amp_len":      best["amp_len"],
        "fwd_mismatch": best["fwd_mismatch"],
        "rev_mismatch": best["rev_mismatch"],
        "primer_name":  primer_pair.name,
        "output_fasta": out_fasta,
        "n_candidates": len(candidates),
    }, "ok"
