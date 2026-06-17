"""
assembly_utils.py
-----------------
Parse FASTA assemblies from multiple assemblers, extract per-contig stats,
and compute standard assembly QC metrics (N50, GC%, etc.).

Supported assemblers (auto-detected from contig header format):
  - SPAdes / metaSPAdes
  - MEGAHIT
  - Flye
  - Hifiasm
  - Velvet
  - Generic (falls back gracefully — coverage from header treated as optional)
"""

import re
import numpy as np
from pathlib import Path
from typing import Optional
from Bio import SeqIO


# ── Regex patterns for contig name formats ──────────────────────────────────

_PATTERNS = {
    # >NODE_1_length_500_cov_45.3 (SPAdes / Velvet)
    "spades": re.compile(r"NODE_\d+_length_(\d+)_cov_([\d.]+)", re.IGNORECASE),
    # >k141_1 flag=1 multi=45.3 len=500 (MEGAHIT)
    "megahit": re.compile(r"multi=([\d.]+)\s+len=(\d+)", re.IGNORECASE),
    # >contig_1 or >edge_1 (Flye)
    "flye": re.compile(r"^(contig|edge)_\d+$", re.IGNORECASE),
    # >ptg000001l or >ctg000001l (Hifiasm)
    "hifiasm": re.compile(r"^(ptg|ctg)\d+[lhp]$", re.IGNORECASE),
}


def detect_assembler(fasta_path: str) -> str:
    """Read the first few headers and infer which assembler produced the FASTA."""
    with open(fasta_path) as f:
        headers_checked = 0
        for line in f:
            if not line.startswith(">"):
                continue
            header = line[1:].strip()
            contig_id = header.split()[0]
            if _PATTERNS["spades"].search(header):
                return "spades"
            if _PATTERNS["megahit"].search(header):
                return "megahit"
            if _PATTERNS["hifiasm"].match(contig_id):
                return "hifiasm"
            if _PATTERNS["flye"].match(contig_id):
                return "flye"
            headers_checked += 1
            if headers_checked >= 5:
                break
    return "unknown"


def _parse_header_coverage(header: str, assembler: str) -> Optional[float]:
    """Extract coverage from contig header if the assembler embeds it."""
    if assembler in ("spades", "velvet"):
        m = _PATTERNS["spades"].search(header)
        return float(m.group(2)) if m else None
    if assembler == "megahit":
        m = _PATTERNS["megahit"].search(header)
        return float(m.group(1)) if m else None
    return None


def _gc_percent(seq: str) -> float:
    seq = seq.upper()
    n = len(seq)
    if n == 0:
        return 0.0
    return (seq.count("G") + seq.count("C")) / n * 100


def _n50(lengths: list) -> tuple:
    """Return (N50, L50) from a list of contig lengths."""
    sorted_len = sorted(lengths, reverse=True)
    total = sum(sorted_len)
    cumsum = 0
    for i, l in enumerate(sorted_len):
        cumsum += l
        if cumsum >= total / 2:
            return l, i + 1
    return 0, 0


def get_assembly_stats(fasta_path: str) -> dict:
    """
    Parse a genome assembly FASTA and return:
      - Summary stats (N50, GC%, contig count, etc.)
      - Per-contig list for downstream use (blob plots, filtering)
    """
    assembler = detect_assembler(fasta_path)

    lengths, gc_list, cov_list = [], [], []
    contigs = []

    for record in SeqIO.parse(fasta_path, "fasta"):
        seq = str(record.seq)
        length = len(seq)
        gc = _gc_percent(seq)
        cov = _parse_header_coverage(record.description, assembler)

        lengths.append(length)
        gc_list.append(gc)
        if cov is not None:
            cov_list.append(cov)

        contigs.append({
            "contig_id": record.id,
            "length": length,
            "gc_percent": round(gc, 2),
            "header_coverage": round(cov, 2) if cov is not None else None,
        })

    n50, l50 = _n50(lengths)
    total_bp = sum(lengths)

    return {
        "assembler": assembler,
        "num_contigs": len(lengths),
        "total_length_bp": total_bp,
        "total_length_mb": round(total_bp / 1e6, 2),
        "n50": n50,
        "l50": l50,
        "largest_contig": max(lengths) if lengths else 0,
        "mean_gc": round(float(np.mean(gc_list)), 2) if gc_list else 0.0,
        "mean_coverage_header": round(float(np.mean(cov_list)), 2) if cov_list else None,
        "contigs": contigs,
    }


def find_assemblies_recursive(root_dir: str,
                               extensions: tuple = (".fasta", ".fa", ".fna")) -> list:
    """Walk a directory tree and return paths to all FASTA files found."""
    found = []
    for path in Path(root_dir).rglob("*"):
        if path.suffix.lower() in extensions and path.is_file():
            found.append(str(path))
    return sorted(found)


def suggest_strain_id(fasta_path: str) -> str:
    """
    Derive a clean strain ID from a filename by stripping common
    assembler/pipeline suffixes.
    """
    name = Path(fasta_path).stem
    for suffix in [
        "_assembly", "_final", "_polished", "_genome",
        "_contigs", "_scaffolds", "_consensus", ".assembly",
    ]:
        name = name.replace(suffix, "")
    return name
