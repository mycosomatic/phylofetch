"""
assembly_utils.py
-----------------
Parse FASTA assemblies from multiple assemblers, extract per-contig stats,
and compute standard assembly QC metrics (N50, GC%, etc.).

Supported assemblers (auto-detected from contig header format, then filename):
  - SPAdes / metaSPAdes
  - MEGAHIT
  - Flye
  - Hifiasm
  - Velvet
  - MaSuRCA, ABySS, Canu, EGAP (detected from filename — these re-header output)
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


# Filename markers for assemblers/pipelines that re-header their output so the
# contig-header patterns above can't see them. Ordered: true assemblers first,
# then pipelines/polishers (lower priority).
_FILENAME_MARKERS: list[tuple[str, str]] = [
    ("spades", "spades"),
    ("masurca", "masurca"),
    ("flye", "flye"),
    ("hifiasm", "hifiasm"),
    ("megahit", "megahit"),
    ("abyss", "abyss"),
    ("canu", "canu"),
    ("velvet", "velvet"),
    ("egap", "egap"),
    ("medaka", "medaka"),
    ("pilon", "pilon"),
    ("racon", "racon"),
    ("polca", "polca"),
]


def _assembler_from_filename(fasta_path: str) -> str:
    """Infer assembler/pipeline from the filename when headers are uninformative."""
    name = Path(fasta_path).name.lower()
    for label, marker in _FILENAME_MARKERS:
        if marker in name:
            return label
    return "unknown"


def detect_assembler(fasta_path: str) -> str:
    """
    Infer which assembler produced the FASTA.

    Contig headers are authoritative (a polished SPAdes assembly keeps its
    ``NODE_..`` headers, so it is still reported as ``spades``). When headers
    are uninformative — e.g. EGAP, MaSuRCA, ABySS which re-header output — fall
    back to recognising the assembler/pipeline name in the filename.
    """
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
    return _assembler_from_filename(fasta_path)


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
        "_final_EGAP_assembly", "_final_polish_assembly", "_best_assembly",
        "_assembly", "_final", "_polished", "_genome", "_EGAP",
        "_contigs", "_scaffolds", "_consensus", ".assembly",
    ]:
        name = name.replace(suffix, "")
    return name


# Distinguishing tokens used to disambiguate strain IDs when several assembly
# versions of one strain collapse to the same suggested ID (e.g. an EGAP final
# and a polish final for the same strain). Ordered most-specific first.
_SOURCE_TAGS: list[str] = [
    "egap", "polish", "best", "consensus",
    "spades", "masurca", "flye", "hifiasm", "megahit",
    "abyss", "canu", "velvet", "medaka", "pilon", "racon", "polca",
    "scaffolds", "contigs",
]


def source_tag_from_filename(fasta_path: str) -> str:
    """
    Return a short tag identifying the assembly version/source from a filename,
    used to disambiguate otherwise-identical suggested strain IDs.
    Returns "" when no known token is present.
    """
    name = Path(fasta_path).name.lower()
    for tag in _SOURCE_TAGS:
        if tag in name:
            return tag
    return ""


def suggest_unique_strain_ids(fasta_paths: list[str]) -> list[str]:
    """
    Derive strain IDs for a batch of FASTAs, guaranteeing uniqueness.

    Clean IDs (from ``suggest_strain_id``) are kept when unique. When several
    files reduce to the same base ID — e.g. ``X_final_EGAP_assembly`` and
    ``X_final_polish_assembly`` both reducing to ``X`` — *every* member of the
    colliding group gets a source tag (``X_egap``, ``X_polish``) so none is
    silently dropped on import. Falls back to the full filename stem, then a
    numeric suffix, if a tag is unavailable or still ambiguous.
    """
    from collections import Counter

    bases = [suggest_strain_id(fp) for fp in fasta_paths]
    base_counts = Counter(bases)

    ids: list[str] = []
    seen: set[str] = set()
    for fp, base in zip(fasta_paths, bases):
        if base_counts[base] > 1:
            tag = source_tag_from_filename(fp)
            sid = f"{base}_{tag}" if tag else Path(fp).stem
        else:
            sid = base
        candidate = sid
        n = 2
        while candidate in seen:
            candidate = f"{sid}_{n}"
            n += 1
        seen.add(candidate)
        ids.append(candidate)
    return ids


# ── QUAST report parsing ─────────────────────────────────────────────────────

# Curated QUAST metrics surfaced in the UI, in display order.
QUAST_DISPLAY_KEYS = [
    "# contigs", "Largest contig", "Total length", "GC (%)",
    "N50", "N75", "N90", "L50", "L75", "L90", "# N's per 100 kbp",
]


def _coerce_number(val: str):
    """Convert a QUAST value string to int or float when possible."""
    try:
        if "." in val or "e" in val.lower():
            return float(val)
        return int(val)
    except (ValueError, AttributeError):
        return val


def parse_quast_report(report_path: str) -> dict:
    """
    Parse a QUAST ``report.tsv`` into ``{metric: value}``.

    report.tsv is tab-separated: metric name in the first column, value in the
    second (single-assembly report). The assembly name row is stored under
    ``assembly_name``. Numeric values are coerced to int/float where possible.
    """
    metrics: dict = {}
    p = Path(report_path)
    if not p.exists():
        return metrics
    for line in p.read_text(encoding="utf-8", errors="replace").splitlines():
        if not line.strip() or "\t" not in line:
            continue
        parts = line.split("\t")
        key = parts[0].strip()
        val = parts[1].strip() if len(parts) > 1 else ""
        if not key:
            continue
        if key.lower() == "assembly":
            metrics["assembly_name"] = val
            continue
        metrics[key] = _coerce_number(val)
    return metrics


def find_quast_report(assembly_path: str) -> Optional[str]:
    """
    Locate the QUAST ``report.tsv`` associated with an assembly FASTA.

    EGAP writes it to a sibling directory ``<stem>_quast/`` (some intermediate
    assemblies use ``<filename>_quast/``). Returns the path to report.tsv if
    found, else None.
    """
    p = Path(assembly_path)
    parent = p.parent
    candidates = [
        parent / f"{p.stem}_quast" / "report.tsv",   # EGAP final: X_quast/report.tsv
        parent / f"{p.name}_quast" / "report.tsv",    # intermediates: X.fasta_quast/report.tsv
    ]
    for c in candidates:
        if c.is_file():
            return str(c)
    # Fallback: any *_quast dir whose name starts with the assembly stem
    for d in sorted(parent.glob("*_quast")):
        if d.is_dir() and d.name.startswith(p.stem):
            rep = d / "report.tsv"
            if rep.is_file():
                return str(rep)
    return None
