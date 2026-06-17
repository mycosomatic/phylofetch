"""
busco_utils.py
--------------
Parse BUSCO v4/v5 and Compleasm output, build occupancy matrices, and export
single-copy ortholog FASTAs for phylogenomic supermatrix construction.

Input directories can come from EGAP, standalone BUSCO runs, or Compleasm.
No upstream pipeline dependency is assumed.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


@dataclass
class BuscoResult:
    sample_id: str
    run_dir: Path
    tool: str          # "busco" | "compleasm"
    lineage: str
    complete: int = 0
    single_copy: int = 0
    duplicated: int = 0
    fragmented: int = 0
    missing: int = 0
    total: int = 0
    completeness_pct: float = 0.0
    # Map busco_id → status (C_S, C_D, F, M)
    gene_status: dict[str, str] = field(default_factory=dict)
    # Map busco_id → Path of single-copy sequence file
    sc_fasta_paths: dict[str, Path] = field(default_factory=dict)


# ── Parsers ──────────────────────────────────────────────────────────────────

def _parse_busco_short_summary(path: Path) -> dict:
    """Parse a BUSCO v4/v5 short_summary*.txt file into a stats dict."""
    stats: dict = {"complete": 0, "single_copy": 0, "duplicated": 0,
                   "fragmented": 0, "missing": 0, "total": 0, "lineage": ""}
    text = path.read_text()
    for line in text.splitlines():
        line = line.strip()
        m = re.search(r"C:(\d+\.\d+)%\[S:(\d+\.\d+)%,D:(\d+\.\d+)%\],F:(\d+\.\d+)%,M:(\d+\.\d+)%,n:(\d+)", line)
        if m:
            stats["total"] = int(m.group(6))
            stats["completeness_pct"] = float(m.group(1))
        if re.match(r"\d+\s+Complete BUSCOs", line):
            stats["complete"] = int(line.split()[0])
        if re.match(r"\d+\s+Complete and single-copy", line):
            stats["single_copy"] = int(line.split()[0])
        if re.match(r"\d+\s+Complete and duplicated", line):
            stats["duplicated"] = int(line.split()[0])
        if re.match(r"\d+\s+Fragmented", line):
            stats["fragmented"] = int(line.split()[0])
        if re.match(r"\d+\s+Missing", line):
            stats["missing"] = int(line.split()[0])
        if re.match(r"\d+\s+Total BUSCO groups", line):
            stats["total"] = int(line.split()[0])
        if "lineage dataset" in line.lower():
            stats["lineage"] = line.split()[-1] if line.split() else ""
    return stats


def _parse_busco_full_table(full_table: Path) -> dict[str, str]:
    """Parse full_table.tsv → {busco_id: status} where status ∈ Complete/Duplicated/Fragmented/Missing."""
    gene_status: dict[str, str] = {}
    if not full_table.exists():
        return gene_status
    for line in full_table.read_text().splitlines():
        if line.startswith("#") or not line.strip():
            continue
        parts = line.split("\t")
        if len(parts) >= 2:
            busco_id = parts[0]
            status   = parts[1]
            # Only record one entry per busco_id; Complete > Fragmented > Missing
            if busco_id not in gene_status or status.startswith("Complete"):
                gene_status[busco_id] = status
    return gene_status


def _collect_sc_fastas(run_dir: Path) -> dict[str, Path]:
    """Locate single_copy_busco_sequences/ FASTA files and index by BUSCO ID."""
    sc_dir = run_dir / "busco_sequences" / "single_copy_busco_sequences"
    if not sc_dir.exists():
        sc_dir = run_dir / "single_copy_busco_sequences"
    if not sc_dir.exists():
        return {}
    paths: dict[str, Path] = {}
    for fp in sc_dir.glob("*.fna"):
        paths[fp.stem] = fp
    for fp in sc_dir.glob("*.faa"):
        paths.setdefault(fp.stem, fp)
    return paths


def _parse_compleasm(run_dir: Path) -> tuple[dict, dict[str, str]]:
    """Parse Compleasm summary.txt and mb_table.txt."""
    stats: dict = {"complete": 0, "single_copy": 0, "duplicated": 0,
                   "fragmented": 0, "missing": 0, "total": 0,
                   "lineage": "", "completeness_pct": 0.0}
    gene_status: dict[str, str] = {}

    summary = run_dir / "summary.txt"
    if summary.exists():
        for line in summary.read_text().splitlines():
            m = re.match(r"S:(\d+\.\d+)%", line)
            if m:
                stats["completeness_pct"] = float(m.group(1))
            for label, key in [("Single", "single_copy"), ("Duplicated", "duplicated"),
                                ("Fragmented", "fragmented"), ("Missing", "missing")]:
                if f"{label}:" in line:
                    n = re.search(r"\d+", line.split(f"{label}:")[1])
                    if n:
                        stats[key] = int(n.group())
        stats["complete"] = stats["single_copy"] + stats["duplicated"]
        stats["total"] = (stats["complete"] + stats["fragmented"] + stats["missing"])

    mb_table = run_dir / "mb_table.txt"
    if not mb_table.exists():
        mb_table = next(run_dir.glob("*_mb_table.txt"), None)
    if mb_table and mb_table.exists():
        for line in mb_table.read_text().splitlines():
            if line.startswith("#") or not line.strip():
                continue
            parts = line.split("\t")
            if len(parts) >= 2:
                gene_status[parts[0]] = parts[1]

    return stats, gene_status


# ── Main scanner ─────────────────────────────────────────────────────────────

def scan_busco_run(run_dir: Path, sample_id: Optional[str] = None) -> Optional[BuscoResult]:
    """
    Detect and parse a BUSCO or Compleasm run directory.
    sample_id defaults to the directory name if not provided.
    """
    run_dir = Path(run_dir).resolve()
    if not run_dir.is_dir():
        return None

    sid = sample_id or run_dir.name

    # Try BUSCO v5 layout first
    short_summaries = list(run_dir.glob("short_summary*.txt")) + \
                      list(run_dir.glob("*/short_summary*.txt"))
    full_tables     = list(run_dir.glob("full_table.tsv")) + \
                      list(run_dir.glob("run_*/full_table.tsv"))

    if short_summaries:
        ss = short_summaries[0]
        stats = _parse_busco_short_summary(ss)
        ft    = full_tables[0] if full_tables else None
        gene_status = _parse_busco_full_table(ft) if ft else {}
        sc_paths = _collect_sc_fastas(ss.parent if ss.parent != run_dir else run_dir)
        return BuscoResult(
            sample_id=sid, run_dir=run_dir, tool="busco",
            lineage=stats.get("lineage", ""),
            complete=stats["complete"], single_copy=stats["single_copy"],
            duplicated=stats["duplicated"], fragmented=stats["fragmented"],
            missing=stats["missing"], total=stats["total"],
            completeness_pct=stats.get("completeness_pct", 0.0),
            gene_status=gene_status, sc_fasta_paths=sc_paths,
        )

    # Try Compleasm layout
    if (run_dir / "summary.txt").exists() or list(run_dir.glob("*_mb_table.txt")):
        stats, gene_status = _parse_compleasm(run_dir)
        return BuscoResult(
            sample_id=sid, run_dir=run_dir, tool="compleasm",
            lineage=stats.get("lineage", ""),
            complete=stats["complete"], single_copy=stats["single_copy"],
            duplicated=stats["duplicated"], fragmented=stats["fragmented"],
            missing=stats["missing"], total=stats["total"],
            completeness_pct=stats.get("completeness_pct", 0.0),
            gene_status=gene_status,
        )

    return None


# ── Occupancy matrix ──────────────────────────────────────────────────────────

def build_occupancy_matrix(results: list[BuscoResult]) -> pd.DataFrame:
    """
    Build a sample × BUSCO occupancy matrix.
    Values: 'S' (single-copy complete), 'D' (duplicated), 'F' (fragmented), 'M' (missing/absent).
    """
    all_ids: set[str] = set()
    for r in results:
        all_ids.update(r.gene_status.keys())

    def _code(status: str) -> str:
        s = status.upper()
        if "SINGLE" in s or s == "COMPLETE":
            return "S"
        if "DUPLIC" in s:
            return "D"
        if "FRAG" in s:
            return "F"
        return "M"

    rows = []
    for r in results:
        row: dict[str, str] = {"sample": r.sample_id}
        for bid in sorted(all_ids):
            raw = r.gene_status.get(bid, "Missing")
            row[bid] = _code(raw)
        rows.append(row)

    df = pd.DataFrame(rows).set_index("sample")
    return df


def filter_single_copy_buscos(
    matrix: pd.DataFrame,
    min_occupancy: float = 0.75,
    exclude_duplicated: bool = True,
) -> list[str]:
    """
    Return BUSCO IDs that are single-copy in at least min_occupancy fraction of samples.
    Optionally exclude any BUSCO that is duplicated in ANY sample.
    """
    n_samples = len(matrix)
    selected: list[str] = []
    for busco_id in matrix.columns:
        col = matrix[busco_id]
        n_single = (col == "S").sum()
        if exclude_duplicated and (col == "D").any():
            continue
        if n_single / n_samples >= min_occupancy:
            selected.append(busco_id)
    return selected


def export_sc_fastas(results: list[BuscoResult],
                     selected_buscos: list[str],
                     output_dir: str) -> dict[str, str]:
    """
    Collect single-copy FASTA sequences for each selected BUSCO across all samples.
    Writes one FASTA per BUSCO containing one sequence per sample.
    Returns {busco_id: fasta_path}.
    """
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)
    written: dict[str, str] = {}

    for busco_id in selected_buscos:
        records: list[SeqRecord] = []
        for r in results:
            fasta_path = r.sc_fasta_paths.get(busco_id)
            if fasta_path and fasta_path.exists():
                for rec in SeqIO.parse(str(fasta_path), "fasta"):
                    rec.id          = r.sample_id
                    rec.description = f"[busco={busco_id}] [sample={r.sample_id}]"
                    records.append(rec)
                    break  # single-copy: take the first (and only) record
        if records:
            dest = str(out / f"{busco_id}.fasta")
            SeqIO.write(records, dest, "fasta")
            written[busco_id] = dest

    return written
