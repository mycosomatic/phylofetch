"""
busco_utils.py
--------------
Parse BUSCO v4/v5 and Compleasm output, build occupancy matrices, export
single-copy ortholog FASTAs, download NCBI genome assemblies, and run
BUSCO/Compleasm on those assemblies.

Genome-comparison workflow:
  1. download_ncbi_genome(accession, output_dir)   → fasta path
  2. run_busco()  / run_compleasm()                → BuscoResult
  3. build_occupancy_matrix(results)               → pd.DataFrame
  4. export_sc_fastas(results, buscos, output_dir) → {busco_id: fasta_path}
"""

from __future__ import annotations

import re
import shutil
import subprocess
import zipfile
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
    tool: str           # "busco" | "compleasm"
    lineage: str
    complete: int = 0
    single_copy: int = 0
    duplicated: int = 0
    fragmented: int = 0
    missing: int = 0
    total: int = 0
    completeness_pct: float = 0.0
    # Map busco_id → status (Complete/Duplicated/Fragmented/Missing)
    gene_status: dict[str, str] = field(default_factory=dict)
    # Map busco_id → Path of single-copy sequence file
    sc_fasta_paths: dict[str, Path] = field(default_factory=dict)


# ── Parsers ──────────────────────────────────────────────────────────────────

def _parse_busco_short_summary(path: Path) -> dict:
    """Parse a BUSCO v4/v5 short_summary*.txt file into a stats dict."""
    stats: dict = {
        "complete": 0, "single_copy": 0, "duplicated": 0,
        "fragmented": 0, "missing": 0, "total": 0,
        "lineage": "", "completeness_pct": 0.0,
    }
    text = path.read_text()
    for line in text.splitlines():
        line = line.strip()
        m = re.search(
            r"C:(\d+\.\d+)%\[S:(\d+\.\d+)%,D:(\d+\.\d+)%\],F:(\d+\.\d+)%,M:(\d+\.\d+)%,n:(\d+)",
            line,
        )
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
    """Parse full_table.tsv → {busco_id: status}."""
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
            if busco_id not in gene_status or status.startswith("Complete"):
                gene_status[busco_id] = status
    return gene_status


def _collect_sc_fastas(run_dir: Path, seq_type: str = "nucleotide") -> dict[str, Path]:
    """Locate single_copy_busco_sequences/ FASTA files indexed by BUSCO ID."""
    sc_dir = run_dir / "busco_sequences" / "single_copy_busco_sequences"
    if not sc_dir.exists():
        sc_dir = run_dir / "single_copy_busco_sequences"
    if not sc_dir.exists():
        return {}
    ext = ".fna" if seq_type == "nucleotide" else ".faa"
    paths: dict[str, Path] = {}
    for fp in sc_dir.glob(f"*{ext}"):
        paths[fp.stem] = fp
    if not paths:
        for fp in sc_dir.glob("*.fna"):
            paths[fp.stem] = fp
        for fp in sc_dir.glob("*.faa"):
            paths.setdefault(fp.stem, fp)
    return paths


def _parse_compleasm(run_dir: Path) -> tuple[dict, dict[str, str]]:
    """Parse Compleasm summary.txt and mb_table.txt."""
    stats: dict = {
        "complete": 0, "single_copy": 0, "duplicated": 0,
        "fragmented": 0, "missing": 0, "total": 0,
        "lineage": "", "completeness_pct": 0.0,
    }
    gene_status: dict[str, str] = {}

    summary = run_dir / "summary.txt"
    if summary.exists():
        for line in summary.read_text().splitlines():
            m = re.match(r"S:(\d+\.\d+)%", line)
            if m:
                stats["completeness_pct"] = float(m.group(1))
            for label, key in [
                ("Single", "single_copy"), ("Duplicated", "duplicated"),
                ("Fragmented", "fragmented"), ("Missing", "missing"),
            ]:
                if f"{label}:" in line:
                    n = re.search(r"\d+", line.split(f"{label}:")[1])
                    if n:
                        stats[key] = int(n.group())
        stats["complete"] = stats["single_copy"] + stats["duplicated"]
        stats["total"] = stats["complete"] + stats["fragmented"] + stats["missing"]

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


# ── Main scanner ──────────────────────────────────────────────────────────────

def scan_busco_run(run_dir: Path | str, sample_id: Optional[str] = None,
                   seq_type: str = "nucleotide") -> Optional[BuscoResult]:
    """
    Detect and parse a BUSCO or Compleasm run directory.
    sample_id defaults to the directory name if not provided.
    """
    run_dir = Path(run_dir).resolve()
    if not run_dir.is_dir():
        return None
    sid = sample_id or run_dir.name

    # BUSCO v5 layout: short_summary*.txt
    short_summaries = (
        list(run_dir.glob("short_summary*.txt")) +
        list(run_dir.glob("*/short_summary*.txt"))
    )
    full_tables = (
        list(run_dir.glob("full_table.tsv")) +
        list(run_dir.glob("run_*/full_table.tsv"))
    )

    if short_summaries:
        ss = short_summaries[0]
        stats = _parse_busco_short_summary(ss)
        ft = full_tables[0] if full_tables else None
        gene_status = _parse_busco_full_table(ft) if ft else {}
        sc_base = ss.parent if ss.parent != run_dir else run_dir
        sc_paths = _collect_sc_fastas(sc_base, seq_type)
        return BuscoResult(
            sample_id=sid, run_dir=run_dir, tool="busco",
            lineage=stats.get("lineage", ""),
            complete=stats["complete"], single_copy=stats["single_copy"],
            duplicated=stats["duplicated"], fragmented=stats["fragmented"],
            missing=stats["missing"], total=stats["total"],
            completeness_pct=stats.get("completeness_pct", 0.0),
            gene_status=gene_status, sc_fasta_paths=sc_paths,
        )

    # Compleasm layout
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


# ── NCBI genome download ──────────────────────────────────────────────────────

def download_ncbi_genome(
    accession: str,
    output_dir: str | Path,
    datasets_bin: str = "datasets",
) -> tuple[int, str, Optional[Path]]:
    """
    Download a genome assembly by GCA/GCF accession using the NCBI Datasets CLI.

    Returns (returncode, log_text, fasta_path).
    fasta_path is None on failure.

    Requires the `datasets` CLI (conda: ncbi-datasets-cli or
    https://www.ncbi.nlm.nih.gov/datasets/docs/v2/download-and-install/).
    """
    if shutil.which(datasets_bin) is None:
        return (
            1,
            f"'{datasets_bin}' not found. Install via: conda install -c conda-forge ncbi-datasets-cli",
            None,
        )

    out = Path(output_dir)
    acc_dir = out / accession
    acc_dir.mkdir(parents=True, exist_ok=True)

    zip_path = acc_dir / "genome.zip"
    log_lines: list[str] = []

    # Download
    cmd = [
        datasets_bin, "download", "genome", "accession", accession,
        "--include", "genome",
        "--filename", str(zip_path),
    ]
    log_lines.append(f"$ {' '.join(cmd)}")
    try:
        result = subprocess.run(
            cmd, capture_output=True, text=True, timeout=600
        )
    except subprocess.TimeoutExpired:
        return 1, f"Download timed out after 10 min for {accession}", None
    except Exception as exc:
        return 1, str(exc), None

    log_lines.append(result.stdout)
    if result.stderr:
        log_lines.append(result.stderr)

    if result.returncode != 0:
        return result.returncode, "\n".join(log_lines), None

    # Unzip
    try:
        with zipfile.ZipFile(zip_path) as zf:
            zf.extractall(acc_dir)
    except Exception as exc:
        log_lines.append(f"Unzip error: {exc}")
        return 1, "\n".join(log_lines), None

    # Find the FASTA (NCBI Datasets v2 layout: ncbi_dataset/data/<ACC>/*.fna)
    fna_files = sorted(acc_dir.rglob("*.fna"))
    if not fna_files:
        fna_files = sorted(acc_dir.rglob("*.fasta"))
    if not fna_files:
        log_lines.append(f"No .fna / .fasta found after unzipping {accession}")
        return 1, "\n".join(log_lines), None

    fasta_path = fna_files[0]
    log_lines.append(f"Assembly FASTA: {fasta_path}")
    return 0, "\n".join(log_lines), fasta_path


# ── BUSCO / Compleasm runners ─────────────────────────────────────────────────

def run_busco(
    assembly_fasta: str | Path,
    lineage: str,
    output_dir: str | Path,
    sample_id: str,
    cpu: int = 4,
    busco_bin: str = "busco",
) -> tuple[int, str, Optional[BuscoResult]]:
    """
    Run BUSCO on an assembly. Returns (returncode, log_text, BuscoResult|None).

    lineage examples: sordariomycota_odb10, ascomycota_odb10, basidiomycota_odb10
    """
    out = Path(output_dir)
    run_out = out / sample_id
    run_out.mkdir(parents=True, exist_ok=True)

    cmd = [
        busco_bin,
        "-i", str(assembly_fasta),
        "-l", lineage,
        "-o", sample_id,
        "--out_path", str(out),
        "-m", "genome",
        "--cpu", str(cpu),
        "--force",
    ]
    log_lines = [f"$ {' '.join(cmd)}"]

    try:
        result = subprocess.run(
            cmd, capture_output=True, text=True, timeout=3600
        )
    except subprocess.TimeoutExpired:
        return 1, f"BUSCO timed out after 60 min for {sample_id}", None
    except FileNotFoundError:
        return 1, f"'{busco_bin}' not found. Install via: conda install -c bioconda busco", None
    except Exception as exc:
        return 1, str(exc), None

    log_lines.append(result.stdout[-5000:] if len(result.stdout) > 5000 else result.stdout)
    if result.stderr:
        log_lines.append(result.stderr[-2000:] if len(result.stderr) > 2000 else result.stderr)

    if result.returncode != 0:
        return result.returncode, "\n".join(log_lines), None

    busco_result = scan_busco_run(run_out, sample_id=sample_id)
    return result.returncode, "\n".join(log_lines), busco_result


def run_compleasm(
    assembly_fasta: str | Path,
    lineage: str,
    output_dir: str | Path,
    sample_id: str,
    cpu: int = 4,
    compleasm_bin: str = "compleasm",
) -> tuple[int, str, Optional[BuscoResult]]:
    """
    Run Compleasm on an assembly. Returns (returncode, log_text, BuscoResult|None).

    lineage examples: sordariomycota, ascomycota, basidiomycota
    (no _odb10 suffix — Compleasm downloads its own databases).
    """
    out = Path(output_dir)
    run_out = out / sample_id
    run_out.mkdir(parents=True, exist_ok=True)

    cmd = [
        compleasm_bin, "run",
        "-a", str(assembly_fasta),
        "-o", str(run_out),
        "-l", lineage,
        "-t", str(cpu),
    ]
    log_lines = [f"$ {' '.join(cmd)}"]

    try:
        result = subprocess.run(
            cmd, capture_output=True, text=True, timeout=3600
        )
    except subprocess.TimeoutExpired:
        return 1, f"Compleasm timed out after 60 min for {sample_id}", None
    except FileNotFoundError:
        return 1, f"'{compleasm_bin}' not found. Install via: pip install compleasm", None
    except Exception as exc:
        return 1, str(exc), None

    log_lines.append(result.stdout[-5000:] if len(result.stdout) > 5000 else result.stdout)
    if result.stderr:
        log_lines.append(result.stderr[-2000:] if len(result.stderr) > 2000 else result.stderr)

    if result.returncode != 0:
        return result.returncode, "\n".join(log_lines), None

    busco_result = scan_busco_run(run_out, sample_id=sample_id)
    return result.returncode, "\n".join(log_lines), busco_result


# ── Occupancy matrix ──────────────────────────────────────────────────────────

def build_occupancy_matrix(results: list[BuscoResult]) -> pd.DataFrame:
    """
    Build a sample × BUSCO occupancy matrix.
    Values: 'S' (single-copy complete), 'D' (duplicated), 'F' (fragmented), 'M' (missing).
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
    """Return BUSCO IDs that are single-copy in at least min_occupancy fraction of samples."""
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


def export_sc_fastas(
    results: list[BuscoResult],
    selected_buscos: list[str],
    output_dir: str | Path,
    seq_type: str = "nucleotide",
) -> dict:
    """
    Collect single-copy FASTA sequences for each selected BUSCO across all samples.
    Writes one FASTA per BUSCO with one sequence per sample.

    Returns {
        "exported_loci": N,
        "n_missing": M,
        "missing_sc_paths": ["sample/busco_id", ...],
        "fasta_paths": {busco_id: str(path)},
    }
    """
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)
    fasta_paths: dict[str, str] = {}
    missing_sc_paths: list[str] = []

    for busco_id in selected_buscos:
        records: list[SeqRecord] = []
        for r in results:
            # Re-collect SC paths with correct seq_type if not already done
            if not r.sc_fasta_paths:
                r.sc_fasta_paths = _collect_sc_fastas(r.run_dir, seq_type)
            fasta_path = r.sc_fasta_paths.get(busco_id)
            if fasta_path and fasta_path.exists():
                for rec in SeqIO.parse(str(fasta_path), "fasta"):
                    rec.id          = r.sample_id
                    rec.description = f"[busco={busco_id}] [sample={r.sample_id}]"
                    records.append(rec)
                    break  # single-copy: take first record
            else:
                missing_sc_paths.append(f"{r.sample_id}/{busco_id}")

        if records:
            dest = str(out / f"{busco_id}.fasta")
            SeqIO.write(records, dest, "fasta")
            fasta_paths[busco_id] = dest

    return {
        "exported_loci": len(fasta_paths),
        "n_missing": len(missing_sc_paths),
        "missing_sc_paths": missing_sc_paths,
        "fasta_paths": fasta_paths,
    }


# ── Common lineage hints ──────────────────────────────────────────────────────

BUSCO_LINEAGE_HINTS: list[str] = [
    "ascomycota_odb10",
    "basidiomycota_odb10",
    "fungi_odb10",
    "sordariomycota_odb10",
    "dothideomycetes_odb10",
    "eurotiomycetes_odb10",
    "leotiomycetes_odb10",
    "hypocreales_odb10",
    "pleosporales_odb10",
    "agaricales_odb10",
    "boletales_odb10",
]

COMPLEASM_LINEAGE_HINTS: list[str] = [
    "ascomycota",
    "basidiomycota",
    "fungi",
    "sordariomycota",
    "dothideomycetes",
    "eurotiomycetes",
    "leotiomycetes",
    "hypocreales",
    "pleosporales",
    "agaricales",
    "boletales",
]
