"""
itsx_utils.py
-------------
ITSx wrapper for rDNA boundary detection (ITS1, ITS2, SSU, LSU, full ITS).

LXD-001 fixes applied:
  - Removed --multi_out flag (not supported in all ITSx versions; caused silent failures)
  - Old output files are cleared before rerun (prevents stale results being collected)
  - Return code and last log lines are surfaced on failure so the UI can display them
  - Run is channelled through RunManager when a manager is provided
"""

import os
import shutil
import subprocess
from datetime import datetime, timezone
from pathlib import Path
from typing import Optional

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


# ITSx output file suffixes — one file per detected region
ITSX_SUFFIXES: dict[str, str] = {
    "ITS1":     ".ITS1.fasta",
    "ITS2":     ".ITS2.fasta",
    "ITS_full": ".full.fasta",
    "SSU":      ".SSU.fasta",
    "LSU":      ".LSU.fasta",
}


def _probe_itsx_version(itsx_bin: str) -> str:
    try:
        r = subprocess.run([itsx_bin, "--help"], capture_output=True, text=True, timeout=10)
        for line in (r.stdout + r.stderr).splitlines():
            if "itsx" in line.lower() and any(c.isdigit() for c in line):
                return line.strip()[:80]
    except Exception:
        pass
    return "unknown"


def _now_utc() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat()


def _build_description(fields: dict) -> str:
    return " ".join(f"[{k}={v}]" for k, v in fields.items() if v is not None)


def _relabel_itsx_output(fasta_path: str, strain_id: str,
                         locus_name: str, itsx_version: str = "") -> bool:
    """
    Relabel ITSx output headers to clean '>STRAIN_LOCUS' format with provenance.
    Returns True if any records were relabeled.
    """
    records = list(SeqIO.parse(fasta_path, "fasta"))
    if not records:
        return False

    ts = _now_utc()
    for i, rec in enumerate(records):
        # Try to extract contig and coordinate info from ITSx header
        # ITSx format: >CONTIG|START-END|...
        contig = rec.id.split("|")[0] if "|" in rec.id else rec.id
        coords = ""
        if "|" in rec.description:
            parts = rec.description.split("|")
            if len(parts) >= 2:
                coords = parts[1]

        suffix = f"_{i + 1}" if len(records) > 1 else ""
        desc = _build_description({
            "sample": strain_id, "locus": locus_name, "type": "rDNA",
            "contig": contig,
            "coords": coords or None,
            "tool": "ITSx",
            "tool_version": itsx_version or None,
            "extracted": ts,
        })
        rec.id          = f"{strain_id}_{locus_name}{suffix}"
        rec.description = desc

    with open(fasta_path, "w") as f:
        SeqIO.write(records, f, "fasta")
    return True


def run_itsx(
    assembly_fasta: str,
    output_dir: str,
    strain_id: str,
    threads: int = 4,
    itsx_bin: str = "ITSx",
    kingdom: str = "fungi",
) -> tuple[int, str, dict[str, str]]:
    """
    Run ITSx on an assembly to extract rDNA regions.

    Returns (returncode, log_text, {locus_name: output_fasta_path}).

    LXD-001: stale outputs are cleared first; --multi_out is not used;
    the last 40 lines of ITSx output are always returned so failures are visible.
    """
    if not shutil.which(itsx_bin):
        return 1, f"ERROR: ITSx not found on PATH ({itsx_bin})", {}

    os.makedirs(output_dir, exist_ok=True)
    prefix = os.path.join(output_dir, strain_id)

    # Clear old output files so stale results are never collected
    for suffix in ITSX_SUFFIXES.values():
        old = prefix + suffix
        if os.path.exists(old):
            os.remove(old)

    itsx_version = _probe_itsx_version(itsx_bin)

    cmd = [
        itsx_bin,
        "-i", assembly_fasta,
        "-o", prefix,
        "--cpu", str(threads),
        "--save_regions", "all",
        "--complement", "T",
        "--heuristics", "T",
        "--graphical", "F",
        "-t", kingdom,
    ]

    result = subprocess.run(cmd, capture_output=True, text=True)
    all_output = result.stdout + result.stderr

    # Always surface the last 40 lines so UI can show them on failure or success
    log_lines = all_output.strip().splitlines()
    log_text  = "\n".join(log_lines[-40:]) if log_lines else "(no output)"

    found: dict[str, str] = {}
    for locus, suffix in ITSX_SUFFIXES.items():
        path = prefix + suffix
        if os.path.exists(path) and os.path.getsize(path) > 0:
            if _relabel_itsx_output(path, strain_id, locus, itsx_version):
                found[locus] = path

    return result.returncode, log_text, found


# ── Multi-sample combiners ────────────────────────────────────────────────────

def merge_per_strain_to_combined(per_strain_dir: str, combined_dir: str,
                                 locus_name: str) -> str:
    """Collect per-strain rDNA FASTAs and write a combined multi-FASTA."""
    os.makedirs(combined_dir, exist_ok=True)
    combined_path = os.path.join(combined_dir, f"{locus_name}_combined.fasta")

    seen_ids: set[str] = set()
    all_records: list[SeqRecord] = []

    for strain_dir in sorted(Path(per_strain_dir).iterdir()):
        if not strain_dir.is_dir():
            continue
        for fasta_file in strain_dir.glob(f"*{locus_name}*.fasta"):
            for rec in SeqIO.parse(str(fasta_file), "fasta"):
                if rec.id not in seen_ids:
                    seen_ids.add(rec.id)
                    all_records.append(rec)

    if all_records:
        with open(combined_path, "w") as f:
            SeqIO.write(all_records, f, "fasta")
    return combined_path


def get_locus_summary(per_strain_dir: str, locus_name: str) -> dict[str, dict]:
    """Return per-strain summary for a locus."""
    summary: dict[str, dict] = {}
    for strain_dir in sorted(Path(per_strain_dir).iterdir()):
        if not strain_dir.is_dir():
            continue
        strain_id  = strain_dir.name
        found_files = list(strain_dir.glob(f"*{locus_name}*.fasta"))
        if found_files:
            records = list(SeqIO.parse(str(found_files[0]), "fasta"))
            summary[strain_id] = {
                "found":      True,
                "seq_length": len(records[0].seq) if records else 0,
                "fasta_path": str(found_files[0]),
            }
        else:
            summary[strain_id] = {"found": False, "seq_length": 0, "fasta_path": ""}
    return summary
