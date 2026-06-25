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


import re as _re

# ITSx output file suffixes — one file per detected region
ITSX_SUFFIXES: dict[str, str] = {
    "ITS1":     ".ITS1.fasta",
    "ITS2":     ".ITS2.fasta",
    "ITS_full": ".full.fasta",
    "SSU":      ".SSU.fasta",
    "LSU":      ".LSU.fasta",
}

# Prefer the rDNA detection on the highest-coverage contig (D-028). The functional rDNA is a
# high-copy tandem array, so it assembles into a contig whose k-mer coverage is many-fold the
# single-copy genomic background. A detection on a low-coverage contig is therefore an *off-array*
# copy — a dispersed/orphan rDNA fragment that may be **RIP-pseudogenized** (the genomic analogue
# of the low-abundance pseudogene ITS amplicons seen in PCR) — or an outright spurious HMM hit on a
# chromosomal contig (which is where ITSx's 60 kb "SSU" / mislabelled regions come from). Within one
# region, keep detections whose contig coverage is ≥ max_cov / RDNA_COV_RATIO; drop the rest.
RDNA_COV_RATIO = 5.0

_COV_RE = _re.compile(r"_cov_([0-9]+(?:\.[0-9]+)?)", _re.IGNORECASE)


def parse_coverage(contig_name: str) -> Optional[float]:
    """k-mer coverage parsed from a SPAdes/pilon-style contig name (``..._cov_45.8_pilon``),
    or ``None`` when the name carries no coverage token (non-SPAdes assembler → no filtering)."""
    m = _COV_RE.search(contig_name or "")
    return float(m.group(1)) if m else None

# ITSx runs HMMER's hmmscan, which aborts on any target sequence > 100 kb
# ("Target sequence length > 100K, over comparison pipeline limit"). Genome-assembly contigs
# are routinely Mb-scale, so ITSx fails outright on assemblies. We chunk long contigs into
# overlapping windows before ITSx: the rDNA cistron (≈5–9 kb; ITS proper < 1 kb) is far smaller
# than the overlap, so every rDNA region stays fully intact in at least one chunk; duplicates
# from the overlap are removed by de-duplicating identical output sequences. (LXD-003)
ITSX_MAX_CONTIG_LEN = 90_000      # safely under hmmscan's 100 kb limit
ITSX_CHUNK_OVERLAP = 20_000       # ≫ longest rDNA subunit, so no rDNA region is split


def chunk_long_contigs(records, max_len: int = ITSX_MAX_CONTIG_LEN,
                       overlap: int = ITSX_CHUNK_OVERLAP):
    """
    Split any contig longer than ``max_len`` into overlapping windows so ITSx/hmmscan can run.
    Contigs ≤ max_len pass through unchanged. Chunk ids are ``<contig>__c<start>`` (the offset
    is preserved for provenance). Returns ``(records_out, was_chunked)``.
    """
    out = []
    changed = False
    step = max(max_len - overlap, 1)
    for rec in records:
        n = len(rec.seq)
        if n <= max_len:
            out.append(rec)
            continue
        changed = True
        start = 0
        while start < n:
            end = min(start + max_len, n)
            sub = rec[start:end]
            sub.id = f"{rec.id}__c{start}"
            sub.name = sub.id
            sub.description = ""
            out.append(sub)
            if end >= n:
                break
            start += step
    return out, changed


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


def _orig_contig(rec_id: str) -> str:
    """The original assembly contig name from an ITSx record id (strip the ``|coords`` field and
    any ``__c<offset>`` chunk suffix added by :func:`chunk_long_contigs`)."""
    contig = rec_id.split("|")[0] if "|" in rec_id else rec_id
    if "__c" in contig:
        contig = contig.rsplit("__c", 1)[0]
    return contig


def _relabel_itsx_output(fasta_path: str, strain_id: str, locus_name: str,
                         itsx_version: str = "", prefer_high_cov: bool = True,
                         cov_ratio: float = RDNA_COV_RATIO) -> dict:
    """
    Relabel ITSx output headers to clean '>STRAIN_LOCUS' format with provenance, de-duplicate
    overlapping-chunk repeats, and (D-028) keep only detections on the highest-coverage contig.

    Returns ``{"kept": n, "dropped": [(contig, cov, length), ...]}`` — ``kept`` is the number of
    records written (0 ⇒ nothing usable for this region); ``dropped`` lists low-coverage off-array
    detections removed (for the run log).
    """
    records = list(SeqIO.parse(fasta_path, "fasta"))
    if not records:
        return {"kept": 0, "dropped": []}

    # De-duplicate identical sequences: overlapping chunks (from chunk_long_contigs) can
    # report the same rDNA region more than once. Keep the first occurrence, preserve order.
    seen_seq: set[str] = set()
    deduped = []
    for rec in records:
        key = str(rec.seq).upper()
        if key in seen_seq:
            continue
        seen_seq.add(key)
        deduped.append(rec)
    records = deduped

    # Coverage filter (D-028): drop detections whose contig coverage is far below the best — these
    # are off-array (possibly RIP'd) copies or spurious chromosomal HMM hits.
    dropped: list[tuple] = []
    if prefer_high_cov and len(records) > 1:
        covs = {id(r): parse_coverage(_orig_contig(r.id)) for r in records}
        known = [c for c in covs.values() if c is not None]
        if known:
            cutoff = max(known) / cov_ratio
            kept_records = []
            for r in records:
                c = covs[id(r)]
                if c is not None and c < cutoff:
                    dropped.append((_orig_contig(r.id), c, len(r.seq)))
                else:
                    kept_records.append(r)
            records = kept_records

    ts = _now_utc()
    for i, rec in enumerate(records):
        contig = _orig_contig(rec.id)
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
            "cov": f"{parse_coverage(contig):.0f}" if parse_coverage(contig) else None,
            "tool": "ITSx",
            "tool_version": itsx_version or None,
            "extracted": ts,
        })
        rec.id          = f"{strain_id}_{locus_name}{suffix}"
        rec.description = desc

    with open(fasta_path, "w") as f:
        SeqIO.write(records, f, "fasta")
    return {"kept": len(records), "dropped": dropped}


def run_itsx(
    assembly_fasta: str,
    output_dir: str,
    strain_id: str,
    threads: int = 4,
    itsx_bin: str = "ITSx",
    kingdom: str = "fungi",
    prefer_high_cov: bool = True,
    cov_ratio: float = RDNA_COV_RATIO,
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

    # Chunk long contigs so ITSx/hmmscan's 100 kb target limit isn't hit on genome
    # assemblies (LXD-003). Contigs ≤ limit pass through unchanged.
    records = list(SeqIO.parse(assembly_fasta, "fasta"))
    chunked, was_chunked = chunk_long_contigs(records)
    itsx_input = assembly_fasta
    if was_chunked:
        itsx_input = prefix + "_itsx_input.fasta"
        SeqIO.write(chunked, itsx_input, "fasta")

    cmd = [
        itsx_bin,
        "-i", itsx_input,
        "-o", prefix,
        "--cpu", str(threads),
        "--save_regions", "all",
        "--complement", "T",
        "--heuristics", "T",
        "--graphical", "F",
        "-t", kingdom,
    ]

    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=3600)
    except subprocess.TimeoutExpired:
        return 1, "ERROR: ITSx timed out after 3600s", {}
    except (FileNotFoundError, OSError) as exc:
        # which() was checked above, but the binary can still fail to exec (perms, broken
        # conda shim, missing hmmer/perl dep). run_itsx must always return its 3-tuple (D-033).
        return 1, f"ERROR: ITSx failed to launch: {exc}", {}
    all_output = result.stdout + result.stderr

    # Always surface the last 40 lines so UI can show them on failure or success
    log_lines = all_output.strip().splitlines()
    log_text  = "\n".join(log_lines[-40:]) if log_lines else "(no output)"

    # ITSx exits 0 even when HMMER aborts on an over-limit sequence — surface it as an error
    # so the caller doesn't read a silent empty result as "no rDNA detected".
    if "over comparison pipeline limit" in all_output or "HMMER run seems to have been" in all_output:
        return 1, ("ERROR: hmmscan hit its 100 kb sequence limit (a contig was too long even "
                   "after chunking).\n" + log_text), {}

    found: dict[str, str] = {}
    drop_notes: list[str] = []
    for locus, suffix in ITSX_SUFFIXES.items():
        path = prefix + suffix
        if os.path.exists(path) and os.path.getsize(path) > 0:
            info = _relabel_itsx_output(path, strain_id, locus, itsx_version,
                                        prefer_high_cov=prefer_high_cov, cov_ratio=cov_ratio)
            if info["kept"]:
                found[locus] = path
            for contig, cov, length in info["dropped"]:
                drop_notes.append(f"  {locus}: dropped off-array copy on {contig} "
                                  f"(cov {cov:.0f}, {length} bp) — low-coverage / possibly RIP'd")
    if drop_notes:
        log_text += ("\n\nrDNA coverage filter (D-028) — kept the high-coverage array, dropped:\n"
                     + "\n".join(drop_notes))

    return result.returncode, log_text, found


# ── rDNA region placement + combine (component-page layout, RM-007 step 4c) ───
#
# The extraction layout is <per_strain>/<strain>/<region>/<region>.fasta and
# <combined>/<region>_combined.fasta. `region` is the user-facing name (ITS / ITS1 / ITS2 /
# LSU / SSU); the full ITS maps to ITSx's "ITS_full" key. These helpers factor the logic the
# monolithic page did inline so the ITSx component page can reuse it (and it is unit-tested).

_RDNA_REGION_TO_ITSX_KEY = {
    "ITS": "ITS_full", "ITS1": "ITS1", "ITS2": "ITS2", "LSU": "LSU", "SSU": "SSU",
}


def place_rdna_regions(found: dict[str, str], strain_out: str,
                       regions: list[str]) -> dict[str, int]:
    """
    Copy the requested rDNA ``regions`` from an ITSx ``found`` map into
    ``<strain_out>/<region>/<region>.fasta``. Returns {region: n_records} (0 when ITSx did
    not detect that region).
    """
    counts: dict[str, int] = {}
    for region in regions:
        key = _RDNA_REGION_TO_ITSX_KEY.get(region, region)
        dest_dir = Path(strain_out) / region
        dest_dir.mkdir(parents=True, exist_ok=True)
        dest = dest_dir / f"{region}.fasta"
        src = found.get(key)
        if src and Path(src).exists():
            shutil.copy(src, dest)
            counts[region] = len(list(SeqIO.parse(str(dest), "fasta")))
        else:
            counts[region] = 0
    return counts


def combine_rdna_regions(per_strain_dir: str, combined_dir: str,
                         regions: list[str]) -> dict[str, tuple[str, int]]:
    """
    Merge ``<strain>/<region>/<region>.fasta`` across strains into
    ``<combined_dir>/<region>_combined.fasta``. Returns {region: (path, n_records)} for each
    region that yielded at least one sequence.
    """
    Path(combined_dir).mkdir(parents=True, exist_ok=True)
    result: dict[str, tuple[str, int]] = {}
    for region in regions:
        recs: list[SeqRecord] = []
        for sd in sorted(Path(per_strain_dir).iterdir()):
            if not sd.is_dir():
                continue
            fp = sd / region / f"{region}.fasta"
            if fp.exists():
                recs.extend(list(SeqIO.parse(str(fp), "fasta")))
        if recs:
            out = Path(combined_dir) / f"{region}_combined.fasta"
            SeqIO.write(recs, str(out), "fasta")
            result[region] = (str(out), len(recs))
    return result


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
