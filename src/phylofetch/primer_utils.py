"""
primer_utils.py
---------------
PCR primer-based locus extraction from genome assemblies ("in-silico PCR").

Strategy: run blastn-short to locate forward and reverse primer binding sites
in assembly contigs with an edit-distance tolerance (max_mismatches that counts
both substitutions and unaligned primer bases), then extract the amplicon
between them on the same contig and strand. Useful when NCBI lacks reference
accessions that map cleanly to the assembly, or when coverage is too sparse for
the reference-CDS BLAST strategy.

Primer orientation conventions (consistent with wetlab usage):
  fwd — 5'->3' sequence that anneals to the minus (antisense) strand
         -> blastn reports a PLUS-strand hit (sstart < send)
  rev — 5'->3' sequence that anneals to the plus (sense) strand
         -> blastn reports a MINUS-strand hit (sstart > send)
When the amplicon itself lies on the minus strand the roles are reversed.

Primer library
--------------
The built-in catalogue ships as packaged data (``phylofetch/data/primers.json``)
and is fully citable: every pair carries a primary-literature ``source``,
``citation`` and ``reference_url``. Users may add their own pairs, which persist
in ``~/.phylofetch/primers.json`` and are merged on top of the built-in set.

All external blastn calls are routed through RunManager when a manager is
supplied, so each primer search is logged with its command and tool version,
and a self-contained ``LOCUS_extraction.log`` is written alongside the FASTA.
"""

from __future__ import annotations

import json
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

# User-editable primer library (merged on top of the packaged catalogue).
USER_PRIMER_PATH = Path.home() / ".phylofetch" / "primers.json"


@dataclass
class PrimerPair:
    name: str
    locus: str
    fwd: str          # 5'->3' forward primer (including IUPAC degenerate codes)
    rev: str          # 5'->3' reverse primer (including IUPAC degenerate codes)
    min_amplicon: int = 100
    max_amplicon: int = 5000
    fwd_name: str = ""
    rev_name: str = ""
    source: str = ""          # short attribution, e.g. "White et al. 1990"
    citation: str = ""        # full citation
    reference_url: str = ""   # DOI or URL
    origin: str = "built-in"  # "built-in" or "user"

    def to_record(self) -> dict:
        """Serialisable dict for the on-disk library (excludes name/origin)."""
        return {
            "locus": self.locus,
            "fwd_name": self.fwd_name,
            "fwd": self.fwd,
            "rev_name": self.rev_name,
            "rev": self.rev,
            "min_amplicon": self.min_amplicon,
            "max_amplicon": self.max_amplicon,
            "source": self.source,
            "citation": self.citation,
            "reference_url": self.reference_url,
        }


# ── Library loading / persistence ─────────────────────────────────────────────

def _builtin_data_path() -> Path:
    return Path(__file__).parent / "data" / "primers.json"


def _load_builtin_raw() -> dict:
    """Load the packaged primers.json, tolerating install or src-tree layout."""
    try:                                   # installed package
        from importlib.resources import files
        text = (files("phylofetch") / "data" / "primers.json").read_text(encoding="utf-8")
    except (ModuleNotFoundError, FileNotFoundError, AttributeError, TypeError):
        text = _builtin_data_path().read_text(encoding="utf-8")
    return json.loads(text)


def _pair_from_record(name: str, d: dict, origin: str) -> PrimerPair:
    return PrimerPair(
        name=name,
        locus=str(d["locus"]).strip().upper(),
        fwd=str(d["fwd"]).strip().upper(),
        rev=str(d["rev"]).strip().upper(),
        min_amplicon=int(d.get("min_amplicon", 100)),
        max_amplicon=int(d.get("max_amplicon", 5000)),
        fwd_name=str(d.get("fwd_name", "")),
        rev_name=str(d.get("rev_name", "")),
        source=str(d.get("source", "")),
        citation=str(d.get("citation", "")),
        reference_url=str(d.get("reference_url", "")),
        origin=origin,
    )


def _records_of(raw: dict) -> dict:
    """Accept both {"primers": {...}} and a flat {name: record} mapping."""
    if isinstance(raw, dict) and isinstance(raw.get("primers"), dict):
        return raw["primers"]
    return raw if isinstance(raw, dict) else {}


def load_builtin_primers() -> dict[str, PrimerPair]:
    """Return the packaged, citable catalogue."""
    out: dict[str, PrimerPair] = {}
    for name, rec in _records_of(_load_builtin_raw()).items():
        if name.startswith("_"):
            continue
        try:
            out[name] = _pair_from_record(name, rec, "built-in")
        except (KeyError, TypeError, ValueError):
            continue
    return out


def load_user_primers(path: Path = USER_PRIMER_PATH) -> dict[str, PrimerPair]:
    """Return the user's saved pairs (empty if none/invalid)."""
    p = Path(path)
    if not p.exists():
        return {}
    try:
        raw = json.loads(p.read_text(encoding="utf-8"))
    except (json.JSONDecodeError, OSError):
        return {}
    out: dict[str, PrimerPair] = {}
    for name, rec in _records_of(raw).items():
        if name.startswith("_"):
            continue
        try:
            out[name] = _pair_from_record(name, rec, "user")
        except (KeyError, TypeError, ValueError):
            continue
    return out


def get_primer_catalogue(include_user: bool = True,
                         user_path: Path = USER_PRIMER_PATH) -> dict[str, PrimerPair]:
    """Built-in catalogue with user pairs merged on top (user wins on name clash)."""
    cat = dict(load_builtin_primers())
    if include_user:
        cat.update(load_user_primers(user_path))
    return cat


def locus_primer_map(catalogue: Optional[dict[str, PrimerPair]] = None) -> dict[str, list[str]]:
    """Group catalogue pair names by locus for UI lookup."""
    cat = catalogue if catalogue is not None else get_primer_catalogue()
    out: dict[str, list[str]] = {}
    for pp in cat.values():
        out.setdefault(pp.locus, []).append(pp.name)
    return out


def save_user_primer(pair: PrimerPair, path: Path = USER_PRIMER_PATH) -> Path:
    """Add/overwrite a pair in the user's library on disk. Returns the path."""
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    records: dict = {}
    if p.exists():
        try:
            records = _records_of(json.loads(p.read_text(encoding="utf-8")))
        except (json.JSONDecodeError, OSError):
            records = {}
    records[pair.name] = pair.to_record()
    p.write_text(
        json.dumps({"schema_version": 1, "primers": records}, indent=2, sort_keys=True),
        encoding="utf-8",
    )
    return p


def delete_user_primer(name: str, path: Path = USER_PRIMER_PATH) -> bool:
    """Remove a pair from the user's library. Returns True if it existed."""
    p = Path(path)
    if not p.exists():
        return False
    try:
        records = _records_of(json.loads(p.read_text(encoding="utf-8")))
    except (json.JSONDecodeError, OSError):
        return False
    if name not in records:
        return False
    del records[name]
    p.write_text(
        json.dumps({"schema_version": 1, "primers": records}, indent=2, sort_keys=True),
        encoding="utf-8",
    )
    return True


# Module-level built-in catalogue (back-compat; user pairs are merged via
# get_primer_catalogue() at call sites that need them).
PRIMER_CATALOGUE: dict[str, PrimerPair] = load_builtin_primers()
LOCUS_PRIMER_MAP: dict[str, list[str]] = locus_primer_map(PRIMER_CATALOGUE)


# ── blastn-short runner ────────────────────────────────────────────────────────

_BLAST_FIELDS = [
    "qseqid", "sseqid", "pident", "length", "mismatch",
    "qstart", "qend", "sstart", "send", "evalue", "bitscore",
]


def _parse_blast_rows(text: str) -> list[dict]:
    rows: list[dict] = []
    for line in text.splitlines():
        parts = line.split("\t")
        if len(parts) < len(_BLAST_FIELDS):
            continue
        try:
            rows.append({
                "qseqid":   parts[0],
                "sseqid":   parts[1],
                "pident":   float(parts[2]),
                "length":   int(parts[3]),
                "mismatch": int(parts[4]),
                "qstart":   int(parts[5]),
                "qend":     int(parts[6]),
                "sstart":   int(parts[7]),
                "send":     int(parts[8]),
                "evalue":   float(parts[9]),
                "bitscore": float(parts[10]),
            })
        except (ValueError, IndexError):
            continue
    return rows


def _run_blastn_short(
    primer_fa: str,
    assembly_fa: str,
    perc_identity: float,
    blastn_bin: str = "blastn",
    manager=None,
    module: str = "loci_extraction",
    action: str = "primer_search",
) -> list[dict]:
    """
    Run blastn-short and return parsed HSP dicts.

    When ``manager`` (a RunManager) is given, the command is executed and logged
    through it (command + blastn version captured in the run folder); otherwise
    it runs directly via subprocess.
    """
    cmd = [
        blastn_bin, "-task", "blastn-short",
        "-query", primer_fa,
        "-subject", assembly_fa,
        "-perc_identity", f"{perc_identity:.1f}",
        "-outfmt", "6 " + " ".join(_BLAST_FIELDS),
        "-word_size", "7",
        "-num_alignments", "1000",
    ]

    if manager is not None:
        rr = manager.run(
            cmd, module=module, action=action,
            tool_version_keys=["blastn"],
            inputs={"primers": primer_fa, "assembly": assembly_fa},
            params={"task": "blastn-short", "perc_identity": round(perc_identity, 1)},
        )
        if rr.returncode != 0:
            return []
        try:
            return _parse_blast_rows(Path(rr.stdout_path).read_text(encoding="utf-8"))
        except OSError:
            return []

    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode != 0:
        return []
    return _parse_blast_rows(proc.stdout)


# ── Amplicon discovery ─────────────────────────────────────────────────────────

def _effective_mismatch(hit: dict, primer_len: int) -> int:
    """
    Edit distance = aligned substitutions + unaligned primer bases.

    blastn's ``mismatch`` only counts within the aligned region, so a primer
    that aligns partially (e.g. a 3' truncation) can report mismatch=0. Counting
    the unaligned tail closes that gap and makes max_mismatches meaningful.
    """
    aligned = int(hit.get("length", primer_len))
    truncation = max(0, primer_len - aligned)
    return int(hit.get("mismatch", 0)) + truncation


def find_primer_amplicons(
    assembly_fasta: str,
    primer_pair: PrimerPair,
    max_mismatches: int = 2,
    blastn_bin: str = "blastn",
    manager=None,
    module: str = "loci_extraction",
    action: str = "primer_search",
) -> list[dict]:
    """
    Search an assembly for amplicons bounded by primer_pair.

    Each primer binding site must have an edit distance <= max_mismatches
    (substitutions + unaligned bases). Returns candidates sorted by ascending
    total edit distance. Each candidate dict:
      contig, strand, amp_start, amp_end, amp_len,
      fwd_mismatch, rev_mismatch, fwd_edit, rev_edit, total_edit,
      fwd_hit, rev_hit
    """
    fwd_len = len(primer_pair.fwd.replace(" ", ""))
    rev_len = len(primer_pair.rev.replace(" ", ""))
    min_len = max(1, min(fwd_len, rev_len))
    # Permissive BLAST pre-filter; the real budget is enforced below on the
    # full edit distance so partial alignments cannot sneak through.
    perc_id = max(60.0, (min_len - max_mismatches - 1) / min_len * 100)

    with tempfile.TemporaryDirectory() as tdir:
        primer_fa = os.path.join(tdir, "primers.fasta")
        with open(primer_fa, "w") as f:
            f.write(f">FWD\n{primer_pair.fwd}\n>REV\n{primer_pair.rev}\n")
        hits = _run_blastn_short(
            primer_fa, assembly_fasta, perc_id, blastn_bin,
            manager=manager, module=module, action=action,
        )

    fwd_by_contig: dict[str, list] = defaultdict(list)
    rev_by_contig: dict[str, list] = defaultdict(list)
    for h in hits:
        if h["qseqid"] == "FWD" and _effective_mismatch(h, fwd_len) <= max_mismatches:
            fwd_by_contig[h["sseqid"]].append(h)
        elif h["qseqid"] == "REV" and _effective_mismatch(h, rev_len) <= max_mismatches:
            rev_by_contig[h["sseqid"]].append(h)

    candidates: list[dict] = []
    for contig in set(fwd_by_contig) & set(rev_by_contig):
        for fh in fwd_by_contig[contig]:
            for rh in rev_by_contig[contig]:
                # + strand: F hit plus (sstart < send), R hit minus (sstart > send)
                if fh["sstart"] < fh["send"] and rh["sstart"] > rh["send"]:
                    amp_start, amp_end, strand = fh["sstart"], rh["sstart"], "+"
                # - strand: F hit minus (sstart > send), R hit plus (sstart < send)
                elif fh["sstart"] > fh["send"] and rh["sstart"] < rh["send"]:
                    amp_start, amp_end, strand = rh["sstart"], fh["sstart"], "-"
                else:
                    continue

                amp_len = amp_end - amp_start + 1
                if not (primer_pair.min_amplicon <= amp_len <= primer_pair.max_amplicon):
                    continue

                fwd_edit = _effective_mismatch(fh, fwd_len)
                rev_edit = _effective_mismatch(rh, rev_len)
                candidates.append({
                    "contig":       contig,
                    "strand":       strand,
                    "amp_start":    amp_start,
                    "amp_end":      amp_end,
                    "amp_len":      amp_len,
                    "fwd_mismatch": fh["mismatch"],
                    "rev_mismatch": rh["mismatch"],
                    "fwd_edit":     fwd_edit,
                    "rev_edit":     rev_edit,
                    "total_edit":   fwd_edit + rev_edit,
                    "fwd_hit":      fh,
                    "rev_hit":      rh,
                })

    candidates.sort(key=lambda c: (c["total_edit"], c["amp_len"]))
    return candidates


# ── Amplicon extraction & logging ──────────────────────────────────────────────

def _now_utc() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def extract_primer_amplicon(
    assembly_fasta: str,
    hit: dict,
    strain_id: str,
    locus_name: str,
    output_dir: str,
    primer_pair: PrimerPair,
) -> Optional[SeqRecord]:
    """Extract and write the amplicon FASTA; return the SeqRecord or None."""
    contig_seq: Optional[Seq] = None
    for rec in SeqIO.parse(assembly_fasta, "fasta"):
        if rec.id == hit["contig"]:
            contig_seq = rec.seq
            break
    if contig_seq is None:
        return None

    # BLAST coords are 1-based, inclusive.
    extracted = contig_seq[hit["amp_start"] - 1 : hit["amp_end"]]
    if hit["strand"] == "-":
        extracted = extracted.reverse_complement()

    description = (
        f"[sample={strain_id}] [locus={locus_name}] [type=amplicon] "
        f"[primer={primer_pair.name}] [primer_source={primer_pair.source or 'n/a'}] "
        f"[contig={hit['contig']}] "
        f"[coords={hit['amp_start']}..{hit['amp_end']}] [strand={hit['strand']}] "
        f"[amplicon_len={hit['amp_len']}bp] "
        f"[fwd_mm={hit.get('fwd_edit', hit['fwd_mismatch'])}] "
        f"[rev_mm={hit.get('rev_edit', hit['rev_mismatch'])}] "
        f"[extracted={_now_utc()}]"
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


def write_primer_extraction_log(
    output_dir: str,
    strain_id: str,
    locus_name: str,
    primer_pair: PrimerPair,
    hit: Optional[dict],
    max_mismatches: int,
    n_candidates: int,
    status: str,
    run_dir: str = "",
) -> str:
    """Self-contained per-locus log: primer, citation, result, run provenance."""
    lines = [
        "# phylofetch primer (in-silico PCR) extraction log",
        f"# Generated: {_now_utc()}",
        "",
        f"[sample]        {strain_id}",
        f"[locus]         {locus_name}",
        f"[type]          amplicon (in-silico PCR)",
        "",
        f"[primer_pair]   {primer_pair.name}",
        f"[fwd]           {primer_pair.fwd_name or 'fwd'}  5'-{primer_pair.fwd}-3'",
        f"[rev]           {primer_pair.rev_name or 'rev'}  5'-{primer_pair.rev}-3'",
        f"[size_window]   {primer_pair.min_amplicon}-{primer_pair.max_amplicon} bp",
        f"[max_mismatch]  {max_mismatches} (substitutions + unaligned bases, per primer)",
        f"[source]        {primer_pair.source or 'unknown'}",
        f"[citation]      {primer_pair.citation or 'n/a'}",
        f"[reference]     {primer_pair.reference_url or 'n/a'}",
        "",
        "[result]",
        f"  status        {status}",
        f"  sites_found   {n_candidates}",
    ]
    if hit is not None:
        lines += [
            f"  contig        {hit['contig']}",
            f"  coords        {hit['amp_start']}..{hit['amp_end']}",
            f"  strand        {hit['strand']}",
            f"  amplicon_len  {hit['amp_len']} bp",
            f"  fwd_edit_dist {hit.get('fwd_edit', hit.get('fwd_mismatch'))}",
            f"  rev_edit_dist {hit.get('rev_edit', hit.get('rev_mismatch'))}",
        ]
    if run_dir:
        lines += ["", f"[run_folder]    {run_dir}",
                  f"[terminal_log]  {os.path.join(run_dir, 'terminal.log')}"]

    Path(output_dir).mkdir(parents=True, exist_ok=True)
    log_path = os.path.join(output_dir, f"{locus_name}_extraction.log")
    Path(log_path).write_text("\n".join(lines) + "\n", encoding="utf-8")
    return log_path


def _result_dict(best: dict, primer_pair: PrimerPair, out_fasta: str,
                 n_candidates: int) -> dict:
    return {
        "contig":       best["contig"],
        "strand":       best["strand"],
        "amp_start":    best["amp_start"],
        "amp_end":      best["amp_end"],
        "amp_len":      best["amp_len"],
        "fwd_mismatch": best["fwd_mismatch"],
        "rev_mismatch": best["rev_mismatch"],
        "fwd_edit":     best.get("fwd_edit", best["fwd_mismatch"]),
        "rev_edit":     best.get("rev_edit", best["rev_mismatch"]),
        "primer_name":  primer_pair.name,
        "output_fasta": out_fasta,
        "n_candidates": n_candidates,
    }


def run_primer_extraction(
    assembly_fasta: str,
    primer_pair: PrimerPair,
    output_dir: str,
    strain_id: str,
    locus_name: str,
    max_mismatches: int = 2,
    blastn_bin: str = "blastn",
    manager=None,
    candidates: Optional[list[dict]] = None,
    chosen_index: int = 0,
) -> tuple[Optional[dict], str]:
    """
    Top-level primer-based locus extraction.

    If ``candidates`` is given (e.g. from a user review/disambiguation step), the
    one at ``chosen_index`` is extracted; otherwise the assembly is searched and
    the lowest-edit-distance candidate is used. When ``manager`` is supplied the
    blastn call is logged through it. Always writes a per-locus extraction log.

    Returns (result_dict, status). status == "ok" on success.
    """
    action = f"primer_{locus_name}_{strain_id}"
    if candidates is None:
        candidates = find_primer_amplicons(
            assembly_fasta, primer_pair,
            max_mismatches=max_mismatches, blastn_bin=blastn_bin,
            manager=manager, action=action,
        )

    if not candidates:
        status = (
            f"No amplicon found with {primer_pair.name} "
            f"(max_mismatches={max_mismatches}, "
            f"size {primer_pair.min_amplicon}-{primer_pair.max_amplicon} bp)"
        )
        write_primer_extraction_log(
            output_dir, strain_id, locus_name, primer_pair,
            None, max_mismatches, 0, status,
        )
        return None, status

    idx = chosen_index if 0 <= chosen_index < len(candidates) else 0
    best = candidates[idx]
    rec = extract_primer_amplicon(
        assembly_fasta, best, strain_id, locus_name, output_dir, primer_pair,
    )
    if rec is None:
        status = f"Contig '{best['contig']}' not found in assembly FASTA"
        write_primer_extraction_log(
            output_dir, strain_id, locus_name, primer_pair,
            None, max_mismatches, len(candidates), status,
        )
        return None, status

    write_primer_extraction_log(
        output_dir, strain_id, locus_name, primer_pair,
        best, max_mismatches, len(candidates), "ok",
    )
    out_fasta = str(Path(output_dir) / f"{locus_name}_amplicon.fasta")
    return _result_dict(best, primer_pair, out_fasta, len(candidates)), "ok"
