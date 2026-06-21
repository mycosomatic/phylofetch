"""
taxon_id_utils.py
-----------------
Provisional taxonomic identification of an assembly from its ITS region (RM-007 step 3).

Pipeline (D-012 / user choice 2026-06-20 = remote NCBI):
  ITSx extracts the ITS region from the assembly  ->  blastn -remote against NCBI `nt`
  ->  ranked organism hits  ->  the user picks one  ->  the assembly's taxon is set with
  taxon_source="its_blast" (see project_manager.set_assembly_taxon).

This is a *suggestion* step, not an authoritative identification: it informs which
references to fetch downstream. Remote BLAST needs internet and the blastn binary, and is
slow (often minutes); the call is logged through RunManager for provenance. Command building
and output parsing are split out as pure functions so they can be unit-tested offline.
"""

import os
import subprocess
from pathlib import Path
from typing import Optional

from Bio import SeqIO

from phylofetch.itsx_utils import run_itsx


# Tabular BLAST fields. NOTE: `sscinames` (the scientific name from a taxid) is resolved
# LOCALLY from a BLAST taxdb even under `-remote` — without taxdb installed it comes back
# "N/A" (verified against blastn 2.x). We therefore also request `stitle` (the subject
# defline, which carries the organism for nt records) and derive the organism from it,
# falling back to `sscinames` only when a real taxdb name is present. `stitle` is kept LAST
# because it contains spaces (never tabs), so it can't corrupt earlier columns.
REMOTE_BLAST_FIELDS = [
    "qseqid", "sacc", "pident", "length", "evalue", "bitscore", "staxids", "sscinames",
    "stitle",
]
_REMOTE_OUTFMT = "6 " + " ".join(REMOTE_BLAST_FIELDS)

# Leading status tags some nt deflines carry before the organism (e.g. "PREDICTED: Homo …").
_TITLE_STATUS_TAGS = ("predicted:", "tpa:", "tpa_exp:", "tpa_inf:", "unverified:",
                      "unverified_org:", "mag:", "recname:")


def _organism_from_title(stitle: str) -> str:
    """
    Best-effort organism (leading binomial) from an nt subject title, e.g.
    'Alternaria alternata strain X 18S …' -> 'Alternaria alternata';
    'Alternaria sp. ALT1 …' -> 'Alternaria sp.'; 'Uncultured fungus …' -> 'Uncultured fungus'.
    Heuristic by design — this is a provisional suggestion the user reviews and picks.
    """
    toks = (stitle or "").split()
    if toks and toks[0].lower() in _TITLE_STATUS_TAGS:   # drop a leading status tag
        toks = toks[1:]
    if not toks:
        return ""
    if len(toks) == 1:
        return toks[0].strip(",;")
    return f"{toks[0]} {toks[1]}".strip(",;")


def _organism_from_hit(sscinames: str, stitle: str) -> str:
    """
    Prefer a real taxdb scientific name; otherwise derive it from the subject title.
    `sscinames` can be a ';'-joined list (subject identical to several taxids), and without
    a local taxdb each element is "N/A" (e.g. "N/A;N/A") — treat a leading N/A as missing.
    """
    first = (sscinames or "").split(";")[0].strip()
    if first and first.upper() != "N/A":
        return first
    return _organism_from_title(stitle)


def build_remote_blast_cmd(query_fasta: str, out_tsv: str, *,
                           blast_bin: str = "blastn",
                           evalue: float = 1e-10,
                           max_target_seqs: int = 20,
                           entrez_query: str = "fungi[ORGN]",
                           task: str = "megablast") -> list[str]:
    """
    Build the `blastn -remote -db nt` argument list (pure; unit-tested).

    Note: `-num_threads` is intentionally omitted — BLAST+ rejects it together with
    `-remote`. `entrez_query` restricts the remote search server-side (e.g. to fungi);
    pass "" to search all of nt.
    """
    cmd = [
        blast_bin,
        "-query", query_fasta,
        "-db", "nt",
        "-remote",
        "-task", task,
        "-evalue", str(evalue),
        "-max_target_seqs", str(max_target_seqs),
        "-outfmt", _REMOTE_OUTFMT,
        "-out", out_tsv,
    ]
    if entrez_query and entrez_query.strip():
        cmd += ["-entrez_query", entrez_query.strip()]
    return cmd


def parse_blast_organism_hits(tsv_text: str) -> list[dict]:
    """Parse tabular (-outfmt REMOTE_BLAST_FIELDS) BLAST output into hit dicts (pure)."""
    hits: list[dict] = []
    for line in tsv_text.splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        parts = line.split("\t")
        if len(parts) < len(REMOTE_BLAST_FIELDS):
            continue
        try:
            sscinames, stitle = parts[7], parts[8]
            hits.append({
                "qseqid":    parts[0],
                "accession": parts[1],
                "pident":    float(parts[2]),
                "length":    int(parts[3]),
                "evalue":    parts[4],
                "bitscore":  float(parts[5]),
                "staxids":   parts[6],
                "sscinames": sscinames.strip(),
                "stitle":    stitle.strip(),
                "organism":  _organism_from_hit(sscinames, stitle),
            })
        except ValueError:
            continue
    return hits


def rank_unique_organisms(hits: list[dict]) -> list[dict]:
    """
    Sort hits by bitscore (desc) and keep the best hit per organism, so the user sees one
    row per candidate taxon rather than many strains of the same species.
    """
    ordered = sorted(hits, key=lambda h: h.get("bitscore", 0.0), reverse=True)
    seen: set[str] = set()
    out: list[dict] = []
    for h in ordered:
        org = (h.get("organism") or "").strip()
        if not org or org.lower() in seen:
            continue
        seen.add(org.lower())
        out.append(h)
    return out


def _pick_longest_record(fasta_path: str):
    """Return the longest SeqRecord in a FASTA (or None) — the single ITS query to BLAST."""
    recs = list(SeqIO.parse(fasta_path, "fasta"))
    return max(recs, key=lambda r: len(r.seq)) if recs else None


def run_remote_blast_its(query_fasta: str, output_dir: str, *,
                         run_manager=None,
                         blast_bin: str = "blastn",
                         evalue: float = 1e-10,
                         max_target_seqs: int = 20,
                         entrez_query: str = "fungi[ORGN]",
                         task: str = "megablast",
                         timeout: int = 600) -> tuple[int, str, list[dict]]:
    """
    Remote-BLAST an ITS query against NCBI `nt`. Returns (returncode, stderr, ranked_hits).

    Logged through RunManager when one is supplied (module="taxon_id"); otherwise a plain
    subprocess. The `-out` TSV is written into ``output_dir`` and parsed regardless of which
    path ran.
    """
    os.makedirs(output_dir, exist_ok=True)
    out_tsv = os.path.join(output_dir, "its_remote_blast.tsv")
    cmd = build_remote_blast_cmd(
        query_fasta, out_tsv, blast_bin=blast_bin, evalue=evalue,
        max_target_seqs=max_target_seqs, entrez_query=entrez_query, task=task,
    )

    if run_manager is not None:
        res = run_manager.run(
            cmd, module="taxon_id", action="its_remote_blast",
            inputs={"query": query_fasta},
            outputs={"tsv": out_tsv},
            params={"db": "nt", "remote": True, "task": task, "evalue": evalue,
                    "max_target_seqs": max_target_seqs, "entrez_query": entrez_query},
            tool_version_keys=[blast_bin], timeout=timeout,
        )
        rc = res.returncode
        stderr = Path(res.stderr_path).read_text() if Path(res.stderr_path).exists() else ""
    else:
        try:
            r = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)
            rc, stderr = r.returncode, r.stderr
        except Exception as exc:                       # noqa: BLE001 — surface to UI
            rc, stderr = 999, f"{type(exc).__name__}: {exc}"

    text = Path(out_tsv).read_text() if os.path.exists(out_tsv) else ""
    hits = rank_unique_organisms(parse_blast_organism_hits(text))
    return rc, stderr, hits


def identify_taxon_from_assembly(assembly_path: str, output_dir: str, strain_id: str, *,
                                 run_manager=None,
                                 itsx_bin: str = "ITSx",
                                 blast_bin: str = "blastn",
                                 kingdom: str = "fungi",
                                 entrez_query: str = "fungi[ORGN]",
                                 evalue: float = 1e-10,
                                 max_target_seqs: int = 20,
                                 threads: int = 4,
                                 timeout: int = 600) -> dict:
    """
    Extract ITS (ITSx) from an assembly and remote-BLAST it for a provisional taxon.

    Returns a result dict:
      {"ok": True,  "its_region": "ITS_full"|"ITS1"|"ITS2", "its_length": int,
       "its_query": str, "hits": [ranked organism hits]}
      {"ok": False, "stage": "itsx"|"blast", "error": str, "its_region": ...|None,
       "hits": [], "itsx_log"?: str}

    The longest ITS record is used as the single query (keeps the remote search fast and the
    suggestion unambiguous). Prefers the full ITS, falling back to ITS1 then ITS2.
    """
    os.makedirs(output_dir, exist_ok=True)
    itsx_dir = os.path.join(output_dir, "itsx")
    rc, log, found = run_itsx(
        assembly_path, itsx_dir, strain_id,
        threads=threads, itsx_bin=itsx_bin, kingdom=kingdom,
    )

    region = next((r for r in ("ITS_full", "ITS1", "ITS2") if found.get(r)), None)
    if region is None:
        return {"ok": False, "stage": "itsx", "its_region": None, "hits": [],
                "error": "ITSx found no ITS region in this assembly.", "itsx_log": log}

    qrec = _pick_longest_record(found[region])
    if qrec is None:
        return {"ok": False, "stage": "itsx", "its_region": region, "hits": [],
                "error": f"ITSx {region} output was empty.", "itsx_log": log}

    query_fasta = os.path.join(output_dir, "id_query.fasta")
    SeqIO.write([qrec], query_fasta, "fasta")

    brc, berr, hits = run_remote_blast_its(
        query_fasta, output_dir, run_manager=run_manager, blast_bin=blast_bin,
        evalue=evalue, max_target_seqs=max_target_seqs, entrez_query=entrez_query,
        timeout=timeout,
    )
    if brc != 0 and not hits:
        return {"ok": False, "stage": "blast", "its_region": region,
                "its_query": str(qrec.seq), "hits": [],
                "error": berr.strip() or "Remote BLAST failed (network or blastn error)."}

    return {"ok": True, "its_region": region, "its_length": len(qrec.seq),
            "its_query": str(qrec.seq), "hits": hits}
