"""
tips_utils.py
-------------
Comparison-taxa ("tree tip") import with automatic locus classification (RM-008 component 3,
D-020).

Tips are the *other taxa* you compare your isolates against in a tree — usually partial
nucleotide barcode amplicons pulled from GenBank. They are stored **separately from extraction
guides** (which are protein) in ``<project>/tips/<locus>/`` so the two never mix in one file
(cf. D-018). The headline feature: paste a *flat, mixed* list of accessions and each is
**auto-classified to its locus** by matching its GenBank title against the locus synonym
catalogue (D-011); only the unmatched ones need manual assignment.

Storage reuses the per-locus library helpers in ``ncbi_utils`` with ``ref_dir`` pointed at the
tips directory, so tips get the same FASTA + metadata sidecar + fetch-log provenance.
"""

import time
from collections import defaultdict
from pathlib import Path
from typing import Optional

from Bio import Entrez

from phylofetch.ncbi_utils import (
    LOCUS_CATALOGUE,
    _require_email,
    fetch_and_store,
)


def project_tips_dir(project_dir) -> Path:
    """Per-project tip store: ``<project>/tips`` (created)."""
    d = Path(project_dir) / "tips"
    d.mkdir(parents=True, exist_ok=True)
    return d


# ── Accession normalization ───────────────────────────────────────────────────

import re as _re

# RefSeq accessions carry a mandatory underscore after the two-letter molecule prefix
# (e.g. NR_135944). Pasted without it (``NR135944``) the bare form still resolves as an NCBI
# *URL*, but Entrez/BLAST do **not** recognise it, so the record silently fails to populate.
# These are the two-letter RefSeq prefixes we repair; GenBank accessions never use them, so
# there is no collision with a real bare GenBank id (D-026).
REFSEQ_PREFIXES = frozenset({
    "AC", "NC", "NG", "NM", "NR", "NT", "NW", "NZ", "XM", "XR",   # nucleotide / RNA
    "AP", "NP", "WP", "YP", "XP",                                  # protein
})

_REFSEQ_BARE = _re.compile(r"^([A-Z]{2})(\d{6,})(\.\d+)?$")


def normalize_accession(acc: str) -> str:
    """
    Repair a RefSeq accession that is missing its underscore (``NR135944`` → ``NR_135944``),
    so Entrez/BLAST can resolve it (D-026). Any other accession — GenBank, or an already
    well-formed RefSeq id — is returned unchanged (only whitespace-stripped). The match is
    anchored on a known RefSeq two-letter prefix + ≥6 digits, so ordinary GenBank ids are
    never rewritten.
    """
    a = (acc or "").strip()
    if not a:
        return a
    m = _REFSEQ_BARE.match(a.upper())
    if m and m.group(1) in REFSEQ_PREFIXES:
        return f"{m.group(1)}_{m.group(2)}{m.group(3) or ''}"
    return a


# ── Locus auto-classification ─────────────────────────────────────────────────


def _locus_terms(catalogue) -> dict[str, list[str]]:
    out: dict[str, list[str]] = {}
    for loc, cat in catalogue.items():
        terms = [cat.get("gene", "")] + list(cat.get("synonyms", []))
        out[loc] = [t.strip().lower() for t in terms if t and t.strip()]
    return out


def classify_locus(text: str, catalogue=LOCUS_CATALOGUE) -> Optional[str]:
    """
    Best-matching locus for a free-text record title (and/or gene qualifiers), by synonym match
    against the catalogue. Returns the locus, or None when nothing matches or two loci tie
    (genuinely ambiguous — e.g. a multi-gene record — so the caller routes it to manual review).
    Each locus is scored by its single strongest (longest) matching synonym, so a specific name
    like "RNA polymerase II second largest subunit" beats a bare symbol.
    """
    blob = (text or "").lower()
    scores: dict[str, int] = {}
    for loc, terms in _locus_terms(catalogue).items():
        best = 0
        for t in terms:
            if _re.search(r"(?<![a-z0-9])" + _re.escape(t) + r"(?![a-z0-9])", blob):
                best = max(best, len(t))
        if best:
            scores[loc] = best
    if not scores:
        return None
    top = max(scores.values())
    winners = [loc for loc, v in scores.items() if v == top]
    return winners[0] if len(winners) == 1 else None


# ── Accession summaries + import ──────────────────────────────────────────────

def _esummary_titles(accessions: list[str], db: str) -> dict[str, str]:
    """
    ``{accession (versioned and base): title}`` for a batch of accessions via esummary. Tries
    one batched call; if that raises (NCBI rejects the *whole* batch when any id is malformed),
    it falls back to per-accession calls so the good ones still resolve and only the bad ones go
    missing (→ caller flags them as "not found", D-026).
    """
    _require_email()
    titles: dict[str, str] = {}

    def _batch(ids: list[str]) -> None:
        time.sleep(0.34)
        handle = Entrez.esummary(db=db, id=",".join(ids))
        recs = Entrez.read(handle)
        handle.close()
        for s in recs:
            av = str(s.get("AccessionVersion", ""))
            t = str(s.get("Title", ""))
            if av:
                titles[av] = t
                titles[av.split(".")[0]] = t

    try:
        _batch(accessions)
    except Exception:                                   # noqa: BLE001 — isolate the offender(s)
        for a in accessions:
            try:
                _batch([a])
            except Exception:                           # noqa: BLE001 — unresolvable id; skip
                pass
    return titles


def lookup_accessions(accessions: list[str], db: str = "nucleotide") -> list[dict]:
    """
    Normalize each pasted accession, pull its GenBank title from NCBI (no sequence fetch yet),
    and suggest a locus — one row per input, preserving order (D-026). Each row is::

        {"input": <as pasted>, "accession": <normalized>, "title": str,
         "found": bool, "locus_guess": <locus or None>}

    ``found=False`` means NCBI returned nothing for that id (typo, wrong database, or a bare
    RefSeq id we could not repair) — the page surfaces these so none fail silently.
    """
    seen: list[tuple[str, str]] = []
    for a in accessions:
        if a and a.strip():
            seen.append((a.strip(), normalize_accession(a)))
    titles = _esummary_titles([n for _, n in seen], db) if seen else {}
    rows: list[dict] = []
    for original, acc in seen:
        title = titles.get(acc) or titles.get(acc.split(".")[0]) or ""
        rows.append({
            "input": original,
            "accession": acc,
            "title": title,
            "found": bool(title),
            "locus_guess": classify_locus(title) if title else None,
        })
    return rows


def import_tips_with_assignments(assignments: dict[str, str], tips_dir,
                                 db: str = "nucleotide",
                                 query: str = "tips import") -> dict:
    """
    Fetch + store each accession under the **locus the user chose for it** (D-026). ``assignments``
    maps accession → locus; an entry with a falsy locus is skipped (the user opted not to import
    it). Accessions are normalized first. Returns
    ``{"assigned": {locus: [acc]}, "errors": [str]}``.
    """
    by_locus: dict[str, list[str]] = defaultdict(list)
    for acc, loc in assignments.items():
        norm = normalize_accession(acc)
        if loc and norm:
            by_locus[loc].append(norm)

    assigned: dict[str, list[str]] = {}
    errors: list[str] = []
    for loc, group in by_locus.items():
        try:
            _added, _skipped, errs = fetch_and_store(group, loc, db=db, query=query,
                                                     ref_dir=tips_dir)
            assigned[loc] = group
            errors.extend(errs)
        except Exception as e:                          # noqa: BLE001 — surface to UI
            errors.append(f"{loc}: {e}")
    return {"assigned": assigned, "errors": errors}


def classify_accessions(accessions: list[str], db: str = "nucleotide") -> dict:
    """
    Look up titles for ``accessions`` and classify each to a locus (no fetching of sequences).
    Returns {"by_locus": {locus: [acc, ...]}, "unassigned": [acc, ...], "titles": {acc: title}}.
    """
    accs = [normalize_accession(a) for a in accessions if a.strip()]
    titles = _esummary_titles(accs, db) if accs else {}
    by_locus: dict[str, list[str]] = defaultdict(list)
    unassigned: list[str] = []
    for a in accs:
        title = titles.get(a) or titles.get(a.split(".")[0]) or ""
        loc = classify_locus(title)
        (by_locus[loc].append(a) if loc else unassigned.append(a))
    return {"by_locus": dict(by_locus), "unassigned": unassigned, "titles": titles}


def import_tip_accessions(accessions: list[str], tips_dir, db: str = "nucleotide",
                          force_locus: Optional[str] = None,
                          query: str = "tips import") -> dict:
    """
    Fetch ``accessions`` and store them as tips in ``tips_dir``, auto-classified to their locus
    (or all to ``force_locus`` when given — e.g. reassigning previously unassigned ones).
    Returns {"assigned": {locus: [acc]}, "unassigned": [acc], "errors": [str]}.
    """
    accs = [normalize_accession(a) for a in accessions if a.strip()]
    if force_locus:
        grouped, unassigned = {force_locus: accs}, []
    else:
        c = classify_accessions(accs, db=db)
        grouped, unassigned = c["by_locus"], c["unassigned"]

    assigned: dict[str, list[str]] = {}
    errors: list[str] = []
    for loc, group in grouped.items():
        try:
            _added, _skipped, errs = fetch_and_store(group, loc, db=db, query=query,
                                                     ref_dir=tips_dir)
            assigned[loc] = group
            errors.extend(errs)
        except Exception as e:                      # noqa: BLE001 — surface to UI
            errors.append(f"{loc}: {e}")
    return {"assigned": assigned, "unassigned": unassigned, "errors": errors}
