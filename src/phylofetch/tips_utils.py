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


# ── Locus auto-classification ─────────────────────────────────────────────────

import re as _re


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
    """{accession (versioned and base): title} for a batch of accessions via one esummary."""
    _require_email()
    time.sleep(0.34)
    handle = Entrez.esummary(db=db, id=",".join(accessions))
    recs = Entrez.read(handle)
    handle.close()
    titles: dict[str, str] = {}
    for s in recs:
        av = str(s.get("AccessionVersion", ""))
        t = str(s.get("Title", ""))
        if av:
            titles[av] = t
            titles[av.split(".")[0]] = t
    return titles


def classify_accessions(accessions: list[str], db: str = "nucleotide") -> dict:
    """
    Look up titles for ``accessions`` and classify each to a locus (no fetching of sequences).
    Returns {"by_locus": {locus: [acc, ...]}, "unassigned": [acc, ...], "titles": {acc: title}}.
    """
    accs = [a.strip() for a in accessions if a.strip()]
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
    accs = [a.strip() for a in accessions if a.strip()]
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
