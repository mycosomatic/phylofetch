"""
protein_guide_utils.py
----------------------
Bundled protein reference guides for protein2genome extraction (D-020 / RM-008).

A *universal core* of full-length RefSeq orthologs (one Ascomycota + one Basidiomycota per
conserved coding marker) ships as packaged data (``phylofetch/data/protein_guides.json``) and is
used as the Exonerate query for the standard gene set — so extracting the standard markers needs
no per-project fetching of guides. These markers are conserved kingdom-wide, so a guide from any
fungus locates the ortholog in any other via protein-level alignment (D-020).

Users can add clade-specific **lineage packs** at ``~/.phylofetch/protein_guides.json`` (same
schema, merged on top) — e.g. Alt a1 for Dothideomycetes.

Guides are EXTRACTION references (protein); they are conceptually and physically separate from
tree tips (nucleotide amplicons) and must never be mixed in one locus file (cf. D-018).
"""

import json
import statistics
from pathlib import Path
from typing import Optional

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

USER_GUIDE_PATH = Path.home() / ".phylofetch" / "protein_guides.json"

# Fraction a fetched reference's length may deviate from the curated bundled-guide length for
# its locus before it is treated as a length outlier (mis-annotated / fused / truncated). The
# bundled universal guides are hand-curated full-length orthologs, so they are the trusted
# expectation. 0.30 keeps honest paralog-length variation (e.g. RPB2 ±10 %) while rejecting the
# gross cases seen in the wild — an 865-aa "partial beta-tubulin" (real ≈ 447) or an 814-aa
# EF1-alpha (real ≈ 460) that out-score the correct guides on raw identity (D-025).
LENGTH_TOLERANCE = 0.30


def _builtin_data_path() -> Path:
    return Path(__file__).parent / "data" / "protein_guides.json"


def _load_builtin_raw() -> dict:
    """Load the packaged protein_guides.json, tolerating install or src-tree layout."""
    try:                                   # installed package
        from importlib.resources import files
        text = (files("phylofetch") / "data" / "protein_guides.json").read_text(encoding="utf-8")
    except (ModuleNotFoundError, FileNotFoundError, AttributeError, TypeError):
        text = _builtin_data_path().read_text(encoding="utf-8")
    return json.loads(text)


def load_protein_guides(include_user: bool = True,
                        user_path: Path = USER_GUIDE_PATH) -> dict[str, list[dict]]:
    """
    Return ``{locus: [guide_record, ...]}`` merging the bundled universal core with any user
    lineage packs (user guides appended per locus; new loci added). Each record carries
    ``accession / organism / clade / protein_name / length / seq``.
    """
    raw = _load_builtin_raw()
    guides: dict[str, list[dict]] = {loc: list(recs)
                                     for loc, recs in raw.get("guides", {}).items()}
    if include_user and Path(user_path).exists():
        try:
            uraw = json.loads(Path(user_path).read_text(encoding="utf-8"))
            for loc, recs in (uraw.get("guides", {}) or {}).items():
                guides.setdefault(loc, [])
                guides[loc].extend(recs)
        except (json.JSONDecodeError, OSError):
            pass                            # a broken user pack must not break the built-ins
    return guides


def guide_loci(include_user: bool = True) -> list[str]:
    """Loci that have at least one protein guide."""
    return sorted(load_protein_guides(include_user=include_user))


def get_guides(locus: str, include_user: bool = True) -> list[dict]:
    """Guide records for one locus (empty list if none)."""
    return load_protein_guides(include_user=include_user).get(locus, [])


def expected_length(locus: str, include_user: bool = True) -> Optional[int]:
    """
    Trusted full-length protein expectation for ``locus`` — the median length of its curated
    bundled guides. Used to flag mis-annotated fetched references (D-025). Returns ``None`` when
    the locus has no bundled guide (no expectation can be formed, so nothing is filtered).
    """
    lens = [g.get("length") or len((g.get("seq") or "").strip())
            for g in get_guides(locus, include_user=include_user)]
    lens = [int(n) for n in lens if n]
    return int(round(statistics.median(lens))) if lens else None


def length_flag(length: int, expected: Optional[int],
                tol: float = LENGTH_TOLERANCE) -> Optional[str]:
    """``"short"`` / ``"long"`` if ``length`` deviates from ``expected`` by more than ``tol``;
    ``None`` if within tolerance or there is no expectation to judge against."""
    if not expected:
        return None
    if length < expected * (1 - tol):
        return "short"
    if length > expected * (1 + tol):
        return "long"
    return None


def filter_records_by_length(records, locus: str, tol: float = LENGTH_TOLERANCE,
                             include_user: bool = True):
    """
    Split protein ``SeqRecord``s into ``(kept, dropped)`` by comparing each to the locus's
    bundled-guide :func:`expected_length` (D-025). ``dropped`` is a list of
    ``(record, length, expected, flag)`` so the caller can report exactly what was excluded and
    why. When the locus has no bundled guide, the expectation is ``None`` and **nothing is
    dropped** (we never silently discard a reference we cannot judge).
    """
    exp = expected_length(locus, include_user=include_user)
    kept, dropped = [], []
    for r in records:
        n = len(r.seq)
        flag = length_flag(n, exp, tol)
        (dropped.append((r, n, exp, flag)) if flag else kept.append(r))
    return kept, dropped


def write_guide_fasta(locus: str, output_path: str, include_user: bool = True,
                      extra_records=None) -> Optional[str]:
    """
    Write the protein guide(s) for ``locus`` to ``output_path`` as FASTA (the Exonerate query).
    Multiple guides (e.g. Ascomycota + Basidiomycota) are written together so Exonerate's
    best-model selection can pick the closest. Returns the path, or None if no guide exists.

    ``extra_records`` (a list of ``SeqRecord``) are appended after the bundled guides — used to
    layer the project's fetched **taxon-closer** protein references on top of the universal core
    (D-023), so Exonerate has both the kingdom-wide floor and a near-relative ortholog and keeps
    whichever scores best. They must be protein (the model is ``protein2genome``); the caller is
    responsible for that. With only ``extra_records`` and no bundled guide for the locus, those
    are written alone.
    """
    seqrecs = []
    for g in get_guides(locus, include_user=include_user):
        seq = (g.get("seq") or "").strip()
        if not seq:
            continue
        seqrecs.append(SeqRecord(
            Seq(seq),
            id=f"{locus}_{g.get('accession', 'guide')}",
            description=(f"[guide] [locus={locus}] [organism={g.get('organism', '')}] "
                        f"[clade={g.get('clade', '')}]"),
        ))
    if extra_records:
        seqrecs.extend(extra_records)
    if not seqrecs:
        return None
    Path(output_path).parent.mkdir(parents=True, exist_ok=True)
    SeqIO.write(seqrecs, output_path, "fasta")
    return output_path
