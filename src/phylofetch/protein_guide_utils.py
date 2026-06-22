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
from pathlib import Path
from typing import Optional

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

USER_GUIDE_PATH = Path.home() / ".phylofetch" / "protein_guides.json"


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


def write_guide_fasta(locus: str, output_path: str, include_user: bool = True) -> Optional[str]:
    """
    Write the protein guide(s) for ``locus`` to ``output_path`` as FASTA (the Exonerate query).
    Multiple guides (e.g. Ascomycota + Basidiomycota) are written together so Exonerate's
    best-model selection can pick the closest. Returns the path, or None if no guide exists.
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
    if not seqrecs:
        return None
    Path(output_path).parent.mkdir(parents=True, exist_ok=True)
    SeqIO.write(seqrecs, output_path, "fasta")
    return output_path
