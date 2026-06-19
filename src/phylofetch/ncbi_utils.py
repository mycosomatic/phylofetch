"""
ncbi_utils.py
-------------
Search and fetch sequences from NCBI via Biopython Entrez.

Genus-agnostic: organism defaults are configurable; no hardcoded taxon names.
Reference library stored at ~/.phylofetch/references/.
"""

import json
import os
import time
from dataclasses import asdict, dataclass, fields
from datetime import datetime, timezone
from pathlib import Path
from typing import Optional

from Bio import Entrez, SeqIO
from Bio.SeqRecord import SeqRecord


_entrez_email: str = ""


def set_email(email: str) -> None:
    global _entrez_email
    _entrez_email = email
    Entrez.email = email


def _require_email() -> None:
    if not _entrez_email:
        raise ValueError(
            "NCBI Entrez email not set. "
            "Enter your email in Project Setup before fetching."
        )


# ── Reference library ────────────────────────────────────────────────────────

REF_DIR = Path.home() / ".phylofetch" / "references"


def locus_ref_fasta(locus_name: str) -> str:
    d = REF_DIR / locus_name
    d.mkdir(parents=True, exist_ok=True)
    return str(d / f"{locus_name}_refs.fasta")


def list_loci() -> list[str]:
    if not REF_DIR.exists():
        return []
    return sorted(
        p.name for p in REF_DIR.iterdir()
        if p.is_dir() and (p / f"{p.name}_refs.fasta").exists()
    )


def load_ref_records(locus_name: str) -> list:
    fasta = locus_ref_fasta(locus_name)
    if not os.path.exists(fasta):
        return []
    return list(SeqIO.parse(fasta, "fasta"))


def accessions_in_library(locus_name: str) -> set:
    accs: set = set()
    for r in load_ref_records(locus_name):
        base = r.id.split(".")[0]
        accs.add(r.id)
        accs.add(base)
    return accs


def count_refs(locus_name: str) -> int:
    return len(load_ref_records(locus_name))


# ── Reference metadata sidecar (type material, voucher, provenance) ──────────
#
# Alongside each <locus>/<locus>_refs.fasta we keep a structured sidecar
# <locus>/<locus>_refs.json mapping accession → RefRecord. The FASTA stays the
# extraction-query input, untouched; the JSON carries the organism / strain /
# voucher / type-material metadata the FASTA-only fetch used to discard. It is
# the join source for informative tip labels and for a publication accession
# table (organism + strain are the cross-locus join key). See DECISIONS.md D-007.

def _now_iso() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat()


@dataclass
class RefRecord:
    """Structured metadata for one reference sequence in a locus library."""
    accession: str
    organism: str = ""
    strain: str = ""
    culture_collection: str = ""
    specimen_voucher: str = ""
    type_material: str = ""    # raw INSDC /type_material qualifier value ("" if absent)
    type_kind: str = ""        # normalised: holotype / ex-holotype / epitype / type strain / …
    is_type: bool = False
    length: int = 0
    db: str = "nucleotide"
    fetched_at: str = ""
    query: str = ""


_REFREC_FIELDS = {f.name for f in fields(RefRecord)}


# INSDC /type_material keywords, longest/most-specific first so e.g. 'isolectotype'
# is matched before 'lectotype'. Ex-type CULTURES (annotated 'culture from <kind>')
# are treated as equal in grade to the nomenclatural type (D-007); the kind is
# preserved as 'ex-<kind>' so the distinction is still visible to the user.
_TYPE_KEYWORDS = [
    "paralectotype", "isolectotype", "isoneotype", "isosyntype", "isoepitype",
    "holotype", "isotype", "syntype", "lectotype", "neotype", "epitype",
    "paratype", "topotype", "ex-type", "type strain",
]


def normalize_type_kind(type_material: str) -> tuple[bool, str]:
    """
    Parse an INSDC /type_material qualifier value into (is_type, type_kind).

    The qualifier is usually '<kind> of <organism>' (e.g. 'holotype of Aspergillus
    flavus') or a culture form (e.g. 'culture from holotype of …'). Any non-empty
    qualifier means the record IS type-grade; unrecognised wording falls back to the
    generic kind 'type'. Returns (False, "") only when the qualifier is absent.
    """
    if not type_material or not type_material.strip():
        return False, ""
    text = type_material.strip().lower()
    is_culture = ("culture from" in text) or text.startswith("ex-")
    for kw in _TYPE_KEYWORDS:
        if kw in text:
            if kw in ("ex-type", "type strain"):
                return True, kw
            if is_culture:
                return True, f"ex-{kw}"      # e.g. ex-holotype
            return True, kw
    return True, "type"


def parse_source_metadata(record: SeqRecord) -> dict:
    """
    Pull organism + voucher/type metadata from a GenBank SeqRecord's source feature.
    Safe on records with no source feature (returns empty strings / is_type False).
    """
    organism = ""
    if getattr(record, "annotations", None):
        organism = str(record.annotations.get("organism", "") or "")

    quals: dict = {}
    for feat in getattr(record, "features", []):
        if getattr(feat, "type", "") == "source":
            quals = feat.qualifiers
            break

    def first(key: str) -> str:
        v = quals.get(key)
        return str(v[0]) if v else ""

    type_material = first("type_material")
    is_type, type_kind = normalize_type_kind(type_material)
    return {
        "organism": organism or first("organism"),
        "strain": first("strain"),
        "culture_collection": first("culture_collection"),
        "specimen_voucher": first("specimen_voucher"),
        "type_material": type_material,
        "type_kind": type_kind,
        "is_type": is_type,
    }


def _refrec_to_dict(r: RefRecord) -> dict:
    d = asdict(r)
    d.pop("accession", None)   # stored as the mapping key
    return d


def locus_ref_meta_path(locus_name: str, ref_dir: Path = REF_DIR) -> Path:
    d = Path(ref_dir) / locus_name
    d.mkdir(parents=True, exist_ok=True)
    return d / f"{locus_name}_refs.json"


def load_ref_meta(locus_name: str, ref_dir: Path = REF_DIR) -> dict[str, RefRecord]:
    """Load the per-locus reference metadata sidecar (empty dict if none/invalid)."""
    path = Path(ref_dir) / locus_name / f"{locus_name}_refs.json"
    if not path.exists():
        return {}
    try:
        raw = json.loads(path.read_text(encoding="utf-8"))
    except (json.JSONDecodeError, OSError):
        return {}
    records = raw.get("records", {}) if isinstance(raw, dict) else {}
    out: dict[str, RefRecord] = {}
    for acc, d in records.items():
        if not isinstance(d, dict):
            continue
        kwargs = {k: v for k, v in d.items() if k in _REFREC_FIELDS and k != "accession"}
        try:
            out[acc] = RefRecord(accession=acc, **kwargs)
        except (TypeError, ValueError):
            continue
    return out


def save_ref_meta(locus_name: str, records: dict[str, RefRecord],
                  ref_dir: Path = REF_DIR) -> str:
    """Persist the per-locus reference metadata sidecar. Returns the path."""
    path = locus_ref_meta_path(locus_name, ref_dir=ref_dir)
    payload = {
        "schema_version": 1,
        "records": {acc: _refrec_to_dict(r) for acc, r in records.items()},
    }
    path.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    return str(path)


def _log_fetch(locus_name: str, db: str, query: str,
               requested: list[str], added: int, skipped: int,
               ref_dir: Path = REF_DIR) -> None:
    """
    Append a self-contained provenance line for one fetch into the locus library.

    NCBI Entrez runs in-process (not a subprocess), so it cannot route through
    RunManager; this mirrors the same fields (query, params, counts, timestamp)
    in a fetch log kept with the reference library.
    """
    d = Path(ref_dir) / locus_name
    d.mkdir(parents=True, exist_ok=True)
    line = {
        "fetched_at": _now_iso(), "db": db, "query": query,
        "requested": list(requested), "n_requested": len(requested),
        "added": added, "skipped": skipped,
    }
    with open(d / f"{locus_name}_fetch_log.jsonl", "a", encoding="utf-8") as f:
        f.write(json.dumps(line, sort_keys=True) + "\n")


def fetch_record_with_meta(
    accession: str, db: str = "nucleotide", query: str = "",
) -> tuple[Optional[SeqRecord], Optional[RefRecord]]:
    """
    Fetch a record as GenBank and return (sequence SeqRecord, RefRecord metadata).

    Uses rettype='gb' (not 'fasta') so the source-feature qualifiers — /type_material,
    /strain, /culture_collection, /specimen_voucher — are available; these are exactly
    what the FASTA-only fetch threw away.
    """
    _require_email()
    time.sleep(0.34)
    handle = Entrez.efetch(db=db, id=accession, rettype="gb", retmode="text")
    records = list(SeqIO.parse(handle, "genbank"))
    handle.close()
    if not records:
        return None, None
    rec = records[0]
    meta = parse_source_metadata(rec)
    refrec = RefRecord(
        accession=rec.id or accession,
        organism=meta["organism"], strain=meta["strain"],
        culture_collection=meta["culture_collection"],
        specimen_voucher=meta["specimen_voucher"],
        type_material=meta["type_material"], type_kind=meta["type_kind"],
        is_type=meta["is_type"], length=len(rec.seq), db=db,
        fetched_at=_now_iso(), query=query,
    )
    return rec, refrec


# ── NCBI search ──────────────────────────────────────────────────────────────

_TYPE_FILTER = "sequence_from_type[filter]"


def _esearch_ids(db: str, term: str, max_results: int) -> list[str]:
    handle = Entrez.esearch(db=db, term=term, retmax=max_results)
    record = Entrez.read(handle)
    handle.close()
    return list(record.get("IdList", []))


def _summary_to_dict(s) -> dict:
    return {
        "uid":       str(s.get("Id", "")),
        "accession": str(s.get("AccessionVersion", "")),
        "title":     str(s.get("Title", "")),
        "organism":  str(s.get("Organism", "")),
        "length":    int(s.get("Length", 0)),
    }


def _mark_type_set(results: list[dict], type_uids: set[str]) -> list[dict]:
    """Flag each result is_type by UID membership (pure; unit-tested)."""
    for r in results:
        r["is_type"] = r.get("uid", "") in type_uids
    return results


def _search_ncbi(db: str, base_query: str, max_results: int,
                 type_mode: str) -> list[dict]:
    """
    Shared search for protein/nucleotide with a type-material mode:
      'all'       — all hits, type-derived ones flagged
      'prefer'    — all hits, type flagged AND sorted to the top
      'type_only' — restrict to sequence_from_type hits

    Type status is resolved cheaply with one extra esearch (the base query AND the
    'sequence from type' filter); accessions in that set are flagged is_type without
    fetching every record. The precise type KIND (holotype vs ex-type) is filled in
    later at fetch time from the GenBank record.
    """
    _require_email()
    if type_mode == "type_only":
        ids = _esearch_ids(db, f"{base_query} AND {_TYPE_FILTER}", max_results)
        type_uids = set(ids)
    else:
        ids = _esearch_ids(db, base_query, max_results)
        if ids:
            time.sleep(0.34)
            type_uids = set(_esearch_ids(db, f"{base_query} AND {_TYPE_FILTER}", max_results))
        else:
            type_uids = set()
    if not ids:
        return []

    time.sleep(0.34)
    handle = Entrez.esummary(db=db, id=",".join(ids))
    summaries = Entrez.read(handle)
    handle.close()

    results = _mark_type_set([_summary_to_dict(s) for s in summaries], type_uids)
    if type_mode == "prefer":
        results.sort(key=lambda r: not r["is_type"])   # type-derived first, stable
    return results


def search_ncbi_protein(gene_name: str, organism: str,
                        max_results: int = 20, type_mode: str = "all") -> list[dict]:
    return _search_ncbi(
        "protein", f"{gene_name}[Gene Name] AND {organism}[Organism]",
        max_results, type_mode,
    )


def search_ncbi_nucleotide(gene_name: str, organism: str,
                           max_results: int = 20, type_mode: str = "all") -> list[dict]:
    return _search_ncbi(
        "nucleotide", f"{gene_name}[Title] AND {organism}[Organism]",
        max_results, type_mode,
    )


# ── NCBI fetch ───────────────────────────────────────────────────────────────

def fetch_protein_by_accession(accession: str) -> Optional[SeqRecord]:
    _require_email()
    time.sleep(0.34)
    handle = Entrez.efetch(db="protein", id=accession, rettype="fasta", retmode="text")
    records = list(SeqIO.parse(handle, "fasta"))
    handle.close()
    return records[0] if records else None


def fetch_nucleotide_by_accession(accession: str) -> Optional[SeqRecord]:
    _require_email()
    time.sleep(0.34)
    handle = Entrez.efetch(db="nucleotide", id=accession, rettype="fasta", retmode="text")
    records = list(SeqIO.parse(handle, "fasta"))
    handle.close()
    return records[0] if records else None


def fetch_and_store(accessions: list[str], locus_name: str,
                    db: str = "protein",
                    custom_label: Optional[str] = None,
                    query: str = "") -> tuple[int, int, list]:
    """
    Fetch accessions, append sequences to the locus FASTA, and record per-accession
    metadata (organism / strain / voucher / type material) in the sidecar JSON.
    ``query`` is stored for provenance. Returns (added, skipped, errors).
    """
    _require_email()
    existing = accessions_in_library(locus_name)
    fasta_path = locus_ref_fasta(locus_name)
    meta = load_ref_meta(locus_name)
    added, skipped, errors = 0, 0, []
    with open(fasta_path, "a") as out_f:
        for acc in accessions:
            base = acc.split(".")[0]
            if acc in existing or base in existing:
                skipped += 1
                continue
            try:
                record, refrec = fetch_record_with_meta(acc, db=db, query=query)
                if record is None:
                    errors.append(f"{acc}: not found")
                    continue
                if custom_label:
                    record.id = custom_label
                    record.description = ""
                SeqIO.write(record, out_f, "fasta")
                if refrec is not None:
                    meta[custom_label or refrec.accession] = refrec
                added += 1
            except Exception as e:
                errors.append(f"{acc}: {e}")
    save_ref_meta(locus_name, meta)
    _log_fetch(locus_name, db, query, accessions, added, skipped)
    return added, skipped, errors


def delete_from_library(locus_name: str, accession: str) -> bool:
    fasta_path = locus_ref_fasta(locus_name)
    if not os.path.exists(fasta_path):
        return False
    records = list(SeqIO.parse(fasta_path, "fasta"))
    before = len(records)
    records = [r for r in records if r.id != accession]
    if len(records) == before:
        return False
    with open(fasta_path, "w") as f:
        SeqIO.write(records, f, "fasta")
    return True


# ── Locus catalogue ───────────────────────────────────────────────────────────
# Standard fungal phylogenetic markers with practical NCBI search terms.
# Protein-coding loci default to nucleotide (CDS) so blastn HSP-as-exon
# extraction works naturally. Users can switch to protein for tblastn.
#
# This catalogue is genus-agnostic. Set organism in the UI per locus.

LOCUS_CATALOGUE: dict[str, dict] = {
    # rDNA — ITSx extracts from assemblies; fetch here for outgroup references
    "ITS":   {"db": "nucleotide", "gene": "internal transcribed spacer",
               "note": "ITSx extracts from assemblies; fetch here for outgroup taxa"},
    "LSU":   {"db": "nucleotide", "gene": "28S large subunit ribosomal RNA",
               "note": "ITSx extracts from assemblies; fetch here for outgroup taxa"},
    "SSU":   {"db": "nucleotide", "gene": "18S small subunit ribosomal RNA",
               "note": "ITSx extracts from assemblies; fetch here for outgroup taxa"},
    # Protein-coding — use CDS for blastn extraction
    "TEF1":  {"db": "nucleotide", "gene": "tef1 complete cds",
               "note": "Translation elongation factor 1-alpha. Use CDS-only refs."},
    "RPB1":  {"db": "nucleotide", "gene": "rpb1 complete cds",
               "note": "RNA pol II largest subunit. Use CDS-only refs."},
    "RPB2":  {"db": "nucleotide", "gene": "rpb2 complete cds",
               "note": "RNA pol II second largest subunit. Use CDS-only refs."},
    "TUB2":  {"db": "nucleotide", "gene": "tub2 complete cds",
               "note": "Beta-tubulin. Has small exons — try 'blastn' task if exons are missed."},
    "GAPDH": {"db": "nucleotide", "gene": "gapdh complete cds",
               "note": "Glyceraldehyde-3-phosphate dehydrogenase. Use CDS-only refs."},
    "CAL":   {"db": "nucleotide", "gene": "calmodulin complete cds",
               "note": "Calmodulin. Use CDS-only refs."},
    "ACT":   {"db": "nucleotide", "gene": "actin complete cds",
               "note": "Actin. Use CDS-only refs."},
    "HIS3":  {"db": "nucleotide", "gene": "histone h3 complete cds",
               "note": "Histone H3. Useful for fine-scale resolution in some groups."},
}
