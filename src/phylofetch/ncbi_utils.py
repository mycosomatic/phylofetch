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
import urllib.error
from dataclasses import asdict, dataclass, fields
from datetime import datetime, timezone
from http.client import HTTPException
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


# ── Entrez transport: throttle, retry, typed errors (D-034) ──────────────────
#
# Naked Entrez calls fail two ways that matter for research correctness: a transient network
# blip (DNS, NCBI 5xx, HTTP 429 rate-limit) either crashes a page or — if a caller swallows it
# — looks like a genuine "0 hits / not found", silently changing which references a study uses.
# Every Entrez call now goes through `_entrez_retry`, which throttles to NCBI's rate limit,
# retries transient failures with backoff, and on final failure raises a typed `NCBIError` the
# UI can distinguish from an honest empty result. An `NCBI_API_KEY` (env or `set_api_key`)
# raises the limit from 3 to ~9 req/s.

_entrez_api_key: str = os.environ.get("NCBI_API_KEY", "")
if _entrez_api_key:
    Entrez.api_key = _entrez_api_key


class NCBIError(RuntimeError):
    """An NCBI Entrez request failed after retries (network / HTTP / parse) — NOT 'no hits'."""


# Bio.Entrez.read raises RuntimeError on an NCBI-side error body; the rest are transport-level.
_TRANSIENT_ERRORS = (
    urllib.error.URLError, urllib.error.HTTPError, HTTPException, OSError, RuntimeError,
)

_MIN_INTERVAL_NO_KEY = 0.34   # ~3 requests/s, the unauthenticated NCBI limit
_MIN_INTERVAL_KEY = 0.11      # ~9 requests/s with an API key
_last_entrez_call = [0.0]


def set_api_key(key: str) -> None:
    """Set the NCBI API key (raises the rate limit); also picked up from $NCBI_API_KEY."""
    global _entrez_api_key
    _entrez_api_key = key or ""
    Entrez.api_key = _entrez_api_key


def _throttle() -> None:
    interval = _MIN_INTERVAL_KEY if getattr(Entrez, "api_key", "") else _MIN_INTERVAL_NO_KEY
    wait = interval - (time.time() - _last_entrez_call[0])
    if wait > 0:
        time.sleep(wait)
    _last_entrez_call[0] = time.time()


def _entrez_retry(thunk, *, what: str, retries: int = 3, base_delay: float = 0.5):
    """
    Run ``thunk`` (which opens a handle and parses it) with rate-throttling and retry-on-
    transient. Returns the thunk's value; raises ``NCBIError`` if it fails every attempt.
    """
    last = None
    for attempt in range(retries):
        _throttle()
        try:
            return thunk()
        except _TRANSIENT_ERRORS as exc:
            last = exc
            if attempt < retries - 1:
                time.sleep(base_delay * (2 ** attempt))
    raise NCBIError(f"NCBI request failed ({what}) after {retries} attempts: {last}") from last


# ── Reference library ────────────────────────────────────────────────────────

# Global (default) reference library. Per-project libraries live under <project>/references
# (D-013); every reference function takes a `ref_dir` defaulting here, so existing callers and
# the global library keep working unchanged until a project-scoped dir is passed.
REF_DIR = Path.home() / ".phylofetch" / "references"


def project_ref_dir(project_dir) -> Path:
    """Per-project reference-library root (D-013): ``<project>/references`` (created)."""
    d = Path(project_dir) / "references"
    d.mkdir(parents=True, exist_ok=True)
    return d


def locus_ref_fasta(locus_name: str, ref_dir: Path = REF_DIR) -> str:
    d = Path(ref_dir) / locus_name
    d.mkdir(parents=True, exist_ok=True)
    return str(d / f"{locus_name}_refs.fasta")


def list_loci(ref_dir: Path = REF_DIR) -> list[str]:
    ref_dir = Path(ref_dir)
    if not ref_dir.exists():
        return []
    return sorted(
        p.name for p in ref_dir.iterdir()
        if p.is_dir() and (p / f"{p.name}_refs.fasta").exists()
    )


def load_ref_records(locus_name: str, ref_dir: Path = REF_DIR) -> list:
    fasta = locus_ref_fasta(locus_name, ref_dir=ref_dir)
    if not os.path.exists(fasta):
        return []
    return list(SeqIO.parse(fasta, "fasta"))


def accessions_in_library(locus_name: str, ref_dir: Path = REF_DIR) -> set:
    accs: set = set()
    for r in load_ref_records(locus_name, ref_dir=ref_dir):
        base = r.id.split(".")[0]
        accs.add(r.id)
        accs.add(base)
    return accs


def count_refs(locus_name: str, ref_dir: Path = REF_DIR) -> int:
    return len(load_ref_records(locus_name, ref_dir=ref_dir))


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
    def _do():
        handle = Entrez.efetch(db=db, id=accession, rettype="gb", retmode="text")
        records = list(SeqIO.parse(handle, "genbank"))
        handle.close()
        return records
    # A transient failure now raises NCBIError (caught per-accession by fetch_and_store and
    # recorded as an error) instead of returning (None, None) → a false "not found" (D-034).
    records = _entrez_retry(_do, what=f"efetch {accession}")
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

# RefSeq genome-annotation restriction (D-024). `srcdb_refseq[prop]` keeps only RefSeq records
# (XP_/NP_ for eukaryotes), i.e. full-length proteins predicted from annotated genome assemblies
# — the "closest annotated genome" the References page wants for extraction guides. Verified live:
# `refseq[filter]` does NOT work for this (0 hits with a [Protein Name] query); `srcdb_refseq[prop]`
# is the property that matches (e.g. Alternaria TEF1: 4172 protein records, only 3 RefSeq — the
# rest are partial-CDS barcode translations, including a 5-aa fragment).
_REFSEQ_FILTER = "srcdb_refseq[prop]"


def _esearch_ids(db: str, term: str, max_results: int) -> list[str]:
    def _do():
        handle = Entrez.esearch(db=db, term=term, retmax=max_results)
        record = Entrez.read(handle)
        handle.close()
        return list(record.get("IdList", []))
    return _entrez_retry(_do, what=f"esearch '{term}'")


def _summary_to_dict(s) -> dict:
    title = str(s.get("Title", ""))
    organism = str(s.get("Organism", ""))
    # Protein esummary doesn't populate Organism, but the title ends with "...[Organism]"
    # (e.g. "elongation factor 1-alpha [Paraphaeosphaeria sp.]") — parse it so the candidate
    # picker can show which taxon each reference comes from (D-024; cf. the D-014 stitle finding).
    if not organism and title.endswith("]") and "[" in title:
        organism = title[title.rfind("[") + 1:-1].strip()
    return {
        "uid":       str(s.get("Id", "")),
        "accession": str(s.get("AccessionVersion", "")),
        "title":     title,
        "organism":  organism,
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
            type_uids = set(_esearch_ids(db, f"{base_query} AND {_TYPE_FILTER}", max_results))
        else:
            type_uids = set()
    if not ids:
        return []

    def _do_summary():
        handle = Entrez.esummary(db=db, id=",".join(ids))
        s = Entrez.read(handle)
        handle.close()
        return s
    summaries = _entrez_retry(_do_summary, what="esummary")

    results = _mark_type_set([_summary_to_dict(s) for s in summaries], type_uids)
    if type_mode == "prefer":
        results.sort(key=lambda r: not r["is_type"])   # type-derived first, stable
    return results


def build_entrez_query(terms, organism: str = "", field: str = "Title") -> str:
    """
    Build an Entrez query that ORs a set of synonym ``terms`` within ``field`` and ANDs
    the organism. ``terms`` may be a single ``str`` or an iterable of ``str``.

    Multi-word terms are wrapped in quotes so NCBI matches them as a phrase; without quotes
    the field tag binds only to the *last* word and the earlier words silently fall back to
    All Fields (one of the ways the old single-string query under-/over-matched). A single
    term skips the surrounding parentheses. Terms are de-duplicated case-insensitively,
    order-preserving.

        build_entrez_query(["gapdh", "glyceraldehyde-3-phosphate dehydrogenase"], "Fungi")
        -> '(gapdh[Title] OR "glyceraldehyde-3-phosphate dehydrogenase"[Title]) AND Fungi[Organism]'

    Intentionally carries NO 'complete cds' / 'partial cds' constraint: fungal phylogenetic
    markers are overwhelmingly deposited as *partial cds* barcode amplicons, so forcing
    'complete cds' silently crushed recall. The relaxed BLAST path and Exonerate both handle
    partial / intron-containing references. See DECISIONS.md D-011.
    """
    if isinstance(terms, str):
        terms = [terms]
    seen: set[str] = set()
    uniq: list[str] = []
    for t in terms:
        t = (t or "").strip()
        key = t.lower()
        if t and key not in seen:
            seen.add(key)
            uniq.append(t)
    if not uniq:
        raise ValueError("build_entrez_query: no non-empty search terms given")

    def fmt(t: str) -> str:
        return f'"{t}"[{field}]' if " " in t else f"{t}[{field}]"

    joined = " OR ".join(fmt(t) for t in uniq)
    if len(uniq) > 1:
        joined = f"({joined})"
    org = (organism or "").strip()
    return f"{joined} AND {org}[Organism]" if org else joined


def taxon_fallbacks(taxon: str) -> list[str]:
    """
    Ordered organism candidates for a reference search: the given taxon first, then its genus
    (the first whitespace token) as a fallback. Lets a search for a novel/unsequenced species
    (e.g. 'Alternaria aff. eureka', absent from NCBI) fall back to the genus 'Alternaria'.
    De-duplicated case-insensitively; empty input → [].
    """
    t = (taxon or "").strip()
    if not t:
        return []
    out = [t]
    genus = t.split()[0]
    if genus and genus.lower() != t.lower():
        out.append(genus)
    return out


def search_ncbi_protein(gene_name, organism: str,
                        max_results: int = 20, type_mode: str = "all",
                        field: str = "Protein Name", refseq_only: bool = False) -> list[dict]:
    query = build_entrez_query(gene_name, organism, field=field)
    if refseq_only:
        query = f"{query} AND {_REFSEQ_FILTER}"     # genome-annotated, full-length (D-024)
    return _search_ncbi("protein", query, max_results, type_mode)


def search_ncbi_nucleotide(gene_name, organism: str,
                           max_results: int = 20, type_mode: str = "all") -> list[dict]:
    return _search_ncbi(
        "nucleotide", build_entrez_query(gene_name, organism, field="Title"),
        max_results, type_mode,
    )


def ncbi_search_count(gene_name, organism: str, db: str = "nucleotide",
                      field: str = "Title", type_mode: str = "all",
                      refseq_only: bool = False) -> int:
    """
    Total number of NCBI matches for a (synonym) query, without fetching the IDs — a cheap
    esearch with ``retmax=0`` returning the ``Count``. Used by the References page to preview
    how many references each locus would yield before fetching. ``type_mode="type_only"``
    counts only sequence-from-type records; ``refseq_only=True`` restricts to RefSeq genome-
    annotated records (D-024); "prefer"/"all" return the full count.
    """
    _require_email()
    query = build_entrez_query(gene_name, organism, field=field)
    if refseq_only:
        query = f"{query} AND {_REFSEQ_FILTER}"
    if type_mode == "type_only":
        query = f"{query} AND {_TYPE_FILTER}"
    def _do():
        handle = Entrez.esearch(db=db, term=query, retmax=0)
        record = Entrez.read(handle)
        handle.close()
        return int(record.get("Count", 0))
    # Raise (not return 0) on a network failure so a transient blip can't masquerade as a
    # genuine "0 references for this locus" — a wrong-but-plausible result (D-034).
    return _entrez_retry(_do, what=f"count '{query}'")


# ── NCBI fetch ───────────────────────────────────────────────────────────────

def fetch_protein_by_accession(accession: str) -> Optional[SeqRecord]:
    _require_email()
    def _do():
        handle = Entrez.efetch(db="protein", id=accession, rettype="fasta", retmode="text")
        records = list(SeqIO.parse(handle, "fasta"))
        handle.close()
        return records
    records = _entrez_retry(_do, what=f"efetch protein {accession}")
    return records[0] if records else None


def fetch_nucleotide_by_accession(accession: str) -> Optional[SeqRecord]:
    _require_email()
    def _do():
        handle = Entrez.efetch(db="nucleotide", id=accession, rettype="fasta", retmode="text")
        records = list(SeqIO.parse(handle, "fasta"))
        handle.close()
        return records
    records = _entrez_retry(_do, what=f"efetch nucleotide {accession}")
    return records[0] if records else None


def fetch_and_store(accessions: list[str], locus_name: str,
                    db: str = "protein",
                    custom_label: Optional[str] = None,
                    query: str = "",
                    ref_dir: Path = REF_DIR) -> tuple[int, int, list]:
    """
    Fetch accessions, append sequences to the locus FASTA, and record per-accession
    metadata (organism / strain / voucher / type material) in the sidecar JSON.
    ``query`` is stored for provenance. ``ref_dir`` selects the reference library root
    (per-project when passed; global by default — D-013). Returns (added, skipped, errors).
    """
    _require_email()
    existing = accessions_in_library(locus_name, ref_dir=ref_dir)
    fasta_path = locus_ref_fasta(locus_name, ref_dir=ref_dir)
    meta = load_ref_meta(locus_name, ref_dir=ref_dir)
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
    save_ref_meta(locus_name, meta, ref_dir=ref_dir)
    _log_fetch(locus_name, db, query, accessions, added, skipped, ref_dir=ref_dir)
    return added, skipped, errors


def delete_from_library(locus_name: str, accession: str,
                        ref_dir: Path = REF_DIR) -> bool:
    fasta_path = locus_ref_fasta(locus_name, ref_dir=ref_dir)
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
#
# Each coding locus carries a canonical `gene` keyword plus `synonyms`: the curated set
# of name variants this marker is deposited under (symbol + spelled-out name + common
# abbreviations). `search_ncbi_nucleotide` / `_protein` OR these together (build_entrez_query)
# so inconsistently-titled records are still recovered. NO entry forces 'complete cds' —
# fungal markers are mostly *partial cds* barcode amplicons (D-011, supersedes the old
# 'complete cds' terms). Set organism in the UI per locus. This catalogue is genus-agnostic.

LOCUS_CATALOGUE: dict[str, dict] = {
    # rDNA — ITSx extracts from assemblies; fetch here for outgroup references
    "ITS":   {"db": "nucleotide", "gene": "internal transcribed spacer",
               "synonyms": ["ITS1", "ITS2", "5.8S ribosomal RNA"],
               "note": "ITSx extracts from assemblies; fetch here for outgroup taxa"},
    "LSU":   {"db": "nucleotide", "gene": "28S large subunit ribosomal RNA",
               "synonyms": ["28S ribosomal RNA", "LSU rRNA", "large subunit ribosomal RNA"],
               "note": "ITSx extracts from assemblies; fetch here for outgroup taxa"},
    "SSU":   {"db": "nucleotide", "gene": "18S small subunit ribosomal RNA",
               "synonyms": ["18S ribosomal RNA", "SSU rRNA", "small subunit ribosomal RNA"],
               "note": "ITSx extracts from assemblies; fetch here for outgroup taxa"},
    # Protein-coding — partial-CDS barcode amplicons are the norm and work fine here.
    "TEF1":  {"db": "nucleotide", "gene": "tef1",
               "synonyms": ["tef-1", "tef1-alpha", "translation elongation factor 1-alpha",
                            "elongation factor 1-alpha", "EF-1alpha"],
               "note": "Translation elongation factor 1-alpha."},
    "RPB1":  {"db": "nucleotide", "gene": "rpb1",
               "synonyms": ["RNA polymerase II largest subunit",
                            "DNA-directed RNA polymerase II subunit RPB1"],
               "note": "RNA pol II largest subunit."},
    "RPB2":  {"db": "nucleotide", "gene": "rpb2",
               "synonyms": ["RNA polymerase II second largest subunit",
                            "DNA-directed RNA polymerase II subunit RPB2"],
               "note": "RNA pol II second largest subunit."},
    "TUB2":  {"db": "nucleotide", "gene": "tub2",
               "synonyms": ["benA", "beta-tubulin", "beta tubulin", "btub"],
               "note": "Beta-tubulin. Has small exons — try 'blastn' task if exons are missed."},
    "GAPDH": {"db": "nucleotide", "gene": "gapdh",
               "synonyms": ["gpd", "gpdh", "gpdA", "glyceraldehyde-3-phosphate dehydrogenase"],
               "note": "Glyceraldehyde-3-phosphate dehydrogenase."},
    "CAL":   {"db": "nucleotide", "gene": "calmodulin",
               "synonyms": ["cmd", "cmdA", "CaM"],
               "note": "Calmodulin."},
    "ACT":   {"db": "nucleotide", "gene": "actin",
               "synonyms": ["act1", "actA"],
               "note": "Actin."},
    "HIS3":  {"db": "nucleotide", "gene": "histone H3",
               "synonyms": ["his3", "hH3"],
               "note": "Histone H3. Useful for fine-scale resolution in some groups."},
}


def locus_search_terms(locus_name: str, user_term: str = "") -> list[str]:
    """
    Ordered search terms for a locus: the user's typed keyword first (if any), then the
    catalogue canonical `gene` and curated `synonyms`. De-duplication is left to
    ``build_entrez_query``. Safe for unknown/custom loci (returns ``[user_term]`` or ``[]``).
    """
    cat = LOCUS_CATALOGUE.get(locus_name, {})
    terms: list[str] = []
    if user_term and user_term.strip():
        terms.append(user_term.strip())
    if cat.get("gene"):
        terms.append(cat["gene"])
    terms.extend(cat.get("synonyms", []))
    return terms
