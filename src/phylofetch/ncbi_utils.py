"""
ncbi_utils.py
-------------
Search and fetch sequences from NCBI via Biopython Entrez.

Genus-agnostic: organism defaults are configurable; no hardcoded taxon names.
Reference library stored at ~/.phylofetch/references/.
"""

import os
import time
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


# ── NCBI search ──────────────────────────────────────────────────────────────

def search_ncbi_protein(gene_name: str, organism: str,
                        max_results: int = 20) -> list[dict]:
    _require_email()
    query = f"{gene_name}[Gene Name] AND {organism}[Organism]"
    handle = Entrez.esearch(db="protein", term=query, retmax=max_results)
    record = Entrez.read(handle)
    handle.close()
    ids = record.get("IdList", [])
    if not ids:
        return []
    time.sleep(0.34)
    handle = Entrez.esummary(db="protein", id=",".join(ids))
    summaries = Entrez.read(handle)
    handle.close()
    return [
        {
            "uid":       str(s.get("Id", "")),
            "accession": str(s.get("AccessionVersion", "")),
            "title":     str(s.get("Title", "")),
            "organism":  str(s.get("Organism", "")),
            "length":    int(s.get("Length", 0)),
        }
        for s in summaries
    ]


def search_ncbi_nucleotide(gene_name: str, organism: str,
                            max_results: int = 20) -> list[dict]:
    _require_email()
    query = f"{gene_name}[Title] AND {organism}[Organism]"
    handle = Entrez.esearch(db="nucleotide", term=query, retmax=max_results)
    record = Entrez.read(handle)
    handle.close()
    ids = record.get("IdList", [])
    if not ids:
        return []
    time.sleep(0.34)
    handle = Entrez.esummary(db="nucleotide", id=",".join(ids))
    summaries = Entrez.read(handle)
    handle.close()
    return [
        {
            "uid":       str(s.get("Id", "")),
            "accession": str(s.get("AccessionVersion", "")),
            "title":     str(s.get("Title", "")),
            "organism":  str(s.get("Organism", "")),
            "length":    int(s.get("Length", 0)),
        }
        for s in summaries
    ]


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
                    custom_label: Optional[str] = None) -> tuple[int, int, list]:
    _require_email()
    existing = accessions_in_library(locus_name)
    fasta_path = locus_ref_fasta(locus_name)
    added, skipped, errors = 0, 0, []
    with open(fasta_path, "a") as out_f:
        for acc in accessions:
            base = acc.split(".")[0]
            if acc in existing or base in existing:
                skipped += 1
                continue
            try:
                if db == "protein":
                    record = fetch_protein_by_accession(acc)
                else:
                    record = fetch_nucleotide_by_accession(acc)
                if record is None:
                    errors.append(f"{acc}: not found")
                    continue
                if custom_label:
                    record.id = custom_label
                    record.description = ""
                SeqIO.write(record, out_f, "fasta")
                added += 1
            except Exception as e:
                errors.append(f"{acc}: {e}")
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
