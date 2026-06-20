"""
Tests for src/phylofetch/ncbi_utils.py — reference metadata layer (RM-001, increment 1).

Network-free: covers INSDC /type_material parsing, GenBank source-feature metadata
extraction, the type-set flagging used by search, and the reference-metadata sidecar
round-trip. Live Entrez calls (search/fetch) are not exercised here.
"""

import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).parent.parent / "src"))
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature
from Bio.SeqRecord import SeqRecord

from phylofetch.ncbi_utils import (
    LOCUS_CATALOGUE,
    RefRecord,
    _mark_type_set,
    build_entrez_query,
    load_ref_meta,
    locus_search_terms,
    normalize_type_kind,
    parse_source_metadata,
    save_ref_meta,
)


# ── /type_material parsing ────────────────────────────────────────────────────

class TestNormalizeTypeKind:
    def test_absent_qualifier(self):
        assert normalize_type_kind("") == (False, "")
        assert normalize_type_kind("   ") == (False, "")

    def test_holotype(self):
        assert normalize_type_kind("holotype of Aspergillus flavus") == (True, "holotype")

    def test_ex_type_culture_is_type_grade(self):
        # 'culture from holotype' = ex-holotype, still type-grade (D-007)
        is_type, kind = normalize_type_kind("culture from holotype of Fusarium oxysporum")
        assert is_type is True
        assert kind == "ex-holotype"

    def test_explicit_ex_type(self):
        assert normalize_type_kind("ex-type") == (True, "ex-type")

    def test_type_strain(self):
        assert normalize_type_kind("type strain of Penicillium chrysogenum") == (True, "type strain")

    def test_more_specific_keyword_wins(self):
        # 'isolectotype' must be matched before the substring 'lectotype'
        assert normalize_type_kind("isolectotype of X") == (True, "isolectotype")

    def test_unrecognised_but_present_is_still_type(self):
        assert normalize_type_kind("type material") == (True, "type")


# ── GenBank source-feature metadata ───────────────────────────────────────────

class TestParseSourceMetadata:
    def _record(self, quals, organism="Aspergillus flavus"):
        rec = SeqRecord(Seq("ACGT" * 20), id="MN123456.1")
        rec.annotations["organism"] = organism
        rec.features.append(SeqFeature(type="source", qualifiers=quals))
        return rec

    def test_full_metadata(self):
        rec = self._record({
            "organism": ["Aspergillus flavus"],
            "strain": ["NRRL 3357"],
            "culture_collection": ["NRRL:3357"],
            "type_material": ["culture from holotype of Aspergillus flavus"],
        })
        meta = parse_source_metadata(rec)
        assert meta["organism"] == "Aspergillus flavus"
        assert meta["strain"] == "NRRL 3357"
        assert meta["culture_collection"] == "NRRL:3357"
        assert meta["is_type"] is True
        assert meta["type_kind"] == "ex-holotype"

    def test_no_source_feature(self):
        rec = SeqRecord(Seq("ACGT" * 10), id="X1.1")
        meta = parse_source_metadata(rec)
        assert meta["is_type"] is False
        assert meta["type_kind"] == ""
        assert meta["strain"] == ""

    def test_non_type_record(self):
        rec = self._record({"organism": ["Cladosporium sp."], "strain": ["ABC1"]})
        meta = parse_source_metadata(rec)
        assert meta["is_type"] is False
        assert meta["strain"] == "ABC1"


# ── search type-flagging ──────────────────────────────────────────────────────

class TestMarkTypeSet:
    def test_flags_by_uid(self):
        results = [{"uid": "1"}, {"uid": "2"}, {"uid": "3"}]
        out = _mark_type_set(results, {"2"})
        assert [r["is_type"] for r in out] == [False, True, False]


# ── Entrez query construction (D-011) ─────────────────────────────────────────

class TestBuildEntrezQuery:
    def test_single_term_with_organism(self):
        assert build_entrez_query("gapdh", "Fungi") == "gapdh[Title] AND Fungi[Organism]"

    def test_single_term_no_organism(self):
        assert build_entrez_query("gapdh", "") == "gapdh[Title]"
        assert build_entrez_query("gapdh", "   ") == "gapdh[Title]"

    def test_multiword_term_is_phrase_quoted(self):
        # the field tag must bind to the whole phrase, not just the last word
        q = build_entrez_query("glyceraldehyde-3-phosphate dehydrogenase", "Fungi")
        assert q == '"glyceraldehyde-3-phosphate dehydrogenase"[Title] AND Fungi[Organism]'

    def test_synonym_or_group_parenthesised(self):
        q = build_entrez_query(
            ["gapdh", "gpd", "glyceraldehyde-3-phosphate dehydrogenase"], "Aspergillus")
        assert q == ('(gapdh[Title] OR gpd[Title] OR '
                     '"glyceraldehyde-3-phosphate dehydrogenase"[Title]) '
                     'AND Aspergillus[Organism]')

    def test_case_insensitive_dedupe_preserves_order(self):
        q = build_entrez_query(["GAPDH", "gapdh", "gpd", "GPD"], "")
        assert q == "(GAPDH[Title] OR gpd[Title])"

    def test_blank_terms_dropped(self):
        assert build_entrez_query(["", "  ", "rpb2"], "Fungi") == "rpb2[Title] AND Fungi[Organism]"

    def test_no_terms_raises(self):
        with pytest.raises(ValueError):
            build_entrez_query(["", "   "], "Fungi")

    def test_protein_field(self):
        assert build_entrez_query("rpb2", "Fungi", field="Gene Name") == \
            "rpb2[Gene Name] AND Fungi[Organism]"

    def test_no_complete_cds_constraint_injected(self):
        # the whole point of D-011: the builder never adds a completeness filter
        q = build_entrez_query(["tef1", "tef-1"], "Fungi")
        assert "complete cds" not in q.lower()


class TestLocusSearchTerms:
    def test_user_term_first_then_catalogue(self):
        terms = locus_search_terms("GAPDH", user_term="gpdA")
        assert terms[0] == "gpdA"
        assert "gapdh" in terms              # canonical gene
        assert "glyceraldehyde-3-phosphate dehydrogenase" in terms  # a synonym

    def test_no_user_term_uses_catalogue_only(self):
        terms = locus_search_terms("RPB2")
        assert terms[0] == "rpb2"
        assert "RNA polymerase II second largest subunit" in terms

    def test_unknown_locus_returns_user_term_only(self):
        assert locus_search_terms("ZZZ", user_term="foo") == ["foo"]
        assert locus_search_terms("ZZZ") == []


class TestCatalogueIsBarcodeFriendly:
    """Regression guards for D-011: no forced completeness; coding loci carry synonyms."""

    CODING = [k for k in LOCUS_CATALOGUE if k not in ("ITS", "LSU", "SSU")]

    def test_no_complete_cds_anywhere(self):
        for name, cat in LOCUS_CATALOGUE.items():
            blob = " ".join([cat.get("gene", ""), *cat.get("synonyms", []),
                             cat.get("note", "")]).lower()
            assert "complete cds" not in blob, f"{name} still forces 'complete cds'"

    def test_coding_loci_have_synonyms(self):
        for name in self.CODING:
            assert LOCUS_CATALOGUE[name].get("synonyms"), f"{name} has no synonyms"


# ── sidecar round-trip ────────────────────────────────────────────────────────

class TestSidecarRoundTrip:
    def test_save_and_load(self, tmp_path):
        records = {
            "MN1.1": RefRecord(
                accession="MN1.1", organism="Aspergillus flavus", strain="NRRL 3357",
                type_material="holotype of Aspergillus flavus", type_kind="holotype",
                is_type=True, length=600, db="nucleotide",
            ),
            "MN2.1": RefRecord(
                accession="MN2.1", organism="Aspergillus oryzae", is_type=False, length=580,
            ),
        }
        save_ref_meta("TEF1", records, ref_dir=tmp_path)
        loaded = load_ref_meta("TEF1", ref_dir=tmp_path)
        assert set(loaded) == {"MN1.1", "MN2.1"}
        assert loaded["MN1.1"].is_type is True
        assert loaded["MN1.1"].type_kind == "holotype"
        assert loaded["MN1.1"].accession == "MN1.1"           # key re-applied on load
        assert loaded["MN2.1"].is_type is False
        assert loaded["MN2.1"].organism == "Aspergillus oryzae"

    def test_load_missing_returns_empty(self, tmp_path):
        assert load_ref_meta("NOPE", ref_dir=tmp_path) == {}
