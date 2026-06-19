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
    RefRecord,
    _mark_type_set,
    load_ref_meta,
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
