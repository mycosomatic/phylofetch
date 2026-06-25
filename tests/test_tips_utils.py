"""
Tests for src/phylofetch/tips_utils.py — tip locus auto-classification (RM-008 c3 / D-020).

Network-free: covers classify_locus on representative GenBank titles and the tips dir helper.
The fetch/import paths (classify_accessions, import_tip_accessions) need Entrez and are not
exercised here.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / "src"))
import phylofetch.tips_utils as tips_utils
from phylofetch.tips_utils import (
    classify_locus,
    import_tips_with_assignments,
    lookup_accessions,
    normalize_accession,
    project_tips_dir,
)


class TestClassifyLocus:
    def test_coding_markers_from_typical_titles(self):
        cases = {
            "Alternaria alternata strain X RNA polymerase II second largest subunit "
            "(rpb2) gene, partial cds": "RPB2",
            "Alternaria sp. glyceraldehyde-3-phosphate dehydrogenase (gpd) gene, partial cds": "GAPDH",
            "Fungus translation elongation factor 1-alpha (tef1) gene, partial cds": "TEF1",
            "Isolate Y beta-tubulin (tub2) gene, partial cds": "TUB2",
            "Isolate Y actin (act) gene, partial cds": "ACT",
            "Isolate Y calmodulin gene, partial cds": "CAL",
            "Isolate Y histone H3 (his3) gene, partial cds": "HIS3",
            "Isolate Y RNA polymerase II largest subunit (rpb1) gene, partial cds": "RPB1",
        }
        for title, loc in cases.items():
            assert classify_locus(title) == loc, f"{title!r} -> {classify_locus(title)} (want {loc})"

    def test_rdna_titles(self):
        assert classify_locus("Fungus internal transcribed spacer 1, 5.8S ribosomal RNA gene, "
                              "and internal transcribed spacer 2, complete sequence") == "ITS"
        assert classify_locus("Fungus 28S large subunit ribosomal RNA gene, partial sequence") == "LSU"
        assert classify_locus("Fungus 18S small subunit ribosomal RNA gene, partial sequence") == "SSU"

    def test_bare_symbol_classifies(self):
        assert classify_locus("Alternaria sp. gpd gene, partial cds") == "GAPDH"
        assert classify_locus("Alternaria sp. partial rpb2 for RNA polymerase") == "RPB2"

    def test_rpb1_not_confused_with_rpb2(self):
        assert classify_locus("RNA polymerase II largest subunit gene") == "RPB1"
        assert classify_locus("RNA polymerase II second largest subunit gene") == "RPB2"

    def test_no_match_returns_none(self):
        assert classify_locus("Homo sapiens beta-globin mRNA") is None
        assert classify_locus("Saccharomyces cerevisiae chromosome IV, complete sequence") is None
        assert classify_locus("") is None


class TestNormalizeAccession:
    """D-026: repair bare RefSeq accessions; leave GenBank ids untouched."""

    def test_refseq_underscore_repaired(self):
        assert normalize_accession("NR135944") == "NR_135944"      # the reported case
        assert normalize_accession("NR135944.1") == "NR_135944.1"  # keeps the version
        assert normalize_accession("nr135944") == "NR_135944"      # case-insensitive
        assert normalize_accession("XM123456") == "XM_123456"      # other RefSeq prefix

    def test_already_correct_and_genbank_untouched(self):
        assert normalize_accession("NR_135944") == "NR_135944"     # already underscored
        assert normalize_accession("KC584456") == "KC584456"       # GenBank 2+6
        assert normalize_accession("KC584456.1") == "KC584456.1"
        assert normalize_accession("MN908947") == "MN908947"
        assert normalize_accession("AB12345") == "AB12345"         # too few digits → as-is
        assert normalize_accession("  NR135944 ") == "NR_135944"   # whitespace stripped


class TestLookupAccessions:
    """D-026: per-accession lookup normalizes, flags not-found, suggests a locus."""

    def test_lookup_rows_found_missing_and_guess(self, monkeypatch):
        # Stub the network: NR_135944 (LSU rRNA) resolves; BADACC1 returns nothing.
        def fake_titles(accs, db):
            assert "NR_135944" in accs                      # normalized before lookup
            return {"NR_135944": "Fungus 28S large subunit ribosomal RNA gene, partial sequence"}
        monkeypatch.setattr(tips_utils, "_esummary_titles", fake_titles)

        rows = lookup_accessions(["NR135944", "BADACC1"])
        by_input = {r["input"]: r for r in rows}
        assert by_input["NR135944"]["accession"] == "NR_135944"   # repaired
        assert by_input["NR135944"]["found"] is True
        assert by_input["NR135944"]["locus_guess"] == "LSU"       # classified from title
        assert by_input["BADACC1"]["found"] is False              # not found → flagged
        assert by_input["BADACC1"]["locus_guess"] is None


class TestImportWithAssignments:
    """D-026: per-accession locus assignment (not just one bulk locus)."""

    def test_each_accession_to_its_own_locus(self, monkeypatch):
        calls = []

        def fake_fetch_and_store(group, locus, db, query, ref_dir):
            calls.append((locus, list(group)))
            return len(group), 0, []
        monkeypatch.setattr(tips_utils, "fetch_and_store", fake_fetch_and_store)

        res = import_tips_with_assignments(
            {"NR135944": "LSU", "KC584456": "TEF1", "DQ677938": ""},   # last one skipped
            tips_dir="/tmp/whatever")
        grouped = {loc: accs for loc, accs in calls}
        assert grouped["LSU"] == ["NR_135944"]      # normalized on the way in
        assert grouped["TEF1"] == ["KC584456"]
        assert "DQ677938" not in str(calls)         # blank-locus row not imported
        assert set(res["assigned"]) == {"LSU", "TEF1"}
        assert res["errors"] == []


def test_project_tips_dir(tmp_path):
    d = project_tips_dir(tmp_path / "proj")
    assert d == (tmp_path / "proj" / "tips") and d.is_dir()
