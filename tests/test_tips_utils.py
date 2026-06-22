"""
Tests for src/phylofetch/tips_utils.py — tip locus auto-classification (RM-008 c3 / D-020).

Network-free: covers classify_locus on representative GenBank titles and the tips dir helper.
The fetch/import paths (classify_accessions, import_tip_accessions) need Entrez and are not
exercised here.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / "src"))
from phylofetch.tips_utils import classify_locus, project_tips_dir


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


def test_project_tips_dir(tmp_path):
    d = project_tips_dir(tmp_path / "proj")
    assert d == (tmp_path / "proj" / "tips") and d.is_dir()
