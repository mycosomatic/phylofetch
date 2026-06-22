"""
Tests for src/phylofetch/protein_guide_utils.py — bundled protein guides (D-020 / RM-008).

Network-free: validates the packaged universal-core guide set (one Ascomycota + one
Basidiomycota full-length ortholog per conserved marker), the loader, FASTA writing, and
user lineage-pack merging.
"""

import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / "src"))
from Bio import SeqIO

from phylofetch.protein_guide_utils import (
    get_guides,
    guide_loci,
    load_protein_guides,
    write_guide_fasta,
)

EXPECTED_LOCI = {"RPB1", "RPB2", "TEF1", "TUB2", "ACT", "GAPDH", "CAL", "HIS3"}
_AA = set("ACDEFGHIKLMNPQRSTVWYXBZJUO")


class TestBundledGuideIntegrity:
    def test_expected_conserved_markers_present(self):
        g = load_protein_guides(include_user=False)
        assert EXPECTED_LOCI <= set(g)

    def test_each_locus_has_asco_and_basidio(self):
        g = load_protein_guides(include_user=False)
        for loc in EXPECTED_LOCI:
            clades = {rec["clade"] for rec in g[loc]}
            assert "Ascomycota" in clades, f"{loc} missing Ascomycota guide"
            assert "Basidiomycota" in clades, f"{loc} missing Basidiomycota guide"

    def test_guides_are_full_length_valid_protein(self):
        g = load_protein_guides(include_user=False)
        for loc, recs in g.items():
            for r in recs:
                seq = r["seq"]
                assert r["accession"] and r["seq"], f"{loc}: missing accession/seq"
                assert len(seq) == r["length"], f"{loc}/{r['accession']}: length mismatch"
                assert len(seq) > 50, f"{loc}/{r['accession']}: implausibly short guide"
                assert "*" not in seq, f"{loc}/{r['accession']}: internal stop in guide"
                assert set(seq.upper()) <= _AA, f"{loc}/{r['accession']}: non-AA characters"

    def test_accessions_are_refseq(self):
        g = load_protein_guides(include_user=False)
        for recs in g.values():
            for r in recs:
                assert r["accession"][:3] in ("XP_", "NP_", "WP_", "YP_")


class TestLoader:
    def test_guide_loci_sorted(self):
        loci = guide_loci(include_user=False)
        assert loci == sorted(loci)
        assert EXPECTED_LOCI <= set(loci)

    def test_get_guides_known_and_unknown(self):
        assert len(get_guides("RPB2", include_user=False)) >= 2
        assert get_guides("NOPE", include_user=False) == []


class TestWriteGuideFasta:
    def test_writes_all_guides_for_locus(self, tmp_path):
        out = tmp_path / "RPB2_guide.fasta"
        path = write_guide_fasta("RPB2", str(out), include_user=False)
        assert path and out.exists()
        recs = list(SeqIO.parse(str(out), "fasta"))
        assert len(recs) == len(get_guides("RPB2", include_user=False))
        assert all(len(r.seq) > 50 for r in recs)
        assert any("clade=Ascomycota" in r.description for r in recs)

    def test_unknown_locus_returns_none(self, tmp_path):
        assert write_guide_fasta("NOPE", str(tmp_path / "x.fasta"), include_user=False) is None


class TestUserLineagePackMerge:
    def test_user_pack_extends_and_adds(self, tmp_path):
        user = tmp_path / "protein_guides.json"
        user.write_text(json.dumps({"guides": {
            "RPB2": [{"accession": "USER1", "organism": "My fungus", "clade": "Ascomycota",
                      "protein_name": "RPB2", "length": 4, "seq": "MEEK"}],
            "ALTA1": [{"accession": "USER2", "organism": "Alternaria sp.", "clade": "Ascomycota",
                       "protein_name": "Alt a 1", "length": 4, "seq": "MKLT"}],
        }}))
        g = load_protein_guides(include_user=True, user_path=user)
        # appended to an existing locus
        assert any(r["accession"] == "USER1" for r in g["RPB2"])
        assert len(g["RPB2"]) == len(get_guides("RPB2", include_user=False)) + 1
        # new lineage-pack locus added
        assert "ALTA1" in g and g["ALTA1"][0]["accession"] == "USER2"

    def test_broken_user_pack_ignored(self, tmp_path):
        user = tmp_path / "protein_guides.json"
        user.write_text("{ not valid json")
        g = load_protein_guides(include_user=True, user_path=user)
        assert EXPECTED_LOCI <= set(g)        # built-ins still load
