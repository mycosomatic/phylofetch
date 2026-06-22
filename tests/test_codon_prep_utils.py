"""
Tests for src/phylofetch/codon_prep_utils.py — RM-008 component 2 (D-022).

Two layers, mirroring tests/test_exonerate_utils.py:
  * Offline — soft-mask / exon-mark rendering and the per-locus orchestration
    (merge isolates + framed tips, strict-QC exclusion, locus gating) run with no
    exonerate binary: the rendering tests build a model from the committed verbatim
    GFF fixtures, and the orchestration tests monkeypatch the single-amplicon
    framing so prepare_codon_locus is exercised deterministically.
  * End-to-end — frame_consistent_amplicon against the real exonerate binary on the
    synthetic 2-exon gene fixtures; skipped where exonerate is absent.
"""

import shutil
import sys
from pathlib import Path

import pytest
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

sys.path.insert(0, str(Path(__file__).parent.parent / "src"))
import phylofetch.codon_prep_utils as cp
from phylofetch.codon_prep_utils import (
    coding_loci_with_tips,
    exon_marked_cds,
    frame_consistent_amplicon,
    prepare_codon_locus,
    _tip_seqrecord,
)
from phylofetch.exonerate_utils import build_result_from_model, parse_exonerate_gff
from phylofetch.ncbi_utils import locus_ref_fasta

FIX = Path(__file__).parent / "fixtures"
HAVE_EXONERATE = shutil.which("exonerate") is not None


def _read1(path: Path) -> str:
    return "".join(ln.strip() for ln in Path(path).read_text().splitlines()
                   if not ln.startswith(">"))


def _model(tag: str):
    contig = _read1(FIX / f"contig_{tag}.fasta")
    res = build_result_from_model(
        parse_exonerate_gff((FIX / f"exo_protein2genome_{tag}.gff").read_text())[0],
        contig, "ST1", "TEF1")
    return contig, res


def _write_tips(tips_dir: Path, locus: str, seqs: dict[str, str]) -> None:
    """Write {acc: seq} as the locus tip library FASTA (no metadata sidecar)."""
    fasta = locus_ref_fasta(locus, ref_dir=tips_dir)
    SeqIO.write([SeqRecord(Seq(s), id=acc, description="") for acc, s in seqs.items()],
                fasta, "fasta")


# ── Boundary-aware rendering (offline) ────────────────────────────────────────
# NB: the genomic soft-masking (exon-upper / intron-lower) now lives in
# exonerate_utils.build_result_from_model so isolates and tips are masked identically — its
# correctness is covered in tests/test_exonerate_utils.py (TestSoftMaskGenomic). Here we only
# confirm the framed-tip record carries that masked sequence through.

class TestExonMarkedCds:
    def test_upper_roundtrips_and_marks_junction(self):
        _contig, res = _model("plus")
        marked = exon_marked_cds(res)
        assert marked.upper() == res["cds_seq"]
        # 2 exons → an upper run followed by a lower run
        assert any(c.islower() for c in marked) and any(c.isupper() for c in marked)


def test_tip_seqrecord_header_is_rich():
    _contig, res = _model("plus")
    rec = _tip_seqrecord(
        {"id": "KC584.1", "organism": "Alternaria alternata", "cds": res["cds_seq"],
         "genomic": res["amplicon_seq"],            # masked at source (build_result_from_model)
         "protein": res["protein"], "n_exons": 2, "n_introns": 1,
         "cds_length": res["cds_length"], "n_internal_stops": 0, "verdict": "PASS"},
        "TEF1", "genomic")
    assert rec.id == "KC584.1"
    for token in ("[tip=yes]", "[locus=TEF1]", "[type=genomic]",
                  "[organism=Alternaria alternata]", "[qc=PASS]"):
        assert token in rec.description


# ── Locus gating: only guided loci with tips ──────────────────────────────────

def test_coding_loci_with_tips_excludes_rdna(tmp_path):
    tips = tmp_path / "tips"
    _write_tips(tips, "TEF1", {"a": "ACGT"})     # coding → has a protein guide
    _write_tips(tips, "ITS", {"b": "ACGT"})      # rDNA  → no guide, must be excluded
    got = coding_loci_with_tips(tips)
    assert "TEF1" in got and "ITS" not in got


# ── prepare_codon_locus orchestration (offline, framing monkeypatched) ────────

def _fake_framer(verdicts: dict[str, bool]):
    """Return a frame_consistent_amplicon stand-in keyed on accession → qc_clean."""
    def _fake(seq, seq_id, guide_fasta, locus, *, organism="", mark_exons=False, **kw):
        clean = verdicts.get(seq_id, True)
        rec = {"id": seq_id, "organism": organism, "cds": "ATGAAA", "genomic": "ATGaaa",
               "protein": "MK", "n_exons": 1, "n_introns": 1, "cds_length": 6,
               "len_mod3": 0, "n_internal_stops": 0 if clean else 1,
               "qc_clean": clean, "verdict": "PASS" if clean else "REVIEW"}
        return rec, ("OK" if clean else "OK (QC: review)")
    return _fake


class TestPrepareCodonLocus:
    def test_no_guide_short_circuits(self, tmp_path):
        tips = tmp_path / "tips"
        _write_tips(tips, "ITS", {"a": "ACGT"})
        out = prepare_codon_locus("ITS", tips, tmp_path / "out")
        assert out["status"] == "no guide" and out["outputs"] == {}

    def test_no_tips_short_circuits(self, tmp_path):
        out = prepare_codon_locus("TEF1", tmp_path / "tips", tmp_path / "out")
        assert out["status"] == "no tips"

    def test_merges_isolates_and_tips_and_writes_three_matrices(self, tmp_path, monkeypatch):
        tips = tmp_path / "tips"
        _write_tips(tips, "TEF1", {"KC1": "ACGTACGT", "KC2": "ACGTACGT"})
        iso = tmp_path / "combined"
        iso.mkdir()
        for prod in ("CDS", "genomic", "protein"):
            SeqIO.write([SeqRecord(Seq("ATGAAA"), id=f"ST1_TEF1_{prod}", description="")],
                        str(iso / f"TEF1_{prod}_combined.fasta"), "fasta")
        monkeypatch.setattr(cp, "frame_consistent_amplicon", _fake_framer({}))

        out = prepare_codon_locus("TEF1", tips, tmp_path / "out",
                                  isolate_combined_dir=iso)
        assert out["n_tips"] == 2 and out["n_framed"] == 2 and out["n_failed"] == 0
        assert out["n_isolates"] == 1 and len(out["rows"]) == 2
        for prod in ("CDS", "genomic", "protein"):
            recs = list(SeqIO.parse(out["outputs"][prod], "fasta"))
            assert len(recs) == 3                      # 1 isolate + 2 tips
            assert recs[0].id == f"ST1_TEF1_{prod}"    # isolate first
            assert {r.id for r in recs[1:]} == {"KC1", "KC2"}

    def test_strict_qc_excludes_flagged_tip(self, tmp_path, monkeypatch):
        tips = tmp_path / "tips"
        _write_tips(tips, "TEF1", {"GOOD": "ACGT", "BAD": "ACGT"})
        monkeypatch.setattr(cp, "frame_consistent_amplicon",
                            _fake_framer({"GOOD": True, "BAD": False}))

        loose = prepare_codon_locus("TEF1", tips, tmp_path / "loose")
        assert loose["n_framed"] == 2 and loose["n_flagged"] == 1
        assert len(list(SeqIO.parse(loose["outputs"]["CDS"], "fasta"))) == 2

        strict = prepare_codon_locus("TEF1", tips, tmp_path / "strict", strict_qc=True)
        assert strict["n_framed"] == 1 and strict["n_flagged"] == 1
        ids = {r.id for r in SeqIO.parse(strict["outputs"]["CDS"], "fasta")}
        assert ids == {"GOOD"} and "BAD" not in ids


# ── End-to-end against the real exonerate binary ──────────────────────────────

@pytest.mark.skipif(not HAVE_EXONERATE, reason="exonerate not installed")
class TestFrameConsistentAmpliconE2E:
    @pytest.mark.parametrize("tag", ["plus", "minus"])
    def test_frames_amplicon_to_clean_cds(self, tag):
        seq = _read1(FIX / f"contig_{tag}.fasta")
        expected_cds = (FIX / "expected_cds.txt").read_text().strip()
        expected_prot = (FIX / "expected_protein.txt").read_text().strip()
        rec, status = frame_consistent_amplicon(
            seq, f"ACC_{tag}", str(FIX / "ref_protein.fasta"), "TEF1",
            organism="Synthetica test")
        assert rec is not None, status
        assert status.startswith("OK")
        assert rec["cds"] == expected_cds              # intron-stripped, frame-pinned
        assert rec["protein"] == expected_prot
        assert rec["qc_clean"] and rec["verdict"] == "PASS"
        assert rec["n_introns"] == 1
        # genomic carries exon/intron case-masking; bases match the amplicon span
        assert any(c.islower() for c in rec["genomic"])
        assert rec["genomic"].upper().replace("-", "") != ""

    def test_unalignable_amplicon_reported_not_framed(self):
        # A non-coding repeat has no ORF matching the protein guide → no model.
        rec, status = frame_consistent_amplicon(
            "TA" * 200, "JUNK", str(FIX / "ref_protein.fasta"), "TEF1")
        assert rec is None
        assert "no gene model" in status or "exonerate failed" in status

    def test_prepare_codon_locus_end_to_end(self, tmp_path, monkeypatch):
        # Inject the synthetic-gene protein as the TEF1 guide so the real exonerate
        # path runs through prepare_codon_locus end to end (the bundled fungal TEF1
        # guide wouldn't match the synthetic contig).
        prot = _read1(FIX / "ref_protein.fasta")
        guides = {"TEF1": [{"accession": "SYN1", "organism": "Synthetica",
                            "clade": "test", "protein_name": "TEF1",
                            "length": len(prot), "seq": prot}]}
        monkeypatch.setattr("phylofetch.protein_guide_utils.load_protein_guides",
                            lambda include_user=True, user_path=None: guides)

        tips = tmp_path / "tips"
        _write_tips(tips, "TEF1", {"ACC_plus": _read1(FIX / "contig_plus.fasta"),
                                   "ACC_minus": _read1(FIX / "contig_minus.fasta")})
        out = prepare_codon_locus("TEF1", tips, tmp_path / "out")
        assert out["status"] == "ok"
        assert out["n_framed"] == 2, out["rows"]
        cds = list(SeqIO.parse(out["outputs"]["CDS"], "fasta"))
        assert {r.id for r in cds} == {"ACC_plus", "ACC_minus"}
        expected_cds = (FIX / "expected_cds.txt").read_text().strip()
        assert all(str(r.seq) == expected_cds for r in cds)
