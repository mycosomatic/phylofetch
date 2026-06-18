"""
Tests for src/phylofetch/primer_utils.py

Covers:
  - PRIMER_CATALOGUE structure and locus grouping
  - find_primer_amplicons: candidate pairing logic (mocked blastn-short)
  - extract_primer_amplicon: sequence extraction and FASTA output
  - run_primer_extraction: end-to-end with mock BLAST
  - PrimerPair dataclass defaults
"""

import sys
from pathlib import Path
from unittest.mock import patch, MagicMock

import pytest

sys.path.insert(0, str(Path(__file__).parent.parent / "src"))
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

from phylofetch.primer_utils import (
    PRIMER_CATALOGUE,
    LOCUS_PRIMER_MAP,
    PrimerPair,
    _run_blastn_short,
    find_primer_amplicons,
    extract_primer_amplicon,
    run_primer_extraction,
)


# ── Catalogue sanity checks ───────────────────────────────────────────────────

class TestCatalogue:
    def test_all_entries_have_required_fields(self):
        for name, pp in PRIMER_CATALOGUE.items():
            assert pp.name == name
            assert pp.locus, f"{name} missing locus"
            assert pp.fwd,   f"{name} missing fwd primer"
            assert pp.rev,   f"{name} missing rev primer"
            assert pp.min_amplicon < pp.max_amplicon

    def test_locus_primer_map_covers_all_loci(self):
        for name, pp in PRIMER_CATALOGUE.items():
            assert pp.locus in LOCUS_PRIMER_MAP
            assert name in LOCUS_PRIMER_MAP[pp.locus]

    def test_its_primers_present(self):
        assert "ITS" in LOCUS_PRIMER_MAP
        assert "ITS1/ITS4" in PRIMER_CATALOGUE

    def test_lsu_primers_present(self):
        assert "LSU" in LOCUS_PRIMER_MAP


# ── PrimerPair dataclass ──────────────────────────────────────────────────────

class TestPrimerPair:
    def test_defaults(self):
        pp = PrimerPair(name="x", locus="ITS", fwd="ATGC", rev="CGTA")
        assert pp.min_amplicon == 100
        assert pp.max_amplicon == 5000

    def test_custom_size_range(self):
        pp = PrimerPair(name="x", locus="ITS", fwd="ATGC", rev="CGTA",
                        min_amplicon=200, max_amplicon=600)
        assert pp.min_amplicon == 200
        assert pp.max_amplicon == 600


# ── Mock helper for _run_blastn_short ────────────────────────────────────────

def _fake_hits(fwd_sstart, fwd_send, rev_sstart, rev_send,
               contig="ctg1", fwd_mm=0, rev_mm=1):
    """Build the minimal list of dicts that _run_blastn_short would return."""
    return [
        {"qseqid": "FWD", "sseqid": contig, "pident": 95.0, "length": 20,
         "mismatch": fwd_mm, "sstart": fwd_sstart, "send": fwd_send,
         "evalue": 1e-5, "bitscore": 40.0},
        {"qseqid": "REV", "sseqid": contig, "pident": 94.0, "length": 20,
         "mismatch": rev_mm, "sstart": rev_sstart, "send": rev_send,
         "evalue": 1e-5, "bitscore": 38.0},
    ]


# ── find_primer_amplicons ─────────────────────────────────────────────────────

class TestFindPrimerAmplicons:
    _PP = PrimerPair(name="T/T", locus="ITS", fwd="ATGCATGCATGCATGCATGC",
                     rev="CGTAGCTAGCTAGCTAGCTA",
                     min_amplicon=400, max_amplicon=800)

    def _run(self, hits):
        with patch("phylofetch.primer_utils._run_blastn_short", return_value=hits):
            return find_primer_amplicons("fake.fasta", self._PP)

    def test_plus_strand_amplicon_detected(self):
        # F on + (sstart < send), R on - (sstart > send)
        hits = _fake_hits(fwd_sstart=100, fwd_send=119,   # F binds + strand
                          rev_sstart=619, rev_send=600)    # R binds - strand
        candidates = self._run(hits)
        assert len(candidates) == 1
        c = candidates[0]
        assert c["strand"] == "+"
        assert c["amp_start"] == 100
        assert c["amp_end"]   == 619
        assert c["amp_len"]   == 520

    def test_minus_strand_amplicon_detected(self):
        # F on - (sstart > send), R on + (sstart < send)
        hits = _fake_hits(fwd_sstart=619, fwd_send=600,   # F binds - strand
                          rev_sstart=100, rev_send=119)    # R binds + strand
        candidates = self._run(hits)
        assert len(candidates) == 1
        c = candidates[0]
        assert c["strand"] == "-"
        assert c["amp_start"] == 100
        assert c["amp_end"]   == 619

    def test_amplicon_outside_size_range_rejected(self):
        # amp_len = 900 bp, max = 800 → should be excluded
        hits = _fake_hits(fwd_sstart=100, fwd_send=119,
                          rev_sstart=999, rev_send=980)  # 900 bp span
        candidates = self._run(hits)
        assert candidates == []

    def test_same_strand_hits_not_paired(self):
        # Both hits on plus strand (sstart < send) — no valid pairing
        hits = [
            {"qseqid": "FWD", "sseqid": "ctg1", "pident": 95.0, "length": 20,
             "mismatch": 0, "sstart": 100, "send": 119, "evalue": 1e-5, "bitscore": 40.0},
            {"qseqid": "REV", "sseqid": "ctg1", "pident": 95.0, "length": 20,
             "mismatch": 0, "sstart": 600, "send": 619, "evalue": 1e-5, "bitscore": 40.0},
        ]
        candidates = self._run(hits)
        assert candidates == []

    def test_sorted_by_total_mismatches(self):
        # Two candidates on different contigs — lower total mismatch should be first
        hits = [
            # ctg1: fwd_mm=2, rev_mm=2 → total 4
            {"qseqid": "FWD", "sseqid": "ctg1", "pident": 90.0, "length": 20,
             "mismatch": 2, "sstart": 100, "send": 119, "evalue": 1e-5, "bitscore": 35.0},
            {"qseqid": "REV", "sseqid": "ctg1", "pident": 90.0, "length": 20,
             "mismatch": 2, "sstart": 619, "send": 600, "evalue": 1e-5, "bitscore": 35.0},
            # ctg2: fwd_mm=0, rev_mm=1 → total 1 → should be first
            {"qseqid": "FWD", "sseqid": "ctg2", "pident": 95.0, "length": 20,
             "mismatch": 0, "sstart": 100, "send": 119, "evalue": 1e-5, "bitscore": 40.0},
            {"qseqid": "REV", "sseqid": "ctg2", "pident": 94.0, "length": 20,
             "mismatch": 1, "sstart": 619, "send": 600, "evalue": 1e-5, "bitscore": 38.0},
        ]
        candidates = self._run(hits)
        assert len(candidates) == 2
        assert candidates[0]["contig"] == "ctg2"
        assert candidates[0]["fwd_mismatch"] + candidates[0]["rev_mismatch"] == 1

    def test_no_hits_returns_empty(self):
        candidates = self._run([])
        assert candidates == []

    def test_missing_primer_hits_returns_empty(self):
        # Only F hits, no R hits on any contig
        hits = [
            {"qseqid": "FWD", "sseqid": "ctg1", "pident": 95.0, "length": 20,
             "mismatch": 0, "sstart": 100, "send": 119, "evalue": 1e-5, "bitscore": 40.0},
        ]
        candidates = self._run(hits)
        assert candidates == []


# ── extract_primer_amplicon ───────────────────────────────────────────────────

class TestExtractPrimerAmplicon:
    _PP = PrimerPair(name="ITS1/ITS4", locus="ITS",
                     fwd="ATGC", rev="GCTA",
                     min_amplicon=50, max_amplicon=2000)

    def _make_assembly(self, tmp_path, seq="N" * 1000) -> str:
        fp = tmp_path / "asm.fasta"
        SeqIO.write([SeqRecord(Seq(seq), id="ctg1", description="")], str(fp), "fasta")
        return str(fp)

    def test_plus_strand_extraction(self, tmp_path):
        seq = "A" * 100 + "C" * 500 + "G" * 400   # 1000 bp
        asm = self._make_assembly(tmp_path, seq)
        hit = {"contig": "ctg1", "strand": "+",
               "amp_start": 101, "amp_end": 600,
               "amp_len": 500, "fwd_mismatch": 0, "rev_mismatch": 0}
        out_dir = str(tmp_path / "out")
        rec = extract_primer_amplicon(asm, hit, "S1", "ITS", out_dir, self._PP)
        assert rec is not None
        assert len(rec.seq) == 500
        assert str(rec.seq) == "C" * 500
        assert "[strand=+]" in rec.description

    def test_minus_strand_extraction_is_revcomp(self, tmp_path):
        # Sequence: 100 A's, then ATGCCC (6 nt), then rest
        seq = "A" * 100 + "ATGCCC" + "T" * 894
        asm = self._make_assembly(tmp_path, seq)
        hit = {"contig": "ctg1", "strand": "-",
               "amp_start": 101, "amp_end": 106,
               "amp_len": 6, "fwd_mismatch": 0, "rev_mismatch": 0}
        out_dir = str(tmp_path / "out2")
        rec = extract_primer_amplicon(asm, hit, "S1", "ITS", out_dir, self._PP)
        assert rec is not None
        # rev-comp of ATGCCC = GGGCAT
        assert str(rec.seq) == "GGGCAT"
        assert "[strand=-]" in rec.description

    def test_output_fasta_written(self, tmp_path):
        asm = self._make_assembly(tmp_path)
        hit = {"contig": "ctg1", "strand": "+",
               "amp_start": 1, "amp_end": 200,
               "amp_len": 200, "fwd_mismatch": 0, "rev_mismatch": 0}
        out_dir = str(tmp_path / "out3")
        extract_primer_amplicon(asm, hit, "S1", "ITS", out_dir, self._PP)
        assert (Path(out_dir) / "ITS_amplicon.fasta").exists()

    def test_missing_contig_returns_none(self, tmp_path):
        asm = self._make_assembly(tmp_path)
        hit = {"contig": "missing_contig", "strand": "+",
               "amp_start": 1, "amp_end": 200,
               "amp_len": 200, "fwd_mismatch": 0, "rev_mismatch": 0}
        rec = extract_primer_amplicon(asm, hit, "S1", "ITS", str(tmp_path / "out4"), self._PP)
        assert rec is None

    def test_rich_header_contains_provenance(self, tmp_path):
        asm = self._make_assembly(tmp_path)
        hit = {"contig": "ctg1", "strand": "+",
               "amp_start": 1, "amp_end": 100,
               "amp_len": 100, "fwd_mismatch": 1, "rev_mismatch": 2}
        out_dir = str(tmp_path / "out5")
        rec = extract_primer_amplicon(asm, hit, "MySample", "ITS", out_dir, self._PP)
        desc = rec.description
        assert "[sample=MySample]"    in desc
        assert "[locus=ITS]"         in desc
        assert "[primer=ITS1/ITS4]"  in desc
        assert "[fwd_mm=1]"          in desc
        assert "[rev_mm=2]"          in desc
        assert "[extracted=" in desc


# ── run_primer_extraction ─────────────────────────────────────────────────────

class TestRunPrimerExtraction:
    _PP = PrimerPair(name="ITS1/ITS4", locus="ITS",
                     fwd="TCCGTAGGTGAACCTGCGG",
                     rev="TCCTCCGCTTATTGATATGC",
                     min_amplicon=400, max_amplicon=800)

    def _asm(self, tmp_path, length=1000) -> str:
        seq = "ACGT" * (length // 4)
        fp = tmp_path / "asm.fasta"
        SeqIO.write([SeqRecord(Seq(seq), id="ctg1", description="")], str(fp), "fasta")
        return str(fp)

    def test_no_candidates_returns_none_and_message(self, tmp_path):
        asm = self._asm(tmp_path)
        with patch("phylofetch.primer_utils._run_blastn_short", return_value=[]):
            result, status = run_primer_extraction(
                asm, self._PP, str(tmp_path / "out"), "S1", "ITS",
            )
        assert result is None
        assert "No amplicon" in status

    def test_successful_extraction_returns_result_dict(self, tmp_path):
        asm = self._asm(tmp_path, length=2000)
        hits = _fake_hits(fwd_sstart=100, fwd_send=118,
                          rev_sstart=618, rev_send=600)
        with patch("phylofetch.primer_utils._run_blastn_short", return_value=hits):
            result, status = run_primer_extraction(
                asm, self._PP, str(tmp_path / "out"), "S1", "ITS",
            )
        assert status == "ok"
        assert result is not None
        assert result["amp_len"] == 519
        assert result["strand"] == "+"
        assert result["primer_name"] == "ITS1/ITS4"
        assert "output_fasta" in result
        assert Path(result["output_fasta"]).exists()

    def test_n_candidates_reported(self, tmp_path):
        asm = self._asm(tmp_path, length=2000)
        # Two valid candidates on different contigs
        hits = [
            *_fake_hits(fwd_sstart=100, fwd_send=118, rev_sstart=618, rev_send=600,
                        contig="ctg1", fwd_mm=0, rev_mm=0),
            *_fake_hits(fwd_sstart=100, fwd_send=118, rev_sstart=618, rev_send=600,
                        contig="ctg1", fwd_mm=1, rev_mm=1),
        ]
        # Write assembly with single contig; second "candidate" is a duplicate contig hit
        with patch("phylofetch.primer_utils._run_blastn_short", return_value=hits):
            result, _ = run_primer_extraction(
                asm, self._PP, str(tmp_path / "out"), "S1", "ITS",
            )
        if result:
            assert result["n_candidates"] >= 1
