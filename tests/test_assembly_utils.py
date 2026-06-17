"""
Tests for src/phylofetch/assembly_utils.py
"""

import sys
from pathlib import Path

import pytest
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

sys.path.insert(0, str(Path(__file__).parent.parent / "src"))
from phylofetch.assembly_utils import get_assembly_stats, suggest_strain_id


def _write_fasta(path: Path, records: list[tuple[str, str]]) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    SeqIO.write(
        [SeqRecord(Seq(seq), id=rid, description="") for rid, seq in records],
        str(path),
        "fasta",
    )
    return path


class TestGetAssemblyStats:
    def test_basic_stats(self, tmp_path):
        # Two contigs: 1000 bp + 500 bp → total 1500 bp
        fp = _write_fasta(
            tmp_path / "asm.fasta",
            [("contig_1", "A" * 1000), ("contig_2", "T" * 500)],
        )
        stats = get_assembly_stats(str(fp))
        assert stats["num_contigs"] == 2
        # total_length_mb is rounded to 2 decimal places
        assert abs(stats["total_length_mb"] - round(1500 / 1e6, 2)) < 1e-4

    def test_n50_two_contigs(self, tmp_path):
        # 1000 + 500 = 1500; half = 750; cumulative at 1000 ≥ 750 → N50 = 1000
        fp = _write_fasta(
            tmp_path / "n50.fasta",
            [("c1", "A" * 1000), ("c2", "A" * 500)],
        )
        stats = get_assembly_stats(str(fp))
        assert stats["n50"] == 1000

    def test_single_contig(self, tmp_path):
        fp = _write_fasta(tmp_path / "single.fasta", [("c1", "ATGC" * 250)])
        stats = get_assembly_stats(str(fp))
        assert stats["num_contigs"] == 1
        assert stats["n50"] == 1000

    def test_gc_percent(self, tmp_path):
        # 50% GC: equal A/T/G/C
        fp = _write_fasta(tmp_path / "gc.fasta", [("c1", "ATGC" * 100)])
        stats = get_assembly_stats(str(fp))
        assert abs(stats["mean_gc"] - 50.0) < 1.0

    def test_high_gc(self, tmp_path):
        fp = _write_fasta(tmp_path / "hgc.fasta", [("c1", "GC" * 500)])
        stats = get_assembly_stats(str(fp))
        assert stats["mean_gc"] > 90.0

    def test_n50_three_contigs(self, tmp_path):
        # N50 algorithm: sort descending [800, 600, 400], cumsum until ≥ total/2 (900)
        # cumsum at 800 = 800 < 900; cumsum at 600 = 1400 ≥ 900 → N50 = 600
        fp = _write_fasta(
            tmp_path / "n50.fasta",
            [("c1", "A" * 800), ("c2", "A" * 600), ("c3", "A" * 400)],
        )
        stats = get_assembly_stats(str(fp))
        # N50 is the first length where cumsum >= total/2
        # total=1800, half=900; cumsum hits 800 first (< 900), then 1400 (≥ 900 at c2)
        assert stats["n50"] == 600


class TestSuggestStrainId:
    def test_fasta_extension_stripped(self):
        sid = suggest_strain_id("/data/samples/CBS_123.fasta")
        assert "fasta" not in sid.lower()

    def test_returns_string(self):
        assert isinstance(suggest_strain_id("/path/to/my_assembly.fa"), str)

    def test_non_empty(self):
        assert suggest_strain_id("/some/path/sample.fna").strip() != ""
