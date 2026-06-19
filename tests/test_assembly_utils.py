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
from phylofetch.assembly_utils import (
    detect_assembler,
    find_quast_report,
    get_assembly_stats,
    parse_quast_report,
    source_tag_from_filename,
    suggest_strain_id,
    suggest_unique_strain_ids,
)


# A trimmed QUAST report.tsv (tab-separated), single-assembly form.
_QUAST_TSV = (
    "Assembly\tNS26-3-C2_final_EGAP_assembly\n"
    "# contigs (>= 0 bp)\t312\n"
    "# contigs\t118\n"
    "Largest contig\t2845300\n"
    "Total length\t34521900\n"
    "GC (%)\t51.23\n"
    "N50\t1204567\n"
    "N90\t245100\n"
    "L50\t9\n"
    "L90\t34\n"
    "# N's per 100 kbp\t12.40\n"
)


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

    def test_egap_suffix_stripped(self):
        # EGAP final assembly filename → clean strain ID
        sid = suggest_strain_id("/data/NS26-3-C2/NS26-3-C2_final_EGAP_assembly.fasta")
        assert sid == "NS26-3-C2"


class TestDetectAssemblerFromFilename:
    def test_egap_detected_from_filename(self, tmp_path):
        # EGAP re-headers its output, so headers won't match — filename must win.
        fp = _write_fasta(
            tmp_path / "NS27-2B-B3_final_EGAP_assembly.fasta",
            [("scaffold_1", "ACGT" * 100)],
        )
        assert detect_assembler(str(fp)) == "egap"

    def test_masurca_detected_from_filename(self, tmp_path):
        fp = _write_fasta(
            tmp_path / "NS26-3-C2_masurca.fasta",
            [("jcf7180000000001", "ACGT" * 100)],
        )
        assert detect_assembler(str(fp)) == "masurca"

    def test_spades_headers_beat_filename(self, tmp_path):
        # A pilon-polished SPAdes assembly keeps NODE_ headers → still "spades".
        fp = _write_fasta(
            tmp_path / "NS26-3-C2_pilon.fasta",
            [("NODE_1_length_5000_cov_42.1", "ACGT" * 100)],
        )
        assert detect_assembler(str(fp)) == "spades"

    def test_unknown_when_no_marker(self, tmp_path):
        fp = _write_fasta(
            tmp_path / "mystery.fasta",
            [("seq1", "ACGT" * 100)],
        )
        assert detect_assembler(str(fp)) == "unknown"


class TestSourceTagFromFilename:
    def test_egap_tag(self):
        assert source_tag_from_filename("/d/X_final_EGAP_assembly.fasta") == "egap"

    def test_polish_tag(self):
        assert source_tag_from_filename("/d/X_final_polish_assembly.fasta") == "polish"

    def test_no_tag(self):
        assert source_tag_from_filename("/d/plain_strain.fasta") == ""


class TestSuggestUniqueStrainIds:
    def test_egap_and_polish_do_not_collide(self):
        # The core bug: both collapse to NS27-2B-B3 → second was dropped on import.
        files = [
            "/d/NS27-2B-B3/NS27-2B-B3_final_EGAP_assembly.fasta",
            "/d/NS27-2B-B3/NS27-2B-B3_final_polish_assembly.fasta",
        ]
        ids = suggest_unique_strain_ids(files)
        assert len(set(ids)) == 2, f"IDs collided: {ids}"
        assert "NS27-2B-B3_egap" in ids
        assert "NS27-2B-B3_polish" in ids

    def test_unique_names_kept_clean(self):
        files = [
            "/d/NS26-3-C2_spades.fasta",
            "/d/NS26-3-C2_racon.fasta",
        ]
        ids = suggest_unique_strain_ids(files)
        # These already differ (suffixes not stripped) → kept verbatim, unique
        assert len(set(ids)) == 2
        assert ids == ["NS26-3-C2_spades", "NS26-3-C2_racon"]

    def test_all_ids_unique_even_with_triple_collision(self):
        files = [
            "/d/S1_final_EGAP_assembly.fasta",
            "/d/S1_final_polish_assembly.fasta",
            "/d/S1_best_assembly.fasta",
        ]
        ids = suggest_unique_strain_ids(files)
        assert len(set(ids)) == 3, f"Collision in {ids}"

    def test_count_preserved(self):
        files = [f"/d/sample_{i}.fasta" for i in range(5)]
        ids = suggest_unique_strain_ids(files)
        assert len(ids) == 5
        assert len(set(ids)) == 5


class TestParseQuastReport:
    def test_parses_metrics(self, tmp_path):
        rep = tmp_path / "report.tsv"
        rep.write_text(_QUAST_TSV)
        m = parse_quast_report(str(rep))
        assert m["assembly_name"] == "NS26-3-C2_final_EGAP_assembly"
        assert m["# contigs"] == 118          # int coercion
        assert m["N50"] == 1204567
        assert m["Total length"] == 34521900
        assert abs(m["GC (%)"] - 51.23) < 1e-6  # float coercion
        assert m["L50"] == 9

    def test_missing_file_returns_empty(self, tmp_path):
        assert parse_quast_report(str(tmp_path / "nope.tsv")) == {}


class TestFindQuastReport:
    def test_finds_egap_sibling_quast_dir(self, tmp_path):
        # .../NS26-3-C2_final_EGAP_assembly.fasta  →  <stem>_quast/report.tsv
        stem = "NS26-3-C2_final_EGAP_assembly"
        asm = tmp_path / f"{stem}.fasta"
        asm.write_text(">c1\nACGT\n")
        quast_dir = tmp_path / f"{stem}_quast"
        quast_dir.mkdir()
        (quast_dir / "report.tsv").write_text(_QUAST_TSV)

        found = find_quast_report(str(asm))
        assert found is not None
        assert Path(found).name == "report.tsv"
        assert parse_quast_report(found)["N50"] == 1204567

    def test_returns_none_when_absent(self, tmp_path):
        asm = tmp_path / "lonely_assembly.fasta"
        asm.write_text(">c1\nACGT\n")
        assert find_quast_report(str(asm)) is None
