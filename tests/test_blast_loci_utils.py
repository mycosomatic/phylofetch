"""
Tests for src/phylofetch/blast_loci_utils.py

Includes regression test for LXD-002: HSPs from different reference
accessions must NOT be stitched into one CDS. The winning reference is
chosen by highest total bitscore from a single qseqid.
"""

import re
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).parent.parent / "src"))
from phylofetch.alignment.concat import _shift_partition_spec
from phylofetch.blast_loci_utils import (
    select_best_locus_group,
    write_codon_partition,
)


# ── LXD-002 regression ────────────────────────────────────────────────────

class TestSelectBestLocusGroup:
    """LXD-002: group by qseqid first; never stitch HSPs across references."""

    def test_single_ref_returns_all_hsps(self, single_ref_hsps):
        hsps, ref_acc = select_best_locus_group(single_ref_hsps)
        assert ref_acc == "ref1"
        assert len(hsps) == 2

    def test_two_refs_best_wins(self, two_ref_hsps):
        """ref2 has higher total bitscore → must win; ref1 HSPs excluded."""
        hsps, ref_acc = select_best_locus_group(two_ref_hsps)
        assert ref_acc == "ref2", f"Expected ref2 to win, got {ref_acc}"
        for h in hsps:
            assert h["qseqid"] == "ref2", (
                f"LXD-002 regression: HSP from '{h['qseqid']}' leaked into result"
            )

    def test_two_refs_higher_bitscore_wins(self):
        """ref1 at 600 total vs ref2 at 400 → ref1 must win."""
        hsps = [
            {"qseqid": "ref1", "sseqid": "ctg1", "qstart": 1, "qend": 300,
             "sstart": 100, "send": 400,
             "bitscore": 600.0, "pident": 95.0, "evalue": 1e-80, "length": 300, "qlen": 500},
            {"qseqid": "ref2", "sseqid": "ctg1", "qstart": 1, "qend": 300,
             "sstart": 100, "send": 400,
             "bitscore": 400.0, "pident": 88.0, "evalue": 1e-60, "length": 300, "qlen": 500},
        ]
        _, ref_acc = select_best_locus_group(hsps)
        assert ref_acc == "ref1"

    def test_empty_input_returns_empty_list(self):
        # Empty input: filtered list is empty → returns [], ""
        hsps, ref_acc = select_best_locus_group([])
        assert hsps == []
        assert ref_acc == ""

    def test_low_pident_filtered_out(self):
        # All HSPs below min_pident threshold → should return empty
        low_pident = [
            {"qseqid": "ref1", "sseqid": "ctg1", "sstart": 100, "send": 400,
             "bitscore": 200.0, "pident": 50.0, "evalue": 1e-20, "length": 300, "qlen": 500},
        ]
        hsps, ref_acc = select_best_locus_group(low_pident, min_pident=70.0)
        assert hsps == []

    def test_returns_tuple(self, single_ref_hsps):
        result = select_best_locus_group(single_ref_hsps)
        assert isinstance(result, tuple) and len(result) == 2

    def test_hsps_sorted_by_qstart(self, single_ref_hsps):
        """Returned HSPs must be ordered by query position (qstart ascending)."""
        hsps, _ = select_best_locus_group(single_ref_hsps)
        qstarts = [h["qstart"] for h in hsps]
        assert qstarts == sorted(qstarts), "HSPs not sorted by qstart"

    def test_strand_consistency(self):
        """HSPs on opposite strands of same contig must NOT be grouped together."""
        # Plus strand (sstart < send) and minus strand (sstart > send)
        mixed = [
            {"qseqid": "ref1", "sseqid": "ctg1", "qstart": 1, "qend": 250,
             "sstart": 1000, "send": 1250, "bitscore": 200.0, "pident": 95.0,
             "evalue": 1e-50, "length": 250, "qlen": 500},
            {"qseqid": "ref1", "sseqid": "ctg1", "qstart": 251, "qend": 500,
             "sstart": 1600, "send": 1350,  # minus strand (sstart > send)
             "bitscore": 200.0, "pident": 95.0, "evalue": 1e-50, "length": 250, "qlen": 500},
        ]
        hsps, _ = select_best_locus_group(mixed)
        if len(hsps) > 0:
            # All returned HSPs must be on the same strand
            first_is_plus = hsps[0]["sstart"] <= hsps[0]["send"]
            for h in hsps:
                h_is_plus = h["sstart"] <= h["send"]
                assert h_is_plus == first_is_plus, "Mixed strands in output"


# ── Codon partition writing ────────────────────────────────────────────────

class TestWriteCodonPartition:
    def test_basic_partition(self, tmp_path):
        out_path = write_codon_partition(
            cds_length=750,
            output_dir=str(tmp_path),
            locus_name="TEF1",
        )
        assert Path(out_path).exists()
        text = Path(out_path).read_text()
        assert "#nexus" in text.lower()
        assert "TEF1" in text
        assert "charset" in text.lower()

    def test_three_codon_positions(self, tmp_path):
        out_path = write_codon_partition(
            cds_length=300,
            output_dir=str(tmp_path),
            locus_name="RPB2",
        )
        text = Path(out_path).read_text()
        charset_count = text.lower().count("charset")
        assert charset_count == 3, f"Expected 3 charsets, got {charset_count}"

    def test_output_file_name(self, tmp_path):
        out_path = write_codon_partition(
            cds_length=600,
            output_dir=str(tmp_path),
            locus_name="ITS",
        )
        assert Path(out_path).name == "ITS_partition.nex"

    def test_positions_in_range(self, tmp_path):
        cds_length = 750
        out_path = write_codon_partition(
            cds_length=cds_length,
            output_dir=str(tmp_path),
            locus_name="ACT",
        )
        text = Path(out_path).read_text()
        # Only check positions on charset lines (skip comments with IQ-TREE flags)
        charset_positions = []
        for line in text.splitlines():
            stripped = line.strip().lower()
            if stripped.startswith("charset"):
                charset_positions.extend(
                    int(m) for m in re.findall(r"\b(\d+)\b", line)
                )
        non_step = [p for p in charset_positions if p > 3]  # exclude step sizes (1,2,3)
        if non_step:
            assert max(non_step) <= cds_length, (
                f"Position {max(non_step)} exceeds CDS length {cds_length}"
            )


# ── Partition spec shifting ────────────────────────────────────────────────

class TestShiftPartitionSpec:
    def test_simple_range(self):
        result = _shift_partition_spec("1-750", 1000)
        assert result == "1001-1750"

    def test_step_range_backslash(self):
        # Partition step notation: 1-750\3 → shift positions, preserve \3 step value
        result = _shift_partition_spec("1-750\\3", 1000)
        assert result == "1001-1750\\3", (
            f"Step value should not be shifted. Got: {result!r}"
        )

    def test_multiple_ranges(self):
        result = _shift_partition_spec("1-250, 251-500", 500)
        assert result == "501-750, 751-1000"

    def test_zero_offset(self):
        spec = "10-300"
        assert _shift_partition_spec(spec, 0) == spec

    def test_preserves_step_notation(self):
        spec = "1-300\\3"
        result = _shift_partition_spec(spec, 100)
        # Position numbers shift but step value (\3) is unchanged
        assert result == "101-400\\3", (
            f"Expected '101-400\\\\3', got {result!r}"
        )
