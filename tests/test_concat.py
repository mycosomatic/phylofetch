"""
Tests for the supermatrix concatenator, focused on the D-032 aligned-length codon
partitioning (the path that actually reaches IQ-TREE from the Alignment Prep page).
"""

import re
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from phylofetch.alignment.concat import concatenate_alignments


def _write_aln(path: Path, seqs: dict[str, str]) -> str:
    recs = [SeqRecord(Seq(s), id=t, description="") for t, s in seqs.items()]
    SeqIO.write(recs, str(path), "fasta")
    return str(path)


def _charsets(nex_path: Path) -> dict[str, str]:
    """Parse a nexus partition file → {charset_name: spec}."""
    out = {}
    for m in re.finditer(r"charset\s+(\S+)\s*=\s*([^;]+);", nex_path.read_text(), re.I):
        out[m.group(1)] = m.group(2).strip()
    return out


class TestAlignedCodonPartitions:
    def test_cds_locus_gets_three_codon_charsets_from_aligned_length(self, tmp_path):
        # Two loci: a CDS (codon) of aligned length 9 and an rDNA (non-codon) of length 4.
        cds = _write_aln(tmp_path / "TEF1_CDS_combined.fasta",
                         {"A": "ATGAAATTT", "B": "ATGAAACTT"})           # len 9
        rdna = _write_aln(tmp_path / "ITS_combined.fasta",
                          {"A": "ACGT", "B": "ACGT"})                    # len 4
        out_fa = tmp_path / "sm.fasta"
        out_nex = tmp_path / "sm.nex"
        stats = concatenate_alignments(
            aligned_fastas=[cds, rdna], output_fasta=out_fa, partition_file=out_nex,
            codon_loci=[True, False],
        )
        assert stats["total_sites"] == 13
        cs = _charsets(out_nex)
        # CDS → three codon-position charsets spanning columns 1..9 with step 3.
        assert cs["TEF1_CDS_combined_pos1"] == r"1-9\3"
        assert cs["TEF1_CDS_combined_pos2"] == r"2-9\3"
        assert cs["TEF1_CDS_combined_pos3"] == r"3-9\3"
        # rDNA → single block, offset past the CDS (cols 10..13), no codon split.
        assert cs["ITS_combined"] == "10-13"
        assert not any("ITS_combined_pos" in k for k in cs)

    def test_offsets_accumulate_across_two_cds_loci(self, tmp_path):
        a = _write_aln(tmp_path / "RPB2_CDS_combined.fasta", {"A": "ATGAAA"})   # len 6
        b = _write_aln(tmp_path / "TUB2_CDS_combined.fasta", {"A": "ATGAAATTT"})  # len 9
        out_nex = tmp_path / "sm.nex"
        concatenate_alignments(
            aligned_fastas=[a, b], output_fasta=tmp_path / "sm.fasta",
            partition_file=out_nex, codon_loci=[True, True],
        )
        cs = _charsets(out_nex)
        assert cs["RPB2_CDS_combined_pos1"] == r"1-6\3"
        # Second locus starts at column 7 (after the 6-bp first locus), ends at 15.
        assert cs["TUB2_CDS_combined_pos1"] == r"7-15\3"
        assert cs["TUB2_CDS_combined_pos2"] == r"8-15\3"
        assert cs["TUB2_CDS_combined_pos3"] == r"9-15\3"

    def test_uses_aligned_not_unaligned_length(self, tmp_path):
        # A gapped CDS alignment: aligned length 12, even though ungapped CDS would be 9.
        # The codon charset MUST span the aligned columns (1..12), not the unaligned 9.
        cds = _write_aln(tmp_path / "ACT_CDS_combined.fasta",
                         {"A": "ATG---AAATTT", "B": "ATGCCCAAATTT"})    # aligned len 12
        out_nex = tmp_path / "sm.nex"
        concatenate_alignments(
            aligned_fastas=[cds], output_fasta=tmp_path / "sm.fasta",
            partition_file=out_nex, codon_loci=[True],
        )
        cs = _charsets(out_nex)
        assert cs["ACT_CDS_combined_pos1"] == r"1-12\3"
        assert cs["ACT_CDS_combined_pos3"] == r"3-12\3"

    def test_no_codon_loci_falls_back_to_one_block_per_locus(self, tmp_path):
        a = _write_aln(tmp_path / "TEF1_CDS_combined.fasta", {"A": "ATGAAA"})
        out_nex = tmp_path / "sm.nex"
        concatenate_alignments(
            aligned_fastas=[a], output_fasta=tmp_path / "sm.fasta",
            partition_file=out_nex, codon_loci=None,
        )
        cs = _charsets(out_nex)
        assert cs == {"TEF1_CDS_combined": "1-6"}

    def test_mismatched_flag_length_is_ignored(self, tmp_path):
        # A wrong-length codon_loci must not corrupt the partition — fall back gracefully.
        a = _write_aln(tmp_path / "TEF1_CDS_combined.fasta", {"A": "ATGAAA"})
        b = _write_aln(tmp_path / "ITS_combined.fasta", {"A": "ACGT"})
        out_nex = tmp_path / "sm.nex"
        concatenate_alignments(
            aligned_fastas=[a, b], output_fasta=tmp_path / "sm.fasta",
            partition_file=out_nex, codon_loci=[True],   # len 1 != 2 loci → ignored
        )
        cs = _charsets(out_nex)
        assert "TEF1_CDS_combined" in cs and "ITS_combined" in cs
        assert not any("_pos1" in k for k in cs)
