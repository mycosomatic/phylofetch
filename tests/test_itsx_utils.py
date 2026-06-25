"""
Tests for src/phylofetch/itsx_utils.py

Covers LXD-001 fixes:
- No --multi_out flag in command
- Stale output files are cleared before rerun
- Log is always populated (last N lines returned)
"""

import sys
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

sys.path.insert(0, str(Path(__file__).parent.parent / "src"))
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from phylofetch.itsx_utils import (
    ITSX_SUFFIXES,
    _probe_itsx_version,
    _relabel_itsx_output,
    chunk_long_contigs,
    combine_rdna_regions,
    parse_coverage,
    place_rdna_regions,
    run_itsx,
)


class TestCoverageFilter:
    """D-028: prefer the high-coverage rDNA array; drop low-coverage off-array (RIP'd/spurious)."""

    def test_parse_coverage(self):
        assert parse_coverage("NODE_1_length_6920269_cov_45.858574_pilon") == 45.858574
        assert parse_coverage("NODE_12_length_22758_cov_1193.220177_pilon") == 1193.220177
        assert parse_coverage("scaffold_no_cov_token") is None
        assert parse_coverage("") is None

    def _write(self, path, recs):
        SeqIO.write([SeqRecord(Seq(s), id=rid, description=rid) for rid, s in recs], str(path), "fasta")

    def test_drops_low_coverage_offarray_copy(self, tmp_path):
        # The real NS26 SSU case: a 60 kb spurious hit on the 45x chromosome + the true 2821 bp
        # region on the 1193x rDNA array. The chromosomal copy must be dropped.
        f = tmp_path / "NS26.SSU.fasta"
        self._write(f, [
            ("NODE_1_length_6920269_cov_45.858574_pilon|1-60084|", "A" * 60084),
            ("NODE_12_length_22758_cov_1193.220177_pilon|1-2821|", "C" * 2821),
        ])
        info = _relabel_itsx_output(str(f), "NS26-3-C2", "SSU")
        assert info["kept"] == 1
        assert len(info["dropped"]) == 1 and info["dropped"][0][1] == 45.858574
        kept = list(SeqIO.parse(str(f), "fasta"))
        assert len(kept) == 1 and len(kept[0].seq) == 2821       # the array copy, not the 60 kb one
        assert "[cov=1193]" in kept[0].description

    def test_single_detection_kept_regardless(self, tmp_path):
        f = tmp_path / "x.SSU.fasta"
        self._write(f, [("NODE_5_length_9000_cov_30.0_pilon|1-1800|", "G" * 1800)])
        info = _relabel_itsx_output(str(f), "S1", "SSU")
        assert info["kept"] == 1 and info["dropped"] == []        # nothing to compare against

    def test_no_coverage_tokens_keeps_all(self, tmp_path):
        f = tmp_path / "y.SSU.fasta"
        self._write(f, [("contigA|1-1800|", "G" * 1800), ("contigB|1-1700|", "T" * 1700)])
        info = _relabel_itsx_output(str(f), "S1", "SSU")
        assert info["kept"] == 2 and info["dropped"] == []        # un-parseable cov → no filtering

    def test_filter_can_be_disabled(self, tmp_path):
        f = tmp_path / "z.SSU.fasta"
        self._write(f, [
            ("NODE_1_length_9_cov_45.0_pilon|1-60000|", "A" * 60000),
            ("NODE_2_length_9_cov_1200.0_pilon|1-2800|", "C" * 2800),
        ])
        info = _relabel_itsx_output(str(f), "S1", "SSU", prefer_high_cov=False)
        assert info["kept"] == 2 and info["dropped"] == []


class TestChunkLongContigs:
    """LXD-003: contigs > hmmscan's 100 kb limit are split into overlapping windows."""

    def _rec(self, rid, n):
        return SeqRecord(Seq("A" * n), id=rid, description="")

    def test_short_contig_passthrough(self):
        out, changed = chunk_long_contigs([self._rec("c1", 1000)])
        assert changed is False
        assert len(out) == 1 and out[0].id == "c1"

    def test_long_contig_split_with_overlap(self):
        out, changed = chunk_long_contigs([self._rec("big", 250_000)],
                                          max_len=90_000, overlap=20_000)
        assert changed is True
        assert all(len(r.seq) <= 90_000 for r in out)          # under hmmscan limit
        starts = [int(r.id.split("__c")[1]) for r in out]
        assert starts == [0, 70_000, 140_000, 210_000]         # step = max_len - overlap
        assert out[-1].id == "big__c210000"
        assert len(out[-1].seq) == 250_000 - 210_000            # last chunk reaches the end

    def test_overlap_keeps_small_features_intact(self):
        # any feature shorter than the overlap is fully contained in at least one chunk
        out, _ = chunk_long_contigs([self._rec("c", 200_000)], max_len=90_000, overlap=20_000)
        # consecutive chunks overlap by exactly `overlap`
        spans = [(int(r.id.split("__c")[1]),
                  int(r.id.split("__c")[1]) + len(r.seq)) for r in out]
        for (s0, e0), (s1, e1) in zip(spans, spans[1:]):
            assert e0 - s1 >= 20_000

    def test_mixed_short_and_long(self):
        out, changed = chunk_long_contigs([self._rec("s", 100), self._rec("L", 120_000)],
                                          max_len=90_000, overlap=20_000)
        assert changed is True
        ids = [r.id for r in out]
        assert "s" in ids and any(i.startswith("L__c") for i in ids)

# run_itsx signature: (assembly_fasta, output_dir, strain_id, threads=4,
#                      itsx_bin="ITSx", kingdom="fungi")
# Returns: (returncode, log_text, found_dict)


def _make_mock_run_result(returncode=0, stdout="ITSx done.", stderr=""):
    mock_result = MagicMock()
    mock_result.returncode = returncode
    mock_result.stdout = stdout
    mock_result.stderr = stderr
    return mock_result


class TestITSxCommand:
    """LXD-001: verify ITSx command construction."""

    @patch("phylofetch.itsx_utils.shutil.which", return_value="/usr/bin/ITSx")
    @patch("phylofetch.itsx_utils.subprocess.run")
    def test_no_multi_out_flag(self, mock_run, mock_which, tmp_path):
        """LXD-001: --multi_out must NOT appear in the ITSx command."""
        fasta = tmp_path / "assembly.fasta"
        fasta.write_text(">contig1\nATGCATGC\n")
        mock_run.return_value = _make_mock_run_result()

        run_itsx(str(fasta), str(tmp_path), "strain1", itsx_bin="ITSx")

        cmd = mock_run.call_args[0][0]
        flat_cmd = " ".join(str(c) for c in cmd)
        assert "--multi_out" not in flat_cmd, (
            f"LXD-001 regression: --multi_out found in command: {flat_cmd}"
        )

    @patch("phylofetch.itsx_utils.shutil.which", return_value="/usr/bin/ITSx")
    @patch("phylofetch.itsx_utils.subprocess.run")
    def test_kingdom_in_command(self, mock_run, mock_which, tmp_path):
        """ITSx should pass kingdom via -t flag."""
        fasta = tmp_path / "assembly.fasta"
        fasta.write_text(">c1\nATGC\n")
        mock_run.return_value = _make_mock_run_result()

        run_itsx(str(fasta), str(tmp_path), "s1", itsx_bin="ITSx", kingdom="fungi")

        cmd = mock_run.call_args[0][0]
        flat_cmd = " ".join(str(c) for c in cmd)
        assert "fungi" in flat_cmd.lower(), (
            f"Kingdom 'fungi' not found in ITSx command: {flat_cmd}"
        )

    @patch("phylofetch.itsx_utils.shutil.which", return_value="/usr/bin/ITSx")
    @patch("phylofetch.itsx_utils.subprocess.run")
    def test_output_dir_in_command(self, mock_run, mock_which, tmp_path):
        """Output prefix (output_dir/strain_id) must be in the command."""
        fasta = tmp_path / "assembly.fasta"
        fasta.write_text(">c1\nATGC\n")
        mock_run.return_value = _make_mock_run_result()

        run_itsx(str(fasta), str(tmp_path), "mystrain", itsx_bin="ITSx")

        cmd = mock_run.call_args[0][0]
        flat_cmd = " ".join(str(c) for c in cmd)
        assert "mystrain" in flat_cmd, (
            f"strain_id 'mystrain' not found as part of prefix in ITSx command: {flat_cmd}"
        )


class TestITSxStaleCleanup:
    """LXD-001: stale output files must be removed before rerun."""

    @patch("phylofetch.itsx_utils.shutil.which", return_value="/usr/bin/ITSx")
    @patch("phylofetch.itsx_utils.subprocess.run")
    def test_stale_files_removed(self, mock_run, mock_which, tmp_path):
        """Old output files matching prefix+suffix must be deleted before ITSx runs."""
        strain_id = "mystrain"
        prefix = str(tmp_path / strain_id)
        fasta = tmp_path / "assembly.fasta"
        fasta.write_text(">c1\nATGC\n")

        # Create stale output files using actual ITSX_SUFFIXES
        stale_files = []
        for suffix in list(ITSX_SUFFIXES.values())[:3]:
            stale = Path(prefix + suffix)
            stale.write_text("stale content")
            stale_files.append(stale)

        mock_run.return_value = _make_mock_run_result()

        run_itsx(str(fasta), str(tmp_path), strain_id, itsx_bin="ITSx")

        for sf in stale_files:
            assert not sf.exists(), (
                f"LXD-001 regression: stale file not removed before rerun: {sf}"
            )


class TestITSxLogReturn:
    """LXD-001: log_text must always be populated."""

    @patch("phylofetch.itsx_utils.shutil.which", return_value="/usr/bin/ITSx")
    @patch("phylofetch.itsx_utils.subprocess.run")
    def test_log_populated_on_success(self, mock_run, mock_which, tmp_path):
        fasta = tmp_path / "assembly.fasta"
        fasta.write_text(">c1\nATGC\n")
        mock_run.return_value = _make_mock_run_result(stdout="ITSx 1.1.3\nDone processing.")

        rc, log_text, found_dict = run_itsx(str(fasta), str(tmp_path), "s1")
        assert isinstance(log_text, str)
        assert len(log_text) > 0, "log_text is empty on success"

    @patch("phylofetch.itsx_utils.shutil.which", return_value="/usr/bin/ITSx")
    @patch("phylofetch.itsx_utils.subprocess.run")
    def test_log_populated_on_failure(self, mock_run, mock_which, tmp_path):
        fasta = tmp_path / "assembly.fasta"
        fasta.write_text(">c1\nATGC\n")
        mock_run.return_value = _make_mock_run_result(
            returncode=1,
            stdout="\n".join(f"line {i}" for i in range(50)),
            stderr="Fatal error occurred",
        )

        rc, log_text, found_dict = run_itsx(str(fasta), str(tmp_path), "s1")
        assert rc == 1
        assert len(log_text) > 0, "log_text is empty on failure"

    @patch("phylofetch.itsx_utils.shutil.which", return_value="/usr/bin/ITSx")
    @patch("phylofetch.itsx_utils.subprocess.run")
    def test_return_tuple_structure(self, mock_run, mock_which, tmp_path):
        fasta = tmp_path / "assembly.fasta"
        fasta.write_text(">c1\nATGC\n")
        mock_run.return_value = _make_mock_run_result(stdout="done")

        result = run_itsx(str(fasta), str(tmp_path), "s1")
        assert isinstance(result, tuple) and len(result) == 3, (
            "run_itsx must return (returncode, log_text, found_dict)"
        )
        rc, log_text, found_dict = result
        assert isinstance(rc, int)
        assert isinstance(log_text, str)
        assert isinstance(found_dict, dict)

    def test_missing_binary_returns_error(self, tmp_path):
        """If ITSx binary is not on PATH, run_itsx returns rc=1 immediately."""
        fasta = tmp_path / "assembly.fasta"
        fasta.write_text(">c1\nATGC\n")

        with patch("phylofetch.itsx_utils.shutil.which", return_value=None):
            rc, log_text, found_dict = run_itsx(
                str(fasta), str(tmp_path), "s1", itsx_bin="ITSx"
            )

        assert rc == 1
        assert "ITSx" in log_text or "not found" in log_text.lower()
        assert found_dict == {}


class TestProbeITSxVersion:
    @patch("phylofetch.itsx_utils.subprocess.run")
    def test_version_returns_string(self, mock_run):
        """_probe_itsx_version always returns a string (may be 'unknown')."""
        mock_result = MagicMock()
        mock_result.returncode = 0
        mock_result.stdout = "ITSx -- Identifies ITS sequences\nVersion 1.1.3"
        mock_result.stderr = ""
        mock_run.return_value = mock_result

        version = _probe_itsx_version("ITSx")
        assert isinstance(version, str)
        assert len(version) > 0

    @patch("phylofetch.itsx_utils.subprocess.run")
    def test_version_returns_string_on_error(self, mock_run):
        """On error, _probe_itsx_version returns 'unknown' (not None)."""
        mock_run.side_effect = FileNotFoundError("ITSx not found")
        version = _probe_itsx_version("ITSx")
        assert isinstance(version, str)
        assert version == "unknown"

    @patch("phylofetch.itsx_utils.subprocess.run")
    def test_version_contains_digit_on_success(self, mock_run):
        """If ITSx responds with a version line, the result should contain digits."""
        mock_result = MagicMock()
        mock_result.returncode = 0
        mock_result.stdout = "ITSx 1.1.3 info line"
        mock_result.stderr = ""
        mock_run.return_value = mock_result

        version = _probe_itsx_version("ITSx")
        # If a version line is parsed, it should contain digits
        if version != "unknown":
            assert any(c.isdigit() for c in version)


class TestRdnaPlacementAndCombine:
    """RM-007 step 4c: per-strain region placement + cross-strain combine."""

    def _fasta(self, path: Path, records):
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text("".join(f">{rid}\n{seq}\n" for rid, seq in records))

    def test_place_maps_full_its_and_counts(self, tmp_path):
        src = tmp_path / "itsx"
        self._fasta(src / "s.full.fasta", [("c1", "ACGT"), ("c2", "TTTT")])
        self._fasta(src / "s.LSU.fasta", [("c1", "GGGG")])
        found = {"ITS_full": str(src / "s.full.fasta"), "LSU": str(src / "s.LSU.fasta")}
        strain_out = tmp_path / "per_strain" / "s"

        counts = place_rdna_regions(found, str(strain_out), ["ITS", "LSU", "SSU"])
        assert counts == {"ITS": 2, "LSU": 1, "SSU": 0}
        assert (strain_out / "ITS" / "ITS.fasta").exists()      # ITS_full -> ITS.fasta
        assert (strain_out / "LSU" / "LSU.fasta").exists()
        assert not (strain_out / "SSU" / "SSU.fasta").exists()  # absent region not written

    def test_combine_merges_present_regions_only(self, tmp_path):
        ps = tmp_path / "per_strain"
        for strain, seq in [("A", "ACGT"), ("B", "TTTT")]:
            self._fasta(ps / strain / "ITS" / "ITS.fasta", [(f"{strain}_ITS", seq)])
        self._fasta(ps / "A" / "SSU" / "SSU.fasta", [("A_SSU", "GG")])

        out = combine_rdna_regions(str(ps), str(tmp_path / "combined"), ["ITS", "SSU", "LSU"])
        assert set(out) == {"ITS", "SSU"}            # LSU had no sequences -> omitted
        its_path, its_n = out["ITS"]
        assert its_n == 2 and its_path.endswith("ITS_combined.fasta")
        assert out["SSU"][1] == 1
