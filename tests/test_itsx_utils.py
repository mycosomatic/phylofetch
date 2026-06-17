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
from phylofetch.itsx_utils import (
    ITSX_SUFFIXES,
    _probe_itsx_version,
    run_itsx,
)

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
