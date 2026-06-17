"""MAFFT multiple sequence alignment wrapper."""

from __future__ import annotations

import shutil
import subprocess
from pathlib import Path
from typing import Optional, TYPE_CHECKING

if TYPE_CHECKING:
    from phylofetch.project_manager import RunManager, RunResult


def run_mafft(
    input_fasta: str,
    output_fasta: str,
    *,
    mode: str = "auto",
    threads: int = 4,
    mafft_bin: str = "mafft",
    extra_args: list[str] | None = None,
    run_manager: Optional["RunManager"] = None,
) -> tuple[int, str]:
    """
    Align sequences with MAFFT.

    mode options:
      'auto'       — mafft --auto (recommended default)
      'localpair'  — L-INS-i (accurate, slower; good for short conserved regions)
      'globalpair' — G-INS-i (for globally alignable sequences)
      'genafpair'  — E-INS-i (for sequences with unalignable regions)
      'retree 2'   — FFT-NS-2 (fast, for many sequences)

    Returns (returncode, alignment_text_or_error).
    """
    if not shutil.which(mafft_bin):
        return 1, f"ERROR: {mafft_bin} not found on PATH"

    Path(output_fasta).parent.mkdir(parents=True, exist_ok=True)

    mode_flags: dict[str, list[str]] = {
        "auto":       ["--auto"],
        "localpair":  ["--localpair", "--maxiterate", "1000"],
        "globalpair": ["--globalpair", "--maxiterate", "1000"],
        "genafpair":  ["--genafpair", "--maxiterate", "1000"],
        "retree 2":   ["--retree", "2"],
    }
    flags = mode_flags.get(mode, ["--auto"])
    cmd = [mafft_bin] + flags + ["--thread", str(threads)]
    if extra_args:
        cmd += extra_args
    cmd += [input_fasta]

    if run_manager is not None:
        result = run_manager.run(
            cmd,
            module="alignment", action=f"mafft_{Path(input_fasta).stem}",
            inputs={"input_fasta": input_fasta},
            outputs={"output_fasta": output_fasta},
            params={"mode": mode, "threads": threads},
            tool_version_keys=[mafft_bin],
        )
        # MAFFT writes alignment to stdout; RunManager captures it in stdout.log
        if result.returncode == 0:
            Path(output_fasta).write_text(Path(result.stdout_path).read_text())
        return result.returncode, Path(result.stderr_path).read_text()
    else:
        r = subprocess.run(cmd, capture_output=True, text=True)
        if r.returncode == 0:
            Path(output_fasta).write_text(r.stdout)
        return r.returncode, r.stderr
