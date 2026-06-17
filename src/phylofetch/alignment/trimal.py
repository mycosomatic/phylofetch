"""trimAl alignment trimming wrapper."""

from __future__ import annotations

import shutil
import subprocess
from pathlib import Path
from typing import Optional, TYPE_CHECKING

if TYPE_CHECKING:
    from phylofetch.project_manager import RunManager


def run_trimal(
    input_fasta: str,
    output_fasta: str,
    *,
    mode: str = "automated1",
    gt: Optional[float] = None,
    cons: Optional[float] = None,
    trimal_bin: str = "trimal",
    run_manager: Optional["RunManager"] = None,
) -> tuple[int, str]:
    """
    Trim a multiple sequence alignment with trimAl.

    mode options:
      'automated1'  — heuristic (default; recommended for phylogenomics)
      'gappyout'    — removes gappy positions
      'strict'      — strict heuristic
      'manual'      — use gt/cons thresholds directly

    Returns (returncode, stderr_or_error).
    """
    if not shutil.which(trimal_bin):
        return 1, f"ERROR: {trimal_bin} not found on PATH"

    Path(output_fasta).parent.mkdir(parents=True, exist_ok=True)

    cmd = [trimal_bin, "-in", input_fasta, "-out", output_fasta, "-fasta"]

    if mode == "manual":
        if gt is not None:
            cmd += ["-gt", str(gt)]
        if cons is not None:
            cmd += ["-cons", str(cons)]
    else:
        cmd += [f"-{mode}"]

    if run_manager is not None:
        result = run_manager.run(
            cmd,
            module="alignment", action=f"trimal_{Path(input_fasta).stem}",
            inputs={"input_fasta": input_fasta},
            outputs={"output_fasta": output_fasta},
            params={"mode": mode, "gt": gt, "cons": cons},
            tool_version_keys=[trimal_bin],
        )
        return result.returncode, Path(result.stderr_path).read_text()
    else:
        r = subprocess.run(cmd, capture_output=True, text=True)
        return r.returncode, r.stderr
