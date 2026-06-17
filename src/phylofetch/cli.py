"""Entry point for the `phylofetch` command."""

import subprocess
import sys
from pathlib import Path


def main() -> None:
    app = Path(__file__).parent.parent.parent / "app.py"
    sys.exit(subprocess.run(["streamlit", "run", str(app)] + sys.argv[1:]).returncode)
