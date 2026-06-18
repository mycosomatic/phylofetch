"""Entry point for the `phylofetch` command."""

import subprocess
import sys
from pathlib import Path


def _find_app() -> Path:
    """
    Locate app.py by walking up from this file.

    Works for both editable installs (src/phylofetch/cli.py → repo root) and
    any other layout where app.py lives 1–4 directories above the package.
    """
    here = Path(__file__).resolve().parent
    for candidate in [here, here.parent, here.parent.parent, here.parent.parent.parent]:
        app = candidate / "app.py"
        if app.exists():
            return app
    raise FileNotFoundError(
        "Cannot find app.py. "
        "If installed from source run: streamlit run app.py"
    )


def main() -> None:
    try:
        app = _find_app()
    except FileNotFoundError as exc:
        print(str(exc), file=sys.stderr)
        sys.exit(1)
    sys.exit(subprocess.run(["streamlit", "run", str(app)] + sys.argv[1:]).returncode)
