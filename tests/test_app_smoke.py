"""
Smoke test for the grouped st.navigation entry point (app.py, D-031/D-032).

The navigation registers every page by a hardcoded on-disk path; a renamed/moved page would
break app startup with nothing else catching it. AppTest runs the entry script headlessly,
which validates that all st.Page paths resolve, the Home page renders, and every page module
imports + runs without raising.
"""

from pathlib import Path

import pytest

pytest.importorskip("streamlit")
from streamlit.testing.v1 import AppTest  # noqa: E402

ROOT = Path(__file__).resolve().parent.parent
APP = str(ROOT / "app.py")

PAGES = [
    "pages/0_Project_Setup.py", "pages/1_Assembly_Manager.py", "pages/2_NCBI_References.py",
    "pages/3_ITSx_rDNA.py", "pages/4_Exonerate.py", "pages/5_Primers.py",
    "pages/6_Workflow.py", "pages/7_Reference_Taxa.py", "pages/8_Codon_Tip_Prep.py",
    "pages/9_Alignment_Prep.py", "pages/10_BUSCO_Phylogenomics.py",
    "pages/11_Tree_Visualization.py",
]


def test_navigation_builds_and_home_renders():
    at = AppTest.from_file(APP, default_timeout=60).run()
    assert not at.exception, f"app.py raised: {at.exception}"
    assert any("phylofetch" in t.value for t in at.title)


@pytest.mark.parametrize("page", PAGES)
def test_each_page_loads(page):
    at = AppTest.from_file(APP, default_timeout=60)
    at.run()
    at.switch_page(page)
    at.run()
    assert not at.exception, f"{page} raised: {at.exception}"
