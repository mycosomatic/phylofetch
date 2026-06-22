"""
pages/6_Workflow.py
-------------------
Workflow / Strategy orchestrator — component page (D-012 / RM-007 step 5).

Pick a named strategy and follow its steps. Each step links to its component page and shows
live status read from the project manifest (workflow.steps), which the extraction pages update
as they run. This is the "chain" that ties the standalone component pages into a coherent,
resumable pipeline without hiding state in st.session_state.
"""

import sys
from pathlib import Path

import streamlit as st

sys.path.insert(0, str(Path(__file__).parent.parent / "src"))
from phylofetch.config import load_config
from phylofetch.project_manager import (
    DEFAULT_PROJECT_DIR,
    load_project_manifest,
    project_output_dir,
    set_workflow_strategy,
)

st.set_page_config(page_title="Workflow", page_icon="🧭", layout="wide")
st.title("🧭 Workflow / Strategy")
st.caption("Pick a strategy and work down its steps. Each links to its component page and shows "
           "live status from the project manifest — the chain that ties the pipeline together.")

cfg = load_config()
project_dir = cfg.get("project_dir", str(DEFAULT_PROJECT_DIR))
manifest = load_project_manifest(project_dir)
wf = manifest["workflow"]
default_taxon = manifest.get("default_taxon", "")

st.caption(f"📁 Project: `{Path(project_dir).name}`"
           + (f"  ·  default taxon: **{default_taxon}**" if default_taxon
              else "  ·  ⚠️ no project taxon set (Assembly Manager)"))

# Named strategies = ordered recipes over the manifest's workflow steps.
STRATEGIES = {
    "Fungal barcodes (ITSx + Exonerate)": ["references", "rDNA", "coding"],
    "Primers only (in-silico PCR)": ["primers"],
    "Everything (rDNA + coding + primers)": ["references", "rDNA", "coding", "primers"],
}
STEP_INFO = {
    "references": ("📚 NCBI References", "pages/2_NCBI_References.py",
                   "Fetch per-locus references for the project taxon."),
    "rDNA": ("🧬 ITSx (rDNA)", "pages/3_ITSx_rDNA.py",
             "Extract ITS / LSU / SSU from assemblies."),
    "coding": ("🧬 Exonerate (coding)", "pages/4_Exonerate.py",
               "Frame-safe CDS for protein-coding loci (or relaxed BLAST amplicon)."),
    "primers": ("🔬 Primers", "pages/5_Primers.py",
                "In-silico PCR amplicons (no references needed)."),
}
BADGE = {"pending": "⬜ pending", "running": "🟡 running", "done": "✅ done",
         "error": "❌ error", "skipped": "➖ skipped"}


def _link(page: str, label: str) -> None:
    """st.page_link, with a graceful caption fallback when the page registry isn't available
    (e.g. headless AppTest bare mode), so the orchestrator never crashes on navigation."""
    try:
        st.page_link(page, label=label)
    except Exception:
        st.caption(f"→ {label}  ·  `{page}`")

# ── Strategy ──────────────────────────────────────────────────────────────────
names = list(STRATEGIES)
current = wf.get("strategy")
idx = names.index(current) if current in names else 0
choice = st.radio("Strategy", names, index=idx, horizontal=True)
if choice != current:
    set_workflow_strategy(project_dir, choice)
    st.rerun()

steps = STRATEGIES[choice]
st.markdown("---")

# ── Checklist ─────────────────────────────────────────────────────────────────
done = 0
for step in steps:
    title, page, desc = STEP_INFO[step]
    s = wf["steps"].get(step, {})
    status = s.get("status", "pending")
    done += status == "done"
    c1, c2, c3 = st.columns([3, 2, 2])
    with c1:
        st.markdown(f"**{title}**")
        st.caption(desc)
    with c2:
        st.write(BADGE.get(status, status))
        if s.get("updated_at"):
            st.caption(f"updated {s['updated_at']}")
        if s.get("notes"):
            st.caption(s["notes"])
    with c3:
        _link(page, "Open →")

st.markdown("---")
st.progress(done / len(steps) if steps else 0.0, text=f"{done}/{len(steps)} steps done")

# ── Outputs + handoff ─────────────────────────────────────────────────────────
combined_dir = project_output_dir(project_dir) / "loci" / "combined"
combos = sorted(combined_dir.glob("*_combined.fasta")) if combined_dir.exists() else []
st.markdown("**Combined outputs**")
if combos:
    for p in combos:
        st.write(f"✅ `{p.name}`")
    _link("pages/7_Alignment_Prep.py", "Continue to Alignment Prep →")
else:
    st.caption("No combined FASTAs yet — run the extraction steps above. They appear in "
               f"`{str(combined_dir).replace(str(Path.home()), '~')}`.")
