"""
pages/2_ITSx_rDNA.py
--------------------
ITSx (rDNA) — component page (D-012 / RM-007 step 4c).

Extract rDNA regions (ITS / ITS1 / ITS2 / LSU / SSU) from this project's assemblies with ITSx.
Standalone (pick assemblies + regions, run, see results) and chainable (per-project outputs +
manifest step state). Outputs land under <project>/results/loci (D-015): per-strain region
FASTAs and cross-strain *_combined.fasta files that feed Alignment Prep.

Interim during the page decomposition: runs alongside the old 2_Loci_Extraction page until the
monolith is retired.
"""

import os
import shutil
import sys
from pathlib import Path

import pandas as pd
import streamlit as st

sys.path.insert(0, str(Path(__file__).parent.parent / "src"))
from phylofetch.config import load_config
from phylofetch.itsx_utils import combine_rdna_regions, place_rdna_regions, run_itsx
from phylofetch.project_manager import (
    DEFAULT_PROJECT_DIR,
    effective_taxon,
    load_assembly_registry,
    load_project_manifest,
    project_output_dir,
    update_step,
)

st.set_page_config(page_title="ITSx (rDNA)", page_icon="🧬", layout="wide")
st.title("🧬 ITSx — rDNA extraction")
st.caption("Extract rDNA regions (ITS / ITS1 / ITS2 / LSU / SSU) from this project's assemblies "
           "with ITSx. Outputs are per-project; the combined FASTAs feed Alignment Prep.")

cfg = load_config()
project_dir = cfg.get("project_dir", str(DEFAULT_PROJECT_DIR))
tp = cfg.get("tool_paths", {})
registry = load_assembly_registry(project_dir)
manifest = load_project_manifest(project_dir)
default_taxon = manifest.get("default_taxon", "")

loci_dir = project_output_dir(project_dir) / "loci"
per_strain_dir = loci_dir / "per_strain"
combined_dir = loci_dir / "combined"

st.caption(f"📁 Project: `{Path(project_dir).name}`  ·  outputs: "
           f"`{str(loci_dir).replace(str(Path.home()), '~')}`")

if not registry:
    st.info("No assemblies registered yet — add them in the Assembly Manager.")
    st.stop()

# ── 1 · Assemblies ────────────────────────────────────────────────────────────
st.subheader("1 · Assemblies")


def _fmt(sid: str) -> str:
    t = effective_taxon(registry[sid], default_taxon)
    return f"{sid}  ·  {t}" if t else sid


sel = st.multiselect("Assemblies to run ITSx on", list(registry),
                     default=list(registry), format_func=_fmt)

# ── 2 · Regions ───────────────────────────────────────────────────────────────
st.subheader("2 · rDNA regions")
REGIONS = ["ITS", "ITS1", "ITS2", "LSU", "SSU"]
DEFAULTS = {"ITS", "LSU", "SSU"}
LABELS = {"ITS": "ITS (full)", "ITS1": "ITS1", "ITS2": "ITS2", "LSU": "LSU", "SSU": "SSU"}
regions = []
for c, r in zip(st.columns(len(REGIONS)), REGIONS):
    with c:
        if st.checkbox(LABELS[r], value=(r in DEFAULTS), key=f"reg_{r}"):
            regions.append(r)

# ── 3 · Settings ──────────────────────────────────────────────────────────────
st.subheader("3 · Settings")
s1, s2, s3 = st.columns(3)
with s1:
    itsx_bin = st.text_input("ITSx binary", value=tp.get("itsx", "ITSx"))
with s2:
    kingdom = st.selectbox(
        "ITSx kingdom",
        ["fungi", "all", "metazoa", "viridiplantae", "tracheophyta", "bacteria"],
        index=0, help="ITSx HMM profile set (-t). 'fungi' for fungal genomes.",
    )
with s3:
    threads = st.number_input("Threads", 1, 64, 4)
    prefer_high_cov = st.checkbox(
        "Prefer high-coverage rDNA array", value=True,
        help="Keep only the rDNA detection on the highest-coverage contig (the functional, "
             "high-copy tandem array). Low-coverage off-array copies — dispersed/orphan rDNA that "
             "may be RIP-pseudogenized, or spurious HMM hits on chromosomal contigs (the source of "
             "giant 60 kb 'SSU' extractions) — are dropped. (D-028)")

itsx_ok = shutil.which(itsx_bin) is not None
if not itsx_ok:
    st.warning(f"⚠️ ITSx not found on PATH (`{itsx_bin}`). "
               "Install with `conda install -c bioconda itsx`.")

# ── 4 · Run ───────────────────────────────────────────────────────────────────
st.subheader("4 · Run")
if not sel:
    st.info("Select at least one assembly.")
elif not regions:
    st.info("Select at least one rDNA region.")

if st.button("🚀 Run ITSx", type="primary", disabled=not (sel and regions and itsx_ok)):
    per_strain_dir.mkdir(parents=True, exist_ok=True)
    combined_dir.mkdir(parents=True, exist_ok=True)
    summary: dict[str, dict] = {}
    prog = st.progress(0.0)
    for i, sid in enumerate(sel):
        assembly = registry[sid].get("assembly_path", "")
        strain_out = per_strain_dir / sid
        if not assembly or not os.path.exists(assembly):
            st.error(f"{sid}: assembly file not found ({assembly})")
            summary[sid] = {r: 0 for r in regions}
            prog.progress((i + 1) / len(sel))
            continue
        with st.spinner(f"[{sid}] ITSx…"):
            rc, log, found = run_itsx(assembly, str(strain_out / "_itsx_tmp"), sid,
                                      threads=int(threads), itsx_bin=itsx_bin, kingdom=kingdom,
                                      prefer_high_cov=prefer_high_cov)
        if rc != 0:
            st.error(f"{sid}: ITSx failed (exit {rc})")
            with st.expander(f"{sid} — ITSx log"):
                st.code(log, language=None)
            summary[sid] = {r: 0 for r in regions}
        else:
            summary[sid] = place_rdna_regions(found, str(strain_out), regions)
            if not found:
                st.warning(f"{sid}: no rDNA regions detected.")
        prog.progress((i + 1) / len(sel))

    combined = combine_rdna_regions(str(per_strain_dir), str(combined_dir), regions)
    update_step(
        project_dir, "rDNA", status="done",
        outputs={r: {"combined": p, "n": n} for r, (p, n) in combined.items()},
        notes=f"strains={len(sel)}; regions={','.join(regions)}; kingdom={kingdom}",
    )
    st.session_state["itsx_summary"] = {"per_strain": summary, "combined": combined}
    st.success("ITSx complete — provenance recorded in the project manifest.")
    st.rerun()

# ── Results ───────────────────────────────────────────────────────────────────
res = st.session_state.get("itsx_summary")
if res:
    st.markdown("---")
    st.markdown("**Per-strain regions recovered**")
    rows = [{"Strain": sid, **counts} for sid, counts in res["per_strain"].items()]
    if rows:
        st.dataframe(pd.DataFrame(rows), width="stretch", hide_index=True)
    st.markdown("**Combined FASTAs**")
    if res["combined"]:
        for r, (p, n) in res["combined"].items():
            st.write(f"✅ {r}_combined.fasta — {n} sequence(s)  ·  "
                     f"`{p.replace(str(Path.home()), '~')}`")
    else:
        st.caption("No combined files produced (no regions detected).")
