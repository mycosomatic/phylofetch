"""
pages/2_Primers.py
------------------
Primers (in-silico PCR) — component page (D-012 / RM-007 step 4e).

Locate fwd+rev primer binding sites in each assembly (blastn-short; IUPAC degenerate primers
expanded to concrete oligos, D-009) and extract the amplicon between them — no NCBI reference
library needed. Pick pairs from the citable built-in catalogue (D-003), your saved library, or
enter custom sequences; define a custom locus for anything not in the catalogue. Per-project
outputs (<project>/results/loci, D-015); provenance to the manifest.

Standalone + chainable; runs alongside the old 2_Loci_Extraction page until the monolith retires.
"""

import os
import sys
from pathlib import Path

import pandas as pd
import streamlit as st
from Bio import SeqIO

sys.path.insert(0, str(Path(__file__).parent.parent / "src"))
from phylofetch.config import load_config
from phylofetch.primer_utils import (
    PrimerPair,
    delete_user_primer,
    find_primer_amplicons_escalating,
    get_primer_catalogue,
    load_user_primers,
    locus_primer_map,
    run_primer_extraction,
    save_user_primer,
)
from phylofetch.project_manager import (
    DEFAULT_PROJECT_DIR,
    RunManager,
    effective_taxon,
    load_assembly_registry,
    load_project_manifest,
    project_output_dir,
    update_step,
)

st.set_page_config(page_title="Primers (in-silico PCR)", page_icon="🔬", layout="wide")
st.title("🔬 Primers — in-silico PCR")
st.caption("Locate primer binding sites in each assembly and extract the amplicon between them. "
           "No NCBI references needed — useful when accessions are missing or map poorly. "
           "IUPAC degenerate primers are handled (D-009).")

cfg = load_config()
project_dir = cfg.get("project_dir", str(DEFAULT_PROJECT_DIR))
tp = cfg.get("tool_paths", {})
blastn_bin = tp.get("blastn", "blastn")
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

primer_cat = get_primer_catalogue()
primer_locus_map = locus_primer_map(primer_cat)

# ── 1 · Assemblies ────────────────────────────────────────────────────────────
st.subheader("1 · Assemblies")


def _fmt(sid: str) -> str:
    t = effective_taxon(registry[sid], default_taxon)
    return f"{sid}  ·  {t}" if t else sid


sel = st.multiselect("Assemblies", list(registry), default=list(registry), format_func=_fmt)

# ── 2 · Loci + primer assignment ──────────────────────────────────────────────
st.subheader("2 · Loci & primer pairs")
max_mm = st.slider("Max primer mismatches (auto-escalates: tries stricter thresholds first)",
                   0, 4, 3,
                   help="Per-primer edit distance (substitutions + unaligned bases). The search "
                        "starts strict and loosens up to this cap only if no amplicon is found — "
                        "so a related species' primers in a divergent genome still bind, without "
                        "surfacing loose off-targets when a clean match exists.")
sel_loci = st.multiselect("Loci with catalogue primers", sorted(primer_locus_map),
                          help="Pick loci to amplify; assign each a pair below. "
                               "Use the custom-locus expander for anything not listed.")

if "primer_assignments" not in st.session_state:
    st.session_state.primer_assignments = {}
assignments: dict[str, PrimerPair] = {}

for locus in sel_loci:
    st.markdown(f"**{locus}**")
    opts = ["— pick —"] + primer_locus_map.get(locus, []) + ["Custom…"]
    pa, pb = st.columns([2, 3])
    with pa:
        src = st.selectbox("Primer pair", opts, key=f"pp_src_{locus}",
                           label_visibility="collapsed")
    with pb:
        if src == "Custom…":
            c1, c2, c3, c4 = st.columns([2, 2, 1, 1])
            fwd = c1.text_input("Fwd 5'→3'", key=f"pp_fwd_{locus}", placeholder="ATGC…")
            rev = c2.text_input("Rev 5'→3'", key=f"pp_rev_{locus}", placeholder="ATGC…")
            mn = c3.number_input("Min bp", value=100, step=50, key=f"pp_min_{locus}")
            mx = c4.number_input("Max bp", value=5000, step=100, key=f"pp_max_{locus}")
            if fwd and rev:
                assignments[locus] = PrimerPair(
                    name=f"custom_{locus}", locus=locus, fwd=fwd.strip().upper(),
                    rev=rev.strip().upper(), min_amplicon=int(mn), max_amplicon=int(mx),
                    origin="user")
            else:
                st.warning(f"Enter both primers for {locus}.")
        elif src != "— pick —":
            pp = primer_cat[src]
            badge = "👤 user" if pp.origin == "user" else "📚 built-in"
            st.caption(f"{badge} · `{pp.fwd}` / `{pp.rev}` · {pp.min_amplicon}–{pp.max_amplicon} bp"
                       + (f" · 📖 {pp.source}" if pp.source else ""))
            assignments[locus] = pp

# ── Custom locus ──────────────────────────────────────────────────────────────
if "primer_custom_loci" not in st.session_state:
    st.session_state.primer_custom_loci = {}
with st.expander("➕ Custom locus — primer pair for a locus not in the catalogue"):
    d1, d2, d3 = st.columns([2, 3, 2])
    with d1:
        cl_name = st.text_input("Locus name", key="cl_name", placeholder="e.g. MCM7, BenA")
    with d2:
        cl_fwd = st.text_input("Fwd 5'→3'", key="cl_fwd", placeholder="ACIMGIGTITC… (IUPAC OK)")
        cl_rev = st.text_input("Rev 5'→3'", key="cl_rev", placeholder="GAYTTDGCIAC…")
    with d3:
        cl_min = st.number_input("Min bp", value=100, step=50, key="cl_min")
        cl_max = st.number_input("Max bp", value=5000, step=100, key="cl_max")
    cl_save = st.checkbox("Also save to my library", key="cl_save")
    if st.button("➕ Add custom locus"):
        name = (cl_name or "").strip().upper().replace(" ", "_")
        if not name:
            st.warning("Enter a locus name.")
        elif not (cl_fwd.strip() and cl_rev.strip()):
            st.warning("Enter both primers.")
        else:
            st.session_state.primer_custom_loci[name] = {
                "locus": name, "fwd": cl_fwd.strip().upper(), "rev": cl_rev.strip().upper(),
                "min_amplicon": int(cl_min), "max_amplicon": int(cl_max), "source": "user-added"}
            if cl_save:
                save_user_primer(PrimerPair(
                    name=f"{name}_custom", locus=name, fwd=cl_fwd.strip().upper(),
                    rev=cl_rev.strip().upper(), min_amplicon=int(cl_min),
                    max_amplicon=int(cl_max), source="user-added", origin="user"))
            st.success(f"Added custom locus '{name}'.")
            st.rerun()
    for cname, rec in list(st.session_state.primer_custom_loci.items()):
        rc1, rc2 = st.columns([5, 1])
        rc1.caption(f"**{cname}** · F:`{rec['fwd']}` R:`{rec['rev']}` · "
                    f"{rec['min_amplicon']}–{rec['max_amplicon']} bp")
        if rc2.button("🗑️ Remove", key=f"cl_del_{cname}"):
            del st.session_state.primer_custom_loci[cname]
            st.rerun()

for cname, rec in st.session_state.primer_custom_loci.items():
    assignments[cname] = PrimerPair(
        name=f"{cname}_custom", locus=cname, fwd=rec["fwd"], rev=rec["rev"],
        min_amplicon=int(rec["min_amplicon"]), max_amplicon=int(rec["max_amplicon"]),
        source=rec.get("source", "user-added"), origin="user")

user_lib = load_user_primers()
if user_lib:
    with st.expander(f"📁 My primer library — {len(user_lib)} saved pair(s)"):
        for uname, upp in user_lib.items():
            uc1, uc2 = st.columns([5, 1])
            uc1.caption(f"**{uname}** · {upp.locus} · F:`{upp.fwd}` R:`{upp.rev}`")
            if uc2.button("🗑️ Delete", key=f"del_user_{uname}"):
                delete_user_primer(uname)
                st.rerun()

# ── 3 · Preview binding sites (optional disambiguation) ───────────────────────
st.subheader("3 · Preview binding sites (optional)")
with st.expander("🔍 Scan & choose binding sites (handle off-target hits)"):
    st.caption("Scan assemblies for every binding site, then pick which amplicon to extract per "
               "strain/locus. Skip this and the lowest-edit-distance site is used automatically.")
    if st.button("Scan binding sites", disabled=not (sel and assignments)):
        scan_mgr = RunManager(project_dir)
        scan: dict[str, list] = {}
        scan_used: dict[str, int] = {}
        jobs = [(s, l) for s in sel for l in assignments]
        sprog = st.progress(0.0)
        for i, (s, l) in enumerate(jobs, 1):
            asm = registry[s].get("assembly_path", "")
            try:
                cands, used = find_primer_amplicons_escalating(
                    asm, assignments[l], start_mismatches=min(2, max_mm), max_mismatches=max_mm,
                    blastn_bin=blastn_bin, manager=scan_mgr, action=f"primer_scan_{l}_{s}")
                scan[f"{s}|{l}"] = cands
                scan_used[f"{s}|{l}"] = used
            except (ValueError, FileNotFoundError) as exc:
                scan[f"{s}|{l}"] = []
                st.warning(f"`{s}` · {l}: {exc}")
            sprog.progress(i / max(len(jobs), 1))
        st.session_state["primer_scan"] = scan
        st.session_state["primer_scan_used"] = scan_used
    scan = st.session_state.get("primer_scan", {})
    scan_used = st.session_state.get("primer_scan_used", {})
    for s in sel:
        for l in assignments:
            cands = scan.get(f"{s}|{l}")
            if cands is None:
                continue
            used = scan_used.get(f"{s}|{l}")
            note = f" · matched at edit distance ≤{used}" if (cands and used is not None) else ""
            st.markdown(f"**`{s}` · {l}** — {len(cands)} site(s){note}")
            if not cands:
                st.caption("⚠️ No binding sites within the size window.")
                continue
            st.dataframe(pd.DataFrame([{
                "#": i, "contig": c["contig"], "strand": c["strand"],
                "coords": f"{c['amp_start']}..{c['amp_end']}", "len": c["amp_len"],
                "fwd_mm": c["fwd_edit"], "rev_mm": c["rev_edit"], "total_mm": c["total_edit"],
            } for i, c in enumerate(cands)]), hide_index=True, width="stretch")
            opts = [f"{i}: {c['contig']} {c['strand']} {c['amp_len']} bp (mm {c['total_edit']})"
                    for i, c in enumerate(cands)]
            choice = st.selectbox(f"Extract which site? ({s} / {l})", opts, key=f"primer_choice_{s}|{l}")
            st.session_state[f"primer_pick_{s}|{l}"] = opts.index(choice)

# ── 4 · Run ───────────────────────────────────────────────────────────────────
st.subheader("4 · Run")
if not sel:
    st.info("Select at least one assembly.")
elif not assignments:
    st.info("Assign a primer pair to at least one locus (or add a custom locus).")

if st.button("🚀 Run in-silico PCR", type="primary", disabled=not (sel and assignments)):
    per_strain_dir.mkdir(parents=True, exist_ok=True)
    combined_dir.mkdir(parents=True, exist_ok=True)
    manager = RunManager(project_dir)
    scan = st.session_state.get("primer_scan", {})
    prog = st.progress(0.0)
    jobs = [(s, l) for s in sel for l in assignments]
    for n, (strain_id, locus) in enumerate(jobs, 1):
        assembly = registry[strain_id].get("assembly_path", "")
        pp = assignments[locus]
        locus_out = per_strain_dir / strain_id / locus
        locus_out.mkdir(parents=True, exist_ok=True)
        key = f"{strain_id}|{locus}"
        with st.spinner(f"[{strain_id}] {locus} ({pp.name})…"):
            result, status = run_primer_extraction(
                assembly_fasta=assembly, primer_pair=pp, output_dir=str(locus_out),
                strain_id=strain_id, locus_name=locus,
                max_mismatches=min(2, max_mm), escalate_to=max_mm,
                blastn_bin=blastn_bin, manager=manager,
                candidates=scan.get(key) or None,
                chosen_index=st.session_state.get(f"primer_pick_{key}", 0))
        if result:
            extra = f" · {result['n_candidates']} site(s)" if result["n_candidates"] > 1 else ""
            st.write(f"✅ {strain_id} · {locus} ({pp.name}): {result['amp_len']} bp on "
                     f"`{result['contig']}` [{result['strand']}] "
                     f"· fwd_mm={result['fwd_edit']} rev_mm={result['rev_edit']}{extra}")
        else:
            st.write(f"⚠️ {strain_id} · {locus}: {status}")
        prog.progress(n / max(len(jobs), 1))

    combined = {}
    for locus in assignments:
        recs = []
        for sd in sorted(per_strain_dir.iterdir()):
            fp = sd / locus / f"{locus}_amplicon.fasta"
            if fp.exists():
                recs.extend(list(SeqIO.parse(str(fp), "fasta")))
        if recs:
            out = combined_dir / f"{locus}_amplicon_combined.fasta"
            SeqIO.write(recs, str(out), "fasta")
            combined[locus] = (str(out), len(recs))
    update_step(project_dir, "primers", status="done",
                outputs={l: {"combined": p, "n": n} for l, (p, n) in combined.items()},
                notes=f"strains={len(sel)}; loci={','.join(assignments)}; max_mm={max_mm}")
    st.session_state["primer_combined"] = combined
    st.success("In-silico PCR complete — provenance recorded in the manifest.")
    st.rerun()

# ── Results ───────────────────────────────────────────────────────────────────
res = st.session_state.get("primer_combined")
if res:
    st.markdown("---")
    st.markdown("**Combined amplicon FASTAs**")
    for locus, (p, n) in res.items():
        st.write(f"✅ {locus}_amplicon_combined.fasta — {n} sequence(s)  ·  "
                 f"`{p.replace(str(Path.home()), '~')}`")
    if not res:
        st.caption("No amplicons recovered.")
