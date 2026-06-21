"""
pages/2_NCBI_References.py
--------------------------
NCBI References — component page (D-012 / RM-007 step 4b).

Build the per-project reference library (D-013): tick loci, preview NCBI hit counts for the
project taxon, then fetch references into <project>/references. Standalone (usable on its own)
and chainable (records provenance + step state in the project manifest). References feed the
BLAST and Exonerate extraction strategies.

Interim during the page decomposition: this lives alongside the old 2_Loci_Extraction page
until the ITSx / Exonerate / Primers component pages exist and the monolith is retired.
"""

import sys
from pathlib import Path

import pandas as pd
import streamlit as st

sys.path.insert(0, str(Path(__file__).parent.parent / "src"))
from phylofetch.config import load_config, save_config
from phylofetch.ncbi_utils import (
    LOCUS_CATALOGUE,
    build_entrez_query,
    count_refs,
    delete_from_library,
    fetch_and_store,
    load_ref_meta,
    load_ref_records,
    locus_search_terms,
    ncbi_search_count,
    project_ref_dir,
    search_ncbi_nucleotide,
    set_email,
)
from phylofetch.project_manager import (
    DEFAULT_PROJECT_DIR,
    load_project_manifest,
    update_step,
)

st.set_page_config(page_title="NCBI References", page_icon="📚", layout="wide")
st.title("📚 NCBI References")
st.caption("Build this project's reference library: pick loci, preview NCBI hit counts for "
           "your taxon, then fetch. References are stored per-project and used by the BLAST "
           "and Exonerate extraction strategies.")

cfg = load_config()
project_dir = cfg.get("project_dir", str(DEFAULT_PROJECT_DIR))
ncbi_email = cfg.get("ncbi_email", "")
if ncbi_email:
    set_email(ncbi_email)

ref_dir = project_ref_dir(project_dir)
manifest = load_project_manifest(project_dir)
default_taxon = manifest.get("default_taxon", "")

st.caption(f"📁 Project: `{Path(project_dir).name}`  ·  library: "
           f"`{str(ref_dir).replace(str(Path.home()), '~')}`")

RDNA = ["ITS", "LSU", "SSU"]
CODING = [k for k in LOCUS_CATALOGUE if k not in RDNA]

# ── NCBI email gate ───────────────────────────────────────────────────────────
with st.expander("🔑 NCBI Entrez email (required for fetching)", expanded=not bool(ncbi_email)):
    new_email = st.text_input("Email", value=ncbi_email, placeholder="you@example.com")
    if st.button("Save email"):
        save_config({"ncbi_email": new_email})
        set_email(new_email)
        st.rerun()
if ncbi_email:
    st.success(f"✅ NCBI email: {ncbi_email}")
else:
    st.error("⚠️ Set NCBI email above before searching.")
    st.stop()

st.markdown("---")

# ── 1 · Taxon ─────────────────────────────────────────────────────────────────
st.subheader("1 · Search taxon")
if not default_taxon:
    st.info("No project default taxon set yet (Assembly Manager → Project taxonomy). "
            "You can type one here for this search.")
taxon = st.text_input(
    "Organism / taxon for the reference search",
    value=default_taxon,
    placeholder="e.g. Alternaria",
    help="Defaults to the project taxon from the Assembly Manager. Broaden it "
         "(species → genus → family) if a locus returns few hits — a good reference set "
         "usually spans a wider clade than your isolates.",
)
type_mode_label = st.radio(
    "Type material", ["Prefer type", "All", "Type only"], horizontal=True,
    help="Prefer type: type-derived records sorted first (recommended). "
         "All: no preference. Type only: restrict to NCBI 'sequence from type material'.",
)
type_mode = {"Prefer type": "prefer", "All": "all", "Type only": "type_only"}[type_mode_label]

st.markdown("---")

# ── 2 · Loci checkboxes ───────────────────────────────────────────────────────
st.subheader("2 · Loci")
bsel, bclr, _ = st.columns([1, 1, 5])
with bsel:
    if st.button("Select all"):
        for loc in LOCUS_CATALOGUE:
            st.session_state[f"chk_{loc}"] = True
        st.rerun()
with bclr:
    if st.button("Clear"):
        for loc in LOCUS_CATALOGUE:
            st.session_state[f"chk_{loc}"] = False
        st.rerun()


def _locus_checkbox(locus: str) -> bool:
    n = count_refs(locus, ref_dir=ref_dir)
    label = f"{locus} · {n} in library" if n else locus
    return st.checkbox(label, key=f"chk_{locus}")


selected: list[str] = []
st.markdown("**rDNA** — ITSx extracts these from assemblies; fetch here mainly for outgroups.")
for c, loc in zip(st.columns(len(RDNA)), RDNA):
    with c:
        if _locus_checkbox(loc):
            selected.append(loc)

st.markdown("**Protein-coding**")
ccols = st.columns(4)
for i, loc in enumerate(CODING):
    with ccols[i % 4]:
        if _locus_checkbox(loc):
            selected.append(loc)

st.markdown("---")

# ── 3 · Preview & fetch ───────────────────────────────────────────────────────
st.subheader("3 · Preview & fetch")
if not selected:
    st.info("Tick one or more loci above.")
elif not taxon.strip():
    st.warning("Enter a taxon above to search.")
else:
    cprev, cmax = st.columns([1, 2])
    with cprev:
        do_preview = st.button("🔍 Preview hit counts", type="primary")
    with cmax:
        max_per = st.slider("Max references to fetch per locus", 1, 50, 15)

    if do_preview:
        counts: dict[str, int] = {}
        prog = st.progress(0.0)
        for i, loc in enumerate(selected):
            try:
                counts[loc] = ncbi_search_count(locus_search_terms(loc), taxon, type_mode=type_mode)
            except Exception as e:                          # noqa: BLE001 — surface to UI
                counts[loc] = -1
                st.warning(f"{loc}: {e}")
            prog.progress((i + 1) / len(selected))
        st.session_state["ref_preview"] = {"taxon": taxon, "type_mode": type_mode, "counts": counts}

    preview = st.session_state.get("ref_preview")
    if preview and preview["taxon"] == taxon and preview["type_mode"] == type_mode:
        counts = preview["counts"]
        rows = [{
            "Locus": loc,
            "NCBI hits": ("error" if counts.get(loc, 0) < 0 else counts.get(loc, "—")),
            "In library": count_refs(loc, ref_dir=ref_dir),
            "Will fetch": (0 if counts.get(loc, 0) <= 0 else min(counts[loc], max_per)),
        } for loc in selected]
        st.dataframe(pd.DataFrame(rows), width="stretch", hide_index=True)

        if st.button("⬇️ Fetch selected loci", type="primary"):
            results: dict[str, tuple] = {}
            prog = st.progress(0.0)
            for i, loc in enumerate(selected):
                n = counts.get(loc, 0)
                if n <= 0:
                    results[loc] = (0, 0, "no hits to fetch")
                    prog.progress((i + 1) / len(selected))
                    continue
                terms = locus_search_terms(loc)
                try:
                    hits = search_ncbi_nucleotide(terms, taxon, max_results=min(n, max_per),
                                                  type_mode=type_mode)
                    accs = [h["accession"] for h in hits if h.get("accession")]
                    resolved = build_entrez_query(terms, taxon, field="Title")
                    added, skipped, errors = fetch_and_store(
                        accs, loc, db="nucleotide",
                        query=f"{resolved} [nucleotide, type_mode={type_mode}]",
                        ref_dir=ref_dir,
                    )
                    results[loc] = (added, skipped, "; ".join(errors))
                except Exception as e:                      # noqa: BLE001 — surface to UI
                    results[loc] = (0, 0, str(e))
                prog.progress((i + 1) / len(selected))

            update_step(
                project_dir, "references", status="done",
                outputs={loc: {"in_library": count_refs(loc, ref_dir=ref_dir)} for loc in selected},
                notes=f"taxon={taxon}; loci={','.join(selected)}; type_mode={type_mode}",
            )
            st.success("Fetch complete — provenance recorded in the project manifest.")
            for loc, (added, skipped, err) in results.items():
                st.write(f"**{loc}**: +{added} added · {skipped} already present"
                         + (f" · ⚠️ {err}" if err else ""))
            st.rerun()

st.markdown("---")

# ── 4 · Review / cull library ─────────────────────────────────────────────────
st.subheader("4 · Review library")
loci_in_lib = [loc for loc in LOCUS_CATALOGUE if count_refs(loc, ref_dir=ref_dir) > 0]
if not loci_in_lib:
    st.caption("No references fetched yet for this project.")
else:
    rloc = st.selectbox("Locus to review", loci_in_lib)
    meta = load_ref_meta(rloc, ref_dir=ref_dir)
    rows = []
    for r in load_ref_records(rloc, ref_dir=ref_dir):
        m = meta.get(r.id) or meta.get(r.id.split(".")[0])
        rows.append({
            "Delete?": False,
            "Accession": r.id,
            "Organism": (m.organism if m else ""),
            "Type": (m.type_kind if (m and m.is_type) else ""),
            "Length": len(r.seq),
        })
    edited = st.data_editor(
        pd.DataFrame(rows), width="stretch", hide_index=True,
        column_config={"Delete?": st.column_config.CheckboxColumn()},
        key=f"review_{rloc}",
    )
    to_del = edited[edited["Delete?"] == True]["Accession"].tolist()
    if to_del and st.button(f"🗑️ Delete {len(to_del)} selected from {rloc}"):
        for acc in to_del:
            delete_from_library(rloc, acc, ref_dir=ref_dir)
        st.success(f"Deleted {len(to_del)} from {rloc}.")
        st.rerun()
