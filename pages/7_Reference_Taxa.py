"""
pages/7_Reference_Taxa.py
-------------------------
Reference Taxa / Tips — component page (D-020 / RM-008 component 3).

Bring in the *comparison taxa* (tree tips) for your loci. Paste a flat, mixed list of GenBank
accessions and each is **auto-classified to its locus** (D-011 synonyms) and stored in the
per-project tips library (`<project>/tips/<locus>/`), kept separate from the protein extraction
guides (cf. D-018). Tips are the nucleotide barcode sequences of other taxa you'll align against
your isolates' extracted loci.

This sits after extraction: extract your standard set first, then layer comparison data on top.
"""

import sys
from pathlib import Path

import pandas as pd
import streamlit as st

sys.path.insert(0, str(Path(__file__).parent.parent / "src"))
from phylofetch.config import load_config
from phylofetch.ncbi_utils import (
    LOCUS_CATALOGUE,
    count_refs,
    delete_from_library,
    load_ref_meta,
    load_ref_records,
    set_email,
)
from phylofetch.project_manager import DEFAULT_PROJECT_DIR
from phylofetch.tips_utils import import_tip_accessions, project_tips_dir

st.set_page_config(page_title="Reference Taxa / Tips", page_icon="🌿", layout="wide")
st.title("🌿 Reference Taxa / Tips")
st.caption("Import comparison taxa (tree tips) by accession — pasted as a flat, mixed list and "
           "auto-sorted to the right locus. Stored per-project, separate from extraction guides; "
           "these align against your isolates' extracted loci.")

cfg = load_config()
project_dir = cfg.get("project_dir", str(DEFAULT_PROJECT_DIR))
ncbi_email = cfg.get("ncbi_email", "")
if ncbi_email:
    set_email(ncbi_email)
tips_dir = project_tips_dir(project_dir)

st.caption(f"📁 Project: `{Path(project_dir).name}`  ·  tips: "
           f"`{str(tips_dir).replace(str(Path.home()), '~')}`")
if not ncbi_email:
    st.error("⚠️ Set your NCBI Entrez email on the **NCBI References** page before importing.")
    st.stop()

# ── 1 · Paste accessions → auto-classify ──────────────────────────────────────
st.subheader("1 · Paste accessions (auto-classified to locus)")
st.caption("One per line or comma/space-separated (e.g. a spreadsheet column, or accessions "
           "copied from a publication or an NCBI BLAST result). Each is fetched, sorted to its "
           "locus by its GenBank title, and stored as a tip.")
pasted = st.text_area("Accessions", height=130, placeholder="MN123456\nKC584456.1\nDQ677938.1 …")

if st.button("⬇️ Import & auto-classify", type="primary"):
    accs = [a for chunk in pasted.replace(",", " ").split() for a in [chunk.strip()] if a]
    if not accs:
        st.warning("Paste at least one accession.")
    else:
        with st.spinner(f"Fetching + classifying {len(accs)} accession(s)…"):
            res = import_tip_accessions(accs, tips_dir)
        st.session_state["tips_import"] = res

res = st.session_state.get("tips_import")
if res:
    assigned = res.get("assigned", {})
    if assigned:
        st.success("Assigned to loci:")
        st.dataframe(pd.DataFrame([{"Locus": loc, "Accessions": ", ".join(accs)}
                                   for loc, accs in assigned.items()]),
                     width="stretch", hide_index=True)
    for err in res.get("errors", []):
        st.warning(err)

    unassigned = res.get("unassigned", [])
    if unassigned:
        st.warning(f"⚠️ {len(unassigned)} accession(s) could not be auto-classified "
                   "(no clear single-locus title match): " + ", ".join(unassigned))
        ua1, ua2 = st.columns([2, 1])
        with ua1:
            man_loc = st.selectbox("Assign these to locus", list(LOCUS_CATALOGUE))
        with ua2:
            st.write(""); st.write("")
            if st.button("Assign unassigned"):
                with st.spinner("Importing…"):
                    r2 = import_tip_accessions(unassigned, tips_dir, force_locus=man_loc)
                st.success(f"Assigned {len(unassigned)} to {man_loc}.")
                st.session_state["tips_import"] = None
                st.rerun()

st.markdown("---")
st.caption("💡 Target-taxa search (≥3 accessions/locus to review) is planned as a second import "
           "mode (RM-008). For now, paste accessions from a recent publication or a web BLAST.")

# ── 2 · Tips library (review / cull) ──────────────────────────────────────────
st.subheader("2 · Tips library")
loci_with_tips = [loc for loc in LOCUS_CATALOGUE if count_refs(loc, ref_dir=tips_dir) > 0]
if not loci_with_tips:
    st.caption("No tips imported yet.")
else:
    st.caption("Per-locus tip counts: " + " · ".join(
        f"**{loc}** {count_refs(loc, ref_dir=tips_dir)}" for loc in loci_with_tips))
    rloc = st.selectbox("Review locus", loci_with_tips)
    meta = load_ref_meta(rloc, ref_dir=tips_dir)
    rows = []
    for r in load_ref_records(rloc, ref_dir=tips_dir):
        m = meta.get(r.id) or meta.get(r.id.split(".")[0])
        acc = r.id
        rows.append({
            "Delete?": False,
            "Accession": acc,
            "Organism": (m.organism if m else ""),
            "Type": (m.type_kind if (m and m.is_type) else ""),
            "Length": len(r.seq),
            "NCBI": f"https://www.ncbi.nlm.nih.gov/nuccore/{acc}",
        })
    edited = st.data_editor(
        pd.DataFrame(rows), width="stretch", hide_index=True, key=f"tips_review_{rloc}",
        column_config={
            "Delete?": st.column_config.CheckboxColumn(),
            "NCBI": st.column_config.LinkColumn("NCBI", display_text="open"),
        },
    )
    to_del = edited[edited["Delete?"] == True]["Accession"].tolist()
    if to_del and st.button(f"🗑️ Delete {len(to_del)} from {rloc}"):
        for acc in to_del:
            delete_from_library(rloc, acc, ref_dir=tips_dir)
        st.success(f"Deleted {len(to_del)} from {rloc}.")
        st.rerun()
