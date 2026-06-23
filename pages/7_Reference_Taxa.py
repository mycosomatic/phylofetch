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
from phylofetch.tips_utils import (
    import_tips_with_assignments,
    lookup_accessions,
    project_tips_dir,
)

SKIP = "— skip —"   # Locus value meaning "don't import this accession"

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

# ── 1 · Paste accessions → look up → assign locus per accession ───────────────
st.subheader("1 · Paste accessions")
st.caption("One per line or comma/space-separated (e.g. a spreadsheet column, or accessions "
           "from a publication or a web BLAST). Each is looked up on NCBI; you then **confirm or "
           "correct the locus for each** before importing. RefSeq accessions pasted without the "
           "underscore (`NR135944` → `NR_135944`) are auto-repaired.")
pasted = st.text_area("Accessions", height=130, placeholder="MN123456\nKC584456.1\nNR_135944 …")

if st.button("🔎 Look up on NCBI", type="primary"):
    accs = [c.strip() for c in pasted.replace(",", " ").split() if c.strip()]
    if not accs:
        st.warning("Paste at least one accession.")
    else:
        with st.spinner(f"Looking up {len(accs)} accession(s)…"):
            st.session_state["tips_lookup"] = lookup_accessions(accs)
        st.session_state.pop("tips_import_result", None)

lookup = st.session_state.get("tips_lookup")
if lookup:
    missing = [r for r in lookup if not r["found"]]
    found = [r for r in lookup if r["found"]]
    if missing:
        st.warning("⚠️ Not found on NCBI (won't be imported — check the accession / database): "
                   + ", ".join(f"`{r['input']}`"
                               + ("" if r["input"] == r["accession"] else f" → `{r['accession']}`")
                               for r in missing))
    if found:
        st.caption(f"Confirm or correct the locus for each, then import. Leave a row on "
                   f"“{SKIP}” to omit it.")
        loc_opts = [SKIP] + list(LOCUS_CATALOGUE)
        df = pd.DataFrame([{
            "Accession": r["accession"],
            "Locus": r["locus_guess"] or SKIP,
            "Title": r["title"],
            "NCBI": f"https://www.ncbi.nlm.nih.gov/nuccore/{r['accession']}",
        } for r in found])
        edited = st.data_editor(
            df, width="stretch", hide_index=True, key="tips_assign",
            column_config={
                "Accession": st.column_config.TextColumn(disabled=True),
                "Locus": st.column_config.SelectboxColumn(options=loc_opts, required=True),
                "Title": st.column_config.TextColumn("Title (NCBI)", disabled=True, width="large"),
                "NCBI": st.column_config.LinkColumn("NCBI", display_text="open"),
            },
        )
        n_guessed = sum(1 for r in found if r["locus_guess"])
        st.caption(f"Auto-classified {n_guessed}/{len(found)} by GenBank title; the rest need a "
                   f"locus picked (or left on “{SKIP}”).")
        if st.button("⬇️ Import with these assignments", type="primary"):
            assignments = {row["Accession"]: ("" if row["Locus"] == SKIP else row["Locus"])
                           for _, row in edited.iterrows()}
            if not any(assignments.values()):
                st.warning("No rows assigned to a locus.")
            else:
                with st.spinner("Fetching + storing tips…"):
                    st.session_state["tips_import_result"] = import_tips_with_assignments(
                        assignments, tips_dir)
                st.session_state.pop("tips_lookup", None)
                st.rerun()

result = st.session_state.get("tips_import_result")
if result:
    assigned = result.get("assigned", {})
    if assigned:
        st.success("Imported tips:")
        st.dataframe(pd.DataFrame([{"Locus": loc, "Accessions": ", ".join(accs)}
                                   for loc, accs in assigned.items()]),
                     width="stretch", hide_index=True)
    for err in result.get("errors", []):
        st.warning(err)
    if not assigned and not result.get("errors"):
        st.info("Nothing imported.")

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
