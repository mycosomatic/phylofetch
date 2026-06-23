"""
pages/2_NCBI_References.py
--------------------------
NCBI References — component page (D-012 / RM-007 step 4b; repurposed D-023).

Fetch **optional, taxon-closer extraction references** into the per-project library
(<project>/references, D-013). Since D-020 made the bundled universal protein guides the default
extraction source, this page is no longer required — its job now is to pull genome-annotated
**protein** orthologs from (or near) the project taxon to *supplement* the bundled Asco/Basidio
core on the Exonerate page's "Bundled + project library (taxon-closer)" mode, and to fetch
**nucleotide** references for the relaxed-BLAST amplicon path. Standalone + chainable (provenance
to the manifest).

Coding loci only: rDNA (ITS/LSU/SSU) is **not** here — ITSx extracts it from assemblies (no refs
needed) and rDNA comparison sequences are imported as tips on the Reference Taxa page (D-023).
"""

import sys
from pathlib import Path

import pandas as pd
import streamlit as st

sys.path.insert(0, str(Path(__file__).parent.parent / "src"))
from phylofetch.blast_loci_utils import detect_fasta_type
from phylofetch.config import load_config, save_config
from phylofetch.ncbi_utils import (
    LOCUS_CATALOGUE,
    build_entrez_query,
    count_refs,
    delete_from_library,
    fetch_and_store,
    load_ref_meta,
    load_ref_records,
    locus_ref_fasta,
    locus_search_terms,
    project_ref_dir,
    search_ncbi_nucleotide,
    search_ncbi_protein,
    set_email,
    taxon_fallbacks,
)
from phylofetch.project_manager import (
    DEFAULT_PROJECT_DIR,
    load_project_manifest,
    update_step,
)

st.set_page_config(page_title="NCBI References", page_icon="📚", layout="wide")
st.title("📚 NCBI References — taxon-closer guides")
st.caption("**Optional.** The bundled universal protein guides (D-020) already extract coding "
           "loci kingdom-wide. Use this page to fetch **taxon-closer protein** orthologs that "
           "*supplement* those guides (Exonerate's 'Bundled + project library' mode), or "
           "**nucleotide** refs for the relaxed-BLAST path. Coding loci only — rDNA is handled by "
           "ITSx (extraction) and the Reference Taxa page (comparison tips).")

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

RDNA = ["ITS", "LSU", "SSU"]      # excluded here (D-023): ITSx extracts rDNA, tips compare it
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
    help="Defaults to the project taxon from the Assembly Manager. If the exact taxon returns "
         "nothing (e.g. a novel 'aff.' species not in NCBI), the search automatically falls "
         "back to the genus.",
)

c1, c2, c3 = st.columns([1.2, 1.2, 1])
with c1:
    ref_type_label = st.radio(
        "Reference type", ["Protein (recommended for coding)", "Nucleotide"],
        help="Protein → Exonerate protein2genome guide (frame-safe; supplements the bundled "
             "universal guides on the Exonerate page). Nucleotide → genomic barcodes for the "
             "relaxed-BLAST amplicon path.",
    )
ref_type = "protein" if ref_type_label.startswith("Protein") else "nucleotide"
with c2:
    if ref_type == "protein":
        refseq_only = st.checkbox(
            "RefSeq genome-annotated only", value=True,
            help="Keep only RefSeq XP_/NP_ proteins — full-length orthologs from annotated "
                 "genomes, not partial-CDS barcode translations (the latter include junk like a "
                 "5-aa 'TEF1' fragment). With this on, a novel species with no RefSeq protein "
                 "correctly falls back to the genus. Uncheck only if your lineage has no RefSeq "
                 "for a marker — then broaden the taxon or use the bundled guides.")
        min_len = st.number_input("Min length (aa)", value=100, min_value=0, max_value=5000, step=25)
    else:
        refseq_only = False
        min_len = st.number_input("Min length (bp)", value=300, min_value=0, max_value=20000, step=50)
with c3:
    preselect = st.number_input(
        "Pre-tick top N / locus", value=2, min_value=1, max_value=10,
        help="Candidates are sorted RefSeq + longest first; the top N are pre-ticked. Exonerate "
             "needs only one good ortholog — adjust the ticks before fetching.")
if ref_type == "protein":
    st.caption("🧬 **Protein** guides layer onto the bundled universal core on the Exonerate page; "
               "RefSeq-only pulls the closest *annotated-genome* ortholog and skips barcode junk.")


def _field_for(db: str) -> str:
    return "Protein Name" if db == "protein" else "Title"


def _existing_ref_type(loc: str):
    """Type ('nucleotide'/'protein') of refs already in this project's locus library, or None."""
    fa = Path(locus_ref_fasta(loc, ref_dir=ref_dir))
    if fa.exists() and fa.stat().st_size > 0:
        return detect_fasta_type(str(fa))
    return None

st.markdown("---")

# ── 2 · Loci checkboxes ───────────────────────────────────────────────────────
st.subheader("2 · Loci")
bsel, bclr, _ = st.columns([1, 1, 5])
with bsel:
    if st.button("Select all"):
        for loc in CODING:
            st.session_state[f"chk_{loc}"] = True
        st.rerun()
with bclr:
    if st.button("Clear"):
        for loc in CODING:
            st.session_state[f"chk_{loc}"] = False
        st.rerun()


def _locus_checkbox(locus: str) -> bool:
    n = count_refs(locus, ref_dir=ref_dir)
    label = f"{locus} · {n} in library" if n else locus
    return st.checkbox(label, key=f"chk_{locus}")


selected: list[str] = []
st.markdown("**Protein-coding loci**")
ccols = st.columns(4)
for i, loc in enumerate(CODING):
    with ccols[i % 4]:
        if _locus_checkbox(loc):
            selected.append(loc)

st.markdown("---")

# ── 3 · Find & pick references ───────────────────────────────────────────────────────
st.subheader("3 · Find & pick references")
_REFSEQ_PREFIXES = ("XP_", "NP_", "WP_", "YP_")


def _is_refseq(acc: str) -> bool:
    return acc[:3] in _REFSEQ_PREFIXES


if not selected:
    st.info("Tick one or more loci above.")
elif not taxon.strip():
    st.warning("Enter a taxon above to search.")
else:
    if st.button("🔍 Find candidates", type="primary"):
        cands: dict[str, list] = {}
        used_tax: dict[str, str] = {}
        prog = st.progress(0.0)
        for i, loc in enumerate(selected):
            terms = locus_search_terms(loc)
            hits, used = [], taxon.strip()
            for cand_tax in taxon_fallbacks(taxon):       # exact taxon → genus
                used = cand_tax
                try:
                    if ref_type == "protein":
                        h = search_ncbi_protein(terms, cand_tax, max_results=40,
                                                field="Protein Name", refseq_only=refseq_only)
                    else:
                        h = search_ncbi_nucleotide(terms, cand_tax, max_results=40)
                except Exception as e:                    # noqa: BLE001 — surface to UI
                    st.warning(f"{loc} ({cand_tax}): {e}")
                    h = []
                h = [x for x in h if x.get("length", 0) >= int(min_len)]
                if h:
                    hits = h
                    break
            hits.sort(key=lambda x: (_is_refseq(x.get("accession", "")), x.get("length", 0)),
                      reverse=True)
            cands[loc], used_tax[loc] = hits, used
            prog.progress((i + 1) / len(selected))
        st.session_state["ref_cands"] = {
            "taxon": taxon, "ref_type": ref_type, "refseq_only": refseq_only,
            "min_len": int(min_len), "preselect": int(preselect),
            "cands": cands, "used_tax": used_tax,
        }

    data = st.session_state.get("ref_cands")
    if data and data["taxon"] == taxon and data["ref_type"] == ref_type:
        cands, used_tax = data["cands"], data["used_tax"]
        ncbi_db = "protein" if ref_type == "protein" else "nuccore"
        rows, empty, mism = [], [], []
        for loc in selected:
            ex = _existing_ref_type(loc)
            is_mism = ex is not None and ex != ref_type
            if is_mism:
                mism.append((loc, ex))
            lst = cands.get(loc, [])
            if not lst:
                empty.append(loc)
            for j, h in enumerate(lst):
                acc = h.get("accession", "")
                rows.append({
                    "Keep?": (not is_mism) and (j < data["preselect"]),
                    "Locus": loc,
                    "Accession": acc,
                    "Organism": h.get("organism", ""),
                    "Length": h.get("length", 0),
                    "RefSeq": "✓" if _is_refseq(acc) else "",
                    "Taxon": used_tax.get(loc, "")
                             + (" (genus)" if used_tax.get(loc, "") != taxon.strip() else ""),
                    "NCBI": f"https://www.ncbi.nlm.nih.gov/{ncbi_db}/{acc}",
                })
        if rows:
            st.caption("Tick the references to fetch — sorted RefSeq + longest first, top "
                       f"{data['preselect']} pre-ticked. One good full-length ortholog per locus "
                       "is plenty for an Exonerate guide.")
            edited = st.data_editor(
                pd.DataFrame(rows), width="stretch", hide_index=True,
                key=f"refed_{ref_type}_{data['refseq_only']}_{data['min_len']}_{len(rows)}",
                column_config={
                    "Keep?": st.column_config.CheckboxColumn(),
                    "Length": st.column_config.NumberColumn(format="%d"),
                    "NCBI": st.column_config.LinkColumn("NCBI", display_text="open"),
                },
                disabled=["Locus", "Accession", "Organism", "Length", "RefSeq", "Taxon", "NCBI"],
            )
            if mism:
                st.warning("⚠️ Type mismatch (left unticked, will skip): "
                           + " · ".join(f"`{loc}` already holds **{ex}**" for loc, ex in mism)
                           + " — clear it in Review library / Manage Data before fetching the "
                           "other type (mixing protein + nucleotide in one file breaks Exonerate).")
            if empty:
                st.caption("No candidates for: " + ", ".join(f"`{l}`" for l in empty)
                           + " — broaden the taxon (genus → family/order), relax RefSeq / min "
                           "length, or just use the bundled guides for these.")

            n_keep = int((edited["Keep?"] == True).sum())
            if st.button(f"⬇️ Fetch {n_keep} ticked reference(s)", type="primary",
                         disabled=n_keep == 0):
                keep = edited[edited["Keep?"] == True]
                results: dict[str, str] = {}
                locs_kept = sorted(set(keep["Locus"]))
                prog = st.progress(0.0)
                for i, loc in enumerate(locs_kept):
                    accs = keep[keep["Locus"] == loc]["Accession"].tolist()
                    ex = _existing_ref_type(loc)
                    if ex is not None and ex != ref_type:
                        results[loc] = f"skipped — library already has {ex} refs"
                    else:
                        resolved = build_entrez_query(locus_search_terms(loc), used_tax.get(loc, ""),
                                                      field=_field_for(ref_type))
                        try:
                            added, skipped, errs = fetch_and_store(
                                accs, loc, db=ref_type, ref_dir=ref_dir,
                                query=f"{resolved} [{ref_type}, taxon={used_tax.get(loc, '')}, "
                                      f"refseq_only={data['refseq_only']}]")
                            results[loc] = (f"+{added} added · {skipped} already present"
                                            + (f" · ⚠️ {'; '.join(errs)}" if errs else ""))
                        except Exception as e:            # noqa: BLE001 — surface to UI
                            results[loc] = f"error: {e}"
                    prog.progress((i + 1) / len(locs_kept))
                update_step(
                    project_dir, "references", status="done",
                    outputs={loc: {"in_library": count_refs(loc, ref_dir=ref_dir),
                                   "db": ref_type, "taxon": used_tax.get(loc, "")}
                             for loc in locs_kept},
                    notes=f"taxon={taxon}; ref_type={ref_type}; refseq_only={data['refseq_only']}; "
                          f"loci={','.join(locs_kept)}")
                st.success("Fetch complete — provenance recorded in the manifest.")
                for loc, msg in results.items():
                    st.write(f"**{loc}**: {msg}")
                st.session_state.pop("ref_cands", None)
                st.rerun()
        else:
            st.warning("No candidates found. Broaden the taxon, relax RefSeq / min length, or just "
                       "use the bundled guides (they already cover coding loci kingdom-wide).")

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
