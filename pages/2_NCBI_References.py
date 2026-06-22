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
    ncbi_search_count,
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

cta, ctb = st.columns(2)
with cta:
    ref_type_label = st.radio(
        "Reference type", ["Protein (recommended for coding)", "Nucleotide"], horizontal=False,
        help="Protein → Exonerate uses protein2genome: a protein query is intron-free and "
             "pins the reading frame, giving clean frame-checked CDS even from a congener "
             "(cross-species works) — this is what supplements the bundled guides. Nucleotide → "
             "genomic barcodes for the relaxed-BLAST amplicon path / GenBank-comparable trees.",
    )
with ctb:
    type_mode_label = st.radio(
        "Type material", ["Prefer type", "All", "Type only"], horizontal=False,
        help="Prefer type: type-derived records sorted first. All: no preference. "
             "Type only: restrict to NCBI 'sequence from type material' (nucleotide barcodes).",
    )
ref_type = "protein" if ref_type_label.startswith("Protein") else "nucleotide"
type_mode = {"Prefer type": "prefer", "All": "all", "Type only": "type_only"}[type_mode_label]
if ref_type == "protein":
    st.caption("🧬 **Protein** references (protein2genome, frame-safe) — these layer onto the "
               "bundled universal guides on the Exonerate page as a taxon-closer ortholog.")


def _locus_db(loc: str) -> str:
    """rRNA has no protein, so rDNA is always nucleotide; coding follows the chosen type."""
    return "nucleotide" if loc in RDNA else ref_type


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
        used_tax: dict[str, str] = {}
        dbs: dict[str, str] = {}
        prog = st.progress(0.0)
        for i, loc in enumerate(selected):
            db = _locus_db(loc)
            fld = _field_for(db)
            terms = locus_search_terms(loc)
            n, used = 0, taxon.strip()
            for cand in taxon_fallbacks(taxon):     # exact taxon → genus
                try:
                    c = ncbi_search_count(terms, cand, db=db, field=fld, type_mode=type_mode)
                except Exception as e:              # noqa: BLE001 — surface to UI
                    c = -1
                    st.warning(f"{loc} ({cand}): {e}")
                used = cand
                if c > 0:
                    n = c
                    break
            counts[loc], used_tax[loc], dbs[loc] = n, used, db
            prog.progress((i + 1) / len(selected))
        st.session_state["ref_preview"] = {
            "taxon": taxon, "type_mode": type_mode, "ref_type": ref_type,
            "counts": counts, "used_tax": used_tax, "dbs": dbs,
        }

    preview = st.session_state.get("ref_preview")
    if (preview and preview["taxon"] == taxon and preview["type_mode"] == type_mode
            and preview.get("ref_type") == ref_type):
        counts = preview["counts"]
        used_tax = preview["used_tax"]
        dbs = preview["dbs"]
        rows, mismatches = [], []
        for loc in selected:
            db = dbs.get(loc, "?")
            ex = _existing_ref_type(loc)
            mism = ex is not None and ex != db
            if mism:
                mismatches.append((loc, ex, db))
            rows.append({
                "Locus": loc,
                "Fetch type": db,
                "Existing": (f"⚠️ {ex}" if mism else (ex or "—")),
                "Taxon used": used_tax.get(loc, taxon) + (" (genus fallback)"
                              if used_tax.get(loc, "") != taxon.strip() else ""),
                "NCBI hits": counts.get(loc, "—"),
                "Will fetch": (0 if (mism or counts.get(loc, 0) <= 0)
                               else min(counts[loc], max_per)),
            })
        st.dataframe(pd.DataFrame(rows), width="stretch", hide_index=True)

        if mismatches:
            st.warning(
                "⚠️ **Reference-type mismatch** — these loci already hold a *different* reference "
                "type and will be **skipped** (mixing protein + nucleotide in one locus file "
                "breaks Exonerate model selection). Clear them first in **Review library** below "
                "or in **Project Setup → Manage Data**:\n\n"
                + "  ·  ".join(f"`{loc}`: have **{ex}**, fetching **{db}**"
                              for loc, ex, db in mismatches)
            )

        if st.button("⬇️ Fetch selected loci", type="primary"):
            results: dict[str, tuple] = {}
            prog = st.progress(0.0)
            for i, loc in enumerate(selected):
                db, used, terms = dbs[loc], used_tax[loc], locus_search_terms(loc)
                ex = _existing_ref_type(loc)
                if ex is not None and ex != db:        # don't mix types in one locus file
                    results[loc] = (0, 0, f"skipped: library already has {ex} refs (not {db}) — "
                                          "clear it first")
                    prog.progress((i + 1) / len(selected))
                    continue
                n = counts.get(loc, 0)
                if n <= 0:
                    results[loc] = (0, 0, "no hits to fetch")
                    prog.progress((i + 1) / len(selected))
                    continue
                fld = _field_for(db)
                want = min(n, max_per)
                try:
                    if db == "protein":
                        # over-fetch then prefer genome-annotated RefSeq (XP_/NP_…) full-length
                        # proteins over short barcode-translation "partial" fragments — these
                        # give complete, frame-clean orthologs for protein2genome.
                        pool = search_ncbi_protein(terms, used, max_results=min(max(want * 3, 30), 100),
                                                   type_mode=type_mode, field=fld)

                        def _rank(h):
                            acc = h.get("accession", "")
                            return (acc[:3] in ("XP_", "NP_", "WP_", "YP_"), h.get("length", 0))

                        pool.sort(key=_rank, reverse=True)
                        hits = pool[:want]
                    else:
                        hits = search_ncbi_nucleotide(terms, used, max_results=want, type_mode=type_mode)
                    accs = [h["accession"] for h in hits if h.get("accession")]
                    resolved = build_entrez_query(terms, used, field=fld)
                    added, skipped, errors = fetch_and_store(
                        accs, loc, db=db,
                        query=f"{resolved} [{db}, taxon={used}, type_mode={type_mode}]",
                        ref_dir=ref_dir,
                    )
                    results[loc] = (added, skipped, "; ".join(errors))
                except Exception as e:                      # noqa: BLE001 — surface to UI
                    results[loc] = (0, 0, str(e))
                prog.progress((i + 1) / len(selected))

            update_step(
                project_dir, "references", status="done",
                outputs={loc: {"in_library": count_refs(loc, ref_dir=ref_dir),
                               "db": dbs.get(loc, ""), "taxon": used_tax.get(loc, "")}
                         for loc in selected},
                notes=f"taxon={taxon}; ref_type={ref_type}; loci={','.join(selected)}; "
                      f"type_mode={type_mode}",
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
