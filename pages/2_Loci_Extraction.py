"""
pages/2_Loci_Extraction.py
--------------------------
Extract phylogenetic marker loci from genome assemblies.

  Tab 1 — Reference Library : search NCBI, fetch and manage reference sequences
  Tab 2 — Run Extraction     : ITSx for rDNA, BLAST for protein-coding loci
  Tab 3 — Results            : recovery summary, alignment viewer, downloads

All commands go through RunManager for full logging (command, stdout/stderr,
tool versions, environment snapshot, command history).
"""

import os
import shutil
import sys
import textwrap
from pathlib import Path

import pandas as pd
import streamlit as st
from Bio import SeqIO

sys.path.insert(0, str(Path(__file__).parent.parent / "src"))
from phylofetch.blast_loci_utils import (
    detect_fasta_type,
    extract_locus,
    merge_per_strain_outputs,
    run_blast_alignment,
)
from phylofetch.config import load_config, save_config
from phylofetch.exonerate_utils import extract_locus_exonerate
from phylofetch.itsx_utils import ITSX_SUFFIXES, run_itsx
from phylofetch.ncbi_utils import (
    LOCUS_CATALOGUE,
    accessions_in_library,
    count_refs,
    delete_from_library,
    fetch_and_store,
    load_ref_meta,
    load_ref_records,
    locus_ref_fasta,
    search_ncbi_nucleotide,
    search_ncbi_protein,
    set_email,
)
from phylofetch.primer_utils import (
    PrimerPair,
    delete_user_primer,
    find_primer_amplicons,
    get_primer_catalogue,
    load_user_primers,
    locus_primer_map,
    run_primer_extraction,
    save_user_primer,
)
from phylofetch.project_manager import DEFAULT_PROJECT_DIR, RunManager

st.set_page_config(page_title="Loci Extraction", page_icon="🧬", layout="wide")
st.title("🧬 Loci Extraction")
st.caption(
    "ITSx for rDNA · Exonerate spliced alignment (frame-safe) for coding loci · "
    "Outputs: CDS, protein, genomic, introns, GFF3, codon partitions, extraction logs"
)

# ── Session / config ──────────────────────────────────────────────────────────
cfg = load_config()
if "assemblies" not in st.session_state:
    st.session_state.assemblies = cfg.get("assemblies", {})
if "tool_paths" not in st.session_state:
    st.session_state.tool_paths = cfg.get("tool_paths", {})

tp = st.session_state.tool_paths
output_base  = tp.get("output_base", str(Path.home() / "phylofetch_output"))
project_dir  = cfg.get("project_dir", str(DEFAULT_PROJECT_DIR))

ncbi_email = cfg.get("ncbi_email", "")
if ncbi_email:
    set_email(ncbi_email)

LOCI_RDNA   = ["ITS", "ITS_full", "ITS1", "ITS2", "LSU", "SSU"]
LOCI_CODING = [k for k in LOCUS_CATALOGUE if k not in ("ITS", "LSU", "SSU")]
ALL_LOCI    = ["ITS", "LSU", "SSU"] + LOCI_CODING

# ── Extraction strategy — page-level control, visible on every tab ────────────
# (Drives the whole workflow; placed above the tabs so it's never hidden.)
ref_strategy = st.radio(
    "Extraction strategy",
    ["BLAST – PCR amplicon refs (relaxed)",
     "Coding loci – Exonerate (frame-safe)",
     "PCR Primers"],
    horizontal=True,
    help=(
        "BLAST (relaxed): NCBI amplicon refs — skips the CDS completeness gate; "
        "genomic amplicon is the product. "
        "Exonerate: tblastn/blastn narrows to the best contig, then Exonerate spliced "
        "alignment (protein2genome / coding2genome) gives a frame-safe CDS with accurate "
        "exon/intron boundaries — preferred for protein-coding loci (D-008). "
        "PCR Primers: locate amplicons directly by primer binding sites (no NCBI refs needed)."
    ),
)
use_exonerate = ref_strategy.startswith("Coding loci")
require_cds = use_exonerate            # frame-safe strategy ⇒ enforce CDS gate on fallback
use_primers = ref_strategy == "PCR Primers"
if use_exonerate:
    st.info(
        "🧬 **Exonerate mode** — for each coding locus, BLAST first narrows to the best "
        "contig, then Exonerate aligns the reference *protein* (`protein2genome`) or "
        "*CDS* (`coding2genome`) across intron boundaries to recover a translatable, "
        "frame-checked CDS. rDNA loci (ITS/LSU/SSU) still go through ITSx. You can also "
        "extract an arbitrary **gene of interest** in the Run Extraction tab."
    )
if use_primers:
    st.info(
        "🔬 **PCR Primer mode** — amplicons are located directly by primer binding sites in "
        "each assembly. No NCBI reference library is needed; assign a primer pair to each "
        "locus in the Run Extraction tab."
    )

tab_refs, tab_run, tab_results = st.tabs(
    ["📚 Reference Library", "▶️ Run Extraction", "📂 Results"]
)


# ════════════════════════════════════════════════════════════════════════════
# TAB 1 — Reference Library
# ════════════════════════════════════════════════════════════════════════════
with tab_refs:
    if use_primers:
        st.info(
            "💡 **PCR Primer mode is active** — no NCBI reference library is needed. "
            "Switch to a BLAST or Exonerate strategy above to manage reference sequences, "
            "or go to the **Run Extraction** tab to assign primer pairs."
        )

    st.caption(
        "Reference sequences are used by the **BLAST** and **Exonerate** strategies. "
        "In **PCR Primer** mode this tab is not needed — assign primer pairs in Run Extraction."
    )

    with st.expander("🔑 NCBI Entrez email (required for fetching)", expanded=not bool(ncbi_email)):
        col_e, col_s = st.columns([3, 1])
        with col_e:
            new_email = st.text_input("Email", value=ncbi_email,
                                      placeholder="you@institution.edu")
        with col_s:
            st.write(""); st.write("")
            if st.button("Save email"):
                save_config({"ncbi_email": new_email})
                set_email(new_email)
                st.success("Saved.")
                st.rerun()

    if ncbi_email:
        st.success(f"✅ NCBI email: {ncbi_email}")
    else:
        st.error("⚠️ Set NCBI email above before searching.")

    st.markdown("---")

    col_locus, col_custom = st.columns([2, 1])
    with col_locus:
        selected_locus = st.selectbox("Select locus", list(LOCUS_CATALOGUE.keys()) + ["Custom…"])
    with col_custom:
        if selected_locus == "Custom…":
            locus_name = st.text_input("Custom locus name", placeholder="ACT1").strip().upper()
        else:
            locus_name = selected_locus

    if not locus_name:
        st.info("Select or enter a locus name above.")
        st.stop()

    cat      = LOCUS_CATALOGUE.get(locus_name, {})
    db_type  = cat.get("db", "nucleotide")
    sug_gene = cat.get("gene", locus_name)
    if cat.get("note"):
        st.info(cat["note"], icon="ℹ️")

    n_refs   = count_refs(locus_name)
    ref_path = locus_ref_fasta(locus_name)
    ref_type = detect_fasta_type(ref_path) if n_refs > 0 and os.path.getsize(ref_path) > 0 else "—"

    m1, m2, m3 = st.columns(3)
    m1.metric("Sequences in library", n_refs)
    m2.metric("Reference type", ref_type)
    m3.metric("Library path", ref_path.replace(str(Path.home()), "~"))

    if n_refs > 0:
        st.caption(
            "🧫 tblastn will be used." if ref_type == "protein"
            else "🧬 blastn will be used. Best with CDS-only references (no introns)."
        )

    st.markdown("---")

    with st.expander("🔍 Search NCBI and fetch sequences", expanded=True):
        col_g, col_o, col_d = st.columns([2, 2, 1])
        with col_g:
            gene_query = st.text_input("Gene / keyword", value=sug_gene)
        with col_o:
            org_query = st.text_input("Organism", value="Fungi",
                                      help="Use broad names for wider searches")
        with col_d:
            db_choice = st.selectbox("Database", ["nucleotide", "protein"],
                                     index=0 if db_type == "nucleotide" else 1)

        if db_choice == "nucleotide":
            st.caption("💡 Search for 'complete cds' entries — these are CDS sequences. "
                       "Avoid PCR amplicons that include introns.")

        max_hits = st.slider("Max results", 5, 50, 20)

        type_mode_label = st.radio(
            "Type material",
            ["Prefer type", "All", "Type only"],
            horizontal=True,
            help=(
                "Prefer type: fetch all hits, flag type-derived ones and sort them first "
                "(use this, then fall back to the best non-type hits when type material is "
                "sparse — common for closely related species). "
                "All: no preference, type still flagged. "
                "Type only: restrict to NCBI 'sequence from type material'."
            ),
        )
        type_mode = {"Prefer type": "prefer", "All": "all",
                     "Type only": "type_only"}[type_mode_label]

        if st.button("🔍 Search NCBI", type="primary"):
            if not ncbi_email:
                st.error("Set NCBI email first.")
            else:
                with st.spinner(f"Searching NCBI {db_choice}…"):
                    try:
                        fn = search_ncbi_protein if db_choice == "protein" else search_ncbi_nucleotide
                        hits = fn(gene_query, org_query, max_results=max_hits, type_mode=type_mode)
                        st.session_state[f"hits_{locus_name}"] = hits
                        st.session_state[f"searched_{locus_name}"] = True
                    except Exception as e:
                        st.error(f"NCBI error: {e}")
                        st.session_state[f"hits_{locus_name}"] = []
                        st.session_state[f"searched_{locus_name}"] = True

        hits    = st.session_state.get(f"hits_{locus_name}", [])
        searched = st.session_state.get(f"searched_{locus_name}", False)

        if searched and not hits:
            st.warning(
                f"0 results for `{gene_query}` in `{org_query}` ({db_choice}). "
                "Try a shorter search term — NCBI [Title] requires every word to appear. "
                "E.g. 'tef1 complete cds' instead of the full gene name."
            )

        if hits:
            st.success(f"{len(hits)} result(s).")
            existing = accessions_in_library(locus_name)
            df_hits  = pd.DataFrame(hits)
            if "is_type" not in df_hits.columns:
                df_hits["is_type"] = False
            df_hits["Type"] = df_hits["is_type"].apply(lambda t: "✓ type" if t else "")
            df_hits["In library"] = df_hits["accession"].apply(
                lambda a: "✓" if a in existing or a.split(".")[0] in existing else ""
            )
            df_hits["Add?"] = df_hits["In library"].apply(lambda x: x == "")

            edited = st.data_editor(
                df_hits[["Add?", "accession", "organism", "Type", "title", "length", "In library"]],
                width="stretch", hide_index=True,
                column_config={
                    "Add?": st.column_config.CheckboxColumn(),
                    "accession": st.column_config.TextColumn(disabled=True),
                    "organism":  st.column_config.TextColumn(disabled=True),
                    "Type":      st.column_config.TextColumn(
                        disabled=True,
                        help="Flagged via NCBI 'sequence from type material'. "
                             "Exact kind (holotype/ex-type/…) is recorded on fetch.",
                    ),
                    "title":     st.column_config.TextColumn(disabled=True),
                    "length":    st.column_config.NumberColumn(disabled=True),
                    "In library": st.column_config.TextColumn(disabled=True),
                },
            )
            col_fl, col_lb = st.columns([1, 2])
            with col_lb:
                use_label    = st.checkbox("Use custom label for sequence ID")
                custom_label = st.text_input("Label", placeholder="Fungus_sp_CBS123") if use_label else ""
            with col_fl:
                to_add = edited[edited["Add?"] == True]["accession"].tolist()
                st.write(f"**{len(to_add)}** selected")
                if st.button("⬇️ Fetch selected"):
                    if not to_add:
                        st.warning("None selected.")
                    else:
                        with st.spinner(f"Fetching {len(to_add)}…"):
                            added, skipped, errors = fetch_and_store(
                                to_add, locus_name, db=db_choice,
                                custom_label=custom_label if use_label and custom_label else None,
                                query=f"{gene_query} AND {org_query} [{db_choice}, type_mode={type_mode}]",
                            )
                        st.success(f"Added: {added} · Already present: {skipped}")
                        for e in errors:
                            st.warning(e)
                        st.rerun()

        st.markdown("---")
        st.markdown("**Add by accession directly:**")
        col_ma, col_md, col_mb = st.columns([3, 1, 1])
        with col_ma:
            manual_accs = st.text_input("Accession(s)", placeholder="NM_001234.1, XM_005678.2",
                                        key=f"macc_{locus_name}")
        with col_md:
            manual_db = st.selectbox("DB", ["nucleotide", "protein"], key=f"mdb_{locus_name}")
        with col_mb:
            st.write(""); st.write("")
            if st.button("Fetch", key=f"mfetch_{locus_name}") and manual_accs:
                accs = [a.strip() for a in manual_accs.split(",") if a.strip()]
                with st.spinner():
                    added, skipped, errors = fetch_and_store(accs, locus_name, db=manual_db)
                st.success(f"Added: {added} · Skipped: {skipped}")
                for e in errors:
                    st.warning(e)
                st.rerun()

    st.markdown("---")
    st.subheader(f"{locus_name} library  ({n_refs} sequences)")
    if n_refs > 0:
        records = load_ref_records(locus_name)
        meta = load_ref_meta(locus_name)

        def _meta_for(rec_id: str):
            return meta.get(rec_id) or meta.get(rec_id.split(".")[0])

        lib_rows = []
        for r in records:
            m = _meta_for(r.id)
            lib_rows.append({
                "ID": r.id,
                "Organism": m.organism if m else "",
                "Strain/voucher": (
                    (m.strain or m.culture_collection or m.specimen_voucher) if m else ""
                ),
                "Type": (m.type_kind if (m and m.is_type) else ""),
                "Length (bp/aa)": len(r.seq),
                "Remove?": False,
            })
        n_with_meta = sum(1 for row in lib_rows if row["Organism"])
        if n_with_meta < len(lib_rows):
            st.caption(
                f"ℹ️ {len(lib_rows) - n_with_meta} sequence(s) predate the metadata sidecar — "
                "re-fetch them to populate Organism / Strain / Type."
            )
        lib_df = st.data_editor(
            pd.DataFrame(lib_rows), width="stretch", hide_index=True,
            column_config={
                "Remove?": st.column_config.CheckboxColumn(),
                "ID": st.column_config.TextColumn(disabled=True),
                "Organism": st.column_config.TextColumn(disabled=True),
                "Strain/voucher": st.column_config.TextColumn(disabled=True),
                "Type": st.column_config.TextColumn(
                    disabled=True, help="Normalised /type_material kind (e.g. holotype, ex-holotype)."
                ),
                "Length (bp/aa)": st.column_config.NumberColumn(disabled=True),
            },
        )
        to_remove = lib_df[lib_df["Remove?"] == True]["ID"].tolist()
        if to_remove and st.button(f"🗑️ Remove {len(to_remove)} sequence(s)"):
            for acc in to_remove:
                delete_from_library(locus_name, acc)
            st.rerun()
        with open(locus_ref_fasta(locus_name), "rb") as f:
            st.download_button(f"⬇️ Download {locus_name}_refs.fasta",
                               data=f, file_name=f"{locus_name}_refs.fasta",
                               mime="text/plain")
    else:
        st.info("No sequences yet. Search NCBI above to add references.")


# ════════════════════════════════════════════════════════════════════════════
# TAB 2 — Run Extraction
# ════════════════════════════════════════════════════════════════════════════
with tab_run:

    if not st.session_state.assemblies:
        st.warning("No assemblies registered. Use the Assembly Manager page to import assemblies first.")
        st.stop()

    col_str, col_loc = st.columns(2)
    with col_str:
        st.markdown("**Assemblies**")
        all_ids = list(st.session_state.assemblies.keys())
        if "sel_strains" not in st.session_state:
            st.session_state["sel_strains"] = all_ids
        ba, bn, _ = st.columns([1, 1, 3])
        with ba:
            if st.button("All", key="btn_strains_all"):
                st.session_state["sel_strains"] = all_ids
                st.rerun()
        with bn:
            if st.button("None", key="btn_strains_none"):
                st.session_state["sel_strains"] = []
                st.rerun()
        sel_strains = st.multiselect("Strains", all_ids, key="sel_strains",
                                     label_visibility="collapsed")
    with col_loc:
        st.markdown("**Loci**")
        # PCR Primer mode does in-silico PCR and needs no NCBI references, so offer
        # every coding locus; BLAST / Exonerate need references, so gate on those.
        if use_primers:
            coding_opts, missing = LOCI_CODING, []
        else:
            coding_opts = [l for l in LOCI_CODING if count_refs(l) > 0]
            missing     = [l for l in LOCI_CODING if count_refs(l) == 0]
        loci_opts = ["ITS", "LSU", "SSU"] + coding_opts
        sel_loci = st.multiselect(
            "Loci", loci_opts, default=loci_opts,
            label_visibility="collapsed",
        )
        if missing:
            st.caption(f"⚠️ No refs yet for: {', '.join(missing)} "
                       "(needed for BLAST / Exonerate, not for PCR primers)")

    with st.expander("⚙️ Settings"):
        ca, cb, cc = st.columns(3)
        with ca:
            threads      = st.slider("Threads", 1, 32, 4)
            blastn_bin   = st.text_input("blastn", value=tp.get("blastn", "blastn"))
            tblastn_bin  = st.text_input("tblastn", value=tp.get("tblastn", "tblastn"))
            itsx_bin     = st.text_input("ITSx", value=tp.get("itsx", "ITSx"))
            exonerate_bin = st.text_input("exonerate", value=tp.get("exonerate", "exonerate"))
            itsx_kingdom = st.selectbox("ITSx kingdom", ["fungi", "all", "metazoa", "viridiplantae"],
                                        help="Restrict ITSx to specific HMM profiles")
        with cb:
            blastn_task = st.selectbox(
                "blastn task",
                ["dc-megablast", "megablast", "blastn"],
                disabled=use_primers or use_exonerate,
                help="dc-megablast: balanced default. blastn: most sensitive (use if exons are missed).",
            )
            min_pident  = st.number_input(
                "Min % identity", value=70, min_value=30, max_value=100,
                disabled=use_primers,
                help="BLAST identity floor (also gates the Exonerate contig-narrowing step).",
            )
            min_cds_pct = st.number_input(
                "Min CDS % of reference", value=50, min_value=10, max_value=100,
                disabled=use_exonerate or use_primers,
                help="CDS-length gate for the BLAST HSP fallback (used only if Exonerate is unavailable).",
            )
            evalue      = st.text_input("E-value", "1e-20", disabled=use_primers)
        with cc:
            st.markdown("**Tool check:**")
            exonerate_ok = shutil.which(exonerate_bin) is not None
            for label, exe in [("blastn", blastn_bin), ("tblastn", tblastn_bin),
                               ("ITSx", itsx_bin), ("exonerate", exonerate_bin)]:
                ok = shutil.which(exe) is not None
                st.write(f"{'✅' if ok else '❌'} {label}")
            if not shutil.which(blastn_bin):
                st.caption("`conda install -c bioconda blast`")
            if not shutil.which(itsx_bin):
                st.caption("`conda install -c bioconda itsx`")
            if not exonerate_ok:
                st.caption("`conda install -c bioconda exonerate`")

        # Exonerate-specific parameters (used by the frame-safe coding strategy).
        if use_exonerate:
            st.markdown("**Exonerate (spliced alignment) options:**")
            ex1, ex2, ex3 = st.columns(3)
            with ex1:
                ex_maxintron = st.number_input(
                    "Max intron (bp)", value=2000, min_value=50, max_value=200000, step=100,
                    help="Fungal introns are short (50–300 bp). Exonerate's own default (200 kb) "
                         "invites spurious giant introns; 2000 is a safe fungal ceiling.",
                )
                ex_minintron = st.number_input("Min intron (bp)", value=20, min_value=4, max_value=200)
            with ex2:
                ex_bestn = st.number_input(
                    "Report top N models", value=1, min_value=1, max_value=20,
                    help=">1 surfaces paralogs / tandem duplicates for review.",
                )
                ex_geneticcode = st.number_input(
                    "Genetic code (NCBI)", value=1, min_value=1, max_value=33,
                    help="1 = standard (fungal nuclear). Change only for organellar/alternative codes.",
                )
            with ex3:
                ex_refine = st.selectbox(
                    "Boundary refinement", ["none", "region", "full"],
                    help="'full' refines exon boundaries (slower, more accurate).",
                )
                ex_strict_qc = st.checkbox(
                    "Strict QC (reject frameshift / internal stop)", value=False,
                    help="When off, imperfect models are still written but flagged "
                         "(consistent with keeping partial / type sequences).",
                )
                ex_narrow = st.checkbox(
                    "BLAST-narrow to best contig (faster)", value=True,
                    help="Off ⇒ run Exonerate against the whole assembly (thorough, slower).",
                )
        else:
            ex_maxintron, ex_minintron, ex_bestn = 2000, 20, 1
            ex_geneticcode, ex_refine, ex_strict_qc, ex_narrow = 1, "none", False, True

    # ── Primer pair assignment (shown only in PCR Primers mode) ───────────────
    if use_primers:
        primer_cat       = get_primer_catalogue()
        primer_locus_map = locus_primer_map(primer_cat)

        st.markdown("#### 🔬 Primer pair assignment")
        st.caption(
            "Assign a primer pair to each selected locus — pick from the built-in "
            "citable catalogue, your saved library, or enter custom sequences. This is "
            "in-silico PCR: each primer pair is located in the assembly and the region "
            "between them is extracted, useful when NCBI lacks references that map cleanly."
        )
        max_mm = st.slider(
            "Max primer mismatches (edit distance per primer: substitutions + unaligned bases)",
            0, 4, 2,
            help="2 = up to 2 mismatches or truncated bases within each primer binding site.",
        )
        if "primer_assignments" not in st.session_state:
            st.session_state.primer_assignments = {}

        primer_assignments: dict[str, PrimerPair] = {}
        for locus in sel_loci:
            with st.container():
                st.markdown(f"**{locus}**")
                catalogue_options = primer_locus_map.get(locus, [])
                source_opts = ["— catalogue —"] + catalogue_options + ["Custom…"]
                prev_src = st.session_state.primer_assignments.get(
                    f"{locus}_src", source_opts[1] if catalogue_options else "Custom…")
                pa, pb = st.columns([2, 3])
                with pa:
                    src = st.selectbox(
                        "Primer pair", source_opts,
                        index=source_opts.index(prev_src) if prev_src in source_opts else 0,
                        key=f"pp_src_{locus}", label_visibility="collapsed",
                    )
                    st.session_state.primer_assignments[f"{locus}_src"] = src
                with pb:
                    if src == "Custom…":
                        c1, c2, c3, c4 = st.columns([2, 2, 1, 1])
                        fwd_seq = c1.text_input("Fwd (5'→3')", key=f"pp_fwd_{locus}", placeholder="ATGCATGC…")
                        rev_seq = c2.text_input("Rev (5'→3')", key=f"pp_rev_{locus}", placeholder="ATGCATGC…")
                        min_amp = c3.number_input("Min bp", value=100, step=50, key=f"pp_min_{locus}")
                        max_amp = c4.number_input("Max bp", value=5000, step=100, key=f"pp_max_{locus}")
                        if fwd_seq and rev_seq:
                            primer_assignments[locus] = PrimerPair(
                                name=f"custom_{locus}", locus=locus,
                                fwd=fwd_seq.strip().upper(), rev=rev_seq.strip().upper(),
                                min_amplicon=int(min_amp), max_amplicon=int(max_amp),
                                origin="user",
                            )
                            sc1, sc2 = st.columns([2, 1])
                            save_name = sc1.text_input(
                                "Save as", value=f"{locus}_custom",
                                key=f"pp_save_name_{locus}", label_visibility="collapsed",
                                placeholder="name for my library",
                            )
                            if sc2.button("💾 Save to my library", key=f"pp_save_{locus}"):
                                if save_name.strip():
                                    save_user_primer(PrimerPair(
                                        name=save_name.strip(), locus=locus,
                                        fwd=fwd_seq.strip().upper(), rev=rev_seq.strip().upper(),
                                        min_amplicon=int(min_amp), max_amplicon=int(max_amp),
                                        source="user-added", origin="user",
                                    ))
                                    st.success(f"Saved '{save_name.strip()}' to ~/.phylofetch/primers.json")
                                    st.rerun()
                                else:
                                    st.warning("Enter a name to save.")
                        else:
                            st.warning(f"Enter both primers for {locus} to enable extraction.")
                    elif src != "— catalogue —":
                        pp = primer_cat[src]
                        badge = "👤 user" if pp.origin == "user" else "📚 built-in"
                        st.caption(
                            f"{badge} · {pp.fwd_name or 'F'}: `{pp.fwd}` · "
                            f"{pp.rev_name or 'R'}: `{pp.rev}` · "
                            f"{pp.min_amplicon}–{pp.max_amplicon} bp"
                        )
                        if pp.source:
                            ref = f" · [ref]({pp.reference_url})" if pp.reference_url else ""
                            st.caption(f"📖 {pp.source}{ref}")
                        primer_assignments[locus] = pp

        # ── Custom locus (primer pair for a locus not in the catalogue) ───────
        if "primer_custom_loci" not in st.session_state:
            st.session_state.primer_custom_loci = {}   # name -> primer record dict
        with st.expander("➕ Custom locus — primer pair for a locus not in the catalogue"):
            st.caption(
                "Define in-silico PCR for any locus by name, even one phylofetch doesn't "
                "ship (e.g. MCM7, BenA). It's extracted just like a catalogue primer and "
                "added to this run; tick the box to also keep it in your saved library."
            )
            cl1, cl2, cl3 = st.columns([2, 3, 2])
            with cl1:
                cl_name = st.text_input("Locus name", key="cl_name",
                                        placeholder="e.g. MCM7, BenA, LROR2")
            with cl2:
                cl_fwd = st.text_input("Fwd primer (5'→3')", key="cl_fwd",
                                       placeholder="ACIMGIGTITC…  (IUPAC degenerate OK)")
                cl_rev = st.text_input("Rev primer (5'→3')", key="cl_rev",
                                       placeholder="GAYTTDGCIAC…")
            with cl3:
                cl_min = st.number_input("Min bp", value=100, step=50, key="cl_min")
                cl_max = st.number_input("Max bp", value=5000, step=100, key="cl_max")
            cl_save = st.checkbox("Also save to my library (~/.phylofetch/primers.json)",
                                  key="cl_save")
            if st.button("➕ Add custom locus", key="cl_add"):
                name = (cl_name or "").strip().upper().replace(" ", "_")
                if not name:
                    st.warning("Enter a locus name.")
                elif not (cl_fwd.strip() and cl_rev.strip()):
                    st.warning("Enter both forward and reverse primers.")
                else:
                    st.session_state.primer_custom_loci[name] = {
                        "locus": name,
                        "fwd": cl_fwd.strip().upper(), "rev": cl_rev.strip().upper(),
                        "min_amplicon": int(cl_min), "max_amplicon": int(cl_max),
                        "source": "user-added",
                    }
                    if cl_save:
                        save_user_primer(PrimerPair(
                            name=f"{name}_custom", locus=name,
                            fwd=cl_fwd.strip().upper(), rev=cl_rev.strip().upper(),
                            min_amplicon=int(cl_min), max_amplicon=int(cl_max),
                            source="user-added", origin="user"))
                        st.success(f"Added '{name}' and saved it to ~/.phylofetch/primers.json")
                    else:
                        st.success(f"Added custom locus '{name}' for this run.")
                    st.rerun()
            for cname, rec in list(st.session_state.primer_custom_loci.items()):
                rcx1, rcx2 = st.columns([5, 1])
                rcx1.caption(f"**{cname}** · F:`{rec['fwd']}` R:`{rec['rev']}` · "
                             f"{rec['min_amplicon']}–{rec['max_amplicon']} bp")
                if rcx2.button("🗑️ Remove", key=f"cl_del_{cname}"):
                    del st.session_state.primer_custom_loci[cname]
                    st.rerun()

        # Merge custom loci into the assignments that actually run (so the preview,
        # run loop and combine step all pick them up like any catalogue locus).
        for cname, rec in st.session_state.primer_custom_loci.items():
            primer_assignments[cname] = PrimerPair(
                name=f"{cname}_custom", locus=cname,
                fwd=rec["fwd"], rev=rec["rev"],
                min_amplicon=int(rec["min_amplicon"]), max_amplicon=int(rec["max_amplicon"]),
                source=rec.get("source", "user-added"), origin="user",
            )

        # Manage saved (user) primers
        user_lib = load_user_primers()
        if user_lib:
            with st.expander(f"📁 My primer library — {len(user_lib)} saved pair(s)"):
                for uname, upp in user_lib.items():
                    uc1, uc2 = st.columns([5, 1])
                    uc1.caption(
                        f"**{uname}** · {upp.locus} · F:`{upp.fwd}` R:`{upp.rev}` · "
                        f"{upp.min_amplicon}–{upp.max_amplicon} bp"
                    )
                    if uc2.button("🗑️ Delete", key=f"del_user_{uname}"):
                        delete_user_primer(uname)
                        st.rerun()

        # Optional: preview binding sites and disambiguate off-target hits
        if primer_assignments:
            with st.expander("🔍 Preview & choose binding sites (handle off-target hits)"):
                st.caption(
                    "Scan assemblies for every primer binding site, then pick which "
                    "amplicon to extract per strain/locus. If you skip this, the "
                    "lowest-edit-distance site is used automatically."
                )
                if st.button("Scan binding sites", key="scan_primers"):
                    scan_mgr = RunManager(project_dir)
                    scan: dict[str, list] = {}
                    jobs = [(s, l) for s in sel_strains for l in primer_assignments]
                    sprog = st.progress(0)
                    for i, (s, l) in enumerate(jobs, 1):
                        asm = st.session_state.assemblies[s]["assembly_path"]
                        try:
                            scan[f"{s}|{l}"] = find_primer_amplicons(
                                asm, primer_assignments[l],
                                max_mismatches=max_mm, blastn_bin=blastn_bin,
                                manager=scan_mgr, action=f"primer_scan_{l}_{s}",
                            )
                        except ValueError as exc:
                            scan[f"{s}|{l}"] = []
                            st.warning(f"`{s}` · {l}: {exc}")
                        sprog.progress(i / max(len(jobs), 1))
                    st.session_state["primer_scan"] = scan

                scan = st.session_state.get("primer_scan", {})
                for s in sel_strains:
                    for l in primer_assignments:
                        key = f"{s}|{l}"
                        cands = scan.get(key)
                        if cands is None:
                            continue
                        st.markdown(f"**`{s}` · {l}** — {len(cands)} site(s)")
                        if not cands:
                            st.caption("⚠️ No binding sites within the size window.")
                            continue
                        st.dataframe(pd.DataFrame([{
                            "#": i, "contig": c["contig"], "strand": c["strand"],
                            "coords": f"{c['amp_start']}..{c['amp_end']}",
                            "len(bp)": c["amp_len"], "fwd_mm": c["fwd_edit"],
                            "rev_mm": c["rev_edit"], "total_mm": c["total_edit"],
                        } for i, c in enumerate(cands)]),
                            hide_index=True, width="stretch")
                        opts = [
                            f"{i}: {c['contig']} {c['strand']} "
                            f"{c['amp_start']}..{c['amp_end']} ({c['amp_len']} bp, mm {c['total_edit']})"
                            for i, c in enumerate(cands)
                        ]
                        choice = st.selectbox(
                            f"Extract which site? ({s} / {l})", opts, index=0,
                            key=f"primer_choice_{key}",
                        )
                        st.session_state[f"primer_pick_{key}"] = opts.index(choice)

    # ── Gene of interest (ad-hoc ortholog → CDS) — Exonerate strategy only ─────
    goi_genes: dict[str, str] = {}
    if use_exonerate:
        goi_dir = os.path.join(output_base, "loci", "_goi_refs")
        if "goi_genes" not in st.session_state:
            st.session_state.goi_genes = {}
        with st.expander("🔎 Gene of interest — extract any ortholog's CDS", expanded=False):
            st.caption(
                "Paste or upload a reference **protein** or **CDS** for any gene of interest "
                "(an ortholog from a related species). Exonerate locates it in each selected "
                "assembly and returns just the CDS + exon model — no NCBI reference library needed."
            )
            gc1, gc2 = st.columns([1, 2])
            with gc1:
                goi_name = st.text_input("Gene name (used for output / tip labels)",
                                         key="goi_name", placeholder="e.g. CAL, mcm7, betaTUB")
                goi_file = st.file_uploader("…or upload FASTA",
                                            type=["fa", "fasta", "faa", "fna", "txt"],
                                            key="goi_file")
            with gc2:
                goi_text = st.text_area("Reference FASTA (protein or nucleotide CDS)",
                                        key="goi_text", height=120,
                                        placeholder=">ref_ortholog\nMSEQ… (protein) or ATG… (CDS)")
            if st.button("➕ Add gene of interest", key="goi_add"):
                name = (goi_name or "").strip().replace(" ", "_")
                content = ""
                if goi_file is not None:
                    content = goi_file.getvalue().decode("utf-8", "replace")
                elif goi_text.strip():
                    content = goi_text
                if not name:
                    st.warning("Enter a gene name.")
                elif not content.strip():
                    st.warning("Paste or upload a reference FASTA.")
                else:
                    if not content.lstrip().startswith(">"):
                        content = f">{name}_ref\n{content.strip()}\n"
                    os.makedirs(goi_dir, exist_ok=True)
                    ref_path = os.path.join(goi_dir, f"{name}.fasta")
                    Path(ref_path).write_text(content)
                    st.session_state.goi_genes[name] = ref_path
                    st.success(f"Added gene of interest '{name}'.")
                    st.rerun()
            for gname, gpath in list(st.session_state.goi_genes.items()):
                rc1, rc2 = st.columns([5, 1])
                exists = os.path.exists(gpath)
                rtype = detect_fasta_type(gpath) if exists else "?"
                rc1.caption(f"**{gname}** · {rtype} ref · `{gpath}`"
                            + ("" if exists else " · ⚠️ file missing"))
                if rc2.button("🗑️ Remove", key=f"goi_del_{gname}"):
                    del st.session_state.goi_genes[gname]
                    st.rerun()
        goi_genes = {n: p for n, p in st.session_state.goi_genes.items() if os.path.exists(p)}

    have_custom_primers = use_primers and bool(st.session_state.get("primer_custom_loci"))
    if not sel_strains or (not sel_loci and not goi_genes and not have_custom_primers):
        st.info("Select strains and loci above (or add a gene of interest / custom primer locus).")
        st.stop()

    rdna_run   = [l for l in sel_loci if l in ("ITS", "LSU", "SSU", "ITS1", "ITS2", "ITS_full")]
    coding_run = [l for l in sel_loci if l not in rdna_run]

    # Primer mode: show assignment coverage warning before run
    if use_primers:
        unassigned = [l for l in sel_loci if l not in primer_assignments]
        assigned   = list(primer_assignments)   # incl. custom loci not in sel_loci
        if unassigned:
            st.warning(f"No primer pair assigned for: {', '.join(unassigned)}. Those loci will be skipped.")
        st.markdown(
            f"**Ready:** {len(sel_strains)} strain(s) · "
            f"PCR Primers: {', '.join(assigned) or '—'} · "
            f"Max mismatches: {max_mm}"
        )
    else:
        coding_engine = "Exonerate" if use_exonerate else "BLAST"
        if use_exonerate and not exonerate_ok:
            st.warning(
                "⚠️ **Exonerate not found on PATH** — coding loci will fall back to the BLAST "
                "HSP-as-exon path, which does **not** validate the reading frame (a 1–2 bp "
                "boundary error can frameshift the CDS silently). Install for frame-safe CDS: "
                "`conda install -c bioconda exonerate`."
            )
            coding_engine = "BLAST HSP fallback ⚠️"
        st.markdown(
            f"**Ready:** {len(sel_strains)} strain(s) · "
            f"rDNA (ITSx): {', '.join(rdna_run) or '—'} · "
            f"CDS ({coding_engine}): {', '.join(coding_run) or '—'}"
            + (f" · Gene(s) of interest: {', '.join(goi_genes)}" if goi_genes else "")
        )

    if st.button("🚀 Run extraction", type="primary"):
        loci_dir       = os.path.join(output_base, "loci")
        per_strain_dir = os.path.join(loci_dir, "per_strain")
        combined_dir   = os.path.join(loci_dir, "combined")
        os.makedirs(per_strain_dir, exist_ok=True)
        os.makedirs(combined_dir,   exist_ok=True)

        evalue_float = 1e-20
        if not use_primers:
            try:
                evalue_float = float(evalue)
            except ValueError:
                st.error(f"Invalid e-value: {evalue}")
                st.stop()

        manager = RunManager(project_dir)

        active_loci = (list(primer_assignments) if use_primers
                       else list(sel_loci) + list(goi_genes))
        total_jobs  = len(sel_strains) * max(len(active_loci), 1)
        prog  = st.progress(0)
        job_n = 0

        # Shared coding-locus extraction: Exonerate (frame-safe) when available,
        # else the BLAST HSP-as-exon fallback. Used for both catalogue loci and
        # ad-hoc genes of interest, so they share one logged code path.
        def _extract_coding(strain_id, assembly, locus, ref_fa, locus_out):
            os.makedirs(locus_out, exist_ok=True)
            if use_exonerate and exonerate_ok:
                _rid, rdir = manager.dry_run(
                    [exonerate_bin, "--model", "(auto)", "--query", ref_fa, "--target", assembly],
                    module="loci_extraction", action=f"exonerate_{locus}_{strain_id}",
                    inputs={"assembly": assembly, "reference": ref_fa},
                    outputs={"locus_dir": locus_out},
                    params={"min_pident": min_pident, "maxintron": int(ex_maxintron),
                            "bestn": int(ex_bestn), "narrow": ex_narrow,
                            "refine": ex_refine, "geneticcode": int(ex_geneticcode)},
                )
                with st.spinner(f"[{strain_id}] Exonerate: {locus}…"):
                    return extract_locus_exonerate(
                        assembly_fasta=assembly, reference_fasta=ref_fa,
                        output_dir=locus_out, strain_id=strain_id, locus_name=locus,
                        exonerate_bin=exonerate_bin, blastn_bin=blastn_bin, tblastn_bin=tblastn_bin,
                        narrow=ex_narrow, minintron=int(ex_minintron), maxintron=int(ex_maxintron),
                        bestn=int(ex_bestn), refine=ex_refine, geneticcode=int(ex_geneticcode),
                        min_pident=float(min_pident), evalue=evalue_float, threads=threads,
                        strict_qc=ex_strict_qc, manager=manager, run_dir=str(rdir),
                    )
            # BLAST HSP-as-exon fallback (no reading-frame guarantee; see D-008).
            cmd_parts = [
                blastn_bin if detect_fasta_type(ref_fa) == "nucleotide" else tblastn_bin,
                "-query", ref_fa, "-subject", assembly, "-evalue", str(evalue_float),
            ]
            rr = manager.dry_run(
                cmd_parts, module="loci_extraction", action=f"{locus}_{strain_id}",
                inputs={"assembly": assembly, "reference": ref_fa},
                outputs={"locus_dir": locus_out},
                params={"evalue": evalue_float, "min_pident": min_pident,
                        "min_cds_pct": min_cds_pct, "blastn_task": blastn_task},
            )
            with st.spinner(f"[{strain_id}] BLAST: {locus}…"):
                return extract_locus(
                    assembly_fasta=assembly, reference_fasta=ref_fa, output_dir=locus_out,
                    strain_id=strain_id, locus_name=locus,
                    min_pident=float(min_pident), min_cds_pct_of_ref=float(min_cds_pct),
                    evalue=evalue_float, blastn_task=blastn_task, threads=threads,
                    blastn_bin=blastn_bin, tblastn_bin=tblastn_bin, run_dir=str(rr[1]),
                    require_complete_cds=require_cds,
                )

        def _report_coding(locus, result, status, goi=False):
            label = f"🔎 {locus}" if goi else locus
            if not result:
                st.write(f"⚠️ {label}: {status}")
                return
            n_introns  = result["n_introns"]
            splices_ok = sum(1 for intr in result["introns"]
                             if intr["splice_5"] == "GT" and intr["splice_3"] == "AG")
            splice_note = f" · GT-AG: {splices_ok}/{n_introns}" if n_introns > 0 else ""
            if result.get("tool") == "exonerate":
                stops, frame = result.get("n_internal_stops", 0), result.get("len_mod3", 0)
                qc = ("✅ frame OK" if (stops == 0 and frame == 0)
                      else f"⚠️ QC review: {stops} internal stop(s), len%3={frame}")
                model_label = result.get("blast_type", "").replace("exonerate:", "").split(":")[0]
                paralog = (f" · {result['n_other_models']} paralog cand(s)"
                           if result.get("n_other_models") else "")
                st.write(
                    f"✅ {label} (Exonerate {model_label or 'protein2genome'}): "
                    f"{result['n_exons']} exon(s), {result['cds_length']} bp CDS, "
                    f"{n_introns} intron(s){splice_note} · {qc} · "
                    f"id {result.get('query_pident', 0):.0f}% cov {result.get('query_coverage', 0):.0f}%"
                    f"{paralog} · ref `{result.get('ref_accession', '?')}`"
                )
            else:
                st.write(
                    f"✅ {label} ({result['blast_type']}): "
                    f"{result['n_exons']} exon(s), {result['cds_length']} bp CDS, "
                    f"{n_introns} intron(s){splice_note} · ref `{result.get('ref_accession', '?')}`"
                )

        for strain_id in sel_strains:
            assembly   = st.session_state.assemblies[strain_id]["assembly_path"]
            strain_out = os.path.join(per_strain_dir, strain_id)
            os.makedirs(strain_out, exist_ok=True)
            st.markdown(f"#### `{strain_id}`")

            # ── PCR Primer strategy ────────────────────────────────────────
            if use_primers:
                scan = st.session_state.get("primer_scan", {})
                for locus, pp in primer_assignments.items():
                    locus_out = os.path.join(strain_out, locus)
                    os.makedirs(locus_out, exist_ok=True)
                    # Use a pre-scanned + user-chosen binding site if available.
                    # Fall back to a fresh search when nothing was scanned.
                    scan_key  = f"{strain_id}|{locus}"
                    chosen    = scan.get(scan_key) or None
                    pick_idx  = st.session_state.get(f"primer_pick_{scan_key}", 0)
                    with st.spinner(f"[{strain_id}] Primer search: {locus} ({pp.name})…"):
                        result, status = run_primer_extraction(
                            assembly_fasta=assembly,
                            primer_pair=pp,
                            output_dir=locus_out,
                            strain_id=strain_id,
                            locus_name=locus,
                            max_mismatches=max_mm,
                            blastn_bin=blastn_bin,
                            manager=manager,
                            candidates=chosen,
                            chosen_index=pick_idx,
                        )
                    if result:
                        mm_note = (
                            f"fwd_mm={result['fwd_edit']}, rev_mm={result['rev_edit']}"
                        )
                        candidates_note = (
                            f" · {result['n_candidates']} site(s) found"
                            if result["n_candidates"] > 1 else ""
                        )
                        st.write(
                            f"✅ {locus} ({pp.name}): {result['amp_len']} bp amplicon "
                            f"on `{result['contig']}` [{result['strand']}] "
                            f"coords {result['amp_start']}..{result['amp_end']} "
                            f"· {mm_note}{candidates_note}"
                        )
                    else:
                        st.write(f"⚠️ {locus}: {status}")
                    job_n += 1; prog.progress(job_n / total_jobs)
                # skip BLAST/ITSx blocks below when in primer mode
                prog.progress(1.0)
                continue

            # ── ITSx for rDNA ─────────────────────────────────────────────
            if rdna_run:
                with st.spinner(f"[{strain_id}] ITSx…"):
                    itsx_tmp = os.path.join(strain_out, "_itsx_tmp")
                    rc, log, found = run_itsx(
                        assembly, itsx_tmp, strain_id,
                        threads=threads, itsx_bin=itsx_bin, kingdom=itsx_kingdom,
                    )

                if rc != 0:
                    st.error(f"ITSx failed (exit {rc}):")
                    st.code(log, language=None)
                elif not found:
                    st.warning("ITSx: no rDNA regions detected.")
                    if log:
                        with st.expander("ITSx log"):
                            st.code(log, language=None)

                for locus in rdna_run:
                    locus_dir_out = os.path.join(strain_out, locus)
                    os.makedirs(locus_dir_out, exist_ok=True)
                    key = "ITS_full" if locus == "ITS" else locus
                    if key in found:
                        import shutil as _sh
                        dest = os.path.join(locus_dir_out, f"{locus}.fasta")
                        _sh.copy(found[key], dest)
                        recs = list(SeqIO.parse(dest, "fasta"))
                        st.write(f"✅ {locus}: {len(recs)} sequence(s)")
                    else:
                        st.write(f"⚠️ {locus}: not detected by ITSx")
                    job_n += 1; prog.progress(job_n / total_jobs)

            # ── Coding loci: Exonerate (frame-safe) or BLAST HSP fallback ──
            for locus in coding_run:
                ref_fa = locus_ref_fasta(locus)
                if not os.path.exists(ref_fa) or os.path.getsize(ref_fa) == 0:
                    st.write(f"⚠️ {locus}: no references in library, skipping")
                    job_n += 1; prog.progress(job_n / total_jobs)
                    continue
                locus_out = os.path.join(strain_out, locus)
                result, status = _extract_coding(strain_id, assembly, locus, ref_fa, locus_out)
                _report_coding(locus, result, status)
                job_n += 1; prog.progress(job_n / total_jobs)

            # ── Gene(s) of interest (Exonerate, ad-hoc ortholog → CDS) ─────
            for gname, gref in goi_genes.items():
                locus_out = os.path.join(strain_out, gname)
                result, status = _extract_coding(strain_id, assembly, gname, gref, locus_out)
                _report_coding(gname, result, status, goi=True)
                job_n += 1; prog.progress(job_n / total_jobs)

        # ── Merge combined files ───────────────────────────────────────────
        st.markdown("**Merging into combined multi-FASTAs…**")
        merge_loci = (list(primer_assignments) if use_primers
                      else list(sel_loci) + list(goi_genes))
        for locus in merge_loci:
            # Primer-based amplicons
            if use_primers:
                recs = []
                for sd in sorted(Path(per_strain_dir).iterdir()):
                    fp = sd / locus / f"{locus}_amplicon.fasta"
                    if fp.exists():
                        recs.extend(list(SeqIO.parse(str(fp), "fasta")))
                if recs:
                    out = os.path.join(combined_dir, f"{locus}_amplicon_combined.fasta")
                    SeqIO.write(recs, out, "fasta")
                    st.write(f"✅ {locus}_amplicon_combined.fasta — {len(recs)} sequence(s)")
            elif locus in ("ITS", "LSU", "SSU", "ITS1", "ITS2", "ITS_full"):
                recs = []
                for sd in sorted(Path(per_strain_dir).iterdir()):
                    fp = sd / locus / f"{locus}.fasta"
                    if fp.exists():
                        recs.extend(list(SeqIO.parse(str(fp), "fasta")))
                if recs:
                    out = os.path.join(combined_dir, f"{locus}_combined.fasta")
                    SeqIO.write(recs, out, "fasta")
                    st.write(f"✅ {locus}_combined.fasta — {len(recs)} sequence(s)")
            else:
                combined = merge_per_strain_outputs(per_strain_dir, combined_dir, locus)
                for kind, path in combined.items():
                    n = sum(1 for _ in SeqIO.parse(path, "fasta"))
                    st.write(f"✅ {locus}_{kind}_combined.fasta — {n} sequence(s)")

        prog.progress(1.0)
        st.success("Extraction complete. See **Results** tab.")


# ════════════════════════════════════════════════════════════════════════════
# TAB 3 — Results
# ════════════════════════════════════════════════════════════════════════════
with tab_results:

    loci_dir       = os.path.join(output_base, "loci")
    per_strain_dir = os.path.join(loci_dir, "per_strain")
    combined_dir   = os.path.join(loci_dir, "combined")

    if not os.path.isdir(loci_dir):
        st.info("No results yet. Run extraction first.")
        st.stop()

    # ── Recovery summary ──────────────────────────────────────────────────────
    st.subheader("Recovery summary")
    st.caption("CDS length (bp) per strain per locus. Green = found · Red = missing")

    strain_dirs = sorted(
        [d for d in Path(per_strain_dir).iterdir() if d.is_dir()],
        key=lambda d: d.name,
    ) if os.path.isdir(per_strain_dir) else []

    if strain_dirs:
        rows = []
        for sd in strain_dirs:
            row = {"Strain": sd.name}
            for locus in ALL_LOCI:
                length = 0
                lsd = sd / locus
                if lsd.is_dir():
                    for fname in [f"{locus}_CDS.fasta", f"{locus}_amplicon.fasta",
                                  f"{locus}.fasta"]:
                        fp = lsd / fname
                        if fp.exists():
                            recs = list(SeqIO.parse(str(fp), "fasta"))
                            length = len(recs[0].seq) if recs else 0
                            break
                row[locus] = length
            rows.append(row)

        df_sum = pd.DataFrame(rows).set_index("Strain")
        df_sum = df_sum.loc[:, (df_sum > 0).any()]

        def _recovery_cell(v) -> str:
            # Pastel background + dark text: stays legible on both light and dark
            # themes (a bare background color renders as white-on-pastel in dark mode).
            if v > 0:
                return "background-color:#c8e6c9; color:#1b5e20; font-weight:600"
            return "background-color:#ffcdd2; color:#b71c1c; font-weight:600"

        st.dataframe(
            df_sum.style.map(_recovery_cell),
            width="stretch",
        )

    st.markdown("---")

    # ── Per-strain results browser ─────────────────────────────────────────────
    st.subheader("Per-strain results")
    strain_names = [d.name for d in strain_dirs]
    if not strain_names:
        st.info("No per-strain results.")
        st.stop()

    browse_strain = st.selectbox("Strain", strain_names, key="results_strain")
    if browse_strain:
        strain_path = Path(per_strain_dir) / browse_strain
        locus_subdirs = sorted(
            [d for d in strain_path.iterdir() if d.is_dir() and not d.name.startswith("_")]
        )

        for locus_dir in locus_subdirs:
            locus = locus_dir.name
            with st.expander(f"**{locus}**", expanded=False):

                # CDS summary
                cds_fp = locus_dir / f"{locus}_CDS.fasta"
                if not cds_fp.exists():
                    cds_fp = locus_dir / f"{locus}_amplicon.fasta"
                if not cds_fp.exists():
                    cds_fp = locus_dir / f"{locus}.fasta"
                if cds_fp.exists():
                    recs = list(SeqIO.parse(str(cds_fp), "fasta"))
                    if recs:
                        r = recs[0]
                        st.markdown(f"**Sequence:** `{r.id}` — {len(r.seq)} bp")
                        # Show provenance from header
                        desc_fields = {}
                        for token in r.description.split("["):
                            token = token.strip().rstrip("]")
                            if "=" in token:
                                k, v = token.split("=", 1)
                                desc_fields[k.strip()] = v.strip()
                        if desc_fields:
                            meta_rows = [{"Field": k, "Value": v}
                                         for k, v in desc_fields.items()]
                            st.dataframe(pd.DataFrame(meta_rows),
                                         width="stretch", hide_index=True)

                        # Frame / QC badge (Exonerate path records internal_stops)
                        stops = desc_fields.get("internal_stops")
                        if stops is not None:
                            try:
                                nstop = int(stops)
                            except ValueError:
                                nstop = -1
                            if nstop == 0:
                                st.success("✅ CDS QC — in-frame, no internal stop codons.")
                            elif nstop > 0:
                                st.warning(
                                    f"⚠️ CDS QC — {nstop} internal stop codon(s); review "
                                    "(possible frameshift, pseudogene, or wrong genetic code)."
                                )

                        # Translated protein (Exonerate path)
                        prot_fp = locus_dir / f"{locus}_protein.fasta"
                        if prot_fp.exists():
                            precs = list(SeqIO.parse(str(prot_fp), "fasta"))
                            if precs:
                                pseq = str(precs[0].seq)
                                with st.expander(
                                    f"🧬 Translated protein — {len(pseq.rstrip('*'))} aa"
                                ):
                                    st.code("\n".join(textwrap.wrap(pseq, 60)), language=None)

                # Exon / intron structure from GFF3
                gff_path = locus_dir / f"{locus}.gff3"
                if gff_path.exists():
                    exon_rows = []
                    with open(gff_path) as f:
                        for line in f:
                            if line.startswith("#"):
                                continue
                            parts = line.strip().split("\t")
                            if len(parts) >= 8 and parts[2] == "CDS":
                                exon_rows.append({
                                    "Start":  int(parts[3]), "End": int(parts[4]),
                                    "Length": int(parts[4]) - int(parts[3]) + 1,
                                    "Strand": parts[6], "Phase": parts[7],
                                    "% Identity": parts[5],
                                })
                    if exon_rows:
                        st.markdown("**Exon structure:**")
                        st.dataframe(pd.DataFrame(exon_rows),
                                     width="stretch", hide_index=True)

                # Intron splice sites
                intron_fp = locus_dir / f"{locus}_introns.fasta"
                if intron_fp.exists():
                    intron_recs = list(SeqIO.parse(str(intron_fp), "fasta"))
                    if intron_recs:
                        st.markdown("**Introns:**")
                        for intr in intron_recs:
                            desc_fields = {}
                            for token in intr.description.split("["):
                                token = token.strip().rstrip("]")
                                if "=" in token:
                                    k, v = token.split("=", 1)
                                    desc_fields[k.strip()] = v.strip()
                            splice_5 = desc_fields.get("splice5", "")
                            splice_3 = desc_fields.get("splice3", "")
                            canonical = splice_5 == "GT" and splice_3 == "AG"
                            icon = "✅" if canonical else "⚠️"
                            st.write(
                                f"{icon} `{intr.id}` — {len(intr.seq)} bp "
                                f"({splice_5}…{splice_3})"
                            )

                # ── Alignment viewer: extracted sequence vs reference ──────
                st.markdown("---")
                st.markdown("**Alignment vs reference query:**")
                blast_aln_txt = locus_dir / "blast_alignment.txt"

                col_aln1, col_aln2 = st.columns([1, 1])
                with col_aln1:
                    if st.button(f"Show BLAST alignment ({locus})",
                                 key=f"btn_aln_{browse_strain}_{locus}"):
                        ref_fa = locus_ref_fasta(locus)
                        if not os.path.exists(ref_fa) or os.path.getsize(ref_fa) == 0:
                            st.warning("No reference library for this locus.")
                        else:
                            assembly_path = (
                                st.session_state.assemblies.get(browse_strain, {})
                                .get("assembly_path", "")
                            )
                            if not assembly_path:
                                st.warning("Assembly path not found in session.")
                            else:
                                blast_bin = tp.get("blastn", "blastn")
                                with st.spinner("Running BLAST alignment…"):
                                    rc, err, aln_path = run_blast_alignment(
                                        ref_fa, assembly_path,
                                        str(locus_dir),
                                        blast_bin=blast_bin,
                                        evalue=1e-10,
                                    )
                                if rc == 0 and Path(aln_path).exists():
                                    st.session_state[f"alntext_{browse_strain}_{locus}"] = \
                                        Path(aln_path).read_text()
                                else:
                                    st.error(f"BLAST failed: {err[:200]}")

                aln_text = st.session_state.get(f"alntext_{browse_strain}_{locus}", "")
                if aln_text or blast_aln_txt.exists():
                    text = aln_text or blast_aln_txt.read_text()
                    with st.expander("📄 Pairwise alignment (query vs assembly)",
                                     expanded=bool(aln_text)):
                        st.code(text[:8000]
                                + ("\n… (truncated)" if len(text) > 8000 else ""),
                                language=None)

                # Extraction log
                log_fp = locus_dir / f"{locus}_extraction.log"
                if log_fp.exists():
                    with st.expander("📋 Extraction log"):
                        st.code(log_fp.read_text(), language=None)

                # File downloads
                st.markdown("**Files:**")
                for fp in sorted(locus_dir.glob("*")):
                    if (fp.is_file() and not fp.name.startswith("_")
                            and fp.suffix in (".fasta", ".gff3", ".nex", ".log")):
                        with open(fp, "rb") as f:
                            st.download_button(
                                f"⬇️ {fp.name}", data=f,
                                file_name=fp.name, mime="text/plain",
                                key=f"dl_{browse_strain}_{locus}_{fp.name}",
                            )

    st.markdown("---")

    # ── Combined FASTA downloads ──────────────────────────────────────────────
    st.subheader("Combined multi-FASTA downloads")
    combined_files = sorted(Path(combined_dir).glob("*_combined.fasta")) \
        if os.path.isdir(combined_dir) else []

    if not combined_files:
        st.info("No combined files yet.")
    else:
        by_locus: dict = {}
        for cf in combined_files:
            stem = cf.stem.replace("_combined", "")
            for kind in ["_CDS", "_protein", "_genomic", "_introns"]:
                if stem.endswith(kind):
                    locus = stem[:-len(kind)]
                    by_locus.setdefault(locus, {})[kind.strip("_")] = cf
                    break
            else:
                by_locus.setdefault(stem, {})["rDNA"] = cf

        for locus, kinds in sorted(by_locus.items()):
            with st.expander(f"**{locus}** — {len(kinds)} file(s)", expanded=True):
                cols = st.columns(len(kinds))
                for col, (kind, fp) in zip(cols, sorted(kinds.items())):
                    n = sum(1 for _ in SeqIO.parse(str(fp), "fasta"))
                    with col:
                        st.markdown(f"**{kind}**")
                        st.caption(f"{n} sequences · {fp.stat().st_size:,} bytes")
                        with open(fp, "rb") as f:
                            st.download_button(
                                f"⬇️ {fp.name}", data=f,
                                file_name=fp.name, mime="text/plain",
                                key=f"dlc_{fp.name}",
                            )

    st.markdown("---")
    with st.expander("📋 Recommended next steps"):
        st.markdown("""
**After extraction — send combined FASTAs to Alignment Prep (page 3):**

```bash
# Manual workflow if needed:
mafft --auto TEF1_CDS_combined.fasta > TEF1_aligned.fasta
trimal -in TEF1_aligned.fasta -out TEF1_trimmed.fasta -automated1
```

Use **CDS** combined FASTAs for phylogenies (intron-free, codon-aware alignment).
Use **genomic** combined FASTAs for comparison to GenBank PCR amplicon entries.
The `_partition.nex` files contain codon-position blocks for IQ-TREE2 `-p` partitioned analysis.
        """)
