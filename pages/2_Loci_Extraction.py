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
from phylofetch.itsx_utils import ITSX_SUFFIXES, run_itsx
from phylofetch.ncbi_utils import (
    LOCUS_CATALOGUE,
    accessions_in_library,
    count_refs,
    delete_from_library,
    fetch_and_store,
    load_ref_records,
    locus_ref_fasta,
    search_ncbi_nucleotide,
    search_ncbi_protein,
    set_email,
)
from phylofetch.project_manager import DEFAULT_PROJECT_DIR, RunManager

st.set_page_config(page_title="Loci Extraction", page_icon="🧬", layout="wide")
st.title("🧬 Loci Extraction")
st.caption(
    "ITSx for rDNA · BLAST HSP-as-exon for coding loci · "
    "Outputs: CDS, genomic, introns, GFF3, codon partitions, extraction logs"
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

tab_refs, tab_run, tab_results = st.tabs(
    ["📚 Reference Library", "▶️ Run Extraction", "📂 Results"]
)


# ════════════════════════════════════════════════════════════════════════════
# TAB 1 — Reference Library
# ════════════════════════════════════════════════════════════════════════════
with tab_refs:

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

        if st.button("🔍 Search NCBI", type="primary"):
            if not ncbi_email:
                st.error("Set NCBI email first.")
            else:
                with st.spinner(f"Searching NCBI {db_choice}…"):
                    try:
                        fn = search_ncbi_protein if db_choice == "protein" else search_ncbi_nucleotide
                        hits = fn(gene_query, org_query, max_results=max_hits)
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
            df_hits["In library"] = df_hits["accession"].apply(
                lambda a: "✓" if a in existing or a.split(".")[0] in existing else ""
            )
            df_hits["Add?"] = df_hits["In library"].apply(lambda x: x == "")

            edited = st.data_editor(
                df_hits[["Add?", "accession", "organism", "title", "length", "In library"]],
                use_container_width=True, hide_index=True,
                column_config={
                    "Add?": st.column_config.CheckboxColumn(),
                    "accession": st.column_config.TextColumn(disabled=True),
                    "organism":  st.column_config.TextColumn(disabled=True),
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
        lib_rows = [{"ID": r.id,
                     "Description": (r.description.replace(r.id, "").strip())[:80],
                     "Length (bp/aa)": len(r.seq),
                     "Remove?": False}
                    for r in records]
        lib_df = st.data_editor(
            pd.DataFrame(lib_rows), use_container_width=True, hide_index=True,
            column_config={
                "Remove?": st.column_config.CheckboxColumn(),
                "ID": st.column_config.TextColumn(disabled=True),
                "Description": st.column_config.TextColumn(disabled=True),
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
        st.warning("No assemblies registered. Use Project Setup → Import Assemblies first.")
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
        avail_coding = [l for l in LOCI_CODING if count_refs(l) > 0]
        missing      = [l for l in LOCI_CODING if count_refs(l) == 0]
        sel_loci = st.multiselect(
            "Loci", ["ITS", "LSU", "SSU"] + avail_coding,
            default=["ITS", "LSU", "SSU"] + avail_coding,
            label_visibility="collapsed",
        )
        if missing:
            st.caption(f"⚠️ No refs yet for: {', '.join(missing)}")

    with st.expander("⚙️ Settings"):
        ca, cb, cc = st.columns(3)
        with ca:
            threads     = st.slider("Threads", 1, 32, 4)
            blastn_bin  = st.text_input("blastn", value=tp.get("blastn", "blastn"))
            tblastn_bin = st.text_input("tblastn", value=tp.get("tblastn", "tblastn"))
            itsx_bin    = st.text_input("ITSx", value=tp.get("itsx", "ITSx"))
            itsx_kingdom = st.selectbox("ITSx kingdom", ["fungi", "all", "metazoa", "viridiplantae"],
                                        help="Restrict ITSx to specific HMM profiles")
        with cb:
            blastn_task = st.selectbox(
                "blastn task",
                ["dc-megablast", "megablast", "blastn"],
                help="dc-megablast: balanced default. blastn: most sensitive (use if exons are missed).",
            )
            min_pident  = st.number_input("Min % identity", value=70, min_value=30, max_value=100)
            min_cds_pct = st.number_input("Min CDS % of reference", value=50, min_value=10, max_value=100)
            evalue      = st.text_input("E-value", "1e-20")
        with cc:
            st.markdown("**Tool check:**")
            for label, exe in [("blastn", blastn_bin), ("tblastn", tblastn_bin), ("ITSx", itsx_bin)]:
                ok = shutil.which(exe) is not None
                st.write(f"{'✅' if ok else '❌'} {label}")
            if not shutil.which(blastn_bin):
                st.caption("`conda install -c bioconda blast`")
            if not shutil.which(itsx_bin):
                st.caption("`conda install -c bioconda itsx`")

    if not sel_strains or not sel_loci:
        st.info("Select strains and loci above.")
        st.stop()

    rdna_run   = [l for l in sel_loci if l in ("ITS", "LSU", "SSU", "ITS1", "ITS2", "ITS_full")]
    coding_run = [l for l in sel_loci if l not in rdna_run]

    st.markdown(
        f"**Ready:** {len(sel_strains)} strain(s) · "
        f"rDNA (ITSx): {', '.join(rdna_run) or '—'} · "
        f"CDS (BLAST): {', '.join(coding_run) or '—'}"
    )

    if st.button("🚀 Run extraction", type="primary"):
        loci_dir       = os.path.join(output_base, "loci")
        per_strain_dir = os.path.join(loci_dir, "per_strain")
        combined_dir   = os.path.join(loci_dir, "combined")
        os.makedirs(per_strain_dir, exist_ok=True)
        os.makedirs(combined_dir,   exist_ok=True)

        try:
            evalue_float = float(evalue)
        except ValueError:
            st.error(f"Invalid e-value: {evalue}")
            st.stop()

        manager = RunManager(project_dir)

        total_jobs = len(sel_strains) * len(sel_loci)
        prog  = st.progress(0)
        job_n = 0

        for strain_id in sel_strains:
            assembly   = st.session_state.assemblies[strain_id]["assembly_path"]
            strain_out = os.path.join(per_strain_dir, strain_id)
            os.makedirs(strain_out, exist_ok=True)
            st.markdown(f"#### `{strain_id}`")

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

            # ── BLAST for coding loci ──────────────────────────────────────
            for locus in coding_run:
                ref_fa = locus_ref_fasta(locus)
                if not os.path.exists(ref_fa) or os.path.getsize(ref_fa) == 0:
                    st.write(f"⚠️ {locus}: no references in library, skipping")
                    job_n += 1; prog.progress(job_n / total_jobs)
                    continue

                locus_out = os.path.join(strain_out, locus)

                # Wrap in RunManager for full logging
                cmd_parts = [
                    blastn_bin if detect_fasta_type(ref_fa) == "nucleotide" else tblastn_bin,
                    "-query", ref_fa, "-subject", assembly,
                    "-evalue", str(evalue_float),
                ]
                rr = manager.dry_run(
                    cmd_parts, module="loci_extraction", action=f"{locus}_{strain_id}",
                    inputs={"assembly": assembly, "reference": ref_fa},
                    outputs={"locus_dir": locus_out},
                    params={"evalue": evalue_float, "min_pident": min_pident,
                            "min_cds_pct": min_cds_pct, "blastn_task": blastn_task},
                )

                with st.spinner(f"[{strain_id}] BLAST: {locus}…"):
                    result, status = extract_locus(
                        assembly_fasta=assembly,
                        reference_fasta=ref_fa,
                        output_dir=locus_out,
                        strain_id=strain_id,
                        locus_name=locus,
                        min_pident=float(min_pident),
                        min_cds_pct_of_ref=float(min_cds_pct),
                        evalue=evalue_float,
                        blastn_task=blastn_task,
                        threads=threads,
                        blastn_bin=blastn_bin,
                        tblastn_bin=tblastn_bin,
                        run_dir=str(rr[1]),
                    )

                if result:
                    n_introns    = result["n_introns"]
                    splices_ok   = sum(
                        1 for intr in result["introns"]
                        if intr["splice_5"] == "GT" and intr["splice_3"] == "AG"
                    )
                    splice_note = (
                        f" · GT-AG splice sites: {splices_ok}/{n_introns}"
                        if n_introns > 0 else ""
                    )
                    st.write(
                        f"✅ {locus} ({result['blast_type']}): "
                        f"{result['n_exons']} exon(s), {result['cds_length']} bp CDS, "
                        f"{n_introns} intron(s){splice_note} · "
                        f"ref: `{result.get('ref_accession', '?')}`"
                    )
                else:
                    st.write(f"⚠️ {locus}: {status}")

                job_n += 1; prog.progress(job_n / total_jobs)

        # ── Merge combined files ───────────────────────────────────────────
        st.markdown("**Merging into combined multi-FASTAs…**")
        for locus in sel_loci:
            if locus in ("ITS", "LSU", "SSU", "ITS1", "ITS2", "ITS_full"):
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
                    for fname in [f"{locus}_CDS.fasta", f"{locus}.fasta"]:
                        fp = lsd / fname
                        if fp.exists():
                            recs = list(SeqIO.parse(str(fp), "fasta"))
                            length = len(recs[0].seq) if recs else 0
                            break
                row[locus] = length
            rows.append(row)

        df_sum = pd.DataFrame(rows).set_index("Strain")
        df_sum = df_sum.loc[:, (df_sum > 0).any()]
        st.dataframe(
            df_sum.style.map(
                lambda v: "background-color:#c8e6c9" if v > 0 else "background-color:#ffcdd2"
            ),
            use_container_width=True,
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
                                         use_container_width=True, hide_index=True)

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
                                     use_container_width=True, hide_index=True)

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
                                 key=f"aln_{browse_strain}_{locus}"):
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
                                    st.session_state[f"aln_{browse_strain}_{locus}"] = \
                                        Path(aln_path).read_text()
                                else:
                                    st.error(f"BLAST failed: {err[:200]}")

                aln_text = st.session_state.get(f"aln_{browse_strain}_{locus}", "")
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
                    if fp.is_file() and fp.suffix in (".fasta", ".gff3", ".nex", ".log"):
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
            for kind in ["_CDS", "_genomic", "_introns"]:
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
