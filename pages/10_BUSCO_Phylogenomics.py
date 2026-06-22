"""
pages/10_BUSCO_Phylogenomics.py
-------------------------------
Genome-comparison BUSCO / Compleasm workflow:

  1. Genome Input  — NCBI accessions (primary) + optional local project assemblies
  2. Download & Analyse — download genomes from NCBI; run BUSCO or Compleasm
  3. Import Existing Runs — scan for already-completed BUSCO/Compleasm directories
  4. Occupancy Matrix — colour-coded heatmap, filtering controls
  5. Export — per-BUSCO SC FASTAs → Alignment Prep
"""

import sys
from pathlib import Path

import pandas as pd
import streamlit as st

sys.path.insert(0, str(Path(__file__).parent.parent / "src"))
from phylofetch.busco_utils import (
    BUSCO_LINEAGE_HINTS,
    COMPLEASM_LINEAGE_HINTS,
    BuscoResult,
    build_occupancy_matrix,
    download_ncbi_genome,
    export_sc_fastas,
    filter_single_copy_buscos,
    run_busco,
    run_compleasm,
    scan_busco_run,
)
from phylofetch.config import load_config
from phylofetch.project_manager import load_assembly_registry

st.set_page_config(
    page_title="BUSCO Phylogenomics", page_icon="🧬", layout="wide"
)
st.title("🧬 BUSCO Phylogenomics")
st.caption(
    "Compare genomes via BUSCO single-copy orthologs. "
    "Add NCBI genome accessions, download assemblies, run BUSCO or Compleasm, "
    "and export a filtered supermatrix for phylogenetic inference."
)

# ── Session state ─────────────────────────────────────────────────────────────

cfg = load_config()
if "project_dir" not in st.session_state:
    st.session_state.project_dir = cfg.get(
        "project_dir",
        str(Path.home() / ".phylofetch" / "projects" / "default"),
    )
if "busco_results" not in st.session_state:
    st.session_state.busco_results = {}   # sample_id → BuscoResult
if "ncbi_genome_queue" not in st.session_state:
    st.session_state.ncbi_genome_queue = []   # list of {accession, label}
if "busco_run_log" not in st.session_state:
    st.session_state.busco_run_log = []   # list of {sample_id, log, status}

_project_dir = Path(st.session_state.project_dir)
_genomes_dir = _project_dir / "ncbi_genomes"
_busco_dir   = _project_dir / "busco_runs"

(
    tab_input,
    tab_run,
    tab_import,
    tab_matrix,
    tab_export,
) = st.tabs([
    "🌐 Genome Input",
    "⚙️ Download & Analyse",
    "📥 Import Existing Runs",
    "📊 Occupancy Matrix",
    "📤 Export",
])

# ══════════════════════════════════════════════════════════════════════════════
# TAB 1 — Genome Input
# ══════════════════════════════════════════════════════════════════════════════
with tab_input:
    st.subheader("Specify comparison genomes")
    st.markdown(
        "Enter NCBI genome accessions (GCA or GCF). "
        "These will be downloaded and analysed with BUSCO/Compleasm in the next tab."
    )

    col_acc, col_local = st.columns([3, 2])

    with col_acc:
        st.markdown("**NCBI genome accessions**")
        raw_accessions = st.text_area(
            "One accession per line (GCA_XXXXXXX.X or GCF_XXXXXXX.X)",
            value="\n".join(
                e["accession"] for e in st.session_state.ncbi_genome_queue
            ),
            height=200,
            key="acc_textarea",
            placeholder=(
                "GCA_000002945.2\n"
                "GCF_000149735.1\n"
                "GCA_001599815.1"
            ),
            help="Paste assembly accessions from NCBI Assembly or Genome pages.",
        )

        add_btn = st.button("➕ Add to queue", type="primary", key="btn_add_acc")
        if add_btn and raw_accessions.strip():
            queue = list(st.session_state.ncbi_genome_queue)
            existing = {e["accession"] for e in queue}
            added = 0
            for line in raw_accessions.strip().splitlines():
                acc = line.strip()
                if not acc:
                    continue
                if not (acc.startswith("GCA_") or acc.startswith("GCF_")):
                    st.warning(f"Skipped (not GCA/GCF): {acc}")
                    continue
                if acc in existing:
                    continue
                queue.append({"accession": acc, "label": acc})
                existing.add(acc)
                added += 1
            st.session_state.ncbi_genome_queue = queue
            if added:
                st.success(f"Added {added} accession(s) to queue.")

    with col_local:
        st.markdown("**Local assemblies (optional)**")
        st.caption(
            "Local project assemblies are not included by default. "
            "Enable below to add them as comparison genomes."
        )
        include_local = st.toggle(
            "Include registered project assemblies",
            value=False,
            key="include_local",
        )
        if include_local:
            try:
                registry = load_assembly_registry(_project_dir)
            except Exception:
                registry = {}
            if not registry:
                st.info("No assemblies registered in this project yet.")
            else:
                local_choices = st.multiselect(
                    "Select assemblies to include",
                    options=list(registry.keys()),
                    default=[],
                    key="local_assembly_sel",
                )
                if local_choices:
                    # Stash selected records under a SEPARATE key — writing the
                    # multiselect's own widget key ("local_assembly_sel") raises
                    # StreamlitAPIException, as it is bound to the live widget.
                    st.session_state["local_assembly_records"] = {
                        sid: registry[sid] for sid in local_choices
                    }
                    st.caption(
                        f"{len(local_choices)} local assembly/assemblies selected."
                    )

    # Current queue table
    st.markdown("---")
    st.markdown(f"**Download queue ({len(st.session_state.ncbi_genome_queue)} accessions)**")

    if st.session_state.ncbi_genome_queue:
        queue_df = pd.DataFrame(st.session_state.ncbi_genome_queue)
        edited_queue = st.data_editor(
            queue_df,
            width="stretch",
            hide_index=True,
            column_config={
                "accession": st.column_config.TextColumn("Accession", disabled=True),
                "label": st.column_config.TextColumn(
                    "Sample label",
                    help="Rename the sample as it appears in the occupancy matrix.",
                ),
            },
        )
        st.session_state.ncbi_genome_queue = edited_queue.to_dict("records")

        if st.button("🗑️ Clear queue", key="btn_clear_queue"):
            st.session_state.ncbi_genome_queue = []
            st.rerun()
    else:
        st.info("Queue is empty — paste accessions above and click Add.")

# ══════════════════════════════════════════════════════════════════════════════
# TAB 2 — Download & Analyse
# ══════════════════════════════════════════════════════════════════════════════
with tab_run:
    st.subheader("Download genomes & run BUSCO / Compleasm")

    if not st.session_state.ncbi_genome_queue and not st.session_state.get("local_assembly_records"):
        st.info("Add genomes in the **Genome Input** tab first.")
    else:
        col_cfg1, col_cfg2 = st.columns(2)

        with col_cfg1:
            tool_choice = st.radio(
                "Ortholog tool",
                options=["BUSCO", "Compleasm"],
                horizontal=True,
                key="busco_tool_choice",
                help=(
                    "BUSCO: robust, widely cited. "
                    "Compleasm: faster, compatible output format."
                ),
            )
            is_busco = tool_choice == "BUSCO"
            lineage_hints = BUSCO_LINEAGE_HINTS if is_busco else COMPLEASM_LINEAGE_HINTS
            lineage = st.selectbox(
                "Lineage",
                options=lineage_hints,
                index=0,
                key="busco_lineage",
                help=(
                    "Choose the lineage database that covers all comparison genomes. "
                    "Use a higher-level lineage (e.g. ascomycota) when mixing orders."
                ),
            )
            custom_lineage = st.text_input(
                "Or enter custom lineage",
                key="busco_lineage_custom",
                placeholder="e.g. pleosporales_odb10",
            )
            if custom_lineage.strip():
                lineage = custom_lineage.strip()

        with col_cfg2:
            cpu_count = st.slider("CPUs per run", min_value=1, max_value=32, value=4, key="busco_cpu")
            genomes_out = st.text_input(
                "Genome download directory",
                value=str(_genomes_dir),
                key="genomes_out_dir",
            )
            busco_out = st.text_input(
                "BUSCO output directory",
                value=str(_busco_dir),
                key="busco_out_dir",
            )

        st.markdown("---")
        st.markdown("**Pending analysis**")

        pending_ncbi = [
            e for e in st.session_state.ncbi_genome_queue
            if e["label"] not in st.session_state.busco_results
        ]
        pending_local = {
            sid: info
            for sid, info in st.session_state.get("local_assembly_records", {}).items()
            if sid not in st.session_state.busco_results
        }

        n_pending = len(pending_ncbi) + len(pending_local)
        if n_pending == 0:
            st.success("All queued genomes already have BUSCO results.")
        else:
            st.write(
                f"**{len(pending_ncbi)}** NCBI genome(s) to download + analyse, "
                f"**{len(pending_local)}** local assembly/assemblies to analyse."
            )

        run_btn = st.button(
            f"▶ Download & Run {tool_choice} ({n_pending} pending)",
            type="primary",
            disabled=n_pending == 0,
            key="btn_run_busco",
        )

        if run_btn:
            busco_results = dict(st.session_state.busco_results)
            run_log: list[dict] = list(st.session_state.busco_run_log)

            # ── NCBI genomes ──────────────────────────────────────────────
            for entry in pending_ncbi:
                acc   = entry["accession"]
                label = entry["label"]

                with st.status(f"⬇️  Downloading {acc} …", expanded=True) as status:
                    rc, log, fasta_path = download_ncbi_genome(
                        accession=acc,
                        output_dir=genomes_out,
                    )
                    if rc != 0 or fasta_path is None:
                        status.update(label=f"❌ Download failed: {acc}", state="error")
                        st.code(log, language=None)
                        run_log.append({"sample_id": label, "stage": "download", "status": "error", "log": log})
                        continue
                    status.update(label=f"✅ Downloaded {acc}", state="complete")

                with st.status(f"🔬 Running {tool_choice} on {label} …", expanded=True) as status:
                    if is_busco:
                        rc2, log2, result = run_busco(
                            assembly_fasta=fasta_path,
                            lineage=lineage,
                            output_dir=busco_out,
                            sample_id=label,
                            cpu=cpu_count,
                        )
                    else:
                        rc2, log2, result = run_compleasm(
                            assembly_fasta=fasta_path,
                            lineage=lineage.replace("_odb10", ""),
                            output_dir=busco_out,
                            sample_id=label,
                            cpu=cpu_count,
                        )

                    if rc2 != 0 or result is None:
                        status.update(label=f"❌ {tool_choice} failed: {label}", state="error")
                        st.code(log2, language=None)
                        run_log.append({"sample_id": label, "stage": tool_choice.lower(), "status": "error", "log": log2})
                    else:
                        busco_results[label] = result
                        status.update(
                            label=(
                                f"✅ {label}: {result.single_copy} S / "
                                f"{result.duplicated} D / {result.fragmented} F / "
                                f"{result.missing} M  "
                                f"({result.completeness_pct:.1f}% complete)"
                            ),
                            state="complete",
                        )
                        run_log.append({"sample_id": label, "stage": tool_choice.lower(), "status": "ok", "log": log2})

            # ── Local assemblies ──────────────────────────────────────────
            for sid, info in pending_local.items():
                asm_path = info.get("assembly_path", "")
                if not asm_path or not Path(asm_path).exists():
                    st.warning(f"Assembly path not found for {sid}: {asm_path}")
                    continue

                with st.status(f"🔬 Running {tool_choice} on local {sid} …", expanded=True) as status:
                    if is_busco:
                        rc2, log2, result = run_busco(
                            assembly_fasta=asm_path,
                            lineage=lineage,
                            output_dir=busco_out,
                            sample_id=sid,
                            cpu=cpu_count,
                        )
                    else:
                        rc2, log2, result = run_compleasm(
                            assembly_fasta=asm_path,
                            lineage=lineage.replace("_odb10", ""),
                            output_dir=busco_out,
                            sample_id=sid,
                            cpu=cpu_count,
                        )

                    if rc2 != 0 or result is None:
                        status.update(label=f"❌ {tool_choice} failed: {sid}", state="error")
                        st.code(log2, language=None)
                        run_log.append({"sample_id": sid, "stage": tool_choice.lower(), "status": "error", "log": log2})
                    else:
                        busco_results[sid] = result
                        status.update(
                            label=(
                                f"✅ {sid}: {result.single_copy} S / "
                                f"{result.duplicated} D / {result.fragmented} F / "
                                f"{result.missing} M  "
                                f"({result.completeness_pct:.1f}% complete)"
                            ),
                            state="complete",
                        )
                        run_log.append({"sample_id": sid, "stage": tool_choice.lower(), "status": "ok", "log": log2})

            st.session_state.busco_results = busco_results
            st.session_state.busco_run_log = run_log

        # Run log expander
        if st.session_state.busco_run_log:
            with st.expander(f"📋 Run log ({len(st.session_state.busco_run_log)} entries)"):
                for entry in reversed(st.session_state.busco_run_log):
                    icon = "✅" if entry["status"] == "ok" else "❌"
                    st.markdown(f"**{icon} {entry['sample_id']}** ({entry['stage']})")
                    st.code(entry["log"][-3000:], language=None)

# ══════════════════════════════════════════════════════════════════════════════
# TAB 3 — Import Existing Runs
# ══════════════════════════════════════════════════════════════════════════════
with tab_import:
    st.subheader("Import already-completed BUSCO / Compleasm runs")
    st.caption(
        "If you have pre-existing BUSCO or Compleasm output directories, "
        "scan a parent folder to register them directly without re-running."
    )

    scan_busco_dir = st.text_input(
        "Directory to scan",
        placeholder="/data/my_busco_runs/",
        key="scan_busco_dir",
    )

    if st.button("🔍 Scan for BUSCO / Compleasm runs", key="btn_scan") and scan_busco_dir:
        if not Path(scan_busco_dir).is_dir():
            st.error(f"Directory not found: {scan_busco_dir}")
        else:
            found_dirs: list[Path] = []
            for summary in Path(scan_busco_dir).rglob("short_summary*.txt"):
                rd = summary.parent
                if rd not in found_dirs:
                    found_dirs.append(rd)
            for summary in Path(scan_busco_dir).rglob("summary.txt"):
                rd = summary.parent
                if rd not in found_dirs:
                    found_dirs.append(rd)

            if not found_dirs:
                st.warning("No BUSCO/Compleasm results found.")
            else:
                st.success(f"Found {len(found_dirs)} run(s).")
                scan_root = Path(scan_busco_dir)
                _skip = {"busco", "compleasm", "qc", "busco_output", "busco_results"}
                rows = []
                for rd in sorted(found_dirs):
                    sample_id = rd.name
                    cand = rd
                    for _ in range(4):
                        cand = cand.parent
                        if cand == scan_root or not cand.name:
                            break
                        if not (cand.name.lower().startswith("run_")
                                or cand.name.lower() in _skip):
                            sample_id = cand.name
                            break
                    rows.append({
                        "Import?": True,
                        "Sample ID": sample_id,
                        "Run directory": str(rd),
                    })
                st.session_state["busco_scan_rows"] = rows

    if "busco_scan_rows" in st.session_state:
        edited = st.data_editor(
            pd.DataFrame(st.session_state["busco_scan_rows"]),
            width="stretch",
            hide_index=True,
            column_config={
                "Import?": st.column_config.CheckboxColumn(),
                "Sample ID": st.column_config.TextColumn(),
                "Run directory": st.column_config.TextColumn(disabled=True),
            },
        )
        to_import = edited[edited["Import?"] == True]
        st.write(f"**{len(to_import)}** selected for import")

        if st.button("⬇️ Import selected", type="primary", key="btn_import_runs") and not to_import.empty:
            imported, errors = 0, []
            busco_results = dict(st.session_state.busco_results)
            for _, row in to_import.iterrows():
                sid  = str(row["Sample ID"]).strip()
                rdir = str(row["Run directory"]).strip()
                if not sid:
                    errors.append(f"Empty sample ID for {rdir}")
                    continue
                if sid in busco_results:
                    st.warning(f"⚠️  {sid} already registered — skipping.")
                    continue
                try:
                    result = scan_busco_run(run_dir=rdir, sample_id=sid)
                    if result is None:
                        errors.append(f"{sid}: could not parse run directory")
                    else:
                        busco_results[sid] = result
                        imported += 1
                except Exception as exc:
                    errors.append(f"{sid}: {exc}")
            st.session_state.busco_results = busco_results
            if imported:
                st.success(f"Imported {imported} run(s).")
            for err in errors:
                st.error(err)
            del st.session_state["busco_scan_rows"]
            st.rerun()

    # Show registered runs
    st.markdown("---")
    st.subheader(f"Registered results ({len(st.session_state.busco_results)})")
    if st.session_state.busco_results:
        rows = []
        for sid, res in st.session_state.busco_results.items():
            rows.append({
                "Sample": sid,
                "Tool": res.tool,
                "Lineage": res.lineage or "—",
                "S (single-copy)": res.single_copy,
                "D (duplicated)": res.duplicated,
                "F (fragmented)": res.fragmented,
                "M (missing)": res.missing,
                "Total": res.total,
                "Complete %": f"{res.completeness_pct:.1f}",
            })
        st.dataframe(pd.DataFrame(rows), width="stretch", hide_index=True)

        to_remove = st.selectbox(
            "Remove a result", [""] + list(st.session_state.busco_results.keys()),
            key="busco_remove_sel",
        )
        if to_remove and st.button(f"🗑️ Remove {to_remove}", key="btn_busco_remove"):
            busco_results = dict(st.session_state.busco_results)
            del busco_results[to_remove]
            st.session_state.busco_results = busco_results
            st.rerun()
    else:
        st.info("No results registered yet.")

# ══════════════════════════════════════════════════════════════════════════════
# TAB 4 — Occupancy Matrix
# ══════════════════════════════════════════════════════════════════════════════
with tab_matrix:
    st.subheader("Occupancy matrix")

    if not st.session_state.busco_results:
        st.info("Add and analyse genomes first (Genome Input → Download & Analyse).")
    else:
        results: list[BuscoResult] = list(st.session_state.busco_results.values())
        matrix = build_occupancy_matrix(results)

        st.caption(
            f"**{matrix.shape[0]} samples × {matrix.shape[1]} BUSCO IDs** — "
            "S=single-copy complete, D=duplicated, F=fragmented, M=missing"
        )

        # Summary metrics
        col_s1, col_s2, col_s3, col_s4 = st.columns(4)
        total_cells = matrix.size
        s_count = (matrix == "S").sum().sum()
        d_count = (matrix == "D").sum().sum()
        f_count = (matrix == "F").sum().sum()
        m_count = (matrix == "M").sum().sum()
        col_s1.metric("Complete (S)", f"{s_count} ({s_count/total_cells*100:.1f}%)")
        col_s2.metric("Duplicated (D)", f"{d_count} ({d_count/total_cells*100:.1f}%)")
        col_s3.metric("Fragmented (F)", f"{f_count} ({f_count/total_cells*100:.1f}%)")
        col_s4.metric("Missing (M)", f"{m_count} ({m_count/total_cells*100:.1f}%)")

        occupancy = (matrix == "S").mean(axis=0) * 100
        st.markdown(
            f"Mean single-copy occupancy across {matrix.shape[1]} BUSCOs: "
            f"**{occupancy.mean():.1f}%**"
        )

        col_filt, col_show = st.columns([1, 3])
        with col_filt:
            st.markdown("**Filter settings**")
            min_occ = st.slider("Min occupancy (%)", 0, 100, 75, 5, key="min_occ_slider")
            excl_dup = st.checkbox("Exclude if any sample duplicated", value=True, key="excl_dup")
            max_missing = st.slider("Max per-sample missingness (%)", 0, 100, 50, 5, key="max_missing")

        with col_show:
            selected_buscos = filter_single_copy_buscos(
                matrix,
                min_occupancy=min_occ / 100,
                exclude_duplicated=excl_dup,
            )
            if max_missing < 100:
                sample_missing = (matrix == "M").mean(axis=1) * 100
                keep_samples = sample_missing[sample_missing <= max_missing].index.tolist()
                sub_matrix = matrix.loc[keep_samples, selected_buscos]
            else:
                keep_samples = matrix.index.tolist()
                sub_matrix = matrix.loc[:, selected_buscos]

            col_b1, col_b2 = st.columns(2)
            col_b1.metric("BUSCOs passing filters", f"{len(selected_buscos)} / {matrix.shape[1]}")
            col_b2.metric("Samples passing filters", f"{len(keep_samples)} / {matrix.shape[0]}")

            if len(selected_buscos) > 0 and len(keep_samples) > 0:
                def _color_cell(val: str) -> str:
                    return {
                        "S": "background-color: #2ecc71; color: #fff",
                        "D": "background-color: #f39c12; color: #fff",
                        "F": "background-color: #e67e22; color: #fff",
                        "M": "background-color: #e74c3c; color: #fff",
                    }.get(str(val), "")

                with st.expander("Show filtered matrix", expanded=False):
                    show_cols = selected_buscos[:100]
                    if len(selected_buscos) > 100:
                        st.caption(f"Showing first 100 of {len(selected_buscos)} BUSCOs.")
                    try:
                        styled = sub_matrix[show_cols].style.map(_color_cell)
                        st.dataframe(styled, width="stretch")
                    except Exception:
                        st.dataframe(sub_matrix[show_cols], width="stretch")

        st.session_state["selected_buscos"]    = selected_buscos
        st.session_state["busco_keep_samples"] = keep_samples

# ══════════════════════════════════════════════════════════════════════════════
# TAB 5 — Export
# ══════════════════════════════════════════════════════════════════════════════
with tab_export:
    st.subheader("Export single-copy BUSCO FASTAs")

    if not st.session_state.busco_results:
        st.info("Add and analyse genomes first.")
    elif "selected_buscos" not in st.session_state or not st.session_state["selected_buscos"]:
        st.info("Set filters in the Occupancy Matrix tab first.")
    else:
        selected_buscos = st.session_state["selected_buscos"]
        keep_samples = st.session_state.get(
            "busco_keep_samples", list(st.session_state.busco_results.keys())
        )

        st.write(
            f"Ready to export **{len(selected_buscos)} BUSCOs** × "
            f"**{len(keep_samples)} samples**."
        )

        col_ex1, col_ex2 = st.columns(2)
        with col_ex1:
            export_out_dir = st.text_input(
                "Export directory",
                value=str(_project_dir / "busco_sc_fastas"),
                key="busco_export_dir",
            )
        with col_ex2:
            seq_type = st.radio(
                "Sequence type",
                options=["nucleotide", "protein"],
                index=0,
                horizontal=True,
                help=(
                    "nucleotide: .fna sequences from BUSCO single_copy_busco_sequences/. "
                    "protein: .faa sequences."
                ),
            )

        if st.button("📤 Export FASTAs", type="primary", key="btn_busco_export"):
            if not export_out_dir:
                st.warning("Specify an export directory.")
            else:
                results_subset = [
                    res
                    for sid, res in st.session_state.busco_results.items()
                    if sid in keep_samples
                ]

                with st.spinner(f"Exporting {len(selected_buscos)} BUSCO loci…"):
                    try:
                        export_stats = export_sc_fastas(
                            results=results_subset,
                            selected_buscos=selected_buscos,
                            output_dir=export_out_dir,
                            seq_type=seq_type,
                        )
                    except Exception as exc:
                        st.error(f"Export error: {exc}")
                        export_stats = None

                if export_stats:
                    st.success(
                        f"Exported {export_stats['exported_loci']} FASTA files to `{export_out_dir}`."
                    )

                    col_e1, col_e2 = st.columns(2)
                    col_e1.metric("Loci exported", export_stats["exported_loci"])
                    col_e2.metric("Missing SC sequences", export_stats["n_missing"])

                    if export_stats.get("missing_sc_paths"):
                        with st.expander(f"⚠️ {export_stats['n_missing']} missing SC sequences"):
                            for item in export_stats["missing_sc_paths"][:50]:
                                st.write(f"- {item}")
                            if len(export_stats["missing_sc_paths"]) > 50:
                                st.caption(
                                    f"… and {len(export_stats['missing_sc_paths']) - 50} more."
                                )

                    st.info(
                        f"Per-BUSCO FASTAs written to `{export_out_dir}`. "
                        "Go to **Alignment Prep** to align and concatenate."
                    )

                    if st.button("→ Send to Alignment Prep", key="btn_send_align"):
                        busco_fastas = sorted(Path(export_out_dir).glob("*.fasta"))
                        st.session_state["trimmed_fastas"] = [str(f) for f in busco_fastas]
                        st.session_state["aligned_fastas"] = []
                        st.success(
                            f"Sent {len(busco_fastas)} BUSCO FASTAs to Alignment Prep. "
                            "Open the Alignment Prep page to continue."
                        )
