"""
pages/4_BUSCO_Phylogenomics.py
-------------------------------
BUSCO / Compleasm occupancy matrix, missingness filtering,
single-copy ortholog FASTA export, and integration with Alignment Prep.
"""

import sys
from pathlib import Path

import pandas as pd
import streamlit as st

sys.path.insert(0, str(Path(__file__).parent.parent / "src"))
from phylofetch.busco_utils import (
    BuscoResult,
    build_occupancy_matrix,
    export_sc_fastas,
    filter_single_copy_buscos,
    scan_busco_run,
)
from phylofetch.config import load_config

st.set_page_config(
    page_title="BUSCO Phylogenomics", page_icon="🧫", layout="wide"
)
st.title("🧫 BUSCO Phylogenomics")
st.caption(
    "Import BUSCO or Compleasm run directories, build an occupancy matrix, "
    "filter single-copy orthologs, and export FASTAs for supermatrix construction."
)

if "project_dir" not in st.session_state:
    st.session_state.project_dir = load_config().get(
        "project_dir",
        str(Path.home() / ".phylofetch" / "projects" / "default"),
    )
if "busco_results" not in st.session_state:
    st.session_state.busco_results = {}  # sample_id → BuscoResult

tab_import, tab_matrix, tab_export = st.tabs(
    ["📥 Import BUSCO Runs", "📊 Occupancy Matrix", "📤 Export"]
)

# ══════════════════════════════════════════════════════════════════════════════
# TAB 1 — Import BUSCO / Compleasm runs
# ══════════════════════════════════════════════════════════════════════════════
with tab_import:
    st.subheader("Register BUSCO / Compleasm runs")
    st.caption(
        "Point to a directory containing one or more BUSCO v4/v5 or Compleasm result folders. "
        "phylofetch auto-detects the format and extracts completeness stats."
    )

    scan_busco_dir = st.text_input(
        "Directory to scan for BUSCO runs",
        placeholder="/data/my_busco_runs/",
        key="scan_busco_dir",
    )

    if st.button("🔍 Scan for BUSCO runs") and scan_busco_dir:
        if not Path(scan_busco_dir).is_dir():
            st.error(f"Directory not found: {scan_busco_dir}")
        else:
            # Auto-detect: look for run_*/short_summary*.txt (BUSCO)
            # or compleasm output directories (results.txt / summary.txt)
            found_dirs: list[Path] = []

            # BUSCO v4/v5: short_summary*.txt in each run folder
            for summary in Path(scan_busco_dir).rglob("short_summary*.txt"):
                run_dir = summary.parent
                if run_dir not in found_dirs:
                    found_dirs.append(run_dir)

            # Compleasm: summary.txt at top-level or in run_*/
            for summary in Path(scan_busco_dir).rglob("summary.txt"):
                run_dir = summary.parent
                if run_dir not in found_dirs:
                    found_dirs.append(run_dir)

            if not found_dirs:
                st.warning(
                    "No BUSCO/Compleasm results found. "
                    "Expected short_summary*.txt (BUSCO) or summary.txt (Compleasm)."
                )
            else:
                st.success(f"Found {len(found_dirs)} run(s).")
                scan_root = Path(scan_busco_dir)
                # Names that indicate an intermediate container dir, not the sample dir
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
                    rows.append(
                        {
                            "Import?": True,
                            "Sample ID": sample_id,
                            "Run directory": str(rd),
                        }
                    )
                st.session_state["busco_scan_rows"] = rows

    if "busco_scan_rows" in st.session_state:
        edited = st.data_editor(
            pd.DataFrame(st.session_state["busco_scan_rows"]),
            use_container_width=True,
            hide_index=True,
            column_config={
                "Import?": st.column_config.CheckboxColumn(),
                "Sample ID": st.column_config.TextColumn(),
                "Run directory": st.column_config.TextColumn(disabled=True),
            },
        )

        to_import = edited[edited["Import?"] == True]
        st.write(f"**{len(to_import)}** selected for import")

        if st.button("⬇️ Import selected BUSCO runs", type="primary") and not to_import.empty:
            imported, errors = 0, []
            busco_results = dict(st.session_state.busco_results)

            for _, row in to_import.iterrows():
                sid = str(row["Sample ID"]).strip()
                rdir = str(row["Run directory"]).strip()

                if not sid:
                    errors.append(f"Empty sample ID for {rdir}")
                    continue
                if sid in busco_results:
                    st.warning(f"⚠️  {sid} already registered — skipping.")
                    continue
                try:
                    result = scan_busco_run(run_dir=rdir, sample_id=sid)
                    busco_results[sid] = result
                    imported += 1
                except Exception as exc:
                    errors.append(f"{sid}: {exc}")

            st.session_state.busco_results = busco_results
            if imported:
                st.success(f"Imported {imported} BUSCO run(s).")
            for err in errors:
                st.error(err)
            if "busco_scan_rows" in st.session_state:
                del st.session_state["busco_scan_rows"]
            st.rerun()

    # Show registered runs
    st.markdown("---")
    st.subheader(f"Registered BUSCO runs ({len(st.session_state.busco_results)})")

    if st.session_state.busco_results:
        rows = []
        for sid, res in st.session_state.busco_results.items():
            rows.append(
                {
                    "Sample ID": sid,
                    "Tool": res.tool,
                    "Lineage": res.lineage or "—",
                    "Complete (S)": res.counts.get("S", "—"),
                    "Duplicate (D)": res.counts.get("D", "—"),
                    "Fragmented (F)": res.counts.get("F", "—"),
                    "Missing (M)": res.counts.get("M", "—"),
                    "Total": res.counts.get("total", "—"),
                    "Run dir": res.run_dir,
                }
            )
        st.dataframe(pd.DataFrame(rows), use_container_width=True, hide_index=True)

        to_remove = st.selectbox(
            "Remove a run", [""] + list(st.session_state.busco_results.keys()),
            key="busco_remove_sel",
        )
        if to_remove and st.button(f"🗑️ Remove {to_remove}", key="btn_busco_remove"):
            busco_results = dict(st.session_state.busco_results)
            del busco_results[to_remove]
            st.session_state.busco_results = busco_results
            st.rerun()
    else:
        st.info("No BUSCO runs registered yet.")

# ══════════════════════════════════════════════════════════════════════════════
# TAB 2 — Occupancy matrix
# ══════════════════════════════════════════════════════════════════════════════
with tab_matrix:
    st.subheader("Occupancy matrix")

    if not st.session_state.busco_results:
        st.info("Import BUSCO runs first.")
    else:
        results: list[BuscoResult] = list(st.session_state.busco_results.values())
        matrix = build_occupancy_matrix(results)

        st.caption(
            f"**{matrix.shape[0]} samples × {matrix.shape[1]} BUSCO IDs** — "
            "S=single-copy complete, D=duplicated, F=fragmented, M=missing"
        )

        # Color-coded display
        def _color_cell(val: str) -> str:
            colors = {
                "S": "background-color: #2ecc71; color: #fff",
                "D": "background-color: #f39c12; color: #fff",
                "F": "background-color: #e67e22; color: #fff",
                "M": "background-color: #e74c3c; color: #fff",
            }
            return colors.get(str(val), "")

        # Show summary stats first
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

        # Occupancy per BUSCO
        occupancy = (matrix == "S").mean(axis=0) * 100
        st.markdown(
            f"Mean single-copy occupancy across {matrix.shape[1]} BUSCOs: "
            f"**{occupancy.mean():.1f}%**"
        )

        col_filt, col_show = st.columns([1, 3])
        with col_filt:
            st.markdown("**Filter settings**")
            min_occ = st.slider(
                "Min occupancy (%)",
                min_value=0,
                max_value=100,
                value=75,
                step=5,
                key="min_occ_slider",
            )
            excl_dup = st.checkbox(
                "Exclude if any sample is Duplicated",
                value=True,
                key="excl_dup",
            )
            max_missing = st.slider(
                "Max per-sample missingness (%)",
                min_value=0,
                max_value=100,
                value=50,
                step=5,
                key="max_missing",
            )

        with col_show:
            selected_buscos = filter_single_copy_buscos(
                matrix,
                min_occupancy=min_occ / 100,
                exclude_duplicated=excl_dup,
            )
            # Apply per-sample missingness filter
            if max_missing < 100:
                sample_missing = (matrix == "M").mean(axis=1) * 100
                keep_samples = sample_missing[sample_missing <= max_missing].index.tolist()
                sub_matrix = matrix.loc[keep_samples, selected_buscos]
            else:
                keep_samples = matrix.index.tolist()
                sub_matrix = matrix.loc[:, selected_buscos]

            st.metric(
                "BUSCOs passing filters",
                f"{len(selected_buscos)} / {matrix.shape[1]}",
            )
            st.metric(
                "Samples passing filters",
                f"{len(keep_samples)} / {matrix.shape[0]}",
            )

            if len(selected_buscos) > 0 and len(keep_samples) > 0:
                with st.expander("Show filtered matrix", expanded=False):
                    show_cols = selected_buscos[:100]  # cap for display
                    if len(selected_buscos) > 100:
                        st.caption(f"Showing first 100 of {len(selected_buscos)} BUSCOs.")
                    try:
                        styled = (
                            sub_matrix[show_cols]
                            .style.applymap(_color_cell)
                        )
                        st.dataframe(styled, use_container_width=True)
                    except Exception:
                        st.dataframe(
                            sub_matrix[show_cols], use_container_width=True
                        )

        st.session_state["selected_buscos"] = selected_buscos
        st.session_state["busco_keep_samples"] = keep_samples

# ══════════════════════════════════════════════════════════════════════════════
# TAB 3 — Export
# ══════════════════════════════════════════════════════════════════════════════
with tab_export:
    st.subheader("Export single-copy BUSCO FASTAs")

    if not st.session_state.busco_results:
        st.info("Import BUSCO runs first.")
    elif "selected_buscos" not in st.session_state or not st.session_state["selected_buscos"]:
        st.info("Set filters in the Occupancy Matrix tab first.")
    else:
        selected_buscos = st.session_state["selected_buscos"]
        keep_samples = st.session_state.get("busco_keep_samples", list(st.session_state.busco_results.keys()))

        st.write(
            f"Ready to export **{len(selected_buscos)} BUSCOs** × "
            f"**{len(keep_samples)} samples**."
        )

        export_out_dir = st.text_input(
            "Export directory",
            value=str(
                Path(st.session_state.project_dir) / "busco_sc_fastas"
            ),
            key="busco_export_dir",
        )

        seq_type = st.radio(
            "Sequence type",
            options=["nucleotide", "protein"],
            index=0,
            horizontal=True,
            help=(
                "nucleotide: single-copy BUSCO sequences from genome (DNA). "
                "protein: translated sequences (AA). "
                "Depends on what your BUSCO run stored."
            ),
        )

        if st.button("📤 Export FASTAs", type="primary", key="btn_busco_export"):
            if not export_out_dir:
                st.warning("Specify an export directory.")
            else:
                out_path = Path(export_out_dir)
                out_path.mkdir(parents=True, exist_ok=True)

                # Filter results to keep_samples
                results_subset = {
                    sid: res
                    for sid, res in st.session_state.busco_results.items()
                    if sid in keep_samples
                }

                with st.spinner(f"Exporting {len(selected_buscos)} BUSCO loci…"):
                    try:
                        export_stats = export_sc_fastas(
                            results=list(results_subset.values()),
                            selected_buscos=selected_buscos,
                            output_dir=out_path,
                            seq_type=seq_type,
                        )
                    except Exception as exc:
                        st.error(f"Export error: {exc}")
                        export_stats = None

                if export_stats:
                    st.success(
                        f"Exported {export_stats['exported_loci']} FASTA files to `{out_path}`."
                    )
                    if export_stats.get("missing_sc_paths"):
                        with st.expander("⚠️ Missing SC sequence files"):
                            for item in export_stats["missing_sc_paths"]:
                                st.write(f"- {item}")

                    col_e1, col_e2 = st.columns(2)
                    col_e1.metric("Loci exported", export_stats["exported_loci"])
                    col_e2.metric(
                        "Loci with missing sequences",
                        export_stats.get("n_missing", 0),
                    )

                    st.info(
                        f"Per-BUSCO FASTAs are in `{out_path}`. "
                        "Go to **Alignment Prep** to align and concatenate them."
                    )

                    if st.button("→ Send to Alignment Prep", key="btn_send_align"):
                        busco_fastas = sorted(out_path.glob("*.fasta"))
                        st.session_state["trimmed_fastas"] = [
                            str(f) for f in busco_fastas
                        ]
                        st.session_state["aligned_fastas"] = []
                        st.success(
                            f"Sent {len(busco_fastas)} BUSCO FASTAs to Alignment Prep. "
                            "Open the Alignment Prep page to continue."
                        )
