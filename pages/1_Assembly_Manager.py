"""
pages/1_Assembly_Manager.py
----------------------------
Register genome assemblies, associate read files, and view QC statistics.
"""

import os
import sys
from pathlib import Path

import pandas as pd
import plotly.express as px
import streamlit as st

sys.path.insert(0, str(Path(__file__).parent.parent / "src"))
from phylofetch.assembly_utils import (
    QUAST_DISPLAY_KEYS,
    find_assemblies_recursive,
    find_quast_report,
    get_assembly_stats,
    parse_quast_report,
    suggest_unique_strain_ids,
)
from phylofetch.config import load_config, save_config
from phylofetch.project_manager import (
    DEFAULT_PROJECT_DIR,
    load_assembly_registry,
    now_iso,
    save_assembly_registry,
)

# ── Page setup ───────────────────────────────────────────────────────────────
st.set_page_config(
    page_title="Assembly Manager", page_icon="📁", layout="wide"
)
st.title("📁 Assembly Manager")
st.caption("Register assemblies, link read files, and explore per-strain statistics.")


def _active_project() -> str:
    return load_config().get("project_dir", str(DEFAULT_PROJECT_DIR))


# ── Session state ─────────────────────────────────────────────────────────────
if "assemblies" not in st.session_state:
    # Prefer the active project's saved registry; fall back to global config.
    st.session_state.assemblies = (
        load_assembly_registry(_active_project())
        or load_config().get("assemblies", {})
    )


def _save():
    # Global cache (other pages read this) + durable per-project registry/manifest.
    save_config({"assemblies": st.session_state.assemblies})
    save_assembly_registry(_active_project(), st.session_state.assemblies)


def _attach_quast(stats: dict, assembly_path: str) -> dict:
    """Auto-discover a QUAST report next to the assembly and merge it in."""
    quast_path = find_quast_report(assembly_path)
    if quast_path:
        stats["quast"] = parse_quast_report(quast_path)
        stats["quast_report"] = quast_path
    return stats


# ── Tabs ──────────────────────────────────────────────────────────────────────
tab_add, tab_list, tab_stats = st.tabs(
    ["🔍 Find / Add Assemblies", "📋 Registered Assemblies", "📊 Assembly Stats"]
)


# ════════════════════════════════════════════════════════
# TAB 1 — Find / Add
# ════════════════════════════════════════════════════════
with tab_add:
    method = st.radio(
        "How would you like to add assemblies?",
        ["Scan a folder recursively", "Enter a path manually"],
        horizontal=True,
    )

    # ── Recursive scan ──────────────────────────────────
    if method == "Scan a folder recursively":
        scan_dir = st.text_input(
            "Root directory to scan",
            placeholder="/home/user/genomes",
        )
        extensions = st.multiselect(
            "Assembly file extensions",
            [".fasta", ".fa", ".fna"],
            default=[".fasta", ".fa", ".fna"],
        )

        if scan_dir and os.path.isdir(scan_dir):
            if st.button("🔍 Scan"):
                found = find_assemblies_recursive(scan_dir, tuple(extensions))
                st.session_state["_scan_results"] = found
        elif scan_dir:
            st.warning("Directory not found.")

        if st.session_state.get("_scan_results"):
            found = st.session_state["_scan_results"]
            st.success(f"Found {len(found)} file(s).")

            # Build dataframe in session state on first appearance after a scan
            if "scan_df" not in st.session_state or st.session_state.get("scan_df_source") != tuple(found):
                unique_ids = suggest_unique_strain_ids(list(found))
                rows = [
                    {
                        "Add?":      True,
                        "Strain ID": sid,
                        "Path":      p,
                    }
                    for p, sid in zip(found, unique_ids)
                ]
                st.session_state["scan_df"]        = pd.DataFrame(rows)
                st.session_state["scan_df_source"] = tuple(found)

            # ── Quick filter presets ──
            qc1, qc2, _ = st.columns([1, 1, 5])
            with qc1:
                if st.button("EGAP", key="qf_egap"):
                    st.session_state["scan_filter"] = "EGAP"
                    st.rerun()
            with qc2:
                if st.button("Clear", key="qf_clear"):
                    st.session_state["scan_filter"] = ""
                    st.rerun()

            # ── Filter ──
            filter_text = st.text_input(
                "Filter paths (case-insensitive substring; leave blank for all)",
                value=st.session_state.get("scan_filter", ""),
                key="scan_filter",
                placeholder="e.g. _best, _curated, polished",
            )

            # Apply filter to a view (not the underlying df, so toggles stay)
            df_full = st.session_state["scan_df"]
            if filter_text:
                mask    = df_full["Path"].str.contains(filter_text, case=False, regex=False)
                df_view = df_full[mask].copy()
            else:
                df_view = df_full.copy()

            # ── Selection buttons ──
            n_total    = len(df_full)
            n_visible  = len(df_view)
            n_selected = int(df_full["Add?"].sum())

            bc1, bc2, bc3, bc4 = st.columns([1, 1, 1, 3])
            with bc1:
                if st.button("✅ Select all"):
                    if filter_text:
                        st.session_state["scan_df"].loc[mask, "Add?"] = True
                    else:
                        st.session_state["scan_df"]["Add?"] = True
                    st.rerun()
            with bc2:
                if st.button("⬜ Select none"):
                    if filter_text:
                        st.session_state["scan_df"].loc[mask, "Add?"] = False
                    else:
                        st.session_state["scan_df"]["Add?"] = False
                    st.rerun()
            with bc3:
                if st.button("🔁 Invert"):
                    if filter_text:
                        st.session_state["scan_df"].loc[mask, "Add?"] = ~st.session_state["scan_df"].loc[mask, "Add?"]
                    else:
                        st.session_state["scan_df"]["Add?"] = ~st.session_state["scan_df"]["Add?"]
                    st.rerun()
            with bc4:
                st.caption(
                    f"**{n_selected}** selected of **{n_total}** total"
                    + (f" · **{n_visible}** match filter" if filter_text else "")
                    + " · selection buttons act on the **filtered view**"
                )

            # ── Common-prefix detection for shorter Path display ──
            try:
                common_prefix = os.path.commonpath(df_view["Path"].tolist()) if len(df_view) > 1 else ""
                # Only collapse if prefix is meaningfully long (>30 chars)
                if len(common_prefix) > 30:
                    df_view["Path (short)"] = df_view["Path"].str.replace(common_prefix + "/", "…/", regex=False)
                    show_short = True
                else:
                    show_short = False
            except ValueError:
                show_short = False

            cols = ["Add?", "Strain ID"]
            if show_short:
                cols += ["Path (short)", "Path"]
                st.caption(f"Paths shown relative to `{common_prefix}/` — hover or scroll for full path")
            else:
                cols += ["Path"]

            edited_view = st.data_editor(
                df_view[cols],
                use_container_width=True,
                hide_index=True,
                column_config={
                    "Add?":         st.column_config.CheckboxColumn(),
                    "Strain ID":    st.column_config.TextColumn(
                        help="Edit to your preferred strain identifier"
                    ),
                    "Path (short)": st.column_config.TextColumn(disabled=True, width="medium"),
                    "Path":         st.column_config.TextColumn(disabled=True, width="large"),
                },
                key="scan_editor",
            )

            # Sync edits back to the master dataframe by Path (stable key)
            for _, row in edited_view.iterrows():
                idx = st.session_state["scan_df"].index[st.session_state["scan_df"]["Path"] == row["Path"]]
                if len(idx):
                    st.session_state["scan_df"].loc[idx, "Add?"]      = row["Add?"]
                    st.session_state["scan_df"].loc[idx, "Strain ID"] = row["Strain ID"]

            # Use the synced master df from here on
            edited = st.session_state["scan_df"]

            # ── Strain ID collision detector ──
            selected_rows = edited[edited["Add?"] == True]
            dup_mask      = selected_rows["Strain ID"].duplicated(keep=False)
            if dup_mask.any():
                dups = selected_rows[dup_mask]["Strain ID"].unique()
                st.warning(
                    f"⚠️ **Duplicate Strain IDs in selection:** {', '.join(dups)}. "
                    "Each strain ID must be unique. Edit the Strain ID column to disambiguate "
                    "before adding."
                )

            if st.button("✅ Add selected"):
                added, skipped = 0, 0
                prog = st.progress(0)
                selected_rows = edited[edited["Add?"] == True]
                total = len(selected_rows)
                for i, (_, row) in enumerate(selected_rows.iterrows()):
                    sid = row["Strain ID"].strip()
                    if not sid or not row["Path"]:
                        skipped += 1
                        continue
                    if sid in st.session_state.assemblies:
                        st.warning(f"'{sid}' already registered — skipping.")
                        skipped += 1
                        continue
                    if not os.path.exists(row["Path"]):
                        st.error(f"File not found: {row['Path']}")
                        skipped += 1
                        continue
                    with st.spinner(f"Computing stats: {sid}"):
                        stats = _attach_quast(get_assembly_stats(row["Path"]), row["Path"])
                    st.session_state.assemblies[sid] = {
                        "strain_id":     sid,
                        "assembly_path": row["Path"],
                        "reads_r1":      "",
                        "reads_r2":      "",
                        "stats":         stats,
                        "registered_at": now_iso(),
                    }
                    added += 1
                    prog.progress((i + 1) / max(total, 1))
                _save()
                st.success(f"Added {added} assembly/assemblies. Skipped {skipped}.")
                st.rerun()

    # ── Manual entry ────────────────────────────────────
    else:
        with st.form("manual_add"):
            col_id, col_path = st.columns([1, 2])
            with col_id:
                strain_id = st.text_input(
                    "Strain ID *", placeholder="ALT_CA_001"
                )
            with col_path:
                assembly_path = st.text_input(
                    "Assembly FASTA path *",
                    placeholder="/data/assemblies/ALT_CA_001/assembly.fasta",
                )

            col_r1, col_r2 = st.columns(2)
            with col_r1:
                r1 = st.text_input(
                    "Illumina R1 reads (optional)",
                    placeholder="/data/reads/ALT_CA_001_R1.fastq.gz",
                )
            with col_r2:
                r2 = st.text_input(
                    "Illumina R2 reads (optional)",
                    placeholder="/data/reads/ALT_CA_001_R2.fastq.gz",
                )

            submitted = st.form_submit_button("Add Assembly")

        if submitted:
            if not strain_id or not assembly_path:
                st.error("Strain ID and assembly path are required.")
            elif not os.path.exists(assembly_path):
                st.error(f"File not found: {assembly_path}")
            elif strain_id in st.session_state.assemblies:
                st.warning(
                    f"'{strain_id}' is already registered. "
                    "Edit it in the Registered Assemblies tab."
                )
            else:
                with st.spinner(f"Computing stats for {strain_id}…"):
                    stats = _attach_quast(get_assembly_stats(assembly_path), assembly_path)
                st.session_state.assemblies[strain_id] = {
                    "strain_id":     strain_id,
                    "assembly_path": assembly_path,
                    "reads_r1":      r1,
                    "reads_r2":      r2,
                    "stats":         stats,
                    "registered_at": now_iso(),
                }
                _save()
                st.success(f"✅ Added {strain_id}")
                st.rerun()


# ════════════════════════════════════════════════════════
# TAB 2 — Registered Assemblies
# ════════════════════════════════════════════════════════
with tab_list:
    if not st.session_state.assemblies:
        st.info("No assemblies registered yet. Use the 'Find / Add' tab.")
    else:
        # Summary table
        rows = []
        for sid, d in st.session_state.assemblies.items():
            s = d.get("stats", {})
            rows.append(
                {
                    "Strain ID":    sid,
                    "Total (Mb)":   s.get("total_length_mb", "—"),
                    "Contigs":      s.get("num_contigs", "—"),
                    "N50 (bp)":     s.get("n50", "—"),
                    "GC (%)":       s.get("mean_gc", "—"),
                    "Assembler":    s.get("assembler", "—"),
                    "QUAST":        "✓" if s.get("quast") else "✗",
                    "R1 linked":    "✓" if d.get("reads_r1") else "✗",
                    "R2 linked":    "✓" if d.get("reads_r2") else "✗",
                }
            )
        st.dataframe(
            pd.DataFrame(rows),
            use_container_width=True,
            hide_index=True,
        )

        st.markdown("---")
        st.subheader("Edit read paths / remove assembly")

        selected = st.selectbox(
            "Select strain to edit",
            list(st.session_state.assemblies.keys()),
        )
        if selected:
            d = st.session_state.assemblies[selected]
            st.caption(f"Assembly: `{d['assembly_path']}`")

            with st.form(f"edit_reads_{selected}"):
                r1 = st.text_input("R1 path", value=d.get("reads_r1", ""))
                r2 = st.text_input("R2 path", value=d.get("reads_r2", ""))
                save_btn = st.form_submit_button("💾 Save read paths")

            if save_btn:
                st.session_state.assemblies[selected]["reads_r1"] = r1
                st.session_state.assemblies[selected]["reads_r2"] = r2
                _save()
                st.success("Read paths saved.")

            if st.button(f"🗑️ Remove '{selected}' from registry"):
                del st.session_state.assemblies[selected]
                _save()
                st.rerun()


# ════════════════════════════════════════════════════════
# TAB 3 — Assembly Stats
# ════════════════════════════════════════════════════════
with tab_stats:
    if not st.session_state.assemblies:
        st.info("No assemblies registered yet.")
    else:
        selected_stat = st.selectbox(
            "Select assembly",
            list(st.session_state.assemblies.keys()),
            key="stats_select",
        )

        if selected_stat:
            d = st.session_state.assemblies[selected_stat]
            s = d.get("stats", {})

            # Key metrics
            m1, m2, m3, m4 = st.columns(4)
            m1.metric("Total Size",     f"{s.get('total_length_mb', '?')} Mb")
            m2.metric("Contigs",         s.get("num_contigs", "?"))
            n50_val = s.get("n50", None)
            m3.metric("N50",             f"{n50_val:,} bp" if n50_val else "?")
            m4.metric("Mean GC",         f"{s.get('mean_gc', '?')} %")

            m5, m6, m7, m8 = st.columns(4)
            m5.metric("L50",             s.get("l50", "?"))
            lg = s.get("largest_contig", None)
            m6.metric("Largest Contig",  f"{lg:,} bp" if lg else "?")
            m7.metric("Assembler",       s.get("assembler", "unknown"))
            hcov = s.get("mean_coverage_header", None)
            m8.metric("Header Cov. (×)", hcov if hcov else "n/a")

            # ── QUAST report (auto-discovered alongside the assembly) ──
            quast = s.get("quast")
            if quast:
                st.markdown("---")
                st.markdown("**QUAST report**")
                st.caption(f"Source: `{s.get('quast_report', '—')}`")
                quast_rows = [
                    {"Metric": k, "Value": quast[k]}
                    for k in QUAST_DISPLAY_KEYS if k in quast
                ]
                if quast_rows:
                    st.dataframe(pd.DataFrame(quast_rows),
                                 use_container_width=True, hide_index=True)
                with st.expander("📄 Full QUAST metrics"):
                    st.dataframe(
                        pd.DataFrame(
                            [{"Metric": k, "Value": v} for k, v in quast.items()]
                        ),
                        use_container_width=True, hide_index=True,
                    )

            contigs = s.get("contigs", [])
            if contigs:
                df_c = pd.DataFrame(contigs)

                st.markdown("---")
                col_a, col_b = st.columns(2)

                with col_a:
                    fig = px.histogram(
                        df_c, x="length", nbins=60, log_y=True,
                        title="Contig Length Distribution",
                        labels={"length": "Contig length (bp)", "count": "Count"},
                        color_discrete_sequence=["steelblue"],
                    )
                    fig.update_layout(showlegend=False)
                    st.plotly_chart(fig, use_container_width=True)

                with col_b:
                    fig2 = px.histogram(
                        df_c, x="gc_percent", nbins=50,
                        title="GC Content Distribution",
                        labels={"gc_percent": "GC content (%)", "count": "Count"},
                        color_discrete_sequence=["seagreen"],
                    )
                    fig2.update_layout(showlegend=False)
                    st.plotly_chart(fig2, use_container_width=True)

                # If header coverage is available, show GC vs coverage scatter
                if df_c["header_coverage"].notna().any():
                    fig3 = px.scatter(
                        df_c.dropna(subset=["header_coverage"]),
                        x="gc_percent", y="header_coverage",
                        size="length",
                        hover_data=["contig_id", "length"],
                        log_y=True,
                        title="GC% vs Coverage (from contig headers)",
                        labels={
                            "gc_percent": "GC content (%)",
                            "header_coverage": "Coverage (×, log scale)",
                        },
                        color_discrete_sequence=["darkorange"],
                        size_max=30,
                    )
                    st.plotly_chart(fig3, use_container_width=True)

                with st.expander("📄 Per-contig table"):
                    st.dataframe(df_c, use_container_width=True, hide_index=True)
