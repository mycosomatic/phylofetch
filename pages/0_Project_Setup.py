"""
pages/0_Project_Setup.py
------------------------
Initialize a phylofetch workspace and register genome assemblies.

Accepts assemblies from any source — EGAP, SPAdes, Flye, Hifiasm, or any
other assembler. BUSCO summaries and QUAST reports are linked optionally.
No upstream pipeline dependency is assumed.
"""

import sys
from pathlib import Path

import pandas as pd
import streamlit as st

sys.path.insert(0, str(Path(__file__).parent.parent / "src"))
from phylofetch.assembly_utils import get_assembly_stats, suggest_strain_id
from phylofetch.config import load_config, save_config
from phylofetch.project_manager import (
    DEFAULT_PROJECT_DIR, DEFAULT_PROJECTS_ROOT, RunManager, check_tools,
    init_project, list_projects, load_assembly_registry, load_json,
    safe_slug, save_assembly_registry,
)

st.set_page_config(page_title="Project Setup", page_icon="⚙️", layout="wide")
st.title("⚙️ Project Setup")

if "assemblies" not in st.session_state:
    st.session_state.assemblies = load_config().get("assemblies", {})
if "tool_paths" not in st.session_state:
    st.session_state.tool_paths = load_config().get("tool_paths", {})
if "project_dir" not in st.session_state:
    st.session_state.project_dir = load_config().get(
        "project_dir", str(DEFAULT_PROJECT_DIR)
    )


def _persist_assemblies(assemblies: dict) -> None:
    """Save the registry to the global cache AND the active project."""
    save_config({"assemblies": assemblies})
    save_assembly_registry(st.session_state.project_dir, assemblies)

tab_proj, tab_import, tab_tools, tab_history = st.tabs(
    ["📂 Project", "📥 Import Assemblies", "🔧 Tool Status", "📋 Command History"]
)

# ══════════════════════════════════════════════════════════════════════════════
# TAB 1 — Project
# ══════════════════════════════════════════════════════════════════════════════
with tab_proj:
    st.subheader("Projects")
    st.caption(
        "A project is a re-openable workspace that durably stores its assembly "
        "registry, file-location manifest, run logs, and command history under "
        f"`{DEFAULT_PROJECTS_ROOT}`."
    )

    existing = list_projects()
    active_dir = st.session_state.project_dir
    active_name = Path(active_dir).name
    st.markdown(f"**Active project:** `{active_name}`  ·  `{active_dir}`")

    col_open, col_new = st.columns(2)

    # ── Open an existing project ──
    with col_open:
        st.markdown("**📂 Open existing**")
        if existing:
            labels = {
                f"{p['name']}  ({p['n_assemblies']} assemblies · {p['created_at'][:10]})": p
                for p in existing
            }
            choice = st.selectbox("Project", list(labels.keys()),
                                  label_visibility="collapsed")
            if st.button("Open project"):
                proj = labels[choice]
                reg = load_assembly_registry(proj["path"])
                st.session_state.project_dir = proj["path"]
                st.session_state.assemblies = reg
                save_config({"project_dir": proj["path"], "assemblies": reg})
                st.success(f"Opened `{proj['name']}` — {len(reg)} assemblies.")
                st.rerun()
        else:
            st.info("No saved projects yet. Create one →")

    # ── Create a new project ──
    with col_new:
        st.markdown("**➕ Create new**")
        new_name = st.text_input("Project name", placeholder="alternaria_2026",
                                 label_visibility="collapsed")
        if st.button("Create project", type="primary"):
            slug = safe_slug(new_name) if new_name.strip() else ""
            if not slug:
                st.error("Enter a project name.")
            else:
                target = DEFAULT_PROJECTS_ROOT / slug
                if target.exists():
                    st.warning(f"Project `{slug}` already exists — open it instead.")
                else:
                    root = init_project(target)
                    st.session_state.project_dir = str(root)
                    st.session_state.assemblies = {}
                    save_config({"project_dir": str(root), "assemblies": {}})
                    st.success(f"Created `{slug}`.")
                    st.rerun()

    with st.expander("⚙️ Advanced: use a custom directory"):
        new_dir = st.text_input("Project directory", value=st.session_state.project_dir)
        if st.button("Initialize / Open path"):
            root = init_project(new_dir)
            reg = load_assembly_registry(root)
            st.session_state.project_dir = str(root)
            st.session_state.assemblies = reg or st.session_state.assemblies
            save_config({"project_dir": str(root), "assemblies": st.session_state.assemblies})
            st.success(f"Workspace ready: `{root}`")
            st.rerun()

    project_dir = st.session_state.project_dir
    manifest_path = Path(project_dir) / "metadata" / "project_manifest.json"
    if manifest_path.exists():
        manifest = load_json(manifest_path, {})
        st.info(
            f"Created: {manifest.get('created_at', '?')}  ·  "
            f"Schema v{manifest.get('schema_version', '?')}",
            icon="📂",
        )
    else:
        st.warning("Workspace not initialized.")

    st.markdown("---")
    st.subheader(f"Registered assemblies ({len(st.session_state.assemblies)})")
    if st.session_state.assemblies:
        rows = []
        for sid, v in st.session_state.assemblies.items():
            s = v.get("stats", {}) if isinstance(v.get("stats"), dict) else {}
            rows.append({
                "ID": sid,
                "Assembly path": v.get("assembly_path", ""),
                "N50": v.get("n50") or s.get("n50", ""),
                "Contigs": v.get("num_contigs") or s.get("num_contigs", ""),
                "Size (Mb)": v.get("total_length_mb") or s.get("total_length_mb", ""),
                "QUAST": "✓" if s.get("quast") else "",
            })
        st.dataframe(pd.DataFrame(rows), use_container_width=True, hide_index=True)

        manifest_tsv = Path(project_dir) / "metadata" / "assembly_manifest.tsv"
        if manifest_tsv.exists():
            with st.expander("📄 Assembly file manifest (all paths)"):
                st.caption(f"`{manifest_tsv}`")
                st.code(manifest_tsv.read_text(), language=None)
                st.download_button("⬇️ Download manifest.tsv",
                                   data=manifest_tsv.read_bytes(),
                                   file_name="assembly_manifest.tsv",
                                   mime="text/tab-separated-values")
    else:
        st.info("No assemblies registered yet. Use the Import tab.")

# ══════════════════════════════════════════════════════════════════════════════
# TAB 2 — Import Assemblies
# ══════════════════════════════════════════════════════════════════════════════
with tab_import:
    st.subheader("Scan a directory for FASTA assemblies")
    st.caption(
        "Point to any directory containing genome assembly FASTAs. "
        "BUSCO/Compleasm and QUAST paths can be linked per-assembly after import."
    )

    scan_dir = st.text_input(
        "Directory to scan",
        placeholder="/data/my_assemblies/",
        key="scan_dir",
    )
    exts = st.multiselect("FASTA extensions", [".fasta", ".fa", ".fna"],
                          default=[".fasta", ".fa", ".fna"])

    if st.button("🔍 Scan for assemblies") and scan_dir:
        if not Path(scan_dir).is_dir():
            st.error(f"Directory not found: {scan_dir}")
        else:
            found_files = []
            for ext in exts:
                found_files.extend(Path(scan_dir).rglob(f"*{ext}"))
            found_files = sorted(set(found_files))

            if not found_files:
                st.warning("No FASTA files found.")
            else:
                st.success(f"Found {len(found_files)} FASTA file(s).")
                rows = []
                for fp in found_files:
                    suggested_id = suggest_strain_id(str(fp))
                    rows.append({
                        "Import?":      True,
                        "Strain ID":    suggested_id,
                        "File":         str(fp),
                        "BUSCO dir":    "",
                        "QUAST report": "",
                    })
                st.session_state["scan_rows"] = rows

    if "scan_rows" in st.session_state:
        edited = st.data_editor(
            pd.DataFrame(st.session_state["scan_rows"]),
            use_container_width=True, hide_index=True,
            column_config={
                "Import?":      st.column_config.CheckboxColumn(),
                "Strain ID":    st.column_config.TextColumn(),
                "File":         st.column_config.TextColumn(disabled=True),
                "BUSCO dir":    st.column_config.TextColumn(
                    help="Optional: path to BUSCO/Compleasm run directory for this sample"),
                "QUAST report": st.column_config.TextColumn(
                    help="Optional: path to QUAST report.tsv for this sample"),
            },
        )

        to_import = edited[edited["Import?"] == True]
        st.write(f"**{len(to_import)}** selected for import")

        if st.button("⬇️ Import selected", type="primary") and not to_import.empty:
            assemblies = dict(st.session_state.assemblies)
            imported, skipped = 0, 0
            for _, row in to_import.iterrows():
                sid = str(row["Strain ID"]).strip()
                fp  = str(row["File"]).strip()
                if not sid or not Path(fp).exists():
                    skipped += 1
                    continue
                if sid in assemblies:
                    st.warning(f"⚠️  {sid} already registered — skipping. Remove first to reimport.")
                    skipped += 1
                    continue
                with st.spinner(f"Computing stats for {sid}…"):
                    stats = get_assembly_stats(fp)
                assemblies[sid] = {
                    "assembly_path":   fp,
                    "assembler":       stats["assembler"],
                    "num_contigs":     stats["num_contigs"],
                    "total_length_mb": stats["total_length_mb"],
                    "n50":             stats["n50"],
                    "gc_percent":      stats["mean_gc"],
                    "busco_dir":       str(row["BUSCO dir"]).strip() or "",
                    "quast_report":    str(row["QUAST report"]).strip() or "",
                }
                imported += 1

            st.session_state.assemblies = assemblies
            _persist_assemblies(assemblies)
            st.success(f"Imported {imported} assemblies. Skipped {skipped}.")
            if "scan_rows" in st.session_state:
                del st.session_state["scan_rows"]
            st.rerun()

    st.markdown("---")
    st.subheader("Remove an assembly")
    if st.session_state.assemblies:
        to_remove = st.selectbox("Select assembly to remove",
                                 [""] + list(st.session_state.assemblies.keys()))
        if to_remove and st.button(f"🗑️ Remove {to_remove}"):
            assemblies = dict(st.session_state.assemblies)
            del assemblies[to_remove]
            st.session_state.assemblies = assemblies
            _persist_assemblies(assemblies)
            st.success(f"Removed {to_remove}.")
            st.rerun()

# ══════════════════════════════════════════════════════════════════════════════
# TAB 3 — Tool Status
# ══════════════════════════════════════════════════════════════════════════════
with tab_tools:
    st.subheader("Tool availability")
    st.caption("Checks PATH for each tool and probes its version string.")

    tp = st.session_state.tool_paths
    core_tools = {
        "ITSx":    tp.get("itsx", "ITSx"),
        "blastn":  tp.get("blastn", "blastn"),
        "tblastn": tp.get("tblastn", "tblastn"),
        "MAFFT":   tp.get("mafft", "mafft"),
        "trimAl":  tp.get("trimal", "trimal"),
    }
    optional_tools = {
        "IQ-TREE2": tp.get("iqtree2", "iqtree2"),
    }

    if st.button("🔍 Check tools"):
        rows = []
        for label, exe in {**core_tools, **optional_tools}.items():
            statuses = check_tools({label: exe})
            s = statuses[0]
            rows.append({
                "Tool":    s.label,
                "Status":  "✅ Found" if s.available else "❌ Not found",
                "Version": s.version or "—",
                "Path":    s.resolved_path or "—",
            })
        # MACSE JAR
        macse_jar = tp.get("macse_jar", "")
        rows.append({
            "Tool":    "MACSE (JAR)",
            "Status":  "✅ Found" if macse_jar and Path(macse_jar).exists() else "❌ Not configured",
            "Version": "—",
            "Path":    macse_jar or "—",
        })
        st.dataframe(pd.DataFrame(rows), use_container_width=True, hide_index=True)

        missing_core = [
            label for label, exe in core_tools.items()
            if not check_tools({label: exe})[0].available
        ]
        if missing_core:
            st.warning(f"Missing core tools: {', '.join(missing_core)}. "
                       f"Install via conda: `conda install -c bioconda blast itsx mafft trimal`")

# ══════════════════════════════════════════════════════════════════════════════
# TAB 4 — Command History
# ══════════════════════════════════════════════════════════════════════════════
with tab_history:
    project_dir = st.session_state.project_dir
    history_tsv = Path(project_dir) / "metadata" / "command_history.tsv"

    st.subheader("Command history")
    st.caption(f"From: `{history_tsv}`")

    if not history_tsv.exists():
        st.info("No command history yet. Run an extraction or alignment to populate this.")
    else:
        import csv
        rows = []
        with open(history_tsv, newline="") as f:
            reader = csv.DictReader(f, delimiter="\t")
            rows = list(reader)

        if not rows:
            st.info("History file exists but is empty.")
        else:
            df = pd.DataFrame(rows)
            df_display = df[["started_at", "module", "action", "returncode", "command"]].copy()
            df_display["returncode"] = df_display["returncode"].astype(str)
            st.dataframe(df_display, use_container_width=True, hide_index=True)

            st.markdown("---")
            st.subheader("Inspect a run folder")
            run_ids = df["run_id"].tolist()
            selected_run = st.selectbox("Select run", [""] + run_ids)
            if selected_run:
                run_dir = Path(project_dir) / "runs" / selected_run
                terminal_log = run_dir / "terminal.log"
                if terminal_log.exists():
                    with st.expander("📄 terminal.log", expanded=True):
                        st.code(terminal_log.read_text(), language=None)
                cmd_json = run_dir / "command.json"
                if cmd_json.exists():
                    with st.expander("📋 command.json"):
                        st.json(load_json(cmd_json, {}))
                env_json = run_dir / "environment.json"
                if env_json.exists():
                    with st.expander("🔬 environment.json (tool versions + env)"):
                        st.json(load_json(env_json, {}))
