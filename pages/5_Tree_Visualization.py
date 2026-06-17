"""
pages/5_Tree_Visualization.py
-------------------------------
IQ-TREE2 runner + tree visualization stub.

Currently working:
  - Run IQ-TREE2 on a supermatrix + partition nexus
  - Display Newick string and raw log

In development (roadmap):
  - Parse ML tree → render with Bio.Phylo / toytree
  - Overlay Bayesian posteriors (BEAST / MrBayes import)
  - Export publication-ready SVG/PDF
"""

import csv
import subprocess
import sys
from pathlib import Path

import streamlit as st

sys.path.insert(0, str(Path(__file__).parent.parent / "src"))
from phylofetch.config import load_config
from phylofetch.project_manager import RunManager, load_json

st.set_page_config(
    page_title="Tree Visualization", page_icon="🌳", layout="wide"
)
st.title("🌳 Tree Visualization")
st.caption("Run IQ-TREE2 and inspect ML trees. Bayesian posterior overlay in development.")

if "tool_paths" not in st.session_state:
    st.session_state.tool_paths = load_config().get("tool_paths", {})
if "project_dir" not in st.session_state:
    st.session_state.project_dir = load_config().get(
        "project_dir",
        str(Path.home() / ".phylofetch" / "projects" / "default"),
    )

tp = st.session_state.tool_paths
project_dir = Path(st.session_state.project_dir)
run_manager = RunManager(project_dir)

tab_iqtree, tab_tree, tab_roadmap = st.tabs(
    ["🧮 IQ-TREE2", "🌿 Tree View", "🗺️ Roadmap"]
)

# ══════════════════════════════════════════════════════════════════════════════
# TAB 1 — IQ-TREE2
# ══════════════════════════════════════════════════════════════════════════════
with tab_iqtree:
    st.subheader("Run IQ-TREE2")

    col_i1, col_i2 = st.columns([2, 1])
    with col_i1:
        supermatrix = st.text_input(
            "Supermatrix FASTA",
            placeholder="/path/to/supermatrix.fasta",
            key="iqtree_fasta",
        )
        partition_nex = st.text_input(
            "Partition nexus file (optional)",
            placeholder="/path/to/supermatrix_partitions.nex",
            key="iqtree_partition",
        )
        iqtree_outdir = st.text_input(
            "Output directory",
            value=str(Path(supermatrix).parent / "iqtree") if supermatrix else "",
            key="iqtree_outdir",
        )

    with col_i2:
        st.markdown("**IQ-TREE2 options**")
        iqtree_model = st.text_input(
            "Model / model testing",
            value="MFP+MERGE",
            help=(
                "MFP+MERGE: ModelFinder + merge partition scheme (recommended). "
                "TEST: full model testing without merging. "
                "Or specify a fixed model e.g. GTR+G"
            ),
        )
        iqtree_bootstrap = st.number_input(
            "Ultrafast bootstrap replicates",
            min_value=0,
            max_value=10000,
            value=1000,
            step=100,
            help="0 to disable bootstrapping.",
        )
        iqtree_threads = st.selectbox(
            "Threads",
            options=["AUTO"] + [str(i) for i in range(1, 33)],
            index=0,
        )
        iqtree_extra = st.text_input(
            "Extra IQ-TREE2 args",
            value="",
            placeholder="e.g. --redo --alrt 1000",
        )

    iqtree_bin = tp.get("iqtree2", "iqtree2")

    if st.button("▶️ Run IQ-TREE2", type="primary", key="btn_iqtree"):
        if not supermatrix or not Path(supermatrix).exists():
            st.error("Supermatrix FASTA not found.")
        else:
            out_path = Path(iqtree_outdir)
            out_path.mkdir(parents=True, exist_ok=True)
            prefix = str(out_path / Path(supermatrix).stem)

            cmd = [iqtree_bin, "-s", supermatrix]
            if partition_nex and Path(partition_nex).exists():
                cmd += ["-p", partition_nex]
            cmd += ["-m", iqtree_model]
            if int(iqtree_bootstrap) > 0:
                cmd += ["-B", str(iqtree_bootstrap)]
            cmd += ["-T", str(iqtree_threads), "--prefix", prefix]
            if iqtree_extra.strip():
                cmd += iqtree_extra.split()

            with st.spinner("IQ-TREE2 running… (this may take a while)"):
                res = run_manager.run(
                    cmd,
                    module="iqtree2",
                    action="ml_tree",
                )

            if res.returncode == 0:
                st.success("IQ-TREE2 complete.")
                treefile = Path(prefix + ".treefile")
                if treefile.exists():
                    tree_text = treefile.read_text().strip()
                    st.session_state["newick"] = tree_text
                    st.session_state["iqtree_log"] = Path(res.stdout_path).read_text()
                    st.markdown("**Newick tree**")
                    st.code(tree_text, language=None)
                    st.download_button(
                        "⬇️ Download .treefile",
                        data=tree_text,
                        file_name=treefile.name,
                        mime="text/plain",
                    )
            else:
                st.error(f"IQ-TREE2 failed (return code {res.returncode}).")

            log_lines = Path(res.stdout_path).read_text() + Path(res.stderr_path).read_text()
            if log_lines.strip():
                with st.expander("IQ-TREE2 log"):
                    st.code(log_lines[-8000:], language=None)

    # Show previously run trees from project history
    st.markdown("---")
    st.subheader("Previous IQ-TREE2 runs")
    history_tsv = project_dir / "metadata" / "command_history.tsv"
    if history_tsv.exists():
        rows = []
        with open(history_tsv, newline="") as f:
            reader = csv.DictReader(f, delimiter="\t")
            for row in reader:
                if row.get("module") == "iqtree2":
                    rows.append(row)
        if rows:
            import pandas as pd
            df = pd.DataFrame(rows)[
                ["started_at", "action", "returncode", "command"]
            ]
            st.dataframe(df, use_container_width=True, hide_index=True)
        else:
            st.info("No IQ-TREE2 runs in history.")
    else:
        st.info("No run history yet.")

# ══════════════════════════════════════════════════════════════════════════════
# TAB 2 — Tree view
# ══════════════════════════════════════════════════════════════════════════════
with tab_tree:
    st.subheader("Tree View")
    st.caption("Basic Newick rendering via Biopython. Interactive rendering coming soon.")

    newick_input = st.text_area(
        "Paste Newick string or run IQ-TREE2 above",
        value=st.session_state.get("newick", ""),
        height=100,
        key="newick_input",
    )

    treefile_path = st.text_input(
        "Or load a .treefile",
        placeholder="/path/to/supermatrix.treefile",
        key="treefile_path",
    )
    if treefile_path and Path(treefile_path).exists():
        newick_input = Path(treefile_path).read_text().strip()
        st.session_state["newick"] = newick_input

    if newick_input and st.button("🌿 Render tree", key="btn_render"):
        try:
            from io import StringIO

            from Bio import Phylo

            tree = Phylo.read(StringIO(newick_input), "newick")
            import matplotlib
            matplotlib.use("Agg")
            import matplotlib.pyplot as plt

            fig, ax = plt.subplots(figsize=(10, max(4, tree.count_terminals() * 0.25)))
            Phylo.draw(tree, axes=ax, do_show=False)
            ax.set_title("ML tree (IQ-TREE2)")
            st.pyplot(fig)
            plt.close(fig)
        except ImportError:
            st.warning(
                "matplotlib is not installed. Showing Newick text only. "
                "Install matplotlib for tree rendering: `pip install matplotlib`"
            )
            st.code(newick_input, language=None)
        except Exception as exc:
            st.error(f"Tree rendering error: {exc}")
            st.code(newick_input, language=None)

    if not newick_input:
        st.info("Run IQ-TREE2 in the IQ-TREE2 tab, or paste a Newick string above.")

# ══════════════════════════════════════════════════════════════════════════════
# TAB 3 — Roadmap
# ══════════════════════════════════════════════════════════════════════════════
with tab_roadmap:
    st.subheader("Development Roadmap")
    st.markdown(
        """
        ### Planned features

        #### Interactive tree rendering
        - Parse ML tree (Newick) with **toytree** or **ETE3**
        - Configurable label display (sample ID, support values, clade colors)
        - Pan/zoom via Plotly or Bokeh
        - Export as SVG / PDF for publication

        #### Bayesian posterior overlay
        - Import **BEAST2** or **MrBayes** posterior tree annotated with node probabilities
        - Map posterior support values onto the ML tree topology
        - Display dual-support nodes: bootstrap (IQ-TREE) + posterior (BEAST/MrBayes)
        - Export labeled tree figure

        #### IQ-TREE2 workflow integration
        - Partition scheme selection from ModelFinder output
        - Per-partition substitution model display
        - Branch support visualization (UFBoot, SH-aLRT)

        #### Concordance factors
        - Site concordance factors (sCF) and gene concordance factors (gCF)
        - Overlay on tree as pie charts per node

        ---
        > *To request a feature or report a bug, open an issue on GitHub.*
        """
    )
