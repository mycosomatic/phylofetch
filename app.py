"""
app.py — phylofetch
===================
Entry point: defines the grouped sidebar navigation (st.navigation, D-031) and the
Home landing page with global tool-path configuration.
Run with:  streamlit run app.py
       or: phylofetch   (after pip install -e .)

Navigation is built explicitly here rather than relying on Streamlit's filename-ordered
`pages/` auto-discovery, so the sidebar reads as the actual pipeline — phase-grouped, with
the Workflow orchestrator lifted to the top as the guided entry point. Pages are registered
with their on-disk relative paths so the Workflow page's `st.page_link(...)` calls keep
resolving. Each page still calls its own `st.set_page_config` (additive in Streamlit ≥1.58),
which sets that page's browser-tab title/icon; the `st.Page` title/icon below set the sidebar
label/icon.
"""

import sys
from pathlib import Path

import streamlit as st

sys.path.insert(0, str(Path(__file__).parent / "src"))
from phylofetch.config import load_config, save_config
from phylofetch.project_manager import check_tools

st.set_page_config(
    page_title="phylofetch",
    page_icon="🌿",
    layout="wide",
    initial_sidebar_state="expanded",
)

if "assemblies" not in st.session_state:
    st.session_state.assemblies = load_config().get("assemblies", {})

if "tool_paths" not in st.session_state:
    st.session_state.tool_paths = load_config().get("tool_paths", {
        "itsx":        "ITSx",
        "blastn":      "blastn",
        "tblastn":     "tblastn",
        "mafft":       "mafft",
        "trimal":      "trimal",
        "macse_jar":   "",
        "iqtree2":     "iqtree2",
        "output_base": str(Path.home() / "phylofetch_output"),
    })


# ── Home (landing) ────────────────────────────────────────────────────────────
def home() -> None:
    st.title("🌿 phylofetch")
    st.caption(
        "Fungal assembly → loci extraction → alignment → BUSCO phylogenomics → tree prep"
    )
    st.markdown("---")

    col1, col2, col3 = st.columns(3)
    with col1:
        n = len(st.session_state.assemblies)
        st.markdown("### 📁 Assembly Manager")
        st.markdown(
            "Register genome assemblies from any source. "
            "View per-strain N50, GC%, contig counts, and assembler metadata."
        )
        st.metric("Assemblies registered", n)

    with col2:
        st.markdown("### 🧩 Extract loci")
        st.markdown(
            "Three strategies: **rDNA** via ITSx (ITS/LSU/SSU), **coding loci** via "
            "Exonerate frame-safe spliced alignment, and **in-silico PCR** via primer "
            "pairs. Outputs CDS, genomic, introns, GFF3, and codon-position partitions."
        )
        st.caption("← Primary workflow")

    with col3:
        st.markdown("### 🔀 Alignment Prep")
        st.markdown(
            "Align per-locus combined FASTAs with MAFFT, trim with trimAl, "
            "optionally use MACSE for codon-aware alignment of CDS loci, "
            "then concatenate into a supermatrix with merged partition file."
        )

    col4, col5, _ = st.columns(3)
    with col4:
        st.markdown("### 🧫 BUSCO Phylogenomics")
        st.markdown(
            "Import BUSCO or Compleasm results. Build occupancy matrix. "
            "Export single-copy orthologs for supermatrix construction."
        )
    with col5:
        st.markdown("### 🌳 Tree Visualization")
        st.markdown(
            "Run IQ-TREE, parse ML trees, and overlay Bayesian posteriors "
            "for publication-ready figures."
        )
        st.caption("← In development")

    st.markdown("---")
    st.info(
        "**Getting started:** Open **Project Setup** to initialize a workspace and "
        "register assemblies, then open **Workflow** to pick an extraction strategy — "
        "it lights up each step in order as you complete it.",
        icon="👈",
    )

    # ── Tool Settings ─────────────────────────────────────────────────────────
    st.markdown("---")
    with st.expander("⚙️ Tool Settings", expanded=False):
        st.markdown(
            "Set paths to external tools. Defaults assume tools are on your PATH "
            "(e.g. via conda). MACSE requires a JAR file path."
        )
        tp = st.session_state.tool_paths
        col_a, col_b = st.columns(2)
        with col_a:
            st.subheader("Executables")
            tp["itsx"]    = st.text_input("ITSx",           value=tp.get("itsx", "ITSx"))
            tp["blastn"]  = st.text_input("blastn (BLAST+)", value=tp.get("blastn", "blastn"))
            tp["tblastn"] = st.text_input("tblastn",         value=tp.get("tblastn", "tblastn"))
            tp["mafft"]   = st.text_input("MAFFT",           value=tp.get("mafft", "mafft"))
            tp["trimal"]  = st.text_input("trimAl",          value=tp.get("trimal", "trimal"))
            tp["iqtree2"] = st.text_input("IQ-TREE",         value=tp.get("iqtree2", "iqtree2"))
        with col_b:
            st.subheader("Paths")
            tp["macse_jar"] = st.text_input(
                "MACSE JAR path",
                value=tp.get("macse_jar", ""),
                placeholder="/opt/macse/macse_v2.07.jar",
                help="Download from https://bioweb.supagro.inra.fr/macse/",
            )
            tp["output_base"] = st.text_input(
                "Base output directory",
                value=tp.get("output_base", str(Path.home() / "phylofetch_output")),
            )

        col_save, col_check = st.columns([1, 2])
        with col_save:
            if st.button("💾 Save settings"):
                st.session_state.tool_paths = tp
                save_config({"tool_paths": tp})
                st.success("Settings saved.")
        with col_check:
            if st.button("🔧 Check tool availability"):
                check_map = {k: v for k, v in tp.items()
                             if k not in ("macse_jar", "output_base") and v}
                statuses = check_tools(check_map)
                for s in statuses:
                    icon = "✅" if s.available else "❌"
                    ver  = f" — {s.version}" if s.available and s.version else ""
                    st.write(f"{icon} {s.label} ({s.executable}){ver}")
                if tp.get("macse_jar"):
                    ok = Path(tp["macse_jar"]).exists()
                    st.write(f"{'✅' if ok else '❌'} MACSE JAR: {tp['macse_jar']}")


# ── Grouped navigation (D-031) ────────────────────────────────────────────────
# Empty-string section ("") renders header-less at the top; remaining keys are phase
# headers. Page paths are the on-disk relative paths so Workflow's st.page_link resolves.
PAGES = "pages"
nav = st.navigation(
    {
        "": [
            st.Page(home, title="Home", icon="🏠", default=True),
            st.Page(f"{PAGES}/6_Workflow.py", title="Workflow", icon="🧭"),
        ],
        "Set up": [
            st.Page(f"{PAGES}/0_Project_Setup.py", title="Project Setup", icon="📁"),
            st.Page(f"{PAGES}/1_Assembly_Manager.py", title="Assembly Manager", icon="🧬"),
        ],
        "References (NCBI)": [
            st.Page(f"{PAGES}/2_NCBI_References.py", title="NCBI References", icon="🌐"),
            st.Page(f"{PAGES}/7_Reference_Taxa.py", title="Reference Taxa", icon="🌳"),
        ],
        "Extract loci": [
            st.Page(f"{PAGES}/3_ITSx_rDNA.py", title="ITSx · rDNA", icon="🧬"),
            st.Page(f"{PAGES}/4_Exonerate.py", title="Exonerate · coding", icon="🧩"),
            st.Page(f"{PAGES}/5_Primers.py", title="Primers · PCR", icon="🔬"),
        ],
        "Tree prep": [
            st.Page(f"{PAGES}/8_Codon_Tip_Prep.py", title="Codon Tip Prep", icon="🔡"),
            st.Page(f"{PAGES}/12_Orthology_Check.py", title="Orthology Check", icon="🔍"),
            st.Page(f"{PAGES}/9_Alignment_Prep.py", title="Alignment Prep", icon="🧷"),
        ],
        "Phylogenomics & tree": [
            st.Page(f"{PAGES}/10_BUSCO_Phylogenomics.py", title="BUSCO Phylogenomics", icon="🧫"),
            st.Page(f"{PAGES}/11_Tree_Visualization.py", title="Tree Visualization", icon="🌲"),
        ],
    },
    expanded=True,
)
nav.run()
