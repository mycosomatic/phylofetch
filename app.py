"""
app.py — phylofetch
===================
Landing page and global tool-path configuration.
Run with:  streamlit run app.py
       or: phylofetch   (after pip install -e .)
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
    st.markdown("### 🧬 Loci Extraction")
    st.markdown(
        "Extract rDNA (ITS/LSU/SSU via ITSx) and protein-coding marker loci "
        "(blastn/tblastn HSP-as-exon). "
        "Outputs CDS, genomic, introns, GFF3, and codon-position partitions."
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
        "Run IQ-TREE2, parse ML trees, and overlay Bayesian posteriors "
        "for publication-ready figures."
    )
    st.caption("← In development")

st.markdown("---")
st.info(
    "**Getting started:** Open **Project Setup** to initialize a workspace "
    "and register assemblies, then run **Loci Extraction**.",
    icon="👈",
)

# ── Tool Settings ─────────────────────────────────────────────────────────────
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
        tp["iqtree2"] = st.text_input("IQ-TREE2",        value=tp.get("iqtree2", "iqtree2"))
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
                from pathlib import Path as _P
                ok = _P(tp["macse_jar"]).exists()
                st.write(f"{'✅' if ok else '❌'} MACSE JAR: {tp['macse_jar']}")
