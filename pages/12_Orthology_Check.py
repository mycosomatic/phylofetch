"""
pages/12_Orthology_Check.py
---------------------------
Orthology / paralog sanity check — component page (D-036).

Aligns each locus's combined matrix (your extracted isolates + the imported amplicon tips) and
surfaces divergence outliers — sequences that don't sit with the rest. A paralog or artifact
doesn't announce itself; it lands on a long branch, far from everything else. The check is
deliberately **source-blind**: it will indict an imported GenBank reference exactly as readily as
one of Exonerate's extractions, because the QC behind a published barcode is usually opaque.

Slots between Codon Tip Prep (which builds the with_tips matrices) and Alignment Prep — run it
before you commit to an alignment so a paralog never silently enters the tree.
"""

import re
import shutil
import sys
from pathlib import Path

import pandas as pd
import streamlit as st
from Bio import SeqIO

sys.path.insert(0, str(Path(__file__).parent.parent / "src"))
from phylofetch.config import load_config
from phylofetch.orthology_utils import is_isolate, orthology_check
from phylofetch.project_manager import (
    DEFAULT_PROJECT_DIR,
    RunManager,
    load_assembly_registry,
    project_output_dir,
)

st.set_page_config(page_title="Orthology Check", page_icon="🔍", layout="wide")
st.title("🔍 Orthology Check")
st.caption("Align each locus's isolates + imported tips and flag divergence outliers — possible "
           "paralogs or artifacts, whether chosen by Exonerate or hiding in a reference amplicon. "
           "Source-blind: a published reference is judged exactly like your own extraction.")

cfg = load_config()
project_dir = cfg.get("project_dir", str(DEFAULT_PROJECT_DIR))
tp = cfg.get("tool_paths", {})
registry = load_assembly_registry(project_dir)
isolate_ids = set(registry)

loci_dir = project_output_dir(project_dir) / "loci"
with_tips_dir = loci_dir / "with_tips"
ortho_dir = project_output_dir(project_dir) / "orthology"

st.caption(f"📁 Project: `{Path(project_dir).name}`  ·  matrices: "
           f"`{str(with_tips_dir).replace(str(Path.home()), '~')}`")

# ── Which loci have a with_tips matrix, and which products? ────────────────────
SUBSTRATES = {  # label → filename suffix
    "Genomic (intron-inclusive)": "genomic",
    "CDS (nucleotide coding)": "CDS",
    "Protein": "protein",
}


def _matrix(locus: str, product: str) -> Path:
    return with_tips_dir / f"{locus}_{product}_combined.fasta"


# Recover the locus by stripping the product suffix from the RIGHT — `split("_")[0]` would mangle
# any locus / gene-of-interest name containing an underscore (e.g. MAT1_1). (D-036 review)
_LOCUS_RE = re.compile(r"^(.+)_(?:genomic|CDS|protein)_combined\.(?:fasta|fa)$")
loci = sorted({m.group(1) for p in with_tips_dir.glob("*_combined.fasta")
               if (m := _LOCUS_RE.match(p.name))}) if with_tips_dir.is_dir() else []

if not loci:
    st.info("No combined matrices with tips found yet. Import comparison taxa on the **Reference "
            "Taxa** page and frame them on **Codon Tip Prep** first — this page checks the "
            "`with_tips/` matrices those steps produce (isolates + tips together).")
    st.stop()

# ── 1 · Locus + substrate ─────────────────────────────────────────────────────
st.subheader("1 · Locus & substrate")
c1, c2 = st.columns(2)
with c1:
    locus = st.selectbox("Locus", loci)
with c2:
    avail = [lab for lab, prod in SUBSTRATES.items() if _matrix(locus, prod).exists()]
    if not avail:
        st.warning(f"No combined matrix for {locus} in `with_tips/`.")
        st.stop()
    substrate_label = st.selectbox(
        "Compare on", avail, index=0,
        help="Genomic (default): every tip is included and it matches your genomic-vs-amplicon "
             "tree; introns add resolution but can obscure deep paralogy. Protein: most sensitive "
             "to paralogs but only the tips that could be codon-framed. CDS: nucleotide coding.")
substrate = SUBSTRATES[substrate_label]
matrix_fa = _matrix(locus, substrate)

# Count via the canonical is_isolate (single pass; no second hand-rolled heuristic to drift).
_ids = [r.id for r in SeqIO.parse(str(matrix_fa), "fasta")] if matrix_fa.exists() else []
n_recs = len(_ids)
n_iso = sum(1 for i in _ids if is_isolate(i, isolate_ids))
n_tip = n_recs - n_iso
st.caption(f"Matrix: `{matrix_fa.name}` · {n_recs} sequences ({n_iso} isolate, {n_tip} tip)")
if n_recs and n_tip == 0:
    st.warning("This matrix has **no tips** — only your isolates. Import comparison taxa on "
               "**Reference Taxa** and frame them on **Codon Tip Prep** to compare against "
               "references; the check still flags isolate-vs-isolate outliers, but the "
               "amplicon-reference comparison this page is for won't happen.")

# ── 2 · Settings ──────────────────────────────────────────────────────────────
st.subheader("2 · Settings")
sa, sb = st.columns(2)
with sa:
    mafft_bin = st.text_input("mafft", value=tp.get("mafft", "mafft"))
    mafft_mode = st.selectbox("MAFFT mode", ["auto", "localpair", "globalpair", "genafpair",
                                             "retree2"], index=0)
with sb:
    threads = st.slider("Threads", 1, 32, 4)
    sensitivity = st.slider(
        "Outlier sensitivity (lower = flags more)", 2.0, 6.0, 3.5, step=0.5,
        help="Robust (median/MAD) threshold on each sequence's median distance to the others. "
             "Lower flags more aggressively; the tree + table are there for your own call.")

mafft_ok = shutil.which(mafft_bin) is not None
st.caption(f"Tool check: {'✅' if mafft_ok else '❌'} mafft")
if not mafft_ok:
    st.error("⚠️ **mafft not found.** This page needs MAFFT to align before comparing. "
             "Install: `conda install -c bioconda mafft`.")

# ── 3 · Run ───────────────────────────────────────────────────────────────────
st.subheader("3 · Run")
if st.button("🔍 Align & check orthology", type="primary", disabled=not mafft_ok):
    manager = RunManager(project_dir)
    with st.spinner(f"Aligning {locus} ({substrate}) and scoring divergence…"):
        report = orthology_check(str(matrix_fa), isolate_ids, mafft_bin=mafft_bin,
                                 mode=mafft_mode, threads=int(threads), manager=manager,
                                 k=float(sensitivity), protein=(substrate == "protein"))
    report["locus"] = locus
    report["substrate"] = substrate
    # Persist newick + a TSV report for provenance / re-use on the Tree page.
    if report.get("newick"):
        ortho_dir.mkdir(parents=True, exist_ok=True)
        (ortho_dir / f"{locus}_{substrate}.nwk").write_text(report["newick"] + "\n")
        pd.DataFrame(report["rows"]).to_csv(
            ortho_dir / f"{locus}_{substrate}_report.tsv", sep="\t", index=False)
    st.session_state["ortho_report"] = report
    st.rerun()

# ── 4 · Results ───────────────────────────────────────────────────────────────
report = st.session_state.get("ortho_report")
if report:
    st.markdown("---")
    if report.get("error"):
        st.error(f"Could not run the check: {report['error']}")
    else:
        st.subheader(f"Results — {report['locus']} ({report['substrate']})")
        for w in report.get("warnings", []):
            st.warning("⚠️ " + w)
        if report.get("n_low_overlap"):
            st.caption(f"ℹ️ {report['n_low_overlap']} sequence(s) overlapped too little to assess "
                       "(short/offset amplicons) — listed but not flagged.")
        nf = report["n_flagged"]
        n_assessed = report.get("n_assessed", report["n_seqs"])
        if nf:
            st.warning(f"⚠️ **{nf} of {report['n_seqs']}** sequence(s) flagged as divergence "
                       "outliers — inspect for paralogy / artifact before aligning for the tree.")
        elif n_assessed < 4:
            # The robust (median/MAD) test needs ≥4 assessable sequences; below that it can't flag
            # anything, so "no outliers" would be misleading — say so (D-036 review).
            st.info(f"ℹ️ Only **{n_assessed}** sequence(s) overlap enough to assess — the robust "
                    "outlier test needs ≥4, so nothing can be flagged here regardless of how "
                    "divergent a sequence is. Read the tree and the distance table directly below.")
        else:
            st.success(f"✅ No divergence outliers among {report['n_seqs']} sequences — nothing "
                       "stands apart from the rest at this sensitivity.")

        df = pd.DataFrame([{
            "": "⚠️" if r["flagged"] else "", "Sequence": r["id"], "Source": r["source"],
            "Organism": r["organism"], "Median dist": r["median_dist"],
            "Nearest": r["nn_id"], "NN dist": r["nn_dist"],
            "Branch len": r["branch_len"], "Note": r["reason"],
        } for r in report["rows"]])
        st.dataframe(df, width="stretch", hide_index=True)
        st.caption("Sorted most-divergent first. **Median dist** = this sequence's median "
                   "p-distance to all others (the paralog signal). **NN dist** = distance to its "
                   "closest relative — large means it stands alone; ~0 with another flagged "
                   "sequence means a shared-artifact cluster.")

        if report.get("ascii_tree"):
            with st.expander("🌳 Neighbour-joining tree (long branch = outlier)", expanded=bool(nf)):
                st.code(report["ascii_tree"])
                st.caption(f"Saved Newick + TSV to "
                           f"`{str(ortho_dir).replace(str(Path.home()), '~')}` — load the `.nwk` "
                           "on the Tree Visualization page for a proper rendering.")

        st.info("This is a QC guide, not a verdict. A flagged **isolate** may be an Exonerate "
                "paralog pick (re-extract / check the contig); a flagged **tip** may be a "
                "mis-annotated or paralogous reference amplicon — drop it, or trace its accession. "
                "Re-run on **Protein** for the most sensitive paralog separation where tips are "
                "framable. Note: the flag assumes **most** sequences are the true ortholog — if a "
                "locus is dominated by paralogs the majority-based call can invert, so trust the "
                "tree over the flag when in doubt.")
