"""
pages/8_Codon_Tip_Prep.py
-------------------------
Codon Tip Prep — component page (RM-008 component 2 / D-022).

The "manual-Mesquite replacement". Comparison tips imported on the Reference Taxa page are raw
GenBank nucleotide barcode amplicons — for fungal coding markers usually *partial-CDS genomic*
records carrying introns, in an arbitrary reading frame, on either strand. This page runs each
coding-locus tip through the **same bundled protein guide** used for extraction (Exonerate
`protein2genome`): introns stripped, reading frame pinned to the guide ORF, strand oriented —
yielding sequences that line up, codon-by-codon, with your isolates' extracted CDS.

Per coding locus it writes three matrices (isolates + tips together) to
`<output>/loci/with_tips/`:
  • `<locus>_CDS_combined.fasta`      intron-stripped, codon-phased CDS  → codon-aware alignment
  • `<locus>_genomic_combined.fasta`  full gene (exons UPPER / introns lower) → nucleotide alignment
  • `<locus>_protein_combined.fasta`  translation                        → AA tree / alignment guide

It RUNS NO ALIGNER and adds no dependency beyond Exonerate (+ blastn for orientation). Alignment +
by-hand curation happen on the next page (Alignment Prep). rDNA tips (ITS/LSU/SSU) are non-coding,
have no guide, and are not handled here — they go straight to MAFFT.

Intron-rich barcodes (D-027): some markers — fungal TEF1 most notably — amplify a largely intronic
region that cannot be codon-framed. Those tips are oriented against the isolate genomic and added
to the `*_genomic_combined.fasta` matrix only (nucleotide-only), so the intron-inclusive tree still
carries them; the CDS/protein matrices stay isolate-only for that locus.
"""

import shutil
import sys
from pathlib import Path

import pandas as pd
import streamlit as st

sys.path.insert(0, str(Path(__file__).parent.parent / "src"))
from phylofetch.codon_prep_utils import coding_loci_with_tips, prepare_codon_locus
from phylofetch.config import load_config
from phylofetch.ncbi_utils import count_refs
from phylofetch.project_manager import (
    DEFAULT_PROJECT_DIR,
    RunManager,
    load_project_manifest,
    project_output_dir,
    update_step,
)
from phylofetch.protein_guide_utils import get_guides
from phylofetch.tips_utils import project_tips_dir

st.set_page_config(page_title="Codon Tip Prep", page_icon="🔡", layout="wide")
st.title("🔡 Codon Tip Prep")
st.caption("Frame comparison tips against the bundled protein guides (Exonerate) so they line up "
           "codon-by-codon with your extracted loci — the by-hand Mesquite step, automated. "
           "Produces codon-ready CDS + full-gene + protein matrices; runs no aligner.")

cfg = load_config()
project_dir = cfg.get("project_dir", str(DEFAULT_PROJECT_DIR))
tp = cfg.get("tool_paths", {})
manifest = load_project_manifest(project_dir)

tips_dir = project_tips_dir(project_dir)
loci_dir = project_output_dir(project_dir) / "loci"
combined_dir = loci_dir / "combined"          # isolate (extraction) outputs
with_tips_dir = loci_dir / "with_tips"        # our merged matrices

st.caption(f"📁 Project: `{Path(project_dir).name}`  ·  tips: "
           f"`{str(tips_dir).replace(str(Path.home()), '~')}`  ·  outputs: "
           f"`{str(with_tips_dir).replace(str(Path.home()), '~')}`")

coding_loci = coding_loci_with_tips(tips_dir)
if not coding_loci:
    st.info("No coding-locus tips to process yet. Import comparison taxa on the **Reference "
            "Taxa / Tips** page first — rDNA tips (ITS/LSU/SSU) are non-coding and skip this step "
            "(they align directly with MAFFT).")
    st.stop()

# ── 1 · Loci ──────────────────────────────────────────────────────────────────
st.subheader("1 · Coding loci with tips")
st.caption("Only loci that have both imported tips and a bundled protein guide are shown.")
sel_loci = st.multiselect(
    "Loci to frame", coding_loci, default=coding_loci,
    format_func=lambda loc: (f"{loc} · {count_refs(loc, ref_dir=tips_dir)} tip(s) · "
                             f"{len(get_guides(loc))} guide(s)"),
)

# ── 2 · Settings ──────────────────────────────────────────────────────────────
st.subheader("2 · Settings")
ca, cb = st.columns(2)
with ca:
    include_isolates = st.checkbox(
        "Merge with my extracted isolate loci", value=True,
        help="Prepend the isolates' CDS / full-gene / protein from the Exonerate extraction "
             f"(`{combined_dir.name}/`) so each matrix holds isolates + tips together.")
    mark_exons = st.checkbox(
        "Mark exon junctions in the CDS (case)", value=False,
        help="Alternate UPPER/lower per exon in the CDS so exon boundaries stay visible in the "
             "codon alignment. Case is inert to aligners; bases unchanged.")
    strict_qc = st.checkbox(
        "Strict QC (exclude frameshift / internal-stop tips)", value=False,
        help="Default: write-and-flag (imperfect tips are kept but flagged REVIEW, per D-008). "
             "Strict: leave them out of the matrices entirely.")
    nt_fallback = st.checkbox(
        "Keep un-framable tips as nucleotide (genomic matrix)", value=True,
        help="Some standard coding barcodes — fungal TEF1 above all — amplify a largely "
             "*intronic* region with too little exon for Exonerate to codon-frame. With this on, "
             "such a tip is oriented (blastn vs your isolate genomic, which also confirms the "
             "locus) and added to the genomic matrix only, flagged nucleotide-only — so the "
             "intron-inclusive tree still includes it. CDS/protein stay isolate-only for that "
             "locus. (D-027)")
with cb:
    exonerate_bin = st.text_input("exonerate", value=tp.get("exonerate", "exonerate"))
    blastn_bin = st.text_input("blastn", value=tp.get("blastn", "blastn"),
                               help="Used to orient nucleotide-only tips against the isolate genomic.")
    maxintron = st.number_input("Max intron (bp)", value=2000, min_value=50, max_value=200000,
                                step=100, help="Fungal introns are short; 2000 is a safe ceiling.")
    minintron = st.number_input("Min intron (bp)", value=20, min_value=4, max_value=200)
    geneticcode = st.number_input("Genetic code (NCBI)", value=1, min_value=1, max_value=33)

exonerate_ok = shutil.which(exonerate_bin) is not None
blastn_ok = shutil.which(blastn_bin) is not None
st.caption(f"Tool check: {'✅' if exonerate_ok else '❌'} exonerate · "
           f"{'✅' if blastn_ok else '❌'} blastn")
if not exonerate_ok:
    st.error("⚠️ **exonerate not found.** This step requires Exonerate (same tool as coding "
             "extraction). Install: `conda install -c bioconda exonerate`.")
if nt_fallback and not blastn_ok:
    st.warning("⚠️ **blastn not found** — the nucleotide fallback needs it to orient tips; "
               "un-framable tips will be reported but not added to the genomic matrix.")

# ── 3 · Run ───────────────────────────────────────────────────────────────────
st.subheader("3 · Run")
if st.button("🚀 Frame tips → codon-ready matrices", type="primary",
             disabled=not (sel_loci and exonerate_ok)):
    manager = RunManager(project_dir)
    with_tips_dir.mkdir(parents=True, exist_ok=True)
    prog = st.progress(0.0)
    summaries, outputs_prov = [], {}
    for i, locus in enumerate(sel_loci):
        with st.spinner(f"Framing {locus} tips…"):
            s = prepare_codon_locus(
                locus, tips_dir, str(with_tips_dir),
                include_isolates=include_isolates, isolate_combined_dir=str(combined_dir),
                exonerate_bin=exonerate_bin, blastn_bin=blastn_bin,
                minintron=int(minintron), maxintron=int(maxintron),
                geneticcode=int(geneticcode), mark_exons=mark_exons, strict_qc=strict_qc,
                nt_fallback=nt_fallback, manager=manager)
        summaries.append(s)
        if s["outputs"]:
            outputs_prov[locus] = s["outputs"]
        prog.progress((i + 1) / len(sel_loci))

    update_step(project_dir, "codon_prep", status="done", outputs=outputs_prov,
                notes=(f"loci={','.join(sel_loci)}; include_isolates={include_isolates}; "
                       f"strict_qc={strict_qc}; nt_fallback={nt_fallback}; "
                       f"framed={sum(x['n_framed'] for x in summaries)}, "
                       f"flagged={sum(x['n_flagged'] for x in summaries)}, "
                       f"nt_only={sum(x.get('n_nt_only', 0) for x in summaries)}, "
                       f"failed={sum(x['n_failed'] for x in summaries)}"))
    st.session_state["codon_prep_summaries"] = summaries
    st.success("Tips framed — provenance recorded in the manifest. Review QC below, then align "
               "(and hand-check) on the **Alignment Prep** page.")
    st.rerun()

# ── 4 · Results / QC ──────────────────────────────────────────────────────────
summaries = st.session_state.get("codon_prep_summaries")
if summaries:
    st.markdown("---")
    st.subheader("Results")
    overview = pd.DataFrame([{
        "Locus": s["locus"], "Tips": s["n_tips"], "Framed (CDS)": s["n_framed"],
        "Flagged (REVIEW)": s["n_flagged"], "Nucleotide-only": s.get("n_nt_only", 0),
        "Could not place": s["n_failed"], "Isolates merged": s["n_isolates"],
        "Status": s["status"],
    } for s in summaries])
    st.dataframe(overview, width="stretch", hide_index=True)
    if any(s.get("n_nt_only", 0) for s in summaries):
        st.caption("🧬 **Nucleotide-only** tips (e.g. intron-rich TEF1 barcodes) could not be "
                   "codon-framed but were oriented and added to the **genomic** matrix only — "
                   "they carry the comparison taxon into the intron-inclusive tree (D-027).")

    for s in summaries:
        if not s["rows"]:
            continue
        with st.expander(f"🔬 {s['locus']} — per-tip QC "
                         f"({s['n_framed']} framed · {s.get('n_nt_only', 0)} nucleotide-only "
                         f"of {s['n_tips']})"):
            df = pd.DataFrame([{
                "Accession": r["id"], "Organism": r["organism"], "Exons": r["n_exons"],
                "Introns": r["n_introns"], "Len (bp)": r["cds_length"],
                "QC": r["verdict"], "Product": r.get("product", ""),
                "Included": "✅" if r["included"] else "—", "Detail": r["status"],
            } for r in s["rows"]])
            st.dataframe(df, width="stretch", hide_index=True)
            if s["outputs"]:
                st.caption("Wrote: " + " · ".join(
                    f"`{Path(p).name}`" for p in s["outputs"].values()))

    st.markdown("---")
    st.info(
        "**Next — align, then check by hand.** These are *codon-phased* but **not codon-"
        "aligned**, and a plain nucleotide aligner does not guarantee reading frame is preserved. "
        "On **Alignment Prep**:\n"
        "- align `*_CDS_combined.fasta` codon-aware (MACSE if you have the JAR), or align the "
        "proteins and curate the CDS to match;\n"
        "- align `*_genomic_combined.fasta` with MAFFT (introns have no codon structure) — the "
        "exon UPPER / intron lower case marks the boundaries as you edit, and move with the gaps. "
        "This matrix also holds any **nucleotide-only** tips (intron-rich barcodes like TEF1); "
        "MAFFT `--adjustdirection` is a good belt-and-braces check on their orientation;\n"
        "- to get an exon/intron track in Geneious/IGV, load `LOCUS_genomic.fasta` and import the "
        "matching `LOCUS_genomic.gff3` (region-relative, lines up with that file) — **not** "
        "`exonerate_raw.txt` (that's the raw tool log, not valid GFF). `LOCUS.gff3` is the "
        "contig-relative version for whole-assembly context.\n\n"
        "Companions (optional, never required): **AliView** / **SeaView** (fast, frame-aware "
        "editors) or **Geneious**. Build a tree from the CDS and from the full gene and compare "
        "topologies — intron and coding signal can differ.")
