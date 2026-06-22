"""
pages/9_Alignment_Prep.py
--------------------------
Per-locus alignment, trimming, and supermatrix concatenation.

Workflow:
  1. Select combined FASTAs from Loci Extraction output
  2. Align with MAFFT (multiple modes)
  3. Trim with trimAl
  4. Optionally run MACSE for codon-aware CDS alignment
  5. Concatenate selected trimmed alignments into supermatrix + partition nexus
"""

import sys
from pathlib import Path

import pandas as pd
import streamlit as st
from Bio import SeqIO

sys.path.insert(0, str(Path(__file__).parent.parent / "src"))
from phylofetch.alignment.concat import concatenate_alignments
from phylofetch.alignment.macse import is_available as macse_available
from phylofetch.alignment.macse import run_macse
from phylofetch.alignment.mafft import run_mafft
from phylofetch.alignment.trimal import run_trimal
from phylofetch.config import load_config
from phylofetch.project_manager import RunManager, load_json, project_output_dir

st.set_page_config(
    page_title="Alignment Prep", page_icon="🔀", layout="wide"
)
st.title("🔀 Alignment Prep")
st.caption(
    "Align per-locus FASTAs with MAFFT, trim with trimAl, "
    "optionally use MACSE for codon-aware CDS alignment, "
    "then concatenate into a supermatrix."
)

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

tab_align, tab_trim, tab_concat, tab_logs = st.tabs(
    ["🧬 Alignment", "✂️ Trimming", "🔗 Concatenation", "📋 Logs"]
)

# ══════════════════════════════════════════════════════════════════════════════
# TAB 1 — Per-locus alignment with MAFFT
# ══════════════════════════════════════════════════════════════════════════════
with tab_align:
    st.subheader("Per-locus alignment (MAFFT)")
    st.caption(
        "Select one or more combined FASTA files (from Loci Extraction). "
        "Each file is aligned independently."
    )

    col_in, col_opts = st.columns([2, 1])
    with col_in:
        align_dir = st.text_input(
            "Directory containing combined FASTAs",
            value=str(project_output_dir(project_dir) / "loci" / "combined"),
            key="align_dir",
            help="Defaults to this project's output dir (set in Project Setup → Manage Data). "
                 "FASTAs are files matching *_combined.fasta or any .fasta files you choose.",
        )
        if align_dir and Path(align_dir).is_dir():
            candidates = sorted(
                list(Path(align_dir).rglob("*_combined.fasta"))
                + list(Path(align_dir).rglob("*_combined.fa"))
            )
            if not candidates:
                candidates = sorted(Path(align_dir).rglob("*.fasta")) + sorted(
                    Path(align_dir).rglob("*.fa")
                )
            if candidates:
                fasta_labels = [str(p) for p in candidates]
                selected_align = st.multiselect(
                    "Select FASTA files to align",
                    options=fasta_labels,
                    key="selected_fastas_align",
                )
            else:
                st.info("No FASTA files found in this directory.")
                selected_align = []
        else:
            selected_align = st.multiselect(
                "Enter FASTA paths manually (one per line) or select a directory above",
                options=[],
                key="selected_fastas_align_manual",
            )
            manual_input = st.text_area(
                "Or paste paths (one per line)",
                key="manual_fasta_paths",
                height=80,
            )
            if manual_input:
                selected_align = [
                    p.strip() for p in manual_input.splitlines() if p.strip()
                ]

    with col_opts:
        st.markdown("**MAFFT options**")
        mafft_mode = st.selectbox(
            "Alignment mode",
            options=["auto", "localpair", "globalpair", "genafpair", "retree2"],
            index=0,
            help=(
                "auto: MAFFT chooses; "
                "localpair: L-INS-i (slow, best for divergent); "
                "globalpair: G-INS-i; "
                "genafpair: E-INS-i (multiple conserved regions); "
                "retree2: FFT-NS-2 (fast, large datasets)"
            ),
        )
        mafft_threads = st.number_input(
            "Threads", min_value=1, max_value=64, value=4, key="mafft_threads"
        )
        mafft_extra = st.text_input(
            "Extra MAFFT args",
            value="",
            placeholder="e.g. --adjustdirectionaccurately",
            key="mafft_extra",
        )
        align_out_dir = st.text_input(
            "Output directory for aligned FASTAs",
            value=str(Path(align_dir).parent / "aligned") if align_dir else "",
            key="align_out_dir",
        )

    if st.button("▶️ Run MAFFT alignment", type="primary", key="btn_mafft"):
        if not selected_align:
            st.warning("Select at least one FASTA file.")
        elif not align_out_dir:
            st.warning("Specify an output directory.")
        else:
            out_path = Path(align_out_dir)
            out_path.mkdir(parents=True, exist_ok=True)
            results = {}
            progress = st.progress(0, text="Aligning…")
            for i, fasta_path in enumerate(selected_align):
                fp = Path(fasta_path)
                out_fp = out_path / (fp.stem + "_aligned.fasta")
                extra = mafft_extra.split() if mafft_extra.strip() else None
                rc, stderr = run_mafft(
                    input_fasta=fp,
                    output_fasta=out_fp,
                    mode=mafft_mode,
                    threads=int(mafft_threads),
                    mafft_bin=tp.get("mafft", "mafft"),
                    extra_args=extra,
                    run_manager=run_manager,
                )
                results[fp.name] = (rc, stderr)
                progress.progress(
                    (i + 1) / len(selected_align),
                    text=f"Aligned {fp.name}",
                )

            progress.empty()
            st.session_state["aligned_fastas"] = [
                str(out_path / (Path(p).stem + "_aligned.fasta"))
                for p in selected_align
            ]

            rows = []
            for fname, (rc, _stderr) in results.items():
                out_fp = out_path / (Path(fname).stem + "_aligned.fasta")
                n_seq, aln_len, gap_pct = 0, 0, 0.0
                if out_fp.exists():
                    seqs = list(SeqIO.parse(str(out_fp), "fasta"))
                    n_seq = len(seqs)
                    if seqs:
                        aln_len = len(seqs[0].seq)
                        total_gaps = sum(str(s.seq).count("-") for s in seqs)
                        total_chars = n_seq * aln_len
                        gap_pct = (
                            round(total_gaps / total_chars * 100, 1)
                            if total_chars
                            else 0.0
                        )
                rows.append(
                    {
                        "Locus": fname,
                        "Status": "✅" if rc == 0 else "❌",
                        "Sequences": n_seq,
                        "Aligned length": aln_len,
                        "Gap %": gap_pct,
                        "Output": str(out_fp),
                    }
                )
            st.dataframe(
                pd.DataFrame(rows), width="stretch", hide_index=True
            )

    if "aligned_fastas" in st.session_state:
        st.info(
            f"✅ {len(st.session_state['aligned_fastas'])} aligned FASTAs ready for trimming."
        )

# ══════════════════════════════════════════════════════════════════════════════
# TAB 2 — Trimming with trimAl
# ══════════════════════════════════════════════════════════════════════════════
with tab_trim:
    st.subheader("Per-locus trimming (trimAl)")

    col_t1, col_t2 = st.columns([2, 1])
    with col_t1:
        trim_dir = st.text_input(
            "Directory containing aligned FASTAs",
            value=(
                str(Path(st.session_state["aligned_fastas"][0]).parent)
                if "aligned_fastas" in st.session_state
                else ""
            ),
            key="trim_dir",
        )
        if trim_dir and Path(trim_dir).is_dir():
            trim_candidates = sorted(Path(trim_dir).rglob("*_aligned.fasta"))
            if not trim_candidates:
                trim_candidates = sorted(Path(trim_dir).rglob("*.fasta"))
            selected_trim = st.multiselect(
                "Select aligned FASTAs to trim",
                options=[str(p) for p in trim_candidates],
                default=(
                    st.session_state.get("aligned_fastas", [])
                    if "aligned_fastas" in st.session_state
                    else []
                ),
                key="selected_fastas_trim",
            )
        else:
            selected_trim = []

        trim_out_dir = st.text_input(
            "Output directory for trimmed FASTAs",
            value=str(Path(trim_dir).parent / "trimmed") if trim_dir else "",
            key="trim_out_dir",
        )

    with col_t2:
        st.markdown("**trimAl options**")
        trim_mode = st.selectbox(
            "Trimming mode",
            options=["automated1", "gappyout", "strict", "manual"],
            index=0,
            help=(
                "automated1: heuristic selection (recommended); "
                "gappyout: remove columns by gap content; "
                "strict: strict gap/conservation threshold; "
                "manual: set gt/cons below"
            ),
        )
        trim_gt = trim_cons = None
        if trim_mode == "manual":
            trim_gt = st.slider(
                "Gap threshold (gt)", min_value=0.0, max_value=1.0, value=0.5, step=0.05
            )
            trim_cons = st.slider(
                "Conservation threshold (cons)",
                min_value=0,
                max_value=100,
                value=60,
            )

    if st.button("▶️ Run trimAl", type="primary", key="btn_trimal"):
        if not selected_trim:
            st.warning("Select at least one aligned FASTA.")
        elif not trim_out_dir:
            st.warning("Specify an output directory.")
        else:
            out_path = Path(trim_out_dir)
            out_path.mkdir(parents=True, exist_ok=True)
            results = {}
            progress = st.progress(0, text="Trimming…")
            for i, fasta_path in enumerate(selected_trim):
                fp = Path(fasta_path)
                stem = fp.stem.replace("_aligned", "")
                out_fp = out_path / (stem + "_trimmed.fasta")
                rc, stderr = run_trimal(
                    input_fasta=fp,
                    output_fasta=out_fp,
                    mode=trim_mode,
                    gt=trim_gt,
                    cons=trim_cons,
                    trimal_bin=tp.get("trimal", "trimal"),
                    run_manager=run_manager,
                )
                results[fp.name] = (rc, stderr)
                progress.progress(
                    (i + 1) / len(selected_trim),
                    text=f"Trimmed {fp.name}",
                )

            progress.empty()
            st.session_state["trimmed_fastas"] = [
                str(out_path / (Path(p).stem.replace("_aligned", "") + "_trimmed.fasta"))
                for p in selected_trim
            ]

            rows = []
            for fname, (rc, _stderr) in results.items():
                stem = Path(fname).stem.replace("_aligned", "")
                out_fp = out_path / (stem + "_trimmed.fasta")
                n_seq, aln_len = 0, 0
                if out_fp.exists():
                    seqs = list(SeqIO.parse(str(out_fp), "fasta"))
                    n_seq = len(seqs)
                    if seqs:
                        aln_len = len(seqs[0].seq)
                rows.append(
                    {
                        "Locus": fname,
                        "Status": "✅" if rc == 0 else "❌",
                        "Sequences": n_seq,
                        "Trimmed length": aln_len,
                        "Output": str(out_fp),
                    }
                )
            st.dataframe(
                pd.DataFrame(rows), width="stretch", hide_index=True
            )

    # MACSE optional section
    st.markdown("---")
    st.subheader("Codon-aware alignment (MACSE — optional)")
    st.caption(
        "MACSE is a Java tool for codon-aware alignment of protein-coding loci. "
        "It handles frameshifts and produces both NT and AA output. "
        "Requires a JAR file (see Tool Settings on the home page)."
    )

    macse_jar = tp.get("macse_jar", "")
    if not macse_jar or not Path(macse_jar).exists():
        st.warning(
            "MACSE JAR not configured. "
            "Set the path in **Tool Settings** on the home page."
        )
    else:
        st.success(f"MACSE JAR: `{macse_jar}`")
        macse_input = st.text_input(
            "Input FASTA for MACSE (CDS sequences)",
            placeholder="/path/to/TEF1_combined.fasta",
            key="macse_input",
        )
        macse_out_dir = st.text_input(
            "Output directory",
            value=str(Path(macse_input).parent / "macse") if macse_input else "",
            key="macse_out_dir",
        )
        if st.button("▶️ Run MACSE", key="btn_macse"):
            if not macse_input or not Path(macse_input).exists():
                st.error("Input FASTA not found.")
            else:
                fp = Path(macse_input)
                out_path = Path(macse_out_dir)
                out_path.mkdir(parents=True, exist_ok=True)
                out_nt = out_path / (fp.stem + "_macse_nt.fasta")
                out_aa = out_path / (fp.stem + "_macse_aa.fasta")
                with st.spinner("Running MACSE…"):
                    rc, stderr = run_macse(
                        input_fasta=fp,
                        output_nt=out_nt,
                        output_aa=out_aa,
                        macse_jar=macse_jar,
                        run_manager=run_manager,
                    )
                if rc == 0:
                    st.success(f"MACSE complete. NT: `{out_nt}` AA: `{out_aa}`")
                else:
                    st.error("MACSE failed.")
                    st.code(stderr or "", language=None)

    if "trimmed_fastas" in st.session_state:
        st.info(
            f"✅ {len(st.session_state['trimmed_fastas'])} trimmed FASTAs ready for concatenation."
        )

# ══════════════════════════════════════════════════════════════════════════════
# TAB 3 — Concatenation
# ══════════════════════════════════════════════════════════════════════════════
with tab_concat:
    st.subheader("Concatenate alignments → supermatrix")
    st.caption(
        "Combine trimmed per-locus alignments into a supermatrix FASTA and "
        "a NEXUS partition file for IQ-TREE2 / RAxML. "
        "Missing taxa in a locus are filled with gap characters."
    )

    col_c1, col_c2 = st.columns([2, 1])
    with col_c1:
        concat_dir = st.text_input(
            "Directory containing trimmed (or aligned) FASTAs",
            value=(
                str(Path(st.session_state["trimmed_fastas"][0]).parent)
                if "trimmed_fastas" in st.session_state
                else ""
            ),
            key="concat_dir",
        )
        if concat_dir and Path(concat_dir).is_dir():
            concat_candidates = sorted(Path(concat_dir).rglob("*_trimmed.fasta"))
            if not concat_candidates:
                concat_candidates = sorted(Path(concat_dir).rglob("*.fasta"))
            selected_concat = st.multiselect(
                "Select FASTAs to concatenate (order matters for partition file)",
                options=[str(p) for p in concat_candidates],
                default=st.session_state.get("trimmed_fastas", []),
                key="selected_fastas_concat",
            )
        else:
            selected_concat = []

        codon_part_dir = st.text_input(
            "Directory with codon partition files (optional)",
            placeholder="/path/to/loci_dir/  (contains *_partition.nex files)",
            key="codon_part_dir",
            help=(
                "If Loci Extraction produced *_partition.nex files (codon-position partitions), "
                "point here to merge them into the output nexus."
            ),
        )

    with col_c2:
        st.markdown("**Concatenation options**")
        missing_char = st.selectbox(
            "Missing-taxon character",
            options=["-", "N", "?"],
            index=0,
            help="Character used to fill positions where a taxon is absent in a locus.",
        )
        concat_out_dir = st.text_input(
            "Output directory",
            value=str(Path(concat_dir).parent / "supermatrix") if concat_dir else "",
            key="concat_out_dir",
        )
        supermatrix_stem = st.text_input(
            "Output file stem",
            value="supermatrix",
            key="supermatrix_stem",
        )

    if selected_concat:
        st.markdown("**Locus preview**")
        rows = []
        for p in selected_concat:
            fp = Path(p)
            n_seq, aln_len = 0, 0
            if fp.exists():
                seqs = list(SeqIO.parse(str(fp), "fasta"))
                n_seq = len(seqs)
                if seqs:
                    aln_len = len(seqs[0].seq)
            rows.append(
                {
                    "Locus": fp.stem,
                    "Sequences": n_seq,
                    "Alignment length": aln_len,
                }
            )
        st.dataframe(pd.DataFrame(rows), width="stretch", hide_index=True)
        total_sites = sum(r["Alignment length"] for r in rows)
        st.metric("Total sites (before concat)", total_sites)

    if st.button("▶️ Concatenate", type="primary", key="btn_concat"):
        if not selected_concat:
            st.warning("Select at least two aligned FASTAs.")
        elif not concat_out_dir:
            st.warning("Specify an output directory.")
        else:
            out_path = Path(concat_out_dir)
            out_path.mkdir(parents=True, exist_ok=True)
            out_fasta = out_path / f"{supermatrix_stem}.fasta"
            out_nex = out_path / f"{supermatrix_stem}_partitions.nex"

            # Look for codon partition files
            codon_parts = None
            if codon_part_dir and Path(codon_part_dir).is_dir():
                codon_parts = []
                for fasta_path in selected_concat:
                    stem = Path(fasta_path).stem.replace("_trimmed", "").replace(
                        "_aligned", ""
                    )
                    nex_candidate = Path(codon_part_dir) / f"{stem}_partition.nex"
                    codon_parts.append(str(nex_candidate))

            with st.spinner("Concatenating…"):
                stats = concatenate_alignments(
                    aligned_fastas=selected_concat,
                    output_fasta=out_fasta,
                    partition_file=out_nex,
                    codon_partition_files=codon_parts,
                    missing_char=missing_char,
                )

            st.success("Supermatrix created.")
            col_m1, col_m2, col_m3 = st.columns(3)
            col_m1.metric("Loci", stats["n_loci"])
            col_m2.metric("Taxa", stats["n_taxa"])
            col_m3.metric("Total sites", stats["total_sites"])

            if stats.get("missing_taxa"):
                with st.expander("⚠️ Missing taxa (gap-padded)"):
                    for locus, taxa in stats["missing_taxa"].items():
                        st.write(f"**{locus}**: {', '.join(taxa)}")

            # Per-locus sites
            st.markdown("**Sites per locus**")
            st.dataframe(
                pd.DataFrame(
                    [
                        {"Locus": k, "Sites": v}
                        for k, v in stats["per_locus_sites"].items()
                    ]
                ),
                width="stretch",
                hide_index=True,
            )

            # Downloads
            col_dl1, col_dl2 = st.columns(2)
            with col_dl1:
                st.download_button(
                    "⬇️ Download supermatrix FASTA",
                    data=out_fasta.read_bytes(),
                    file_name=out_fasta.name,
                    mime="application/octet-stream",
                    key="dl_supermatrix",
                )
            with col_dl2:
                st.download_button(
                    "⬇️ Download partition nexus",
                    data=out_nex.read_bytes(),
                    file_name=out_nex.name,
                    mime="text/plain",
                    key="dl_partitions",
                )

            # Show partition file content
            with st.expander("📋 Partition file"):
                st.code(out_nex.read_text(), language=None)

# ══════════════════════════════════════════════════════════════════════════════
# TAB 4 — Logs
# ══════════════════════════════════════════════════════════════════════════════
with tab_logs:
    st.subheader("Alignment run history")
    history_tsv = project_dir / "metadata" / "command_history.tsv"

    if not history_tsv.exists():
        st.info("No run history yet. Run an alignment to populate this.")
    else:
        import csv

        rows = []
        with open(history_tsv, newline="") as f:
            reader = csv.DictReader(f, delimiter="\t")
            for row in reader:
                if row.get("module") in ("mafft", "trimal", "macse", "concat"):
                    rows.append(row)

        if not rows:
            st.info("No alignment runs in history yet.")
        else:
            df = pd.DataFrame(rows)
            df_display = df[
                ["started_at", "module", "action", "returncode", "command"]
            ].copy()
            df_display["returncode"] = df_display["returncode"].astype(str)
            st.dataframe(df_display, width="stretch", hide_index=True)

            st.markdown("---")
            run_ids = df["run_id"].tolist()
            selected_run = st.selectbox("Inspect run", [""] + run_ids, key="align_run_sel")
            if selected_run:
                run_dir = project_dir / "runs" / selected_run
                for fname, label in [
                    ("terminal.log", "📄 terminal.log"),
                    ("command.json", "📋 command.json"),
                    ("environment.json", "🔬 environment.json"),
                ]:
                    fpath = run_dir / fname
                    if fpath.exists():
                        with st.expander(label):
                            if fname.endswith(".json"):
                                st.json(load_json(fpath, {}))
                            else:
                                st.code(fpath.read_text(), language=None)
