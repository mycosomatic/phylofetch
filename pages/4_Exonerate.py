"""
pages/2_Exonerate.py
--------------------
Exonerate (coding loci) — component page (D-012 / RM-007 step 4d).

Frame-safe protein-coding locus extraction (D-008): for each assembly × locus, BLAST narrows
to the best contig, then Exonerate spliced alignment (protein2genome / coding2genome) recovers
a translatable, frame-checked CDS. Two reference sources: the project's per-locus library
(NCBI References page) and an ad-hoc **gene of interest** (paste/upload any ortholog). When
`exonerate` is absent it falls back to the BLAST HSP-as-exon path with a frame-safety warning
(D-008). Outputs are per-project (<project>/results/loci, D-015); provenance to the manifest.

Standalone + chainable; runs alongside the old 2_Loci_Extraction page until the monolith retires.
"""

import os
import shutil
import sys
from pathlib import Path

import pandas as pd
import streamlit as st

sys.path.insert(0, str(Path(__file__).parent.parent / "src"))
from phylofetch.blast_loci_utils import (
    detect_fasta_type,
    extract_locus,
    merge_per_strain_outputs,
)
from phylofetch.config import load_config
from phylofetch.exonerate_utils import extract_locus_exonerate
from phylofetch.ncbi_utils import (
    LOCUS_CATALOGUE,
    count_refs,
    locus_ref_fasta,
    project_ref_dir,
)
from phylofetch.project_manager import (
    DEFAULT_PROJECT_DIR,
    RunManager,
    effective_taxon,
    load_assembly_registry,
    load_project_manifest,
    update_step,
)

st.set_page_config(page_title="Exonerate (coding)", page_icon="🧬", layout="wide")
st.title("🧬 Exonerate — coding loci")
st.caption("Frame-safe CDS extraction: BLAST narrows to the best contig, then Exonerate refines "
           "exon/intron boundaries into a translatable, frame-checked CDS. A relaxed BLAST "
           "amplicon mode (genomic, no frame guarantee) is also available. Use the project's "
           "reference library or paste an ad-hoc gene of interest.")

cfg = load_config()
project_dir = cfg.get("project_dir", str(DEFAULT_PROJECT_DIR))
tp = cfg.get("tool_paths", {})
registry = load_assembly_registry(project_dir)
manifest = load_project_manifest(project_dir)
default_taxon = manifest.get("default_taxon", "")

ref_dir = project_ref_dir(project_dir)
loci_dir = Path(project_dir) / "results" / "loci"
per_strain_dir = loci_dir / "per_strain"
combined_dir = loci_dir / "combined"
goi_dir = loci_dir / "_goi_refs"

st.caption(f"📁 Project: `{Path(project_dir).name}`  ·  outputs: "
           f"`{str(loci_dir).replace(str(Path.home()), '~')}`")

if not registry:
    st.info("No assemblies registered yet — add them in the Assembly Manager.")
    st.stop()

CODING = [k for k in LOCUS_CATALOGUE if k not in ("ITS", "LSU", "SSU")]

# ── 1 · Assemblies ────────────────────────────────────────────────────────────
st.subheader("1 · Assemblies")


def _fmt(sid: str) -> str:
    t = effective_taxon(registry[sid], default_taxon)
    return f"{sid}  ·  {t}" if t else sid


sel = st.multiselect("Assemblies", list(registry), default=list(registry), format_func=_fmt)

# ── 2 · Coding loci (project library) ─────────────────────────────────────────
st.subheader("2 · Coding loci (from this project's reference library)")
avail = [loc for loc in CODING if count_refs(loc, ref_dir=ref_dir) > 0]
if not avail:
    st.caption("No coding-locus references in this project yet — fetch them on the "
               "**NCBI References** page first.")
sel_loci = st.multiselect(
    "Loci", avail, default=avail,
    format_func=lambda loc: f"{loc} · {count_refs(loc, ref_dir=ref_dir)} refs",
)

# ── 3 · Gene of interest (ad-hoc ortholog) ────────────────────────────────────
st.subheader("3 · Gene of interest (optional)")
if "goi_genes" not in st.session_state:
    st.session_state.goi_genes = {}
with st.expander("🔎 Extract any ortholog's CDS — paste/upload a reference protein or CDS"):
    st.caption("Exonerate locates the ortholog in each selected assembly and returns just the "
               "CDS + exon model — no NCBI reference library needed.")
    g1, g2 = st.columns([1, 2])
    with g1:
        goi_name = st.text_input("Gene name (for outputs / tip labels)", key="goi_name",
                                 placeholder="e.g. mcm7, betaTUB")
        goi_file = st.file_uploader("…or upload FASTA",
                                    type=["fa", "fasta", "faa", "fna", "txt"], key="goi_file")
    with g2:
        goi_text = st.text_area("Reference FASTA (protein or nucleotide CDS)", key="goi_text",
                                height=110, placeholder=">ref\nMSEQ… (protein) or ATG… (CDS)")
    if st.button("➕ Add gene of interest"):
        name = (goi_name or "").strip().replace(" ", "_")
        content = ""
        if goi_file is not None:
            content = goi_file.getvalue().decode("utf-8", "replace")
        elif goi_text.strip():
            content = goi_text
        if not name:
            st.warning("Enter a gene name.")
        elif not content.strip():
            st.warning("Paste or upload a reference FASTA.")
        else:
            if not content.lstrip().startswith(">"):
                content = f">{name}_ref\n{content.strip()}\n"
            goi_dir.mkdir(parents=True, exist_ok=True)
            ref_path = goi_dir / f"{name}.fasta"
            ref_path.write_text(content)
            st.session_state.goi_genes[name] = str(ref_path)
            st.success(f"Added gene of interest '{name}'.")
            st.rerun()
    for gname, gpath in list(st.session_state.goi_genes.items()):
        rc1, rc2 = st.columns([5, 1])
        exists = os.path.exists(gpath)
        rtype = detect_fasta_type(gpath) if exists else "?"
        rc1.caption(f"**{gname}** · {rtype} ref · `{gpath}`"
                    + ("" if exists else " · ⚠️ file missing"))
        if rc2.button("🗑️ Remove", key=f"goi_del_{gname}"):
            del st.session_state.goi_genes[gname]
            st.rerun()
goi_genes = {n: p for n, p in st.session_state.goi_genes.items() if os.path.exists(p)}

# ── 4 · Settings ──────────────────────────────────────────────────────────────
st.subheader("4 · Settings")
coding_mode = st.radio(
    "Coding extraction mode",
    ["Exonerate (frame-safe)", "BLAST amplicon (relaxed, genomic)"],
    horizontal=True,
    help="Exonerate: spliced alignment → a translatable, frame-checked CDS (recommended). "
         "BLAST amplicon (relaxed): HSP-as-exon genomic amplicon with NO reading-frame "
         "guarantee — the former 'PCR amplicon refs (relaxed)' strategy; use with "
         "amplicon-style references where introns/partials are expected (D-008).",
)
relaxed_blast = coding_mode.startswith("BLAST amplicon")
sa, sb, sc = st.columns(3)
with sa:
    threads = st.slider("Threads", 1, 32, 4)
    exonerate_bin = st.text_input("exonerate", value=tp.get("exonerate", "exonerate"))
    blastn_bin = st.text_input("blastn", value=tp.get("blastn", "blastn"))
    tblastn_bin = st.text_input("tblastn", value=tp.get("tblastn", "tblastn"))
with sb:
    min_pident = st.number_input("Min % identity", value=70, min_value=30, max_value=100,
                                 help="BLAST identity floor (also gates Exonerate contig-narrowing).")
    evalue = st.text_input("E-value", "1e-20")
    ex_maxintron = st.number_input("Max intron (bp)", value=2000, min_value=50, max_value=200000,
                                   step=100, help="Fungal introns are short; 2000 is a safe ceiling.")
    ex_minintron = st.number_input("Min intron (bp)", value=20, min_value=4, max_value=200)
with sc:
    ex_bestn = st.number_input("Report top N models", value=1, min_value=1, max_value=20,
                               help=">1 surfaces paralogs / tandem duplicates.")
    ex_geneticcode = st.number_input("Genetic code (NCBI)", value=1, min_value=1, max_value=33)
    ex_refine = st.selectbox("Boundary refinement", ["none", "region", "full"])
    ex_narrow = st.checkbox("BLAST-narrow to best contig (faster)", value=True)
    ex_strict_qc = st.checkbox("Strict QC (reject frameshift / internal stop)", value=False)

exonerate_ok = shutil.which(exonerate_bin) is not None
tcheck = " · ".join(f"{'✅' if shutil.which(b) else '❌'} {n}"
                    for n, b in [("exonerate", exonerate_bin), ("blastn", blastn_bin),
                                 ("tblastn", tblastn_bin)])
st.caption("Tool check: " + tcheck)
min_cds_pct = 50
if relaxed_blast:
    st.info("🧬 **Relaxed BLAST amplicon mode** — genomic amplicon via HSP-as-exon with "
            "`require_complete_cds=False`; **no reading-frame guarantee** (D-008). For "
            "amplicon-style references; Exonerate is skipped even if installed.")
elif not exonerate_ok:
    st.warning("⚠️ **exonerate not found** — coding loci will fall back to the BLAST HSP-as-exon "
               "path, which does **not** validate the reading frame (a 1–2 bp boundary error can "
               "frameshift the CDS silently). Install: `conda install -c bioconda exonerate`.")
    min_cds_pct = st.number_input("Min CDS % of reference (BLAST fallback gate)", value=50,
                                  min_value=10, max_value=100)

try:
    evalue_float = float(evalue)
except ValueError:
    evalue_float = 1e-20
    st.error(f"Invalid e-value '{evalue}', using 1e-20.")

# ── 5 · Run ───────────────────────────────────────────────────────────────────
st.subheader("5 · Run")
targets = list(sel_loci) + list(goi_genes)
if not sel:
    st.info("Select at least one assembly.")
elif not targets:
    st.info("Select a coding locus or add a gene of interest.")

blastn_ok = shutil.which(blastn_bin) is not None
run_disabled = not (sel and targets and (exonerate_ok or blastn_ok))


def _report(locus, result, status, goi=False):
    label = f"🔎 {locus}" if goi else locus
    if not result:
        st.write(f"⚠️ {label}: {status}")
        return
    n_introns = result["n_introns"]
    splices_ok = sum(1 for intr in result["introns"]
                     if intr["splice_5"] == "GT" and intr["splice_3"] == "AG")
    splice_note = f" · GT-AG {splices_ok}/{n_introns}" if n_introns > 0 else ""
    if result.get("tool") == "exonerate":
        stops, frame = result.get("n_internal_stops", 0), result.get("len_mod3", 0)
        qc = ("✅ frame OK" if (stops == 0 and frame == 0)
              else f"⚠️ QC review: {stops} stop(s), len%3={frame}")
        model = result.get("blast_type", "").replace("exonerate:", "").split(":")[0]
        paralog = (f" · {result['n_other_models']} paralog cand(s)"
                   if result.get("n_other_models") else "")
        st.write(f"✅ {label} (Exonerate {model or 'protein2genome'}): {result['n_exons']} exon(s), "
                 f"{result['cds_length']} bp CDS, {n_introns} intron(s){splice_note} · {qc} · "
                 f"id {result.get('query_pident', 0):.0f}% cov {result.get('query_coverage', 0):.0f}%"
                 f"{paralog} · ref `{result.get('ref_accession', '?')}`")
    else:
        st.write(f"✅ {label} ({result['blast_type']}): {result['n_exons']} exon(s), "
                 f"{result['cds_length']} bp CDS, {n_introns} intron(s){splice_note} · "
                 f"ref `{result.get('ref_accession', '?')}`")


if st.button("🚀 Run extraction", type="primary", disabled=run_disabled):
    per_strain_dir.mkdir(parents=True, exist_ok=True)
    combined_dir.mkdir(parents=True, exist_ok=True)
    manager = RunManager(project_dir)
    prog = st.progress(0.0)
    jobs = [(s, loc, False) for s in sel for loc in sel_loci] + \
           [(s, g, True) for s in sel for g in goi_genes]
    total = max(len(jobs), 1)

    def _extract_coding(strain_id, assembly, locus, ref_fa, locus_out):
        os.makedirs(locus_out, exist_ok=True)
        params = {"min_pident": min_pident, "maxintron": int(ex_maxintron),
                  "bestn": int(ex_bestn), "narrow": ex_narrow, "refine": ex_refine,
                  "geneticcode": int(ex_geneticcode), "evalue": evalue_float}
        # Exonerate frame-safe path (skipped in relaxed mode even if exonerate is installed).
        if exonerate_ok and not relaxed_blast:
            _rid, rdir = manager.dry_run(
                [exonerate_bin, "--model", "(auto)", "--query", ref_fa, "--target", assembly],
                module="exonerate", action=f"exonerate_{locus}_{strain_id}",
                inputs={"assembly": assembly, "reference": ref_fa},
                outputs={"locus_dir": locus_out}, params=params)
            return extract_locus_exonerate(
                assembly_fasta=assembly, reference_fasta=ref_fa, output_dir=locus_out,
                strain_id=strain_id, locus_name=locus, exonerate_bin=exonerate_bin,
                blastn_bin=blastn_bin, tblastn_bin=tblastn_bin, narrow=ex_narrow,
                minintron=int(ex_minintron), maxintron=int(ex_maxintron), bestn=int(ex_bestn),
                refine=ex_refine, geneticcode=int(ex_geneticcode), min_pident=float(min_pident),
                evalue=evalue_float, threads=int(threads), strict_qc=ex_strict_qc,
                manager=manager, run_dir=str(rdir))
        # BLAST HSP-as-exon path: relaxed genomic amplicon (require_complete_cds=False, the
        # former 'PCR amplicon refs (relaxed)' strategy) or the frame-safe-intended fallback
        # (require_complete_cds=True) when exonerate is absent. Neither guarantees frame (D-008).
        require_cds = not relaxed_blast
        tag = "blast_relaxed" if relaxed_blast else "blast_fallback"
        rr = manager.dry_run(
            [blastn_bin if detect_fasta_type(ref_fa) == "nucleotide" else tblastn_bin,
             "-query", ref_fa, "-subject", assembly, "-evalue", str(evalue_float)],
            module="exonerate", action=f"{tag}_{locus}_{strain_id}",
            inputs={"assembly": assembly, "reference": ref_fa},
            outputs={"locus_dir": locus_out}, params={**params, "require_complete_cds": require_cds})
        return extract_locus(
            assembly_fasta=assembly, reference_fasta=ref_fa, output_dir=locus_out,
            strain_id=strain_id, locus_name=locus, min_pident=float(min_pident),
            min_cds_pct_of_ref=float(min_cds_pct), evalue=evalue_float,
            blastn_task="dc-megablast", threads=int(threads), blastn_bin=blastn_bin,
            tblastn_bin=tblastn_bin, run_dir=str(rr[1]), require_complete_cds=require_cds)

    done = 0
    for strain_id in sel:
        assembly = registry[strain_id].get("assembly_path", "")
        if not assembly or not os.path.exists(assembly):
            st.error(f"{strain_id}: assembly file not found ({assembly})")
            done += len([j for j in jobs if j[0] == strain_id]); prog.progress(done / total)
            continue
        strain_out = per_strain_dir / strain_id
        st.markdown(f"**{strain_id}**" + (f" · {effective_taxon(registry[strain_id], default_taxon)}"
                                          if effective_taxon(registry[strain_id], default_taxon) else ""))
        for locus in sel_loci:
            ref_fa = locus_ref_fasta(locus, ref_dir=ref_dir)
            with st.spinner(f"[{strain_id}] {locus}…"):
                result, status = _extract_coding(strain_id, assembly, locus, ref_fa,
                                                 str(strain_out / locus))
            _report(locus, result, status)
            done += 1; prog.progress(done / total)
        for gname, gref in goi_genes.items():
            with st.spinner(f"[{strain_id}] {gname} (GOI)…"):
                result, status = _extract_coding(strain_id, assembly, gname, gref,
                                                 str(strain_out / gname))
            _report(gname, result, status, goi=True)
            done += 1; prog.progress(done / total)

    combined = {}
    for locus in targets:
        merged = merge_per_strain_outputs(str(per_strain_dir), str(combined_dir), locus)
        if merged:
            combined[locus] = merged
    update_step(project_dir, "coding", status="done",
                outputs={loc: list(m) for loc, m in combined.items()},
                notes=f"strains={len(sel)}; loci={','.join(sel_loci) or '—'}; "
                      f"goi={','.join(goi_genes) or '—'}; "
                      f"engine={'exonerate' if exonerate_ok else 'blast-fallback'}")
    st.session_state["exonerate_combined"] = {loc: {k: v for k, v in m.items()}
                                              for loc, m in combined.items()}
    st.success("Coding extraction complete — provenance recorded in the manifest.")
    st.rerun()

# ── Results ───────────────────────────────────────────────────────────────────
res = st.session_state.get("exonerate_combined")
if res:
    st.markdown("---")
    st.markdown("**Combined FASTAs**")
    rows = []
    for loc, kinds in res.items():
        for kind, path in kinds.items():
            try:
                from Bio import SeqIO
                n = len(list(SeqIO.parse(path, "fasta")))
            except Exception:
                n = "?"
            rows.append({"Locus": loc, "Product": kind, "Sequences": n,
                         "Path": path.replace(str(Path.home()), "~")})
    if rows:
        st.dataframe(pd.DataFrame(rows), width="stretch", hide_index=True)
    else:
        st.caption("No combined outputs produced.")
