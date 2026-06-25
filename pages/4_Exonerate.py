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
import tempfile
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
from phylofetch.exonerate_utils import (
    extract_locus_exonerate,
    resolve_guide_path,
    scan_flagged_cds,
)
from phylofetch.ncbi_utils import (
    LOCUS_CATALOGUE,
    count_refs,
    load_ref_records,
    locus_ref_fasta,
    project_ref_dir,
)
from phylofetch.project_manager import (
    DEFAULT_PROJECT_DIR,
    RunManager,
    effective_taxon,
    load_assembly_registry,
    load_project_manifest,
    project_output_dir,
    update_step,
)
from phylofetch.protein_guide_utils import (
    LENGTH_TOLERANCE,
    filter_records_by_length,
    get_guides,
    guide_loci,
    write_guide_fasta,
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
loci_dir = project_output_dir(project_dir) / "loci"
per_strain_dir = loci_dir / "per_strain"
combined_dir = loci_dir / "combined"
goi_dir = loci_dir / "_goi_refs"

st.caption(f"📁 Project: `{Path(project_dir).name}`  ·  outputs: "
           f"`{str(loci_dir).replace(str(Path.home()), '~')}`")

if not registry:
    st.info("No assemblies registered yet — add them in the Assembly Manager.")
    st.stop()

CODING = [k for k in LOCUS_CATALOGUE if k not in ("ITS", "LSU", "SSU")]


def _project_protein_refs(loc: str) -> list:
    """Project-library records for `loc`, but only when that library is PROTEIN — so they can
    join the protein2genome guide set. A nucleotide library returns [] (mixing models is wrong);
    that case is handled by the 'Project reference library' mode instead."""
    fa = locus_ref_fasta(loc, ref_dir=ref_dir)
    if os.path.exists(fa) and os.path.getsize(fa) > 0 and detect_fasta_type(fa) == "protein":
        return load_ref_records(loc, ref_dir=ref_dir)
    return []

# ── 1 · Assemblies ────────────────────────────────────────────────────────────
st.subheader("1 · Assemblies")


def _fmt(sid: str) -> str:
    t = effective_taxon(registry[sid], default_taxon)
    return f"{sid}  ·  {t}" if t else sid


sel = st.multiselect("Assemblies", list(registry), default=list(registry), format_func=_fmt)

# ── 2 · Coding loci ───────────────────────────────────────────────────────────
st.subheader("2 · Coding loci")
ref_source = st.radio(
    "Reference source",
    ["Bundled protein guides (recommended — no fetching)",
     "Bundled guides + project library (taxon-closer)",
     "Project reference library"],
    help="Bundled: universal cross-fungal protein guides shipped with phylofetch (D-020) → "
         "Exonerate protein2genome, works kingdom-wide, no NCBI fetching. "
         "Bundled + project library: the universal guides PLUS the project's fetched "
         "target-taxon protein refs in one query — Exonerate keeps whichever scores best "
         "(D-023). Project library: only the protein/nucleotide refs you fetched on the NCBI "
         "References page.",
)
ref_mode = ("both" if ref_source.startswith("Bundled guides +")
            else "bundled" if ref_source.startswith("Bundled")
            else "library")
keep_flagged = False  # D-025 length-filter opt-out (set in the bundled/both branch)

if ref_mode in ("bundled", "both"):
    avail = [loc for loc in CODING if loc in set(guide_loci())]

    def _fmt_guide(loc: str) -> str:
        base = f"{loc} · {len(get_guides(loc))} guide(s)"
        if ref_mode == "both":
            kept, dropped = filter_records_by_length(_project_protein_refs(loc), loc)
            tail = f" + {len(kept)} proj protein" if kept else ""
            tail += f" · ⚠️ {len(dropped)} length-flagged" if dropped else ""
            return base + tail
        return base

    sel_loci = st.multiselect("Loci (bundled guides)", avail, default=avail,
                              format_func=_fmt_guide)
    if ref_mode == "both":
        st.caption("ℹ️ Only **protein** project references are layered onto the bundled guides "
                   "(protein2genome); a locus whose library is nucleotide just uses the bundled "
                   "guides here. Fetch taxon-closer proteins on the **NCBI References** page.")
        st.caption(f"🛡️ **Length sanity filter (D-025):** project refs whose length deviates "
                   f">{LENGTH_TOLERANCE:.0%} from the curated bundled-guide length for the locus "
                   "are dropped from the Exonerate query (they're usually mis-annotated / fused "
                   "GenBank models that out-score the correct guides on raw identity). The "
                   "bundled guides always remain, so the locus is never left without a query.")
        keep_flagged = st.checkbox("Keep length-flagged project refs anyway (not recommended)",
                                   value=False,
                                   help="Override the D-025 filter and feed every fetched protein "
                                        "ref to Exonerate regardless of length.")
else:
    avail = [loc for loc in CODING if count_refs(loc, ref_dir=ref_dir) > 0]
    if not avail:
        st.caption("No coding-locus references in this project yet — fetch them on the "
                   "**NCBI References** page, or switch to bundled guides above.")
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
    ex_effort = st.selectbox(
        "Refinement effort", ["Balanced — escalate to region (recommended)",
                              "Fast — single pass, no refinement",
                              "Thorough — escalate to full (slow)"],
        help="How hard Exonerate works to fix a frameshifted / internal-stop CDS by re-optimizing "
             "splice boundaries (D-035). **Balanced**: none → `--refine region` (fixes the common "
             "splice misplacement, fast). **Fast**: one pass, no rescue. **Thorough**: also tries "
             "`--refine full` — exhaustive DP over the *whole* narrowed contig, minutes per strain "
             "on a long multi-intron gene like RPB2, and it rarely beats region. Prefer Balanced, "
             "then run full on just the flagged sequences via the deep-refine pass at the bottom.")
    if ex_effort.startswith("Fast"):
        ex_refine, ex_escalate, ex_ceiling = "none", False, "none"
    elif ex_effort.startswith("Thorough"):
        ex_refine, ex_escalate, ex_ceiling = "none", True, "full"
    else:
        ex_refine, ex_escalate, ex_ceiling = "none", True, "region"
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
        refined = (f" · 🔧 refined:{result['refine_used']}" if result.get("refine_escalated") else "")
        st.write(f"✅ {label} (Exonerate {model or 'protein2genome'}): {result['n_exons']} exon(s), "
                 f"{result['cds_length']} bp CDS, {n_introns} intron(s){splice_note} · {qc}{refined} · "
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
                  "escalate_ceiling": ex_ceiling, "geneticcode": int(ex_geneticcode),
                  "evalue": evalue_float}
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
                refine=ex_refine, escalate_refine=ex_escalate, escalate_ceiling=ex_ceiling,
                geneticcode=int(ex_geneticcode), min_pident=float(min_pident),
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
            tblastn_bin=tblastn_bin, run_dir=str(rr[1]), require_complete_cds=require_cds,
            manager=manager)

    flagged_notice: dict = {}   # locus -> [(ref_id, length, expected, flag)] dropped by D-025
    flagged_shown: set = set()

    def _guide_fa_for(locus):
        """Build the Exonerate query FASTA for `locus`, applying the D-025 bundled-length sanity
        filter to fetched project refs (unless the user opted out). Dropped refs are recorded in
        flagged_notice for a post-build notice. The bundled guides always remain, so 'both' mode
        is never left without a query even if every project ref is dropped."""
        from Bio import SeqIO as _SeqIO
        guide_path = str(Path(project_dir) / "scratch" / "guides" / f"{locus}_guide.fasta")
        if ref_mode == "library":
            recs = _project_protein_refs(locus)        # protein-only; [] when library is nt
            if not recs:                               # nucleotide library → length n/a, as-is
                return locus_ref_fasta(locus, ref_dir=ref_dir)
            kept = recs
            if not keep_flagged:
                kept, dropped = filter_records_by_length(recs, locus)
                if dropped:
                    flagged_notice[locus] = [(r.id, n, e, f) for r, n, e, f in dropped]
            if not kept:
                return None
            Path(guide_path).parent.mkdir(parents=True, exist_ok=True)
            _SeqIO.write(kept, guide_path, "fasta")    # library-only: no bundled guides added
            return guide_path
        extra = None
        if ref_mode == "both":
            recs = _project_protein_refs(locus)
            if keep_flagged:
                extra = recs or None
            else:
                kept, dropped = filter_records_by_length(recs, locus)
                if dropped:
                    flagged_notice[locus] = [(r.id, n, e, f) for r, n, e, f in dropped]
                extra = kept or None
        return write_guide_fasta(locus, guide_path, extra_records=extra)

    qc_outcomes: list = []   # (sample, locus, verdict, detail) — for a loud, persistent QC summary

    def _qc_record(sample, locus, result, status):
        if result is None:
            verdict = "DROPPED" if "strict QC" in status else "FAILED"
        elif result.get("n_internal_stops", 0) or result.get("len_mod3", 0):
            verdict = "REVIEW"
        else:
            verdict = "PASS"
        qc_outcomes.append((sample, locus, verdict, status))

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
            ref_fa = _guide_fa_for(locus)
            if locus in flagged_notice and locus not in flagged_shown:
                flagged_shown.add(locus)
                drops = "; ".join(f"`{rid}` ({n} aa vs ~{e}, {flag})"
                                  for rid, n, e, flag in flagged_notice[locus])
                st.caption(f"🛡️ {locus}: dropped length-flagged project ref(s) — {drops}")
            if not ref_fa or not os.path.exists(ref_fa):
                _report(locus, None, "no reference available for this locus")
                done += 1; prog.progress(done / total)
                continue
            with st.spinner(f"[{strain_id}] {locus}…"):
                result, status = _extract_coding(strain_id, assembly, locus, ref_fa,
                                                 str(strain_out / locus))
            _report(locus, result, status)
            _qc_record(strain_id, locus, result, status)
            done += 1; prog.progress(done / total)
        for gname, gref in goi_genes.items():
            with st.spinner(f"[{strain_id}] {gname} (GOI)…"):
                result, status = _extract_coding(strain_id, assembly, gname, gref,
                                                 str(strain_out / gname))
            _report(gname, result, status, goi=True)
            _qc_record(strain_id, gname, result, status)
            done += 1; prog.progress(done / total)

    combined = {}
    for locus in targets:
        merged = merge_per_strain_outputs(str(per_strain_dir), str(combined_dir), locus)
        if merged:
            combined[locus] = merged
    n_dropped = sum(len(v) for v in flagged_notice.values())
    update_step(project_dir, "coding", status="done",
                outputs={loc: list(m) for loc, m in combined.items()},
                notes=f"strains={len(sel)}; loci={','.join(sel_loci) or '—'}; "
                      f"goi={','.join(goi_genes) or '—'}; "
                      f"ref_source={ref_mode}; "
                      f"engine={'exonerate' if exonerate_ok else 'blast-fallback'}"
                      + ("" if not relaxed_blast else "; mode=relaxed_blast")
                      + (f"; length_filter={'off' if keep_flagged else 'on'}")
                      + (f"; dropped_refs={n_dropped}" if n_dropped else ""))
    st.session_state["exonerate_combined"] = {loc: {k: v for k, v in m.items()}
                                              for loc, m in combined.items()}
    st.session_state["exonerate_qc"] = qc_outcomes
    st.success("Coding extraction complete — provenance recorded in the manifest.")
    st.rerun()

# ── Results ───────────────────────────────────────────────────────────────────
qc = st.session_state.get("exonerate_qc")
if qc:
    from collections import Counter
    vc = Counter(v for _, _, v, _ in qc)
    st.markdown("---")
    st.markdown(f"**QC:** {vc.get('PASS', 0)} clean · {vc.get('REVIEW', 0)} review · "
                f"{vc.get('DROPPED', 0)} dropped · {vc.get('FAILED', 0)} failed")
    attn = [r for r in qc if r[2] in ("REVIEW", "DROPPED", "FAILED")]
    if attn:
        # Loud + persistent: flagged / dropped sequences never just vanish from the combined files.
        with st.expander(f"⚠️ {len(attn)} sequence(s) need attention — not in the clean set",
                         expanded=True):
            st.dataframe(pd.DataFrame([{"Sample": s, "Locus": l, "Verdict": v, "Detail": d}
                                       for s, l, v, d in attn]),
                         width="stretch", hide_index=True)
            if vc.get("DROPPED"):
                st.warning("**Dropped** rows were excluded from the combined FASTAs by **Strict QC** "
                           "(frameshift / internal stop). Uncheck Strict QC to keep them written-and-"
                           "flagged, or check the assembly. Note: a splice-boundary frameshift is "
                           "often recovered automatically by boundary-refinement escalation (D-030) — "
                           "if these persist, the genomic DNA may still be clean (verify the gene).")

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

# ── Deep refinement (optional, slow; re-runnable any time) ─────────────────────
# A targeted `--refine full` rescue pass on just the flagged CDS, skipped by the default
# Balanced/Fast effort. Stateless (scans disk), so it works in a later session too — e.g. if a
# tree/alignment looks off, come back and deep-refine the still-flagged sequences (D-035).
st.markdown("---")
st.subheader("🔧 Deep refinement (optional, slow)")
st.caption("Run `--refine full` on just the **flagged** coding CDS (frameshift / internal stop) "
           "from a previous extraction — the exhaustive boundary pass the default effort skips. "
           "Safe to run any time (including a later session); it reuses the guides saved during "
           "extraction and **keeps the new CDS only if it is strictly cleaner** (D-035).")


def _guide_for(locus: str):
    """Guide for a deep-refine re-run; a gene of interest (its own `_goi_refs/` file exists) takes
    precedence over a same-named catalogue guide (D-036 review)."""
    return resolve_guide_path(locus, project_dir, goi_dir, ref_dir,
                              is_goi=(goi_dir / f"{locus}.fasta").exists())


flagged = scan_flagged_cds(str(per_strain_dir), geneticcode=int(ex_geneticcode))
if not flagged:
    st.caption("✅ No flagged coding CDS in this project's outputs — nothing to deep-refine.")
elif not exonerate_ok:
    st.error("⚠️ exonerate not found — deep refinement requires it.")
else:
    st.warning("⚠️ Deep refinement **re-extracts using the current Settings above** (intron sizes, "
               "min identity, narrowing) plus the saved guide — not necessarily the exact "
               "parameters that produced the original. It keeps the new CDS **only if strictly "
               "cleaner**; otherwise the original is left untouched.")
    st.dataframe(pd.DataFrame([{
        "Strain": r["strain"], "Locus": r["locus"], "CDS bp": r["cds_length"],
        "len%3": r["len_mod3"], "Internal stops": r["n_internal_stops"],
        "Refine already run": r["refine_used"] or "—",
        "Guide on disk": "✅" if _guide_for(r["locus"]) else "❌ re-extract first",
    } for r in flagged]), width="stretch", hide_index=True)

    rescuable = [r for r in flagged if _guide_for(r["locus"])]
    if len(rescuable) < len(flagged):
        miss = sorted({r["locus"] for r in flagged if not _guide_for(r["locus"])})
        st.warning(f"No saved guide for: {', '.join(miss)} — re-run extraction for these first.")

    if rescuable and st.button(f"🔧 Run --refine full on {len(rescuable)} flagged sequence(s)",
                               type="primary"):
        manager = RunManager(project_dir)
        prog = st.progress(0.0)
        loci_touched, n_fixed, n_kept = set(), 0, 0
        for i, r in enumerate(rescuable):
            strain, locus = r["strain"], r["locus"]
            asm = registry.get(strain, {}).get("assembly_path", "")
            if not asm or not os.path.exists(asm):
                st.write(f"⚠️ {strain}/{locus}: assembly file not found — skipped")
                prog.progress((i + 1) / len(rescuable))
                continue
            locus_dir = Path(r["locus_dir"])
            before = (r["n_internal_stops"], r["len_mod3"])
            before_key = (before[1] == 0, -before[0])          # (frame_ok, -stops): higher=cleaner
            # Snapshot the original so a worse / failed re-run (e.g. drifted settings) can't
            # clobber a better result — keep the new CDS only if strictly cleaner (D-035 review).
            backup_root = Path(tempfile.mkdtemp(prefix="deeprefine_"))
            backup = backup_root / "bak"
            shutil.copytree(locus_dir, backup)
            try:
                with st.spinner(f"Deep-refining {strain}/{locus} (--refine full)…"):
                    result, status = extract_locus_exonerate(
                        assembly_fasta=asm, reference_fasta=_guide_for(locus),
                        output_dir=str(locus_dir), strain_id=strain, locus_name=locus,
                        exonerate_bin=exonerate_bin, blastn_bin=blastn_bin, tblastn_bin=tblastn_bin,
                        narrow=ex_narrow, minintron=int(ex_minintron), maxintron=int(ex_maxintron),
                        bestn=int(ex_bestn), refine="none", escalate_refine=True,
                        escalate_ceiling="full", geneticcode=int(r.get("geneticcode", ex_geneticcode)),
                        min_pident=float(min_pident), evalue=evalue_float, threads=int(threads),
                        strict_qc=False, manager=manager)
            except Exception as exc:  # noqa: BLE001 — never let one strain abort the batch
                result, status = None, f"error: {exc}"

            keep = bool(result) and (result.get("len_mod3", 1) == 0,
                                     -result.get("n_internal_stops", 99)) > before_key
            if keep:
                after = (result.get("n_internal_stops", 0), result.get("len_mod3", 0))
                loci_touched.add(locus)
                if after == (0, 0):
                    n_fixed += 1
                st.write(f"✅ improved {strain}/{locus}: stops {before[0]}→{after[0]}, "
                         f"len%3 {before[1]}→{after[1]} · refine={result.get('refine_used', '?')}")
                shutil.rmtree(backup_root, ignore_errors=True)
            else:
                # Not cleaner (or failed) → restore the original, untouched.
                shutil.rmtree(locus_dir, ignore_errors=True)
                shutil.move(str(backup), str(locus_dir))
                shutil.rmtree(backup_root, ignore_errors=True)
                n_kept += 1
                detail = status if not result else f"no improvement (stops {before[0]}, len%3 {before[1]})"
                st.write(f"↩️ {strain}/{locus}: {detail} — kept original")
            prog.progress((i + 1) / len(rescuable))

        for locus in sorted(loci_touched):
            merge_per_strain_outputs(str(per_strain_dir), str(combined_dir), locus)
        if loci_touched:
            update_step(project_dir, "coding", status="done",
                        notes=(f"deep_refine (--refine full on flagged): {n_fixed} now clean, "
                               f"{len(loci_touched)} locus dir(s) updated, {n_kept} kept as-is"))
        st.success(f"Deep refinement complete — {n_fixed}/{len(rescuable)} now frame-clean, "
                   f"{len(loci_touched)} updated, {n_kept} kept unchanged (not cleaner).")
        st.rerun()
