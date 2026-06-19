# phylofetch — Changelog

> **Append-only.** Do not delete past entries. Newest at the top. This is the "what actually
> changed" record. Rationale lives in `DECISIONS.md`; roadmap in `PLANNING.md`.

## 2026-06-19

- **Exonerate frame-safe CDS + gene-of-interest extraction (D-008, addresses RM-002).**
  New module `src/phylofetch/exonerate_utils.py`: `run_exonerate` (RunManager-logged),
  `parse_exonerate_gff` (verbatim Exonerate 2.4.0 GFF + `--ryo %tcs` parsing, multi-result
  aware), `select_best_model` (paralog disambiguation), `validate_cds` (reading-frame /
  internal-stop QC), `build_result_from_model` (rebuilds the shared result-dict shape from
  GFF coords; plus/minus-strand correct), `write_exonerate_fastas` (CDS / protein / genomic
  / introns) and `extract_locus_exonerate` (hybrid: tblastn/blastn narrows to the best
  contig, then `protein2genome`/`coding2genome` refines; whole-assembly fallback). Reuses
  `write_gff3` / `write_codon_partition` from `blast_loci_utils`. Validated end-to-end
  against the installed exonerate binary on both strands and the nucleotide-CDS model.
- **Loci Extraction page (`pages/2_Loci_Extraction.py`).** The strict coding strategy is now
  **"Coding loci – Exonerate (frame-safe)"**: BLAST narrows, Exonerate refines, with a
  visible BLAST-HSP fallback + frame-safety warning when `exonerate` is absent. Added
  Exonerate settings (max/min intron, `--bestn` paralogs, refine, genetic code, strict QC,
  narrow toggle), an `exonerate` tool-check row, a **"Gene of interest"** input
  (paste/upload an ortholog → CDS), and Results-tab surfacing of the translated protein +
  a CDS frame/stop QC badge. Catalogue loci and genes-of-interest share one logged code path.
- **`blast_loci_utils.py`.** `extract_from_hsps` documented as the **fallback** for coding
  CDS (Exonerate preferred; D-008), still primary for the relaxed PCR-amplicon strategy.
  `merge_per_strain_outputs` now also merges the translated `*_protein.fasta`.
- **`environment.yml`.** Added `bioconda::exonerate=2.4.0` (not a pip package; no
  `pyproject.toml` runtime-dependency change).
- **Tests.** `tests/test_exonerate_utils.py` — 23 tests (offline parser/QC/build against
  committed verbatim fixtures in `tests/fixtures/`, plus an exonerate-binary-guarded
  end-to-end class). Full suite **139 passing** (was 116).

- **UI legibility — Recovery summary table (`pages/2_Loci_Extraction.py`).** The per-strain ×
  per-locus CDS-length matrix styled each cell with a pastel background
  (`#c8e6c9` found / `#ffcdd2` missing) but set **no text color**, so under Streamlit's dark
  theme the near-white default text rendered white-on-pale-green and the length numbers were
  barely readable. Added an explicit dark text color matching each background
  (`#1b5e20` on green, `#b71c1c` on red) plus `font-weight:600`, via a named `_recovery_cell`
  styler. No semantics changed (green=found / red=missing). Cosmetic only.

## 2026-06-18

- **Proactive runtime-bug sweep (all pages).** A read-only audit for the two
  Streamlit-only bug classes that had been surfacing during review found three more, now
  fixed: (1) **Arrow mixed-type columns** in the Project-Setup registry table
  (`pages/0_Project_Setup.py`, N50/Contigs/Size/GC) and the Assembly-Manager summary table
  (`pages/1_Assembly_Manager.py`, Total/Contigs/N50/GC) — numeric values mix with the "—"
  fallback in one column; the affected columns are now stringified before render. (2) A
  **widget-key/session-state collision** in BUSCO Phylogenomics
  (`pages/4_BUSCO_Phylogenomics.py`): the local-assembly `st.multiselect(key="local_assembly_sel")`
  was overwritten by `st.session_state["local_assembly_sel"] = {…}` (also a list-vs-dict type
  clash with tab 2's reader); the selected-records dict now lives under a separate
  `local_assembly_records` key, with both readers updated.
- **Bug fix — QUAST table Arrow crash (`pages/1_Assembly_Manager.py`).** The QUAST metrics
  table built a `Value` column holding mixed Python types (int N50, float GC%, string fields),
  which pyArrow could not serialise (`ArrowTypeError: Expected bytes, got a 'int' object`).
  Stringify the values (`str(...)`) so the column is uniformly typed; values are displayed as
  text anyway. (The page-2 header-metadata `Value` column is already uniform string — built
  from `token.split("=")` — so it was left as-is.)
- **Bug fix — Results tab crash (`pages/2_Loci_Extraction.py`).** The "Show BLAST alignment"
  button used the same key (`aln_<strain>_<locus>`) as the `session_state` slot the handler
  then wrote to, which Streamlit forbids (`StreamlitAPIException: cannot be modified after the
  widget … is instantiated`) — it crashed the Results tab whenever an extracted locus was
  browsed. Split into distinct keys: button `btn_aln_…`, stored text `alntext_…`. Pre-existing
  bug, unrelated to the RM-001 work.
- **Housekeeping — Streamlit deprecation sweep.** Replaced all 28 `use_container_width=True`
  (deprecated, removal slated post-2025-12-31) with `width="stretch"` across the 6 pages;
  only `st.dataframe` / `st.data_editor` / `st.plotly_chart` were affected (the new literal
  `width` API landed in Streamlit 1.43), so the dependency floor was raised `1.28 → 1.43`.
- **RM-001 increment 2 — type material surfaced in the Reference Library UI
  (`pages/2_Loci_Extraction.py`).** Search gained a **Type material** toggle
  (Prefer type / All / Type only) wired to `type_mode`; search results show a **Type**
  flag column; the locus library table now joins the sidecar to show
  **Organism / Strain-voucher / Type (kind)** per reference, with a notice for pre-sidecar
  sequences. Fetches pass the search `query` through for provenance. Non-breaking; full
  suite 116 passing. (Reviewable in-app; tip-merge is increment 3.)
- **RM-001 increment 1 — reference metadata layer (`ncbi_utils.py`).** Added a structured
  `RefRecord` + per-locus sidecar JSON (`<locus>/<locus>_refs.json`) alongside the existing
  refs FASTA, holding organism / strain / culture-collection / voucher / `type_material` /
  normalised `type_kind` / `is_type` / length / fetched_at / query. Fetching now uses
  GenBank (`rettype="gb"`) and parses the source-feature qualifiers the FASTA-only fetch
  discarded; `normalize_type_kind()` decodes the INSDC `/type_material` value (ex-type
  cultures treated as type-grade per D-007). Search gained a `type_mode`
  (`all` / `prefer` / `type_only`) using NCBI's `sequence_from_type[filter]`, flagging hits
  via one cheap extra esearch. Each fetch appends a provenance line
  (`<locus>_fetch_log.jsonl`). 13 new network-free tests (`tests/test_ncbi_utils.py`); full
  suite 116 passing. Non-breaking to the UI — no page changes yet (tip-merge + UI surfacing
  are later increments). Closes the fetch half of the review's provenance gap.
- **Strategy/implementation review.** Conducted an in-depth, skeptical review of the whole
  project against its stated goals (barcoding → delimitation; BUSCO phylogenomics; type-
  material comparison). Recorded a roadmap (RM-001..RM-005) and a risk register in
  `PLANNING.md`, and two decisions: D-005 (inference boundary — thin runners + coalescent +
  ANI in-app, Bayesian/delimitation ported out) and D-006 (type-material acquisition is the
  first build priority). No code changed in this step.
- Added the working-agreement / session protocol to `CLAUDE.md` (read first every session).
- Created living documents: `PLANNING.md`, `DECISIONS.md`, `CHANGELOG.md`. (See D-001.)
- Merged browser/cloud session work (commit `4a03968`) into `dev`: **citable, user-extensible
  PCR primer library** — packaged `data/primers.json` (14 fully-cited fungal pairs), user
  library at `~/.phylofetch/primers.json`, edit-distance matching hardening, binding-site
  disambiguation UI, and RunManager provenance logging. (See D-003.)
- **PRM-001 data fix:** corrected the corrupted built-in `ACT-512F` primer to its canonical
  20-nt sequence (Carbone & Kohn 1999), with regression test. (See D-004.)

## Earlier (from git history, pre-changelog)

- Fixed assembly stats display, added BUSCO integration, surfaced primer mode.
- Consolidated import scanners; fixed assembly schema inconsistencies.
- Fixed EGAP assembly detection; added bulk-select to import scanner.
- Redesigned BUSCO Phylogenomics for NCBI-first genome comparison workflow.
- Added PCR primer-based locus extraction strategy.
