# phylofetch — Changelog

> **Append-only.** Do not delete past entries. Newest at the top. This is the "what actually
> changed" record. Rationale lives in `DECISIONS.md`; roadmap in `PLANNING.md`.

## 2026-06-20

- **ITSx (rDNA) component page + per-project extraction outputs (D-015 / RM-007 step 4c).**
  New standalone `pages/2_ITSx_rDNA.py`: pick assemblies + rDNA regions (ITS/ITS1/ITS2/LSU/SSU)
  → run ITSx → per-strain region FASTAs + cross-strain `*_combined.fasta`, with provenance to
  `workflow.steps.rDNA`. **D-015:** extraction outputs are now **per-project** under
  `<project>/results/loci/{per_strain,combined}` (not the shared `output_base`), avoiding
  cross-project collisions; Alignment Prep reads any directory (`rglob`) so it stays compatible.
  New tested helpers `itsx_utils.place_rdna_regions` / `combine_rdna_regions` factor the
  region-placement + combine logic the monolith did inline. Tests: 2 new → **210 passing**;
  page render verified headless via AppTest (7 assemblies, region checkboxes, run button).
- **NCBI References component page (D-012/D-013 / RM-007 step 4b).** New standalone
  `pages/2_NCBI_References.py`: tick loci as **checkboxes** → **preview NCBI hit counts** for
  the project taxon → **fetch** top-N (type-preferred) into the **per-project** library
  (`<project>/references`), with a **review/cull** table. The search organism **defaults to the
  project taxon** from the manifest (editable; hint to broaden species→genus→family). Fetch
  provenance is recorded in the manifest (`workflow.steps.references`). New
  `ncbi_utils.ncbi_search_count` (cheap `esearch retmax=0`) backs the preview. Runs alongside
  the old `2_Loci_Extraction` page during the transition (sidebar shows both, by design).
  Tests: 2 new (mocked `ncbi_search_count`) → **208 passing**; page render verified headless
  via `AppTest` on the real project (taxon guard + preview/fetch controls appear correctly).
- **Per-project reference libraries — data layer (D-013 / RM-007 step 4a).** `ncbi_utils.py`:
  the reference functions (`locus_ref_fasta`, `list_loci`, `load_ref_records`,
  `accessions_in_library`, `count_refs`, `fetch_and_store`, `delete_from_library`) now take a
  `ref_dir` parameter defaulting to the global library, plus a new
  `project_ref_dir(project_dir)` → `<project>/references`. Backward-compatible (the global
  default is unchanged); per-project storage takes effect once the References component page
  passes the project-scoped dir (step 4b). Tests: 6 new (`TestPerProjectRefDir`, incl.
  two-project isolation + default-is-global guard). **206 passing** (was 200).
- **ITS-based provisional taxon ID (D-014 / RM-007 step 3).** New
  `src/phylofetch/taxon_id_utils.py` and an "Identify taxon by ITS (NCBI remote BLAST)" control
  in the Assembly Manager: ITSx pulls the ITS region, `blastn -remote -db nt` (optionally
  fungi-restricted) finds the closest organisms, and the user picks one → sets the assembly
  taxon with `taxon_source="its_blast"`.
  - **Organism resolution finding:** `sscinames` is resolved from a *local* BLAST taxdb even
    under `-remote`, so it returns `N/A`/`N/A;N/A` without taxdb installed. The organism is
    therefore derived from the subject title (`stitle`), preferring a real `sscinames` when
    present; ";"-joined `N/A` handled. Avoids forcing a large taxdb download.
  - Remote call routed through RunManager (provenance); pure command-build / parse / rank split
    out for testing. ITSx absence and BLAST errors surface cleanly in the UI.
  - Tests: 15 network-free (`tests/test_taxon_id_utils.py`); **live-verified** against NCBI
    (Alternaria ITS → *Alternaria alternata* 100 %). **200 passing** (was 185).
- **Assembly taxonomy UI (RM-007 step 2).** `pages/1_Assembly_Manager.py`: a project-level
  **default taxon** control (manifest-backed via `set_default_taxon`, shown in an expander that
  opens until set) and a **per-assembly taxon override** in the Manage-assembly section
  (session-mutate + existing `_save()` path), plus **Taxon / Source** columns in the
  registered-assemblies summary showing the effective taxon (`effective_taxon`) and whether it
  is a manual override or inherited from the project default. New records carry
  `taxon`/`taxon_source`. No new `src` logic — reuses the increment-1 manifest helpers.
  Verified end-to-end on a temp project (default + override round-trip, incl. the TSV);
  185 tests still passing (UI not unit-tested).
- **Workflow architecture — manifest schema (D-012 / RM-007 step 1).** First increment of the
  component-page + manifest-chained workflow; pure data layer, no UI yet.
  - `project_manager.py`: project manifest bumped to **schema v2** — adds a project-level
    `default_taxon` and a `workflow` block (`strategy`, `loci`, and per-step `steps` carrying
    `status`/`updated_at`/`outputs`/`notes` over `references`/`rDNA`/`coding`/`primers`/`combine`).
    v1 manifests are read tolerantly (`_ensure_manifest_defaults` fills missing keys in memory;
    a write upgrades the file and preserves any unknown/extra steps).
  - Per-assembly **taxonomy**: `taxon` + `taxon_source` (`manual` | `its_blast`) on registry
    records (backfilled on load/migrate), plus new `assembly_manifest.tsv` columns.
  - New helpers: `load_project_manifest` / `save_project_manifest`, `set_default_taxon`,
    `get_workflow`, `set_workflow_strategy`, `set_workflow_loci`, `update_step`,
    `effective_taxon`, `set_assembly_taxon`. Verified read-tolerant against existing on-disk
    projects (no files modified on load).
  - Tests: 16 new in `tests/test_project_manager.py` (schema/upgrade, workflow-state helpers,
    taxonomy). **185 passing** (was 169).
- **Synonym OR-group NCBI search; dropped forced "complete cds" (D-011, addresses RM-001
  risk-register item).** Reference fetching returned ≈0 hits whenever `"complete cds"` was in
  the query, because fungal markers are deposited mostly as *partial-cds* barcode amplicons.
  - `ncbi_utils.py`: removed the hardcoded `"complete cds"` from every coding-locus
    `LOCUS_CATALOGUE` entry; each coding locus now carries a canonical `gene` keyword plus a
    curated `synonyms` list (TEF1/RPB1/RPB2/TUB2/GAPDH/CAL/ACT/HIS3 + rDNA variants). New
    `build_entrez_query(terms, organism, field)` ORs the terms within the field, phrase-quotes
    multi-word terms, dedupes case-insensitively, and ANDs the organism. `search_ncbi_nucleotide`
    / `search_ncbi_protein` now accept a `str` **or** `list[str]`; new `locus_search_terms`
    helper assembles `[user keyword] + gene + synonyms`. Notes corrected ("Use CDS-only refs"
    removed).
  - `pages/2_Loci_Extraction.py`: per-locus **"Also search N known synonym(s)"** checkbox
    (default on), the **resolved Entrez query shown verbatim** under the search box, fetch
    provenance now stores that resolved query, and the misleading "search 'complete cds' /
    avoid PCR amplicons / best with CDS-only references" captions corrected globally
    (partial-CDS amplicons are fine; Exonerate resolves introns; relaxed BLAST expects
    amplicons).
  - Tests: `tests/test_ncbi_utils.py` — `TestBuildEntrezQuery`, `TestLocusSearchTerms`,
    `TestCatalogueIsBarcodeFriendly` (regression-guards that no catalogue term forces
    'complete cds' and every coding locus has synonyms). **169 passing** (was 155).

## 2026-06-19

- **Degenerate-primer handling for in-silico PCR (D-009) + Matheny RPB1/RPB2 primers
  (D-010).** Established empirically (blastn 2.14.1) that NCBI BLAST does *not* resolve IUPAC
  degenerate codes: a degenerate primer both fails to **seed** (`blastn-short` needs an exact
  7-bp word; e.g. `fRPB2-5F`'s longest concrete run is 5 → silently 0 hits) and is **mis-scored**
  (a biologically-compatible `Y`-over-`C` counts as a mismatch, same as an incompatible
  `R`-over-`C`), so degenerate barcodes were silently under-/non-recovered.
  - `primer_utils.py`: new `expand_degenerate_primer` / `degeneracy_count` /
    `IUPAC_CODES` / `MAX_PRIMER_EXPANSION`. `find_primer_amplicons` now expands fwd+rev into
    concrete oligos, searches each as `FWD_*`/`REV_*` through the existing `blastn-short` +
    RunManager path, buckets by prefix, and **collapses variant duplicates** to one candidate
    per amplicon (lowest edit-distance kept). Over-cap (>8192 variants) raises `ValueError`;
    `run_primer_extraction` reports it as a "primer too degenerate" status and the per-locus
    log records fwd/rev variant counts. Fixes the already-shipped degenerate primers too.
  - `pages/2_Loci_Extraction.py`: the binding-site preview scan catches the over-cap case
    and shows a per-strain `st.warning` instead of crashing.
  - `data/primers.json`: **two new cited pairs** — RPB1 `RPB1-Af`/`RPB1-Cr` (Stiller & Hall
    1997 + Matheny et al. 2002; 400–2000 bp) and RPB2 `fRPB2-5F`/`bRPB2-11R1` (Liu et al.
    1999 + Matheny et al. 2007; 1200–4000 bp), pairing the outermost fwd/rev per gene for the
    largest in-silico span (not bound by wet-lab PCR). Fills the previously primer-less
    `RPB1` locus. Catalogue now **16 pairs**.
  - Tests: `tests/test_primer_utils.py` — expansion correctness, FWD/REV bucketing + dedup,
    over-cap status, the two new pairs, and a `blastn`-guarded end-to-end proof that the
    degenerate RPB2 pair now finds an amplicon the raw primer cannot seed. **155 passing**
    (was 139).
- **Custom-locus primers in PCR Primer mode (`pages/2_Loci_Extraction.py`).** Added an
  "➕ Custom locus" section so in-silico PCR can target a locus not in the catalogue
  (name + fwd/rev + size window; IUPAC degenerate codes allowed), optionally saved to the
  user library. Custom loci merge into `primer_assignments`, so the binding-site preview,
  run loop and combine step pick them up like any catalogue locus; the strain/loci guard
  and the "Ready" summary account for them.
- **Bug fix — PCR Primer mode hid reference-less loci (`pages/2_Loci_Extraction.py`).** The
  Loci selector gated coding loci on `count_refs(l) > 0` for *all* strategies, so loci with
  no NCBI references (RPB2, TUB2, ACT, CAL, GAPDH, HIS3) were unselectable — even in PCR
  Primer mode, which does in-silico PCR and needs no references at all. Now the reference
  gate applies only to the BLAST / Exonerate strategies; primer mode offers every coding
  locus (built-in primers exist for RPB2/TUB2/ACT/CAL/GAPDH/HIS3/TEF1, and Custom… covers
  the rest). The "No refs yet" caption now notes refs are needed for BLAST/Exonerate, not
  primers.
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
