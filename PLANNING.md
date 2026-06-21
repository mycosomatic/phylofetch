# phylofetch — Planning

> **Append-only.** Do not delete past entries. When a plan item is superseded, comment it
> out with `<!-- -->` and add a dated pointer to what replaced it (see `DECISIONS.md`).
> Review this file and `DECISIONS.md` at the start of every session.

## Purpose

Overarching roadmap and goals for phylofetch. This is the "where are we going and why"
document. Day-to-day completed work is recorded in `CHANGELOG.md`; the reasoning behind
choices lives in `DECISIONS.md`.

## Project goal

A Streamlit-based, pip-installable bioinformatics app for fungal genome work:

> assembly → loci extraction → alignment/trimming/concatenation → BUSCO phylogenomics →
> ML/Bayesian tree prep

Designed to stand alone and to plug into a future umbrella Streamlit app as a dependency.

## Roadmap / open threads

<!-- Add roadmap items below. Use a short id + date so DECISIONS.md can reference them. -->

<!-- Superseded by RM-001..RM-005 below (2026-06-18, post strategy review).
- _(2026-06-18)_ No roadmap items recorded yet. First entry to be added when we scope our
  next piece of work together.
-->

These items came out of the 2026-06-18 strategy/implementation review (see "Review
findings — risk register" below and decisions D-005/D-006). Priority order reflects the
user's decisions in that session.

- **RM-001 (2026-06-18) — Type-material-aware reference acquisition. _[FIRST PRIORITY]_**
  Make "compare against type material" a first-class, operational step rather than something
  done by hand. NCBI-first: `sequence_from_type[Filter]` for barcode references; capture
  GenBank voucher/strain/`/type_material` metadata; flag `type_material` /
  `refseq_category` on downloaded assemblies; surface type status in reference tables and on
  tip labels; log the exact Entrez query via RunManager. (See D-006.)
- **RM-002 (2026-06-18) — CDS extraction QC.** Add translation / frame / internal-stop /
  splice-site validation to `blast_loci_utils`; use the GT-AG signal as a gate, not just a
  log line; add paralog awareness; route its BLAST call through RunManager. Prevents silent
  corruption of protein-coding loci.
  - _(2026-06-19) — substantially addressed by **D-008**:_ `exonerate_utils.py` replaces
    HSP-as-exon stitching for coding loci with Exonerate spliced alignment (frame-safe CDS,
    explicit GT-AG introns), adds translation / `len % 3` / internal-stop QC + `%tcs`
    cross-check, and `--bestn` paralog awareness. Exonerate calls route through RunManager.
    The HSP path is demoted to a documented fallback. _Remaining RM-002 follow-ups:_ apply
    the same frame/stop QC to the relaxed BLAST-amplicon path, and route the Exonerate
    contig-narrowing BLAST through RunManager too (currently only the Exonerate call is
    logged).
- **RM-003 (2026-06-18) — Species-delimitation workflow.** Treat per-locus alignments as
  first-class: per-locus gene trees + a concordance/conflict view (GCPSR), instead of
  defaulting to the concatenated supermatrix for the novelty goal.
- **RM-004 (2026-06-18) — Phylogenomics breadth.** Add thin ASTRAL (coalescent species
  tree, robust to ILS) and ANI (fastANI/pyANI genome distance) runners alongside the
  existing IQ-TREE2 concatenation path. (See D-005.)
- **RM-005 (2026-06-18) — Headless reproducibility.** A config-driven `phylofetch run
  project.yaml` over the existing functions, input-file hashing, and pinned tool versions in
  `environment.yml`, so the whole pipeline is re-runnable off the browser session.
- **RM-006 (2026-06-20) — rDNA sequence-variant consensus from Illumina reads. _[FUTURE,
  deferred by user]_** Extracting rDNA (ITS/LSU/SSU) from an *assembly* yields a single
  collapsed copy that does not represent the intragenomic **sequence variants** of the
  multi-copy rDNA tandem array (and the region may not assemble cleanly from short reads
  alone, though hybrid assemblies often resolve it). Planned step: map **Illumina reads only**
  to the rDNA region and emit a *single* consensus carrying **IUPAC ambiguity codes** at
  variable sites — i.e. represent the variation rather than phasing out individual rDNA
  variants, *unless* a meaningful phasing approach emerges. This introduces a **new input
  type** (no fastq/mapping/variant-calling code exists today; the pipeline is currently
  assembly-only). Intended as a QC / marker-reliability signal, **not** as multiple tree
  tips. Terminology fixed to "sequence variants" (user preference). Not to be built now;
  recorded per the user's "tuck into future plans" instruction (2026-06-20).
- **RM-007 (2026-06-20) — Component-page + manifest-chained extraction workflow.**
  Implements **D-012**. Phased so each increment ships behind tests and the old page keeps
  working until parity:
  1. **Manifest schema** — extend the project manifest + assembly registry with per-assembly
     **taxonomy** (project default + override) and a **workflow/step state** block (selected
     loci, chosen strategy, per-step status/outputs); add load/save/update helpers + tests.
  2. **Assembly taxonomy UI** — set project-default taxon + per-assembly override in the
     Assembly Manager; persist to the manifest.
  3. **ITS→BLAST provisional-ID feeder** — extract ITS (ITSx) → BLAST → suggest closest taxon
     → write back to the assembly's taxonomy (the "pick from ITS BLAST result" path).
  4. **Component pages** — split the monolith into standalone *NCBI References*, *ITSx*,
     *Exonerate*, *Primers* pages, each reading/writing the manifest but runnable ad-hoc.
  5. **Workflow / Strategy page** — orchestrator with named strategies and a stepwise,
     manifest-driven checklist that chains the component pages.
  6. **Reference provenance into manifest** — depends on **D-013** (references global vs
     per-project).
  Old monolithic `pages/2_Loci_Extraction.py` retained/redirected until the new pages reach
  feature parity (working-agreement caution on removing code).
  - _(2026-06-20) Step 1 (manifest schema) **done**:_ project manifest v2 (`default_taxon` +
    `workflow`/step-state block), per-assembly `taxon`/`taxon_source`, and the helper API
    (`load_project_manifest`, `update_step`, `set_assembly_taxon`, `effective_taxon`, …) landed
    in `project_manager.py` with 16 tests. Backward-compatible / read-tolerant.
  - _(2026-06-20) Step 2 (Assembly Manager taxonomy UI) **done**:_ project-default taxon
    control + per-assembly override + Taxon/Source summary columns in
    `pages/1_Assembly_Manager.py`, reusing the step-1 helpers.
  - _(2026-06-20) Step 3 (ITS→BLAST provisional-ID feeder) **done** (D-014):_ remote NCBI
    `blastn` of the ITSx ITS region with a ranked-organism picker in the Assembly Manager;
    chosen taxon written with `taxon_source=its_blast`. New `src/phylofetch/taxon_id_utils.py`
    + 15 tests; live-verified against NCBI.
  - _(2026-06-20) **D-013 resolved = per-project.** Step 4a (ncbi_utils `ref_dir` threading +
    `project_ref_dir` + 6 tests) **done**, backward-compatible (global default unchanged).
    Next: step 4b (NCBI References component page) using the project-scoped dir, then the
    ITSx / Exonerate / Primers component pages, then step 5 (Workflow orchestrator). The old
    monolithic `pages/2_Loci_Extraction.py` stays until the new pages reach parity._
  - _(2026-06-20) Step 4b (NCBI References component page) **done**:_ new standalone
    `pages/2_NCBI_References.py` — loci checkboxes → preview hit counts (`ncbi_search_count`)
    → fetch into the per-project library → review/cull; organism defaults to the project
    taxon; provenance to `workflow.steps.references`. Render verified headless (AppTest).
    Runs alongside the monolith (interim duplicate "2" sidebar entry, by user choice).
    **Deferred:** an "import from global library" button (only needed if old global refs
    exist). Next: step 4c (ITSx page), 4d (Exonerate), 4e (Primers), then step 5 (Workflow).

### Review findings — risk register (2026-06-18)

Recorded so we don't lose them; each maps to a roadmap item above.

- **(corruption risk) HSP-as-exon CDS stitching has no frame validation.**
  `blast_loci_utils.extract_from_hsps` concatenates BLAST HSPs as exons with no translate /
  stop-codon / `len % 3` check; a 1-2 bp HSP-boundary error frameshifts the CDS silently.
  Codon partitions assume frame starts at base 1. Splice sites are computed but not gated.
  No paralog screening at `min_pident=70`. → RM-002.
- **(provenance gap) Scientific core bypasses RunManager.** `blast_loci_utils.run_blast` is a
  raw subprocess; all `ncbi_utils` fetching is unlogged (no record of which query produced
  which reference set, when). → RM-001 (fetch logging) + RM-002 (BLAST logging).
- **(fragility) Pipeline state lives in `st.session_state`** (busco_results, queues,
  selections) and isn't persisted until export; closing the tab loses in-progress work.
  → RM-005.
- **(reference quality) `{gene}[Title]` NCBI search is fragile** — biases toward "complete
  cds" titles, misses inconsistently-titled barcode records. → RM-001.
  - _(2026-06-20) — addressed by **D-011**:_ dropped the forced `"complete cds"` from the
    coding-locus catalogue, added curated per-locus **synonym OR-groups** + a phrase-quoting
    query builder (`build_entrez_query`), surfaced the resolved Entrez query in the UI, and
    corrected the self-contradictory "use CDS-only / avoid amplicons" guidance (partial-CDS
    barcode amplicons are the norm and work fine — Exonerate resolves introns, relaxed BLAST
    expects amplicons). Recall of partial-cds barcodes and synonym-titled records restored.
- **(method) Nucleotide BUSCO default saturates at depth;** consider amino-acid / codon-aware
  for deeper comparisons. Concatenation-only ignores ILS (high near speciation). → RM-004.
- **(coverage) No tests for `concat.py`, `busco_utils.py`, `ncbi_utils.py`** — the
  partition-offset math and occupancy filtering are silent-error-prone.

## Done (high-level milestones)

See `CHANGELOG.md` for the detailed, dated record. Recent shipped work (from git history):

- BUSCO integration and assembly-stats display; primer mode surfaced in UI.
- Import scanners consolidated; assembly schema inconsistencies fixed.
- PCR primer-based locus extraction strategy added.
