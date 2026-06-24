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
  - _(2026-06-20) Step 4c (ITSx rDNA component page) **done** (D-015):_ new standalone
    `pages/2_ITSx_rDNA.py` (assemblies + regions → ITSx → per-strain + combined), provenance to
    `workflow.steps.rDNA`. **D-015:** extraction outputs are per-project under
    `<project>/results/loci`. Tested helpers `place_rdna_regions` / `combine_rdna_regions`;
    render verified headless. Next: step 4d (Exonerate component page — the coding-loci path,
    the most involved), then 4e (Primers), then step 5 (Workflow orchestrator), then retire the
    monolith + finalize page numbering.
  - _(2026-06-20) Step 4d (Exonerate coding-loci component page) **done**:_ new standalone
    `pages/2_Exonerate.py` — assemblies × coding loci (per-project library) and/or gene-of-interest
    → BLAST-narrow + Exonerate frame-safe CDS (D-008), BLAST HSP fallback when exonerate absent;
    per-project outputs + combined products; provenance to `workflow.steps.coding`. Reuses tested
    extraction functions (no new src logic); render verified headless. Next: step 4e (Primers
    component page), then step 5 (Workflow orchestrator), then retire monolith + finalize numbering.
  - _(2026-06-20) Step 4e (Primers in-silico PCR component page) **done**:_ new standalone
    `pages/2_Primers.py` — primer assignment (catalogue/library/custom) + custom loci + binding-site
    preview/disambiguation → amplicon extraction (D-009 degenerate handling); per-project outputs +
    combined; provenance to `workflow.steps.primers`. Reuses tested `primer_utils`; render verified.
    **All four extraction component pages (4b–4e) are now done.** Next: step 5 (Workflow
    orchestrator — named strategies + manifest-driven checklist), then retire the monolithic
    `pages/2_Loci_Extraction.py` and finalize page numbering.
  - _(2026-06-20) Step 5 + retirement **done** (D-016) → **RM-007 COMPLETE**:_ added
    `pages/6_Workflow.py` (named-strategy selector + manifest-driven checklist linking the
    component pages); ported the relaxed BLAST-amplicon strategy into the Exonerate page as a
    mode toggle (no capability lost); retired the monolith and renumbered to the final layout
    (2 References · 3 ITSx · 4 Exonerate · 5 Primers · 6 Workflow · 7 Alignment · 8 BUSCO ·
    9 Tree). CLAUDE.md updated; all 10 pages render-verified; 210 tests pass. The
    component-page + manifest-chained workflow (D-012) is now the app's structure.

- **RM-008 (2026-06-21) — Standardized "setup" phase (extraction-first workflow).**
  Establish a turnkey, reproducible front end that runs before any tree-building, designed to
  generalize well beyond Alternaria (see D-020). Components, in build order:
  1. **Bundled protein guide set** — ships with the package (like `primers.json`): a *universal
     core* of conserved coding markers (RPB1/RPB2/TEF1/TUB2/ACT/GAPDH/CAL/HIS3) with ≥1
     Ascomycota + ≥1 Basidiomycota full-length protein per locus, plus *swappable lineage packs*
     for clade-specific markers (e.g. Alt a1, endoPG for Dothideomycetes). Drives `protein2genome`
     extraction with no per-project fetching of guides. _[FIRST — extraction foundation]_
  2. **Exonerate amplicon-tip processing → frame-consistent CDS → MACSE** (the manual-Mesquite
     replacement): run reference amplicon tips through the same protein guide to get CDS in a
     common frame, then codon-aware align.
  3. **Reference Taxa / Tips page** (after extraction): import comparison sequences for target
     taxa via (a) paste accessions, (b) target-taxa search (≥3/locus, NCBI-linked review),
     (c) paste accessions from a web BLAST; **auto-classify each accession to its locus** via the
     synonym catalogue; stored **separately from extraction guides** (`<project>/tips/`).
  4. **Availability matrix** (locus × taxon coverage) — missing data is expected/surfaced, not
     fatal.
  5. **Tree-set page** — choose the final locus set for the matrix (separate from tips import).
  6. **GenBank submission export** (final step) — FASTA + feature table from the Exonerate GFF3.
  All steps record status/outputs in the project manifest so a re-opened project knows what's done.
  - _(2026-06-21) Component 1 (bundled protein guide set) **done** (D-020):_
    `data/protein_guides.json` (8 conserved markers × Asco+Basidio full-length RefSeq guides) +
    `protein_guide_utils.py` (loader + user lineage packs) + Exonerate "Reference source" toggle
    (bundled default, no fetching). Verified: bundled cross-genus guides → clean 0-stop CDS from
    a real Alternaria assembly. 10 tests. Next: component 2 (Exonerate-process amplicon tips →
    MACSE), then the Reference Taxa / Tips page.
  - _(2026-06-21) Component 3 (Reference Taxa / Tips page) **done** (D-020):_ new
    `pages/7_Reference_Taxa.py` + `tips_utils.py` — paste a mixed accession list →
    `classify_locus` auto-sorts each to its locus (D-011 synonyms) → stored in a separate
    per-project tips store (`<project>/tips/`); unassigned → manual assign; per-locus review/cull.
    Analysis pages renumbered (Alignment 8 · BUSCO 9 · Tree 10). 6 classifier tests.
    **Deferred:** target-taxa search mode (≥3/locus review) + the availability matrix + tree-set
    page. Next: component 2 (tips → CDS → MACSE) ties tips into the alignment.
  - _(2026-06-22) Component 2 (Exonerate amplicon-tip processing → frame-consistent CDS) **done**
    (D-022):_ new single-purpose page `pages/8_Codon_Tip_Prep.py` + `codon_prep_utils.py` run each
    coding-locus tip through the bundled protein guide (Exonerate `protein2genome`) and **merge the
    framed tips with the user's extracted isolate loci** into three per-locus matrices in
    `<output>/loci/with_tips/`: CDS (intron-stripped, codon-phased), genomic (full gene, exons
    UPPER / introns lower for boundary visibility), protein. References are treated identically to
    the user's own sequences (raw full gene **and** intron-stripped CDS). **Runs no aligner, no new
    dependency** beyond Exonerate; alignment + by-hand codon-structure curation stay on Alignment
    Prep (MACSE/AliView/Geneious are optional companions, never assumed — see D-022). rDNA tips
    bypass to MAFFT; write-and-flag QC; provenance to `workflow.steps.codon_prep`. Analysis pages
    renumbered (**Alignment 9 · BUSCO 10 · Tree 11**). 16 tests → **266 passing**. _Refined from the
    original "→ MACSE" wording per the user:_ the page stops at codon-ready matrices and never
    assumes MACSE (every alignment is hand-checked); back-translation was considered and dropped as
    a built-in. **Still deferred (RM-008):** components 4 (availability matrix), 5 (tree-set page),
    6 (GenBank submission export), and the target-taxa tips search mode.

### Architecture refinements (not numbered RM items)

- _(2026-06-23) **Codon Tip Prep: nucleotide fallback for intron-rich barcodes** (D-027)._ The
  standard fungal TEF1 barcode amplifies a largely **intronic** region (matches the genomic gene but
  not the CDS), so `protein2genome` can't codon-frame it — 0/25 TEF1 tips framed vs 24/24 RPB2. Such
  tips are now oriented (blastn vs the isolate genomic) and added to the **genomic** matrix only,
  flagged nucleotide-only, so the intron-inclusive tree still carries them; CDS/protein stay
  isolate-only for that locus. This is a partial, locus-appropriate realization of the
  orthology-vs-tips idea below (the blastn orient step also confirms the tip belongs to the locus).
- _(2026-06-23) **Reference-quality QC: guide-length filter** (D-025)._ A real coding run exposed
  that fetched taxon-closer Alternaria RefSeq proteins can be **mis-annotated over-long** (an 865-aa
  "partial β-tubulin", 812–814-aa EF1-α) and, being high-identity, out-score the correct guides in
  "both" mode — yielding ~2× over-long / frameshifted CDS that the frame/stop QC can't catch when
  in-frame. Fix: filter project refs by length against the curated bundled-guide expectation
  (`expected_length` ± 30 %), drop-and-flag outliers, opt-out retained. The bundled guides remain the
  trusted full-length anchor.
- _(2026-06-23) **Tip import: per-accession assignment + accession normalization** (D-026)._ Paste →
  look up each on NCBI → assign a locus **per accession** (not one bulk locus) → import; bare RefSeq
  accessions (`NR135944` → `NR_135944`) are auto-repaired and accessions that don't resolve are
  warned, not silently dropped.
- _(2026-06-23, **planned — orthology sanity check vs tips**)._ Reliable *full-CDS* references for
  the standard coding markers are scarce in Alternaria (even the reference genome lacks clean
  annotations for several), and most GenBank entries are **partial PCR amplicons** — which for some
  species could even amplify a **paralog**. So a curated full-length guide alone can't confirm we
  extracted the *intended* ortholog. Planned step: **align each isolate's extracted CDS against the
  pulled tips** for that locus (the partial amplicons) and report agreement, so the user can confirm
  orthology / catch paralog mismatches before committing to a tree. The aim is a defensible
  phylogeny *against antecedent published data* whether or not that data is perfectly annotated.
  Pairs naturally with the deferred RM-008 target-taxa tips search mode. (Tracked here so it isn't
  lost; not yet built.)
  - _(2026-06-24) — **BUILT (D-036)**:_ the **Orthology Check** page (`12_Orthology_Check.py`,
    `orthology_utils.py`) aligns each locus's isolates + tips (MAFFT), builds a p-distance NJ tree,
    and flags divergence outliers **source-blind** (a reference amplicon is judged like an extraction),
    so a paralog — Exonerate's *or* a published reference's — is caught before the tree. Default
    substrate genomic; CDS/protein selectable. **Still deferred:** wiring it as a manifest
    `workflow.steps` entry + Workflow-page checklist item (currently a standalone page, no orphan
    step) — do this when it joins a named strategy.

- _(2026-06-22) **NCBI References page repurposed** (D-023)._ With bundled universal protein guides
  the default extraction source (D-020), the References page is now optional: it fetches
  **taxon-closer protein** orthologs that *supplement* the bundled core (new Exonerate "Bundled +
  project library" mode) and nucleotide refs for the relaxed-BLAST path. **rDNA removed** from it
  (ITSx extracts rDNA; tips compare it). Clarifies the three sequence sources: bundled guides ·
  taxon-closer guide supplement · comparison tips. (Closer guides aid alignment confidence, *not*
  intron finding — Exonerate locates introns in the target; intron-structure hand-checks route
  through tips, D-022.)

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
  - _(2026-06-23) — partially addressed:_ `concat.py` codon-partition path now tested
    (`test_concat.py`, D-032) and `ncbi_utils` transport/retry tested (`test_ncbi_retry.py`,
    D-034). `busco_utils.py` occupancy filtering still untested.

## Known issues — 2026-06-23 architecture audit (Tier-2/3 deferred)

A four-reviewer audit after the D-025→D-031 run fixed the high-impact items (D-032/D-033/D-034).
These confirmed-but-lower-severity findings were deliberately deferred and must not be lost:

- **(robustness) "Tool failed" conflated with "no result"** in `_run_blastn_short` (and the
  primer/relabel paths): a blastn crash returns `[]`, indistinguishable from a real "no binding
  sites" — a false-negative locus reads as a true biological negative. Distinguish rc≠0 from an
  empty parse. (audit M-4)
- **(robustness) Shared fixed-name scratch files** (`blast_hsps.tsv`, `_exonerate_target.fasta`)
  in the per-locus `output_dir` aren't namespaced/cleared before reuse — the LXD-001 stale-file
  class; only ITSx clears stale outputs. Namespace via `tempfile` or clear at runner entry. (H-5)
- **(robustness) `parse_exonerate_gff` `%tcs` truncates on IUPAC codes** → a spurious QC
  cross-check "mismatch" (the CDS itself is correct). Loosen the accept-set; treat empty `tcs`
  as "cross-check unavailable". (H-4)
- **(robustness) `_orig_contig` `__c` strip isn't anchored** — a contig name containing `__c`
  silently disables its coverage filter; fix `re.sub(r"__c\d+$", "", contig)`. (src H-2)
- **(clarity) nt-fallback row status** can't distinguish "no orientation reference available"
  from "wrong locus / contamination" when `include_isolates=False`; and **library-mode** skips a
  locus entirely if the length filter drops all its refs. Make both reasons explicit. (src M-1/M-2)
- **(housekeeping) Dead/unwired config:** `output_base` (written in Tool Settings, never read),
  orphaned manifest step `combine` (never written/read), and `exonerate`/`busco`/`compleasm`/
  `datasets`/`iqtree3` paths not persistable via any UI. (audit L-1/L-2/L-3/L-4)
- **(test depth) Multi-intron coordinate coverage:** `write_region_gff3` and the codon/soft-mask
  paths are tested only on a 2-exon/1-intron synthetic gene; real RPB2/TEF1/RPB1 have several
  introns. Add a ≥3-exon fixture. Also: `filter_records_by_length(include_user=True)` (the
  production path) and `_relabel_itsx_output` dedup/mixed-cov are untested.
- **(method, pre-existing) ITSx still over-extends LSU to the contig end** on the tandem array
  (separate from the D-028 coverage filter) — a length-cap concern, not built.

### Deferred from the D-035/D-036 pre-commit review (2026-06-24, D-037)
- **(collision) Gene-of-interest outputs share `per_strain/<strain>/<name>/` with catalogue loci.**
  A GOI named the same as a catalogue locus (e.g. "RPB2") overwrites that locus's output dir at
  extraction time. `resolve_guide_path` now prefers the GOI guide, but the real fix is to namespace
  GOI outputs (e.g. `GOI_<name>/`). Not built.
- **(robustness) `scan_flagged_cds` validates only `recs[0]`** of a CDS FASTA — a frameshifted
  second record (from `--bestn>1`) is missed. Normal output is single-record; revisit if multi-model
  CDS files are ever written.
- **(strictness) `extract_locus_exonerate` silently coerces an unknown `refine`/`escalate_ceiling`**
  to a default rather than raising — fine for the UI (valid values only), but a programmatic caller
  typo passes silently. Make it raise.
- **(coverage) The page-4 deep-refine loop body and page-12 result rendering are smoke-only.** The
  testable pieces were extracted (`resolve_guide_path`, `analyze_alignment`, `orthology_check`); the
  Streamlit glue (snapshot/restore, re-merge, table rendering) is exercised only by the import smoke
  test.

### Deferred from the second robustness review (2026-06-24, D-038)
- **(durability) Deep-refine snapshot lives in `/tmp` and restore is `rmtree`+`move`.** Python
  exceptions are caught and the original is restored, but a *hard* kill (OS/power) in the
  rmtree→move window, or mid-rewrite, could lose the locus dir (only copy in `/tmp`, cleared on
  reboot). Accepted for now — the whole deep-refine pass is explicitly re-runnable and the original
  is regenerable by re-extracting. If hardened: snapshot to a sibling on the same filesystem and
  swap with `os.replace` (atomic rename) instead of cross-fs `move`.
- **(contrived) `_dedupe_ids` can re-collide.** A matrix containing both two `MH1.1` *and* a literal
  `MH1.1__dup2` would synthesize a second `MH1.1__dup2`. Loop the suffix until unseen if it ever
  matters; not worth the code today.
- **(cosmetic) Page-12 isolate/tip caption is computed pre-dedup** (a fresh `SeqIO.parse`), so with a
  double-imported accession the caption count can be one off from the deduped results table. The
  dedup warning already surfaces the slip; left as-is.
- **(low risk) Nucleotide p-distance treats RNA `U` vs DNA `T` as a mismatch.** Matrices here are
  DNA, so no real-world impact; note it if RNA tips are ever fetched.

## Done (high-level milestones)

See `CHANGELOG.md` for the detailed, dated record. Recent shipped work (from git history):

- BUSCO integration and assembly-stats display; primer mode surfaced in UI.
- Import scanners consolidated; assembly schema inconsistencies fixed.
- PCR primer-based locus extraction strategy added.
