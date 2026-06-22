# phylofetch — Decision History

> **Append-only.** Do not delete past entries. When a decision is reversed, **comment out**
> the old entry with `<!-- -->` and add a dated pointer to the new decision that supersedes
> it. Review this file and `PLANNING.md` at the start of every session.

Format for each entry:

```
### D-NNN (YYYY-MM-DD) — short title
- **Decision:** what we decided.
- **Why:** the reasoning, including scientific justification where relevant.
- **Alternatives considered:** what we rejected and why.
- **Status:** active | superseded by D-XXX (YYYY-MM-DD).
```

---

### D-001 (2026-06-18) — Working agreement and living documents
- **Decision:** Adopt a session protocol recorded in `CLAUDE.md` (read first, every
  session) plus three append-only living documents: `PLANNING.md`, `DECISIONS.md`,
  `CHANGELOG.md`. Two branches only (`main`, `dev`); superseded entries are commented out
  with a dated pointer rather than deleted.
- **Why:** The work supports active scientific research and must be transparent, traceable,
  and repeatable. A durable file-based record survives across sessions where in-context
  memory does not, and lets the user audit for backtracking.
- **Alternatives considered:** (a) Single combined log file — rejected, harder to scan as it
  grows. (b) Relying on session memory / CLAUDE.md alone — rejected, in-session context can
  be summarized and lose detail.
- **Status:** active.

### D-002 (pre-2026-06-18) — Pre-existing technical decisions
- **Decision:** Carry forward decisions already documented in `CLAUDE.md`: config at
  `~/.phylofetch/config.json`; genus-agnostic (no Alternaria defaults); NCBI bracket-style
  rich FASTA headers; `RunManager` logs every external tool call with versions/timestamps;
  per-locus extraction logs; src layout; MACSE optional with graceful fallback. Critical
  bug fixes LXD-001 (ITSx) and LXD-002 (BLAST HSP grouping) remain in force with regression
  tests.
- **Why:** These were established before this decision log existed; recording them here so
  the log is the single source of truth going forward.
- **Alternatives considered:** n/a — retroactive capture.
- **Status:** active. See `CLAUDE.md` for full detail.

### D-003 (2026-06-18) — Citable, user-extensible PCR primer library
- **Decision:** Convert the hardcoded primer catalogue into a packaged data file
  (`src/phylofetch/data/primers.json`, shipped via `[tool.setuptools.package-data]`). Every
  built-in pair carries `source`, `citation`, and `reference_url` (primary literature; rDNA
  cross-checked against the UNITE primer table). Users may add custom pairs that persist in
  `~/.phylofetch/primers.json` and merge on top of the built-in catalogue (user wins on name
  clash). Primer searches route through `RunManager`; each `LOCUS_extraction.log` records the
  primer, its citation, the command, and the chosen binding site.
- **Why:** Phylogenetic results must be reproducible and defensible. A primer choice is a
  scientific claim; shipping it with its citation makes every extraction traceable to the
  literature, and provenance logging makes runs repeatable. Made by the browser/cloud session
  (commit `4a03968`, Claude Sonnet 4.6); recorded here so our decision log stays authoritative.
- **Alternatives considered:** (a) Keep primers hardcoded in `primer_utils.py` — rejected, not
  citable or user-extensible. (b) Citations in code comments only — rejected, not surfaced to
  the user or the per-locus logs.
- **Status:** active. See `CLAUDE.md` → "Primer library" for detail.

### D-004 / PRM-001 (2026-06-18) — Data correction: ACT-512F primer sequence
- **Decision:** The built-in `ACT-512F` primer was stored as a corrupted 35-nt sequence (a
  concatenation). Corrected to the canonical 20-nt sequence `ATGTGCAAGGCCGGTTTCGC`
  (Carbone & Kohn 1999). A regression test guards against reintroduction.
- **Why:** A wrong primer sequence silently corrupts in-silico PCR results — exactly the kind
  of well-established-knowledge contradiction our working agreement guards against. The
  canonical sequence is the literature standard. Made in commit `4a03968`.
- **Alternatives considered:** n/a — the stored sequence was simply incorrect.
- **Status:** active. Identifier `PRM-001` is also used in `CLAUDE.md` and the test suite.

### D-005 (2026-06-18) — Inference boundary: what lives "under one roof" vs. ported out
- **Decision:** phylofetch owns the *data-assembly layer* and *thin runners* for inference
  that does not require trust-establishing diagnostics: locus extraction, NCBI/assembly
  acquisition, BUSCO/Compleasm, alignment/trim/concat/partitioning, IQ-TREE2 (ML), and —
  to be added — thin **ASTRAL** (coalescent species tree) and **ANI** (fastANI/pyANI)
  runners. Bayesian inference & divergence dating (**BEAST2 / MrBayes**), rigorous species
  delimitation (**BPP / ASAP / GMYC / bPTP**), and publication-grade tree figures
  (**iTOL / FigTree / ggtree**) are **ported out**: phylofetch generates their inputs but
  does not run or diagnose them. In-app tree rendering stays quick-look only.
- **Why:** Tools whose output cannot be trusted without diagnostics (MCMC convergence/ESS,
  delimitation model adequacy) create *false confidence* when thinly wrapped in Streamlit,
  and reproducing those diagnostics in-app is a maintenance tarpit. The genuine, hard value
  is assembling the dataset (extraction + acquisition + matrix building); inference of that
  kind is a commodity best left to the established, audited tools. Coalescent + ANI are
  cheap, deterministic-enough additions that directly serve the two project goals (ILS-aware
  isolate relationships; fast same-species check) and reuse data already on disk.
- **Alternatives considered:** (a) Maximal in-app inference incl. Bayesian/dating/delimitation
  — rejected (false-confidence + tarpit). (b) Thin ML runner only (no coalescent/ANI) —
  rejected as under-powered for the recent-divergence regime the project targets.
- **Status:** active. Chosen by the user in the 2026-06-18 review session.

### D-006 (2026-06-18) — Next-work priority: type-material-aware reference acquisition first
- **Decision:** Of the threads surfaced in the 2026-06-18 strategy review, build
  **type-material support (RM-001)** first, NCBI-first. The remaining roadmap items
  (RM-002 CDS QC, RM-003 delimitation workflow, RM-004 coalescent/ANI, RM-005 headless)
  follow; ordering may be revisited.
- **Why:** "Prioritize comparing with type/reference material" is an explicit project goal
  and is currently *unsupported* in code — NCBI fetching ignores the
  `sequence_from_type[Filter]`, discards GenBank voucher/strain/type metadata, and surfaces
  no type status. It is also the highest-value-per-effort gap: small Entrez/metadata changes
  unlock a stated objective.
- **Alternatives considered:** (a) CDS-extraction QC first (RM-002) — recommended by the
  assistant as the top correctness risk, but the user prioritized the type-material goal;
  RM-002 remains the next correctness item. (b) Delimitation or reproducibility first —
  deferred.
- **Status:** active. First-pass scope confirmed in D-007.

### D-007 (2026-06-18) — RM-001 first-pass scope: references as comparison tips + queries, NCBI-only
- **Decision:** The type-material first pass will (1) make fetched references **first-class
  comparison tips** merged into the per-locus combined FASTAs — *new* plumbing, since today
  the combine step only walks per-strain assembly-derived outputs and never merges
  references — **and** keep them usable as BLAST extraction queries; (2) be **NCBI-only**
  (nucleotide type filter + GenBank voucher/strain/`/type_material` metadata + assembly
  type/representative flags); UNITE and RefSeq Targeted Loci are deferred. Behaviour is
  **prefer-and-flag**: type-grade sequences sorted first, **falling back to best non-type
  hits by % identity then coverage** when type material is sparse/absent (common for
  closely-related species). Ex-type cultures are treated as **equal in grade** to
  holotype/syntype/isotype/neotype/epitype (the *kind* is recorded and displayed, but does
  not change ranking). "Validity" of an ex-type is a **human judgement**: we surface
  `/type_material` kind + strain + culture-collection voucher for the user to cull; we do not
  auto-trust the tag. Type tips are **force-included even when partial** (labelled partial),
  consistent with "types are useful in any analysis."
- **Why:** Comparison against type material only works if type sequences are *tips in the
  tree*; the filter is the easy part, the tip plumbing is the real value. NCBI-only covers
  most of the goal with the smallest surface; UNITE/RefSeq use a different
  identifier/curation-based flow and can be added later without rework.
- **Alternatives considered:** (a) References as extraction queries only — rejected, doesn't
  place isolates against types. (b) Type-only hard filter — rejected, discards unflagged
  type-derived records and leaves thin coverage in poorly-sampled genera. (c) Include
  UNITE/RefSeq-TL now — deferred to keep the first pass small.
- **Status:** active. Implementation plan to be approved before coding.

### D-008 (2026-06-19) — Exonerate spliced alignment for frame-safe CDS + gene-of-interest extraction
- **Decision:** Adopt **Exonerate** for protein-coding locus extraction via a new module
  `src/phylofetch/exonerate_utils.py`. The coding strategy is a **hybrid pipeline**:
  tblastn/blastn first narrows to the single best contig (reusing
  `select_best_locus_group`), then Exonerate refines the gene on that contig with an
  explicit intron/splice model — `protein2genome` for protein references, `coding2genome`
  for nucleotide CDS references (model auto-selected from `detect_fasta_type`). Coordinates
  stay in contig space (no offset math). When BLAST narrowing finds nothing, Exonerate runs
  against the whole assembly. Two query sources are supported: the existing per-locus
  Reference Library, **and** a free "gene of interest" input (paste/upload an arbitrary
  ortholog → get just the CDS + exon model). Outputs mirror the BLAST path plus a translated
  `*_protein.fasta`, and add CDS QC (reading-frame `len % 3`, internal-stop count, GT-AG
  splice tally, %identity/coverage, `%tcs` cross-check). The `extract_from_hsps` HSP-as-exon
  path is **demoted to a documented fallback** for coding CDS (used only when `exonerate`
  is not on PATH, with a visible frame-safety warning) and remains the *primary* path for
  the relaxed PCR-amplicon strategy. By default imperfect models are written-and-flagged,
  not dropped (a `strict_qc` toggle rejects internal-stop / frameshift CDS); consistent with
  D-007 keeping partial/type sequences.
- **Why:** The HSP-as-exon stitching had no reading-frame check — a 1-2 bp HSP-boundary
  error silently frameshifts the CDS (PLANNING.md → RM-002 risk register). Exonerate's
  spliced alignment is the established tool for accurate exon/intron boundaries and a
  translatable CDS, directly addressing RM-002, and it unlocks the user's "locate any gene
  of interest from an ortholog" goal. It is deterministic gene-finding in the data-assembly
  layer, so it fits the D-005 inference boundary (no trust-requiring diagnostics). Demoting
  rather than deleting the HSP path honours the working agreement's caution about removing
  code and keeps a graceful fallback (cf. the MACSE pattern). Validated against the
  installed Exonerate 2.4.0 on plus/minus strands and `coding2genome`.
- **Alternatives considered:** (a) New parallel 4th strategy leaving the BLAST CDS path
  intact — rejected by the user in favour of actually upgrading the coding path. (b) Replace
  the coding path with whole-genome Exonerate (no BLAST narrowing) — rejected as too slow on
  30-40 Mb assemblies × many loci × many strains. (c) **miniprot** (modern protein-to-genome
  aligner) instead of Exonerate — noted as a possible future alternative, but the user chose
  Exonerate after discussion and it is the long-established standard for this task. (d) Hard
  reading-frame gate that drops imperfect CDS — rejected as default (would silently discard
  partial/type sequences, contra D-007); offered as opt-in `strict_qc`.
- **Status:** active. Implemented 2026-06-19; 23 network-free + binary-guarded tests
  (`tests/test_exonerate_utils.py`), full suite 139 passing. `extract_from_hsps` retained
  as documented fallback per the working agreement.

### D-009 (2026-06-19) — Degenerate-primer expansion for the in-silico PCR (BLAST) search
- **Decision:** Make the primer-based locus search correctly handle IUPAC degenerate
  bases by **expanding each primer into the full set of concrete oligonucleotides it
  represents** (`expand_degenerate_primer` / `degeneracy_count` in `primer_utils.py`),
  writing every variant as a separate `FWD_*` / `REV_*` query, searching them together
  through the existing `blastn-short` + RunManager path, then collapsing variant
  duplicates back to one candidate per amplicon (lowest edit-distance pairing kept). A
  non-degenerate primer expands to a single variant, so concrete and degenerate primers
  share one code path. A `MAX_PRIMER_EXPANSION = 8192` cap guards pathological input
  (e.g. a primer full of N's): over-cap raises `ValueError`, surfaced by
  `run_primer_extraction` as a "primer too degenerate to search" status and by the
  binding-site preview as a per-strain warning (never a silent empty result). The
  per-locus extraction log records the fwd/rev variant counts.
- **Why:** Empirically verified (blastn 2.14.1) that NCBI BLAST does **not** resolve IUPAC
  codes, and degenerate primers fail two ways: **(1) seeding** — `blastn-short` requires an
  exact word match (`word_size 7`), so a primer with no 7-bp run of plain ACGT seeds
  *nothing* and silently returns no hits (e.g. the existing `fRPB2-5F`, longest concrete
  run = 5 → 0 hits at word_size 7); **(2) scoring** — a degenerate base is counted as a
  mismatch even when biologically compatible (`Y` over a template `C` scores identically to
  an incompatible `R` over `C`), so every degenerate position eats the `max_mismatches`
  budget and drags the hit under the `perc_identity` floor. The net effect was silent
  under-/non-recovery of exactly the protein-coding barcodes that are *most* degenerate
  (RPB1/RPB2/TEF1). Expansion fixes both at once and, because each variant is concrete,
  makes the residual mismatch count reflect *real* biological mismatch against the
  best-matching variant — and it repairs the already-shipped degenerate primers too. This
  was the prerequisite for adding the heavily-degenerate Matheny-lab RPB1/RPB2 primers
  (some with 2 N's + multiple R's); adding them without it would have shipped cited primers
  that silently don't work, contra the working agreement's scientific-soundness rule.
- **Alternatives considered:** (a) Lower the BLAST `word_size` and inflate `max_mismatches`
  by the degeneracy — rejected: fixes seeding but not scoring, and a biologically-compatible
  degenerate position should not count as a mismatch at all. (b) Swap the engine for EMBOSS
  `primersearch`, which understands IUPAC natively — rejected for now: new external
  dependency + a larger rewrite that discards the existing edit-distance / provenance
  machinery; revisitable. (c) Custom IUPAC-aware Python sliding-window matcher — rejected:
  reimplements alignment, loses the blastn/RunManager provenance, slower on large
  assemblies. (d) Add the Matheny primers and only document the limitation — rejected by the
  user in favour of fixing the engine.
- **Status:** active. Implemented 2026-06-19; expansion + bucketing + dedup + over-cap
  handling with mocked and `blastn`-guarded end-to-end tests in `tests/test_primer_utils.py`
  (51 primer tests; full suite 155 passing, was 139). Paired with the catalogue additions
  below.

### D-010 (2026-06-19) — Matheny-lab RPB1/RPB2 primers added as largest-span in-silico pairs
- **Decision:** Add two built-in, fully-cited pairs from the Matheny lab primer compilation
  (https://wordpress.clarku.edu/polypeet/datasets/primer-information/): **RPB1**
  `RPB1-Af` (Stiller & Hall 1997) / `RPB1-Cr` (Matheny et al. 2002), and **RPB2**
  `fRPB2-5F` (Liu et al. 1999; domain 5) / `bRPB2-11R1` (Matheny et al. 2007; domain 11).
  The Matheny page lists primers as a *menu* for two nested PCRs, not fixed pairs; because
  in-silico extraction from an assembly is **not constrained by wet-lab PCR**, we pair the
  **outermost forward and outermost reverse** of each gene to span the largest region, with
  deliberately wide amplicon windows (RPB1 400–2000 bp; RPB2 1200–4000 bp) sized to include
  introns rather than the wet-lab product. `reference_url` points at the Matheny
  compilation page (the verifiable source the user supplied) with full primary-literature
  citations in the `citation` field, rather than a hand-typed DOI we could not confirm.
  This also fills a real gap: `RPB1` previously had **no** primer pair in the catalogue.
- **Why:** RPB1 and RPB2 are core protein-coding phylogenetic markers for Agaricales /
  Basidiomycota (and broadly across fungi); the Matheny set is the standard source. Pairing
  outermost-to-outermost maximises recovered sequence for tree-building, which is the goal
  in silico, and the binding-site preview/disambiguation already lets the user cull any
  off-target pairing that a wide window admits. Depends on D-009 — these primers are heavily
  degenerate and would silently fail without expansion.
- **Alternatives considered:** (a) Encode the faithful nested wet-lab pairs (domains 5–7 and
  7–11 separately) — rejected per the user's steer to default to the largest extraction.
  (b) `reference_url` as a per-primer DOI — partially rejected to avoid shipping an
  unverified DOI; the page URL is exact provenance and the citations carry the papers.
  (c) Add the full menu of individual primers — deferred; the two largest-span pairs cover
  the stated need and the per-locus "Custom…" entry covers any other combination.
- **Status:** active. Catalogue now 16 built-in pairs; sequences regression-tested
  (`TestMathenyRPBPairs`). Candidate tie-break order (smallest-first) left unchanged pending
  user confirmation — see open question in the session note.

### D-011 (2026-06-20) — Synonym OR-group NCBI search; drop forced "complete cds"
- **Decision:** Rebuild the reference-search query for coding loci. (1) **Remove the
  hardcoded `"complete cds"`** from every `LOCUS_CATALOGUE` coding entry (was
  `"gapdh complete cds"`, `"tef1 complete cds"`, …). (2) Give each coding locus a canonical
  `gene` keyword plus a curated `synonyms` list (symbol + spelled-out name + common
  abbreviations), and add a query builder `build_entrez_query(terms, organism, field)` that
  ORs the terms within the field, phrase-quotes multi-word terms, de-duplicates
  case-insensitively, and ANDs `organism[Organism]`. `search_ncbi_nucleotide` /
  `search_ncbi_protein` now accept a `str` **or** a `list[str]` and route through the
  builder; a `locus_search_terms(locus, user_term)` helper assembles `[user keyword] + gene
  + synonyms`. (3) In the UI: a per-locus **"Also search N known synonym(s)"** checkbox
  (default on), the **resolved Entrez query shown verbatim** under the search box
  (transparency), the stored fetch-provenance `query` now records that resolved query, and
  the misleading captions/notes ("Search for 'complete cds' entries / Avoid PCR amplicons /
  Use CDS-only refs / Best with CDS-only references") are corrected globally to: partial-CDS
  barcode amplicons are the norm for fungal markers and work fine; complete-CDS refs are
  optional. The per-locus search field stays editable, so the user can still hand-type
  `complete cds` if they ever want it. Synonym sets used (title field): TEF1 = tef1 / tef-1 /
  tef1-alpha / translation elongation factor 1-alpha / elongation factor 1-alpha / EF-1alpha;
  RPB1 = rpb1 / RNA polymerase II largest subunit / DNA-directed RNA polymerase II subunit
  RPB1; RPB2 = rpb2 / RNA polymerase II second largest subunit / …RPB2; TUB2 = tub2 / benA /
  beta-tubulin / beta tubulin / btub; GAPDH = gapdh / gpd / gpdh / gpdA /
  glyceraldehyde-3-phosphate dehydrogenase; CAL = calmodulin / cmd / cmdA / CaM; ACT = actin
  / act1 / actA; HIS3 = histone H3 / his3 / hH3; plus rDNA variants for ITS/LSU/SSU.
- **Why:** Empirically the user got ≈0 hits whenever `"complete cds"` was in the query.
  Fungal phylogenetic markers are deposited overwhelmingly as **partial-cds barcode
  amplicons** ("…gene, partial cds"), so forcing `complete cds` excludes exactly the records
  wanted for closely-related-taxon phylogenetics — the recall-crushing fragility already in
  the PLANNING.md risk register (RM-001: "`{gene}[Title]` search biases toward 'complete cds'
  titles, misses inconsistently-titled barcode records"). The original rationale for clean
  complete-CDS refs (so naïve blastn HSP-as-exon stitching wouldn't trip over introns) is
  largely obsolete after **D-008**: Exonerate resolves introns/splice sites, and the relaxed
  BLAST path *expects* genomic amplicons (`require_complete_cds=False`). The catalogue/caption
  guidance ("Use CDS-only refs", "Avoid PCR amplicons") was also **self-contradictory** with
  the relaxed-amplicon strategy the user had selected. Synonym OR-groups fix the second half
  of the fragility (gene-name inconsistency: gpd vs gpdh vs the spelled-out name), and
  showing the resolved query keeps the change transparent/traceable per the working
  agreement.
- **Alternatives considered:** (a) **Just drop "complete cds"** (single-keyword query) —
  rejected by the user as not robust to name variation. (b) **Synonyms + keep "complete cds"
  as an opt-in toggle** — declined; the editable keyword box already lets a user re-add it,
  so a dedicated toggle was unnecessary UI. (c) **Strategy-aware Reference-Library guidance**
  (different advice per extraction strategy) — the user chose to fix the wording globally
  instead, simpler and avoids coupling the library tab to the strategy radio. (d) Use
  `[Gene]`/`[All Fields]` instead of `[Title]` for nucleotide — deferred; `[Title]` with
  phrase-quoted synonyms is predictable and high-recall, and `[All Fields]` adds noise.
  (e) Auto-trust no completeness filter at all even for Exonerate clean queries — acceptable
  because a partial CDS still works as an Exonerate/BLAST query.
- **Status:** active. Implemented 2026-06-20; 14 new network-free tests
  (`TestBuildEntrezQuery`, `TestLocusSearchTerms`, `TestCatalogueIsBarcodeFriendly` in
  `tests/test_ncbi_utils.py`), full suite **169 passing** (was 155). Synonym sets are
  curated defaults — open to per-locus adjustment by the user (mycology domain call).
  Addresses the RM-001 "complete cds" risk-register item.

### D-012 (2026-06-20) — Workflow architecture: standalone component pages chained by a disk-backed project manifest
- **Decision:** Reorganise the loci-extraction workflow into (a) standalone, single-purpose
  **component pages** — *NCBI References*, *ITSx (rDNA)*, *Exonerate (coding / gene-fishing)*,
  *Primers (in-silico PCR)* — each usable on its own, **plus** (b) a **Workflow / Strategy**
  page that chains them in a logical, stepwise order for a chosen *strategy* (a named recipe,
  e.g. "fungal barcodes" = ITSx for rDNA + BLAST→Exonerate for coding + type-material refs).
  Shared state is a **disk-backed project manifest** (extending the existing
  `metadata/project_manifest.json` + assembly registry), **not** `st.session_state`: it
  records assemblies, per-assembly **taxonomy** (project-level default + per-assembly
  override), selected loci, the chosen strategy, and per-step status/outputs. Every page
  reads/writes the manifest, which becomes the single source of truth and the basis for
  headless re-runs (RM-005). Taxonomy is set **manually (closest taxon)** or **picked from the
  ITS-extraction BLAST result**. The rDNA-vs-coding split is retained (ITSx for rDNA;
  BLAST→Exonerate for coding, D-008). The existing standalone "gene of interest" Exonerate
  input is promoted to the first-class Exonerate page.
- **Why:** The current monolithic page 2 couples reference fetching, four extraction
  strategies, and combining behind one strategy radio, and chains them through
  `st.session_state` — the exact fragility in the risk register ("pipeline state lives in
  session_state, lost on reload"). The user wants components usable both standalone (e.g. fish
  a gene out of one assembly with Exonerate) and chained into a repeatable phylogenetic
  pipeline. A manifest-backed design gives both, makes runs transparent/traceable/repeatable
  (the working-agreement mandate), and directly enables RM-005 headless runs. Most extraction
  *logic* already lives in `src/` and is reused unchanged — this is largely a UI/orchestration
  reorganisation, not a science rewrite. Chosen by the user (2026-06-20): "Component pages +
  Workflow" layout and "Disk-backed project manifest" chaining.
- **Alternatives considered:** (a) Hybrid "Extraction Tools" page (tools as tabs) + Workflow
  page — rejected in favour of fully separate component pages. (b) Single stepwise wizard in
  one page — rejected; doesn't give first-class standalone tools. (c) `st.session_state`
  chaining persisted only on export — rejected as the known reload-fragility (RM-005).
- **Status:** active (architecture). Implementation is phased via **RM-007** (PLANNING.md);
  per the working agreement the phased plan is to be approved before coding, and the old
  monolithic page is retained until the new pages reach parity. **Open sub-decision:**
  reference library **global vs per-project** (global cache + per-project selection/provenance
  is the leading option) — to be resolved as **D-013** before the References-page increment.

### D-014 (2026-06-20) — ITS-based provisional taxon ID via remote NCBI BLAST
- **Decision:** The Assembly Manager's "Identify taxon by ITS" feeder (RM-007 step 3) extracts
  the ITS region with **ITSx**, then identifies the closest taxon by **`blastn -remote -db nt`**
  against NCBI (optionally restricted with `entrez_query="fungi[ORGN]"`), surfacing ranked
  candidate organisms for the user to pick; the chosen organism is written to the assembly with
  `taxon_source="its_blast"`. New module `taxon_id_utils.py` provides pure command-building,
  output parsing, and organism ranking, a RunManager-logged remote-BLAST call, and the
  ITSx→BLAST orchestration. Organism names come from the subject title (`stitle`), preferring
  `sscinames` only when a real local taxdb name is present.
- **Why:** The user chose **remote NCBI** (2026-06-20) for zero-setup, broad, cold-start ID —
  no local database is needed *before* references are fetched, which resolves the
  taxonomy↔reference chicken-and-egg. Verified live (blastn 2.x): an *Alternaria* ITS resolves
  to *Alternaria alternata* at 100 %. **Implementation finding:** `sscinames` is resolved from a
  *local* BLAST taxdb even under `-remote`, so without taxdb installed (the default) it returns
  `N/A` / `N/A;N/A`; deriving the organism from `stitle` avoids forcing a large taxdb download
  while staying accurate for nt barcode deflines. The step is provisional by design — the user
  reviews ranked hits and picks; we never auto-set the taxon.
- **Alternatives considered:** (a) Local DB the user configures (UNITE) — faster/offline/
  reproducible but needs DB setup; deferred. (b) Reuse the project's fetched ITS refs — not a
  cold-start path (refs are fetched later). (c) Hybrid local-then-remote — more code; deferred.
  (d) Resolve `staxids`→name via Entrez taxonomy — authoritative but adds an email-dependent
  extra call; `stitle` chosen as dependency-free and sufficient for a reviewed suggestion.
- **Status:** active. Implemented 2026-06-20; `taxon_id_utils.py` + 15 network-free tests
  (`tests/test_taxon_id_utils.py`); live remote-BLAST verified. Full suite 200 passing.

### D-013 (2026-06-20) — Reference libraries are per-project
- **Decision:** Per-locus reference libraries move from the global `~/.phylofetch/references/`
  to **per-project** storage under `<project>/references/<locus>/`, with provenance recorded
  in the project manifest (`workflow.steps.references`). The `ncbi_utils` reference functions
  gain a `ref_dir` parameter (defaulting to the global path for backward compatibility), and a
  `project_ref_dir(project_dir)` helper resolves `<project>/references`. The new References
  component page and the combine step use the project-scoped directory. (Resolves the open
  sub-decision flagged in D-012.)
- **Why:** Chosen by the user (2026-06-20) for self-contained, reproducible, taxonomy-tailored
  projects, consistent with the D-012 manifest backbone and RM-005 headless re-runs. The global
  library predated the per-project manifest and kept no record of which references a given
  analysis used. Re-fetching the same accession across projects (the cost) is acceptable for a
  reproducible research workflow.
- **Alternatives considered:** (a) Hybrid global download-cache + per-project selection — more
  moving parts; deferred and revisitable if duplicate downloads become a pain. (b) Keep global
  — zero refactor but weakest reproducibility and not taxonomy-tailored; rejected.
- **Status:** active. Implemented as RM-007 step 4 — **4a** (this commit): `ref_dir` threading
  through `ncbi_utils` + `project_ref_dir` + tests, backward-compatible (global default
  unchanged); **4b**: the NCBI References component page consumes the project-scoped dir.

### D-015 (2026-06-20) — Extraction outputs are per-project
- **Decision:** Extraction outputs from the **new component pages** (per-strain + combined
  locus FASTAs, logs) go under the project at **`<project>/results/loci/{per_strain,combined}`**,
  not the shared global `output_base`. Alignment Prep already reads a user-chosen directory
  (`rglob *_combined.fasta`), so it stays compatible by pointing at the project results dir
  (its default will be repointed there when the monolith is retired). The legacy monolith keeps
  using `output_base` until retirement.
- **Why:** Consistent with per-project references (D-013) and the manifest model (D-012), and it
  removes a real reproducibility hazard: `output_base` is a single config path shared by every
  project, so two projects' extraction outputs would collide/overwrite. The project already has
  a `results/` subdir (`init_project`), so no new structure is needed.
- **Alternatives considered:** Keep the global `output_base` — zero downstream change but
  collision-prone and inconsistent with per-project refs; rejected.
- **Status:** active. Applied to the ITSx component page (RM-007 step 4c) onward; the
  Exonerate/Primers pages follow the same convention; the monolith is unchanged until retired.

### D-016 (2026-06-20) — Decomposition complete: Workflow orchestrator + monolith retired
- **Decision:** Finished RM-007. (1) Added the **Workflow / Strategy orchestrator**
  (`pages/6_Workflow.py`): a named-strategy selector ("Fungal barcodes (ITSx + Exonerate)",
  "Primers only", "Everything") over a **manifest-driven checklist** that shows each step's live
  status (from `workflow.steps`) and links to its component page — the chain that ties the
  standalone pages together without hiding state in `st.session_state`. (2) **Ported the relaxed
  BLAST-amplicon strategy** into the Exonerate page as a mode toggle ("Exonerate (frame-safe)" |
  "BLAST amplicon (relaxed, genomic)", `require_complete_cds=False`) so no capability was lost
  on retirement (user choice). (3) **Retired the monolithic `pages/2_Loci_Extraction.py`** and
  renumbered to the final layout: 2 NCBI References · 3 ITSx (rDNA) · 4 Exonerate · 5 Primers ·
  6 Workflow · 7 Alignment Prep · 8 BUSCO · 9 Tree.
- **Why:** D-012 chose full decomposition; with all four extraction components (4b–4e) built and
  render-verified, the monolith was redundant. The user confirmed retire-now (2026-06-20) and
  chose to **port** the relaxed BLAST strategy rather than drop it. Deletion (not a stub) is clean
  because the monolith's logic lives entirely in `src/` (`extract_locus`,
  `extract_locus_exonerate`, `run_itsx`, `run_primer_extraction`, `merge_per_strain_outputs`) and
  the component pages, and git history preserves the file. This satisfies the working agreement's
  "document abandoned/removed code (what + why)" rule.
- **Alternatives considered:** (a) Keep a legacy stub page — rejected (clutter; components fully
  cover it). (b) Drop the relaxed BLAST strategy — rejected by the user (regression).
  (c) `st.page_link` without a fallback — hardened with a caption fallback so the orchestrator
  never crashes if the page registry is unavailable (e.g. headless AppTest bare mode).
- **Status:** active. **RM-007 complete.** CLAUDE.md repository-layout + extraction-strategy
  sections updated. All 10 pages render-verified (AppTest); full suite 210 passing.

### LXD-003 (2026-06-20) — ITSx fails on genome assemblies (hmmscan 100 kb limit)
- **Fix:** `run_itsx` now **chunks contigs > 90 kb into overlapping 20 kb-overlap windows**
  before ITSx (`chunk_long_contigs`), de-duplicates identical output regions (overlap dups),
  and surfaces the hmmscan over-limit signature as a real error (rc=1) instead of a silent
  empty result. Contigs ≤ limit pass through unchanged.
- **Why:** ITSx runs HMMER `hmmscan`, which aborts on any target sequence > 100 kb
  ("Target sequence length > 100K, over comparison pipeline limit"). Genome-assembly contigs
  are routinely Mb-scale, so ITSx **aborted and returned no rDNA** (exit 0 + empty → read as
  "no regions detected"). This affected the old monolith too — a pre-existing limitation, not a
  decomposition regression. Chunking with overlap ≫ the rDNA cistron (≈5–9 kb; ITS < 1 kb)
  keeps every rDNA region fully intact in at least one chunk; dedup removes the overlap
  duplicates. Verified on a real assembly (NS26-3-C2): now yields ITS_full 481 bp + ITS1/ITS2/
  SSU/LSU in ~74 s. Found because the user reported ITSx retrieving nothing while primer scan
  found ITS (blastn has no length limit).
- **Status:** active. `chunk_long_contigs` + 4 tests (`TestChunkLongContigs`); full suite 214.
  Benefits the ITSx component page automatically (it calls `run_itsx`).

### D-017 (2026-06-20) — Protein references for coding loci + taxon fallback (fixes Exonerate stops)
- **Decision:** The NCBI References page gains a **Reference type** toggle — **Protein** (default
  for coding loci) vs **Nucleotide** — and a **taxon fallback** (`taxon_fallbacks`: exact taxon →
  genus when the species returns nothing). For coding loci, Protein fetches genome-annotated
  proteins (`search_ncbi_protein`, db=protein, `[Protein Name]`, **RefSeq/full-length preferred**)
  so Exonerate auto-runs **protein2genome**. rDNA loci stay nucleotide (rRNA has no protein, ITSx
  path). `ncbi_search_count`/`search_ncbi_protein` carry a `field` param; the preview shows the
  per-locus DB and which taxon level was used.
- **Why:** Diagnosed from the user's real data that Exonerate **internal stops in barcoding genes
  were a reference artifact, not a bad assembly**: the D-011 synonym search fetches *partial-cds
  barcodes* which for fungi are **genomic (intron-containing)**, and those self-translate with
  stops (CAL ref 10, RPB2 ref 14); used as `coding2genome` queries (which assume an intron-free
  CDS) they mis-place exon/intron boundaries → frameshifts → high-identity-with-many-stops.
  **BUSCO is "fantastic"** for these assemblies (the definitive quality check), and a control
  proved it: the *same* NS26-3-C2 assembly, ACT locus, went from **15 stops (nucleotide barcode
  ref) → 0 stops, 100% id, clean 5-exon CDS** with a full-length **protein** ref from
  *A. dauci* (a different species — cross-genus protein refs work, the canonical Exonerate use).
  Protein + protein2genome is intron-immune and frame-pinned. The genus fallback handles the
  user's target *Alternaria* aff. *eureka* (a novel species absent from NCBI → falls back to
  *Alternaria*; verified live: RPB2 4564 / ACT 612 protein hits at genus level).
- **Alternatives considered:** (a) Nucleotide-only + a "this ref translates with stops" guard —
  insufficient: warns but still can't produce a clean CDS from a genomic barcode; rejected as the
  fix (may add the guard later). (b) Require complete-cds nucleotide refs — scarce (the reason
  D-011 dropped that filter). (c) `[Gene Name]`/`[All Fields]` for protein — `[Protein Name]` is
  more precise for "proteins named X" and matches the descriptive synonyms.
- **Status:** active. `taxon_fallbacks` + `search_ncbi_protein` field param + References-page
  changes; `TestTaxonFallbacks` (+4) → 218 tests. Render + live genus-fallback/protein-count
  verified. Nucleotide references remain right for the relaxed-BLAST / primer / GenBank-barcode
  paths; protein is the default only for the Exonerate (frame-safe CDS) coding path.

### D-018 (2026-06-20) — In-app project/cache management + mixed-reference-type guard
- **Decision:** (1) Add a **Manage Data** tab to Project Setup to inspect and clear per-project
  caches (references / results / run-logs / workflow state), clear the assembly registry, clear
  the legacy global reference cache, and **delete a project** (guarded). (2) Add a guard on the
  References page: when the reference type being fetched (protein/nucleotide) differs from what a
  locus already holds, the locus is **skipped with a warning** (shown in the preview and at
  fetch) rather than mixing types in one `*_refs.fasta`.
- **Why:** The user is starting projects from scratch (new assemblies replacing old) and the
  protein-reference fix (D-017) requires clearing the old *nucleotide* coding refs first — there
  was no in-app way to clear caches (only create/open). Without the guard, fetching protein into
  a locus that already holds nucleotide refs silently produces a mixed file, and
  `detect_fasta_type` (which samples the first sequence) would then pick the wrong Exonerate
  model — a quiet corruption. Skip+warn (vs hard block) was the user's call ("add the warn").
- **Alternatives considered:** (a) Hard-block the mismatched fetch — heavier; the user chose a
  warning. (b) Auto-separate protein/nucleotide into sub-files per locus — larger change;
  deferred. (c) No project deletion in-app — rejected; needed for "start fresh".
- **Status:** active. `project_manager` helpers (`project_data_summary`, `clear_project_data`,
  `reset_workflow`, `delete_project`, `clear_global_reference_cache` — all guarded; delete refuses
  non-projects and protected dirs) + 8 tests; Project Setup "Manage Data" tab; References
  skip+warn on type mismatch. Full suite 226. Verified the guard fires on the real project
  (coding loci are nucleotide → would skip a protein fetch).

### D-019 (2026-06-21) — Escalating edit-distance for in-silico PCR
- **Decision:** The primer search escalates the **per-primer edit-distance threshold**
  (strict → loose) up to a cap, rather than using a single fixed threshold.
  `find_primer_amplicons_escalating(start_mismatches, max_mismatches)` tries each threshold from
  start to cap and returns the first that yields an amplicon plus the threshold used;
  `run_primer_extraction` gains `escalate_to`; the Primers page slider is the cap (default 3,
  escalates from 2) and the preview reports the matched edit distance. Two existing behaviours
  were confirmed (user questions) and left as-is: the search checks **both strands/orientations**
  (fwd+rev required on opposite strands pointing inward; amplicon on either strand) and
  **expands IUPAC degenerate primers** into all concrete oligos before BLAST (D-009).
- **Why:** A congener's primers in a divergent target (the user's novel *Alternaria* aff.
  *eureka*) can need >2 mismatches to bind; the fixed default of 2 missed them and forced manual
  slider-bumping. Strict-first escalation recovers divergent sites automatically while not
  surfacing loose off-targets when a clean match exists.
- **Alternatives considered:** (a) Single permissive search at the cap — surfaces loose
  off-targets even when a strict match exists; escalation is strict-first. (b) One BLAST at the
  cap + re-filter per threshold (more efficient for heavily-degenerate primers) — deferred; the
  simple loop re-runs BLAST per step but only escalates when early thresholds fail (cheap in the
  common case).
- **Status:** active. `find_primer_amplicons_escalating` + `escalate_to` + page wiring;
  `TestEscalatingSearch` (+4) → 230 tests.

### D-020 (2026-06-21) — Standardized, extraction-first "setup" phase
- **Decision:** The front of the workflow is **extraction-first and standardized**: import
  assemblies → extract a *standard gene set* from every assembly (rDNA via ITSx; conserved
  coding markers via protein-guided Exonerate; anonymous markers via primers), *independent of
  whether comparison references exist* — these are the user's own data and GenBank submissions.
  Only afterward does a dedicated **Reference Taxa / Tips** stage bring in comparison sequences
  and assess per-locus availability for tree-building. Key design choices:
  (1) **Bundled protein guides** ship with the package as a *universal core* (≥1 Ascomycota +
  ≥1 Basidiomycota full-length protein per conserved marker — these markers are conserved
  kingdom-wide, so `protein2genome` extracts them anywhere) **plus swappable lineage packs** for
  clade-specific markers (Alt a1, OPA10-2, … are *not* universal). (2) **Extraction guides
  (protein) are stored and conceived separately from tree tips (nucleotide amplicons)** — never
  mixed in one locus file (cf. D-018). (3) The tips importer **auto-classifies pasted accessions
  to loci** via the synonym catalogue (D-011). (4) Tree-locus selection is its **own page**,
  separate from tips import. (5) **GenBank submission export** is a final, optional step.
  (6) Everything is tracked in the **project manifest** so a re-opened project is self-describing.
- **Why:** The user (and the package's broader users) want a reproducible, genus-agnostic setup
  that yields a consistent marker set regardless of reference availability, and a clear record of
  what's been done across sessions. Decoupling "what I extracted" from "what I can compare
  against" matches real practice (missing reference data per taxon/locus is normal) and makes the
  extractions standalone-valuable (GenBank deposits). Bundled universal guides remove per-project
  guide-fetching; lineage packs preserve flexibility beyond Alternaria. Backed by published
  precedent (UnFATE, Syst. Biol. 2025; Feau et al. 2011) — see `docs/method_references.md`.
- **Alternatives considered:** (a) Fetch guides per project per taxon — rejected (repeated work,
  no standard set). (b) One universal guide per locus — kept option open, but ≥1 Asco + ≥1 Basidio
  improves robustness on divergent targets. (c) Tips + guides in one reference store — rejected
  (protein/nucleotide mix breaks model selection, D-018). (d) Reference-availability-first (fetch
  refs before extracting) — rejected in favour of extraction-first.
- **Status:** active. Roadmap = **RM-008**; building component 1 (bundled protein guides) first.
