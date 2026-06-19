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
