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
