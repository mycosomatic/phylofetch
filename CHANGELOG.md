# phylofetch — Changelog

> **Append-only.** Do not delete past entries. Newest at the top. This is the "what actually
> changed" record. Rationale lives in `DECISIONS.md`; roadmap in `PLANNING.md`.

## 2026-06-18

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
