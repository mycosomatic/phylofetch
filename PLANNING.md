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

- _(2026-06-18)_ No roadmap items recorded yet. First entry to be added when we scope our
  next piece of work together.

## Done (high-level milestones)

See `CHANGELOG.md` for the detailed, dated record. Recent shipped work (from git history):

- BUSCO integration and assembly-stats display; primer mode surfaced in UI.
- Import scanners consolidated; assembly schema inconsistencies fixed.
- PCR primer-based locus extraction strategy added.
