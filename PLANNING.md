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
- **RM-003 (2026-06-18) — Species-delimitation workflow.** Treat per-locus alignments as
  first-class: per-locus gene trees + a concordance/conflict view (GCPSR), instead of
  defaulting to the concatenated supermatrix for the novelty goal.
- **RM-004 (2026-06-18) — Phylogenomics breadth.** Add thin ASTRAL (coalescent species
  tree, robust to ILS) and ANI (fastANI/pyANI genome distance) runners alongside the
  existing IQ-TREE2 concatenation path. (See D-005.)
- **RM-005 (2026-06-18) — Headless reproducibility.** A config-driven `phylofetch run
  project.yaml` over the existing functions, input-file hashing, and pinned tool versions in
  `environment.yml`, so the whole pipeline is re-runnable off the browser session.

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
- **(method) Nucleotide BUSCO default saturates at depth;** consider amino-acid / codon-aware
  for deeper comparisons. Concatenation-only ignores ILS (high near speciation). → RM-004.
- **(coverage) No tests for `concat.py`, `busco_utils.py`, `ncbi_utils.py`** — the
  partition-offset math and occupancy filtering are silent-error-prone.

## Done (high-level milestones)

See `CHANGELOG.md` for the detailed, dated record. Recent shipped work (from git history):

- BUSCO integration and assembly-stats display; primer mode surfaced in UI.
- Import scanners consolidated; assembly schema inconsistencies fixed.
- PCR primer-based locus extraction strategy added.
