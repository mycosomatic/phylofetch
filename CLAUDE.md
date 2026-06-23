# phylofetch — Claude Code Notes

## Working agreement (READ FIRST, EVERY SESSION)

This project is used for **active scientific research**. Everything must be rock-solid,
transparent, traceable, and repeatable. These rules override default behavior.

**At the start of every session, before making any change:**
1. Read `PLANNING.md` and `DECISIONS.md` to confirm we are not backtracking on or
   contradicting a prior decision. If the requested work conflicts with a recorded
   decision, stop and flag it before proceeding.

**How we work:**
- **Two branches only.** `main` is the release branch; `dev` is where all work happens.
  Never create per-session or feature branches unless the user explicitly asks. Open PRs
  from `dev` into `main`.
- **Ask often.** The user has basic Python knowledge but intermediate–advanced mycology /
  fungal genomics knowledge. When intent is unclear, ask rather than assume. Work slowly
  and carefully.
- **Explain changes clearly**, including the bioinformatics reasoning, not just the code.
  Prefer over-explaining tools/concepts to assuming familiarity.
- **Scientific soundness is non-negotiable.** If a request appears to contradict
  well-established contemporary knowledge (mycology, genomics, phylogenetics, statistics),
  do not silently comply — explain the concern and make the user argue for it before
  proceeding.
- **Abandoned/unused code** must be documented in `DECISIONS.md` (what, why abandoned).
  Do not resurrect or delete it without asking first.

**Living documents (append-only; never delete past entries):**
- `PLANNING.md` — overarching roadmap and goals.
- `DECISIONS.md` — decision history with rationale and alternatives considered.
- `CHANGELOG.md` — what actually changed, session by session.
- When we reverse a past decision, **comment out** the old entry (`<!-- -->`) and add a
  dated pointer to the new decision that supersedes it. Keep the full history visible.

## Project overview

phylofetch is a Streamlit-based bioinformatics app for fungal genome assembly processing:
**assembly → loci extraction → alignment/trimming/concatenation → BUSCO phylogenomics → ML/Bayesian tree prep**

It is designed as a standalone pip-installable package that can also plug in as a dependency of a future umbrella Streamlit app.

## Repository layout

```
phylofetch/
├── src/phylofetch/        # importable package (src layout)
│   ├── alignment/         # MAFFT, trimAl, MACSE, concat wrappers
│   ├── assembly_utils.py  # N50, GC%, assembler detection
│   ├── blast_loci_utils.py
│   ├── busco_utils.py
│   ├── config.py          # ~/.phylofetch/config.json
│   ├── exonerate_utils.py # spliced protein/CDS→genome, frame-safe CDS (D-008)
│   ├── itsx_utils.py
│   ├── ncbi_utils.py
│   ├── primer_utils.py
│   ├── project_manager.py # RunManager, tool version probing, project manifest
│   ├── protein_guide_utils.py # bundled protein guides → protein2genome (D-020 / RM-008 c1)
│   ├── codon_prep_utils.py # tips → frame-consistent CDS + full-gene + protein (D-022 / RM-008 c2)
│   ├── taxon_id_utils.py  # ITS→remote BLAST provisional taxon ID (D-014)
│   └── tips_utils.py      # comparison-tip import + locus auto-classification (D-020 / RM-008 c3)
├── pages/                 # Streamlit multi-page app (component-page workflow, D-012)
│   ├── 0_Project_Setup.py
│   ├── 1_Assembly_Manager.py    # + per-assembly taxonomy + ITS→BLAST provisional ID (D-014)
│   ├── 2_NCBI_References.py     # optional taxon-closer guides, coding-only, RefSeq candidate picker (D-024, D-023)
│   ├── 3_ITSx_rDNA.py          # rDNA extraction (per-project outputs, D-015)
│   ├── 4_Exonerate.py          # coding loci: Exonerate frame-safe | relaxed BLAST amplicon | gene-of-interest
│   ├── 5_Primers.py            # in-silico PCR (degenerate-aware D-009, edit-distance escalation D-019)
│   ├── 6_Workflow.py           # strategy orchestrator: manifest-driven checklist (D-012)
│   ├── 7_Reference_Taxa.py     # tree tips: paste accessions → auto-classify to locus (D-020)
│   ├── 8_Codon_Tip_Prep.py     # frame coding tips → codon-ready CDS+gene+protein matrices (D-022)
│   ├── 9_Alignment_Prep.py
│   ├── 10_BUSCO_Phylogenomics.py
│   └── 11_Tree_Visualization.py
│   # NB: the old monolithic 2_Loci_Extraction.py was retired 2026-06-20 (D-016);
│   # its logic lives in src/ (extract_locus*, run_itsx, run_primer_extraction).
├── tests/
│   ├── conftest.py
│   ├── test_assembly_utils.py
│   ├── test_blast_loci_utils.py  # MUST cover LXD-002 regression
│   ├── test_itsx_utils.py        # MUST cover LXD-001 fixes
│   ├── test_exonerate_utils.py   # parser/QC/build offline + binary-guarded E2E (D-008)
│   ├── test_codon_prep_utils.py  # soft-mask/merge offline + binary-guarded E2E (D-022)
│   └── fixtures/                 # verbatim exonerate 2.4.0 output + synthetic gene
├── app.py                 # streamlit run app.py
├── pyproject.toml
└── environment.yml
```

## Development commands

```bash
# Install in editable mode
pip install -e .

# Launch the app
phylofetch
# or
streamlit run app.py

# Run tests
pytest tests/ -v

# Run tests with coverage
pytest tests/ -v --cov=phylofetch --cov-report=term-missing
```

## Critical bugs (already fixed)

### LXD-001 — ITSx wrapper (itsx_utils.py)
- `--multi_out T` is not a valid ITSx flag → removed
- Stale `<prefix>.*` files must be deleted before rerun
- `run_itsx()` returns `(returncode, log_text, found_dict)` — log always populated

### LXD-002 — BLAST HSP grouping (blast_loci_utils.py)
- HSPs from different reference accessions (`qseqid`) must NOT be stitched
- `select_best_locus_group()` groups by `qseqid` first, picks highest-bitscore reference
- Returns `(hsps, ref_accession)` tuple (not just list)

### LXD-003 — ITSx on genome assemblies (itsx_utils.py)
- HMMER `hmmscan` aborts on sequences > 100 kb; genome contigs are Mb-scale, so ITSx returned
  no rDNA (exit 0 + empty). `run_itsx` now `chunk_long_contigs()` (> 90 kb → overlapping 20 kb
  windows) before ITSx, de-dupes identical output regions, and returns rc=1 on the over-limit
  signature. Overlap ≫ rDNA cistron, so no region is split. (2026-06-20)

## Extraction strategies (now component pages, D-012/RM-007)

As of 2026-06-20 the single "Run Extraction tab" was decomposed into standalone **component
pages** chained by the **project manifest** (`metadata/project_manifest.json` → `workflow.steps`;
see `project_manager`). Each extraction strategy is its own page; the **Workflow** page
(`6_Workflow.py`) picks a named strategy and links the steps with live status. References are
per-project (`<project>/references`, D-013) and extraction outputs are per-project
(`<project>/results/loci/{per_strain,combined}`, D-015).

1. **rDNA — ITSx** (`3_ITSx_rDNA.py`, `itsx_utils.py`). Extract ITS/ITS1/ITS2/LSU/SSU from
   assemblies; no Exonerate. (`place_rdna_regions` / `combine_rdna_regions`.)
2. **Coding loci — Exonerate (frame-safe)** (`4_Exonerate.py`, `exonerate_utils.py`).
   tblastn/blastn narrows to the best contig, then Exonerate spliced alignment
   (`protein2genome`/`coding2genome`, auto-picked) yields a translatable, frame-checked CDS.
   A **relaxed BLAST amplicon** mode (`extract_locus`, `require_complete_cds=False`, genomic
   amplicon, no frame guarantee — the former "PCR amplicon refs (relaxed)" strategy) is a
   selectable toggle; the BLAST HSP path is also the automatic fallback when `exonerate` is not
   on PATH (with a frame-safety warning). Also extracts an arbitrary **gene of interest**
   (paste/upload ortholog). (D-008, addresses RM-002.) **Reference source** (D-020/D-023):
   *Bundled guides* (universal Asco/Basidio protein core, default, no fetching) · *Bundled +
   project library (taxon-closer)* (bundled guides plus the project's fetched protein refs in one
   `protein2genome` query, best model wins) · *Project reference library* (only the fetched refs).
3. **PCR Primers (in-silico PCR)** (`5_Primers.py`, `primer_utils.py`). Locate fwd+rev primer
   binding sites with `blastn-short` (IUPAC degenerate primers expanded to concrete oligos
   first, D-009), pair on the same contig opposite strands, extract the amplicon between. No
   NCBI reference library needed.

### Exonerate library (`exonerate_utils.py`) — D-008
- **Hybrid pipeline**: `extract_locus_exonerate` narrows via `select_best_locus_group` to the single best contig (coords stay in contig space → no offset math), then runs Exonerate on that contig; `narrow=False` / no BLAST hit ⇒ whole-assembly run.
- **Models**: `protein2genome` (protein query) / `coding2genome` (nucleotide CDS query), via `MODEL_FOR_QUERYTYPE` keyed on `detect_fasta_type`.
- **Parsing**: `parse_exonerate_gff` reads verbatim `--showtargetgff yes` GFF (one result per `START…END OF GFF DUMP` block) + a `--ryo` line whose `%tcs` is the authoritative spliced CDS (cross-checked against the coord-rebuilt CDS). Exonerate emits introns + splice5/splice3 explicitly (no gap inference). `--bestn N` surfaces paralogs (`select_best_model` returns `(best, others)`).
- **QC (RM-002)**: `validate_cds` reports reading-frame (`len % 3`), internal-stop count, translation; the per-locus log records a PASS/REVIEW verdict and GT-AG tally. Default = write-and-flag; `strict_qc=True` rejects internal-stop / frameshift CDS (D-007: don't silently drop partial/type seqs).
- **Outputs**: `LOCUS_CDS/protein/genomic/introns.fasta`, `LOCUS.gff3`, `LOCUS_partition.nex`, `LOCUS_extraction.log` (GFF3 + partition writers reused from `blast_loci_utils`). `LOCUS_genomic.fasta` is **soft-masked** — exons UPPERCASE / introns lowercase via `soft_mask_genomic` (D-022), the same annotation the Codon Tip Prep page applies to reference tips.
- **Gene of interest**: the Run-Extraction page also accepts an arbitrary pasted/uploaded ortholog (protein or CDS) → returns just the CDS + exon model, sharing the same logged code path as catalogue loci.
- **Provenance**: Exonerate runs route through RunManager (`run_exonerate`, `tool_version_keys=["exonerate"]`). Install via `conda install -c bioconda exonerate` (also pinned in `environment.yml`).

### Primer library (`primer_utils.py` + `data/primers.json`)
- **Citable, packaged catalogue**: `src/phylofetch/data/primers.json` ships via `[tool.setuptools.package-data]`. Every pair carries `source`, `citation`, `reference_url` (primary literature; rDNA cross-checked vs UNITE). 16 pairs across ITS/SSU/LSU/TEF1/RPB1/RPB2/TUB2/ACT/CAL/HIS3/GAPDH (RPB1 `RPB1-Af/RPB1-Cr` + RPB2 `fRPB2-5F/bRPB2-11R1` are Matheny-lab pairs, D-010).
- **Degenerate (IUPAC) base handling (D-009)**: BLAST does *not* resolve IUPAC codes — a degenerate primer neither seeds (`blastn-short` needs a 7-bp exact word) nor scores compatibly (`Y` over `C` counts as a mismatch like an incompatible `R`). `find_primer_amplicons` therefore **expands** each primer into all concrete oligos (`expand_degenerate_primer`, `MAX_PRIMER_EXPANSION=8192`), searches each as `FWD_*`/`REV_*`, and collapses variant duplicates to one candidate per amplicon. Concrete primers → 1 variant (same code path). Over-cap → `ValueError` surfaced as a "too degenerate" status, never a silent empty result.
- **User library**: custom pairs persist in `~/.phylofetch/primers.json`; merged on top of built-in via `get_primer_catalogue()` (user wins on name clash). `save_user_primer` / `load_user_primers` / `delete_user_primer`.
- **In-app entry**: PCR Primer mode offers all loci (no reference gate, unlike BLAST/Exonerate), a per-locus "Custom…" pair entry, and an "➕ Custom locus" section for a locus *not* in the catalogue (custom loci merge into `primer_assignments`, so preview/run/combine treat them like any other).
- **Edit-distance matching**: `max_mismatches` = substitutions + unaligned primer bases (`_effective_mismatch`), so partial/truncated primer alignments can't slip through.
- **Disambiguation**: `find_primer_amplicons()` returns all candidate sites sorted by edit distance; the UI "Preview & choose binding sites" lets the user override the auto-pick (off-target handling).
- **Logged + provenanced**: primer search routes through RunManager; `LOCUS_extraction.log` records primer, citation, command, and chosen site.
- **PRM-001 (data fix)**: built-in `ACT-512F` was a corrupted 35-nt sequence → corrected to canonical `ATGTGCAAGGCCGGTTTCGC` (Carbone & Kohn 1999). Regression-tested.

## Key design decisions

- **Config path**: `~/.phylofetch/config.json` (not `~/.alternaria_toolkit`)
- **Genus-agnostic**: no Alternaria-specific defaults anywhere
- **Rich FASTA headers**: NCBI bracket-style `[key=value]` for all outputs
- **RunManager**: every external tool call logged with command + tool versions + timestamps
- **Per-locus extraction logs**: `LOCUS_extraction.log` alongside FASTAs
- **Frame-safe coding CDS**: Exonerate spliced alignment is the preferred coding-locus path (D-008); `extract_from_hsps` HSP-as-exon is a *documented fallback* for coding CDS (still primary for relaxed PCR-amplicon), used only when `exonerate` is absent, with a UI frame-safety warning
- **Alignment viewer**: BLAST `-outfmt 0` pairwise text in `st.code()` (no extra deps)
- **Exonerate / ITSx / MACSE**: optional external tools; graceful fallback if not on PATH
- **src layout**: `sys.path.insert(0, .../src)` in pages for dev mode; package importable after `pip install -e .`

## Plug-in architecture

The umbrella app can depend on phylofetch:
```
phylofetch @ git+https://github.com/mycosomatic/phylofetch
```
Then import directly:
```python
from phylofetch.blast_loci_utils import extract_locus
from phylofetch.project_manager import RunManager
```

## Branch

Active development branch: `dev`

All Claude Code sessions should develop on `dev` and open PRs into `main` when ready to merge. Do not create per-session branches.
