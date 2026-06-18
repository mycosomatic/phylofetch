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
│   ├── itsx_utils.py
│   ├── ncbi_utils.py
│   └── project_manager.py # RunManager, tool version probing
├── pages/                 # Streamlit multi-page app
│   ├── 0_Project_Setup.py
│   ├── 1_Assembly_Manager.py
│   ├── 2_Loci_Extraction.py
│   ├── 3_Alignment_Prep.py
│   ├── 4_BUSCO_Phylogenomics.py
│   └── 5_Tree_Visualization.py
├── tests/
│   ├── conftest.py
│   ├── test_assembly_utils.py
│   ├── test_blast_loci_utils.py  # MUST cover LXD-002 regression
│   └── test_itsx_utils.py        # MUST cover LXD-001 fixes
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

## Extraction strategies (Run Extraction tab)

Three selectable strategies drive loci extraction:

1. **BLAST – PCR amplicon refs (relaxed)** — NCBI amplicon refs, skips the CDS completeness gate (`require_complete_cds=False`).
2. **BLAST – CDS / protein (strict)** — curated CDS/protein refs; enforces `min_cds_pct_of_ref`. ITSx handles rDNA.
3. **PCR Primers (in-silico PCR)** — `primer_utils.py`. Locate fwd+rev primer binding sites with `blastn-short`, pair them on the same contig with opposite strands, extract the amplicon between. No NCBI reference library needed — useful when accessions are missing or map poorly.

### Primer library (`primer_utils.py` + `data/primers.json`)
- **Citable, packaged catalogue**: `src/phylofetch/data/primers.json` ships via `[tool.setuptools.package-data]`. Every pair carries `source`, `citation`, `reference_url` (primary literature; rDNA cross-checked vs UNITE). 14 pairs across ITS/SSU/LSU/TEF1/RPB2/TUB2/ACT/CAL/HIS3/GAPDH.
- **User library**: custom pairs persist in `~/.phylofetch/primers.json`; merged on top of built-in via `get_primer_catalogue()` (user wins on name clash). `save_user_primer` / `load_user_primers` / `delete_user_primer`.
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
- **Alignment viewer**: BLAST `-outfmt 0` pairwise text in `st.code()` (no extra deps)
- **MACSE**: optional Java JAR; graceful fallback if not configured
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
