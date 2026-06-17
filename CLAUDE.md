# phylofetch — Claude Code Notes

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
