# phylofetch

Fungal genome assembly → phylogenetic loci extraction → alignment/concatenation → BUSCO phylogenomics → ML tree prep.

A Streamlit app for reproducible, logging-first phylogenetics workflows. Accepts assemblies from any assembler (EGAP, SPAdes, Flye, Hifiasm, etc.).

## Features

- **Assembly Manager** — register genome assemblies; compute N50, GC%, contig stats
- **Loci Extraction** — extract rDNA (ITS/LSU/SSU via ITSx) and protein-coding loci (blastn/tblastn); CDS, genomic, intron FASTA + GFF3 + codon-position partitions
- **Alignment Prep** — MAFFT alignment, trimAl trimming, optional MACSE codon-aware alignment, supermatrix concatenation with merged partition nexus
- **BUSCO Phylogenomics** — import BUSCO v4/v5 or Compleasm results; occupancy matrix; single-copy FASTA export
- **Tree Visualization** — IQ-TREE2 integration; ML tree rendering; Bayesian posterior overlay (in development)

All steps are logged with full command, tool versions, timestamps, and environment snapshots for reproducibility.

## Installation

```bash
# Recommended: conda environment
conda env create -f environment.yml
conda activate phylofetch

# Install package in editable mode
pip install -e .
```

**Required tools** (via conda/bioconda):
```bash
conda install -c bioconda itsx blast mafft trimal iqtree
```

**Optional:**
- MACSE (Java JAR): download from https://bioweb.supagro.inra.fr/macse/

## Usage

```bash
phylofetch
# or
streamlit run app.py
```

## Workflow

1. **Project Setup** — initialize workspace, register assemblies
2. **Loci Extraction** — download reference sequences from NCBI, run ITSx and/or BLAST extraction
3. **Alignment Prep** — align per-locus FASTAs, trim, concatenate into supermatrix
4. **BUSCO Phylogenomics** — import BUSCO runs, filter single-copy orthologs, export for supermatrix
5. **Tree Visualization** — run IQ-TREE2, inspect ML tree

## Reproducibility

Every external tool call produces a run folder in `~/.phylofetch/projects/default/runs/` containing:
- `command.json` — exact command, working directory, environment variables
- `stdout.log` / `stderr.log` — captured output
- `terminal.log` — combined timestamped output
- `environment.json` — Python version, conda env, tool versions

FASTA headers include sample, contig, coordinates, tool, tool version, and extraction timestamp.

## Plug-in architecture

phylofetch is designed to be used as a pip-installable dependency by an umbrella Streamlit app:

```python
from phylofetch.blast_loci_utils import extract_locus
from phylofetch.project_manager import RunManager
from phylofetch.busco_utils import build_occupancy_matrix
```

## License

MIT
