# phylofetch — Planning

## Vision

A focused, installable Streamlit app that takes a fungal genome assembly and produces analysis-ready phylogenetic datasets:
- Extracted loci (rDNA + protein-coding markers) with full provenance
- Aligned, trimmed, concatenated supermatrices with partition files
- BUSCO single-copy ortholog matrices for phylogenomics
- IQ-TREE2 integration with Bayesian posterior overlay (publication trees)

phylofetch plugs into a future umbrella app (Fungal Genomics Toolkit) as a pip dependency.

---

## Milestones

### M1 — Core infrastructure ✅
- [x] `src/phylofetch/` package with `pyproject.toml`
- [x] `RunManager` with tool version capture and environment snapshot
- [x] `config.py` at `~/.phylofetch/config.json`
- [x] `assembly_utils.py` — N50, GC%, assembler detection
- [x] `ncbi_utils.py` — genus-agnostic NCBI search/fetch
- [x] `blast_loci_utils.py` — HSP-as-exon extraction + LXD-002 fix
- [x] `itsx_utils.py` — ITSx wrapper + LXD-001 fix
- [x] Rich FASTA headers (NCBI bracket-style)
- [x] Per-locus extraction logs

### M2 — Streamlit pages ✅
- [x] `app.py` — landing page + tool settings
- [x] `pages/0_Project_Setup.py` — generic FASTA import, tool check, command history
- [x] `pages/1_Assembly_Manager.py` — assembly stats viewer
- [x] `pages/2_Loci_Extraction.py` — ITSx + BLAST extraction, alignment viewer
- [x] `pages/3_Alignment_Prep.py` — MAFFT, trimAl, MACSE, concatenation
- [x] `pages/4_BUSCO_Phylogenomics.py` — occupancy matrix, SC FASTA export
- [x] `pages/5_Tree_Visualization.py` — IQ-TREE2 runner + tree view stub

### M3 — Alignment module ✅
- [x] `alignment/mafft.py`
- [x] `alignment/trimal.py`
- [x] `alignment/macse.py` (optional)
- [x] `alignment/concat.py` — supermatrix + merged partition nexus
- [x] `busco_utils.py` — occupancy matrix, SC FASTA export

### M4 — Tests ✅
- [x] `tests/conftest.py` — shared fixtures
- [x] `tests/test_assembly_utils.py`
- [x] `tests/test_blast_loci_utils.py` — LXD-002 regression
- [x] `tests/test_itsx_utils.py` — LXD-001 regression

### M5 — Documentation ✅
- [x] `README.md`
- [x] `CLAUDE.md`
- [x] `docs/decisions.md`
- [x] `docs/planning.md`
- [x] `environment.yml`

---

## Next (M6+)

### M6 — Tree visualization (page 5, full)
- [ ] toytree or ETE3 interactive rendering
- [ ] Bayesian posterior overlay (BEAST2 / MrBayes import)
- [ ] Publication SVG/PDF export
- [ ] Node annotation: bootstrap + posterior dual-support

### M7 — Alignment viewer upgrades
- [ ] pymsaviz MSA preview in Alignment Prep
- [ ] Option B: Plotly colored alignment grid for extracted locus vs reference

### M8 — Umbrella app integration
- [ ] Publish to PyPI or GitHub releases
- [ ] Umbrella app `pyproject.toml` adds `phylofetch` as dependency
- [ ] Shared RunManager project workspace across apps

---

## Known issues / backlog

- UI-001: Verify no Streamlit deprecation warnings at runtime (use_container_width, etc.)
- MACSE licensing: guide users to JAR download in UI
- Compleasm: verify `summary.txt` parsing covers all field names
- NCBI rate limits: add optional API key field to Project Setup
- Large assemblies: ITSx memory usage on scaffolded genomes >100Mb
