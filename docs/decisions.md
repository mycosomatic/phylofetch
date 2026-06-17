# phylofetch — Architecture Decisions

## D01 — Single-process streamlit app
**Decision**: Run as a single Streamlit process, no background job queue.  
**Rationale**: Keeps deployment simple; tool runs are fast enough (minutes, not hours) for interactive use. Long-running IQ-TREE jobs show a spinner with live log tailing.

## D02 — src layout for installable package
**Decision**: `src/phylofetch/` with `pyproject.toml` entry point.  
**Rationale**: Clean separation between package code and app code. Enables `pip install -e .` for development and `pip install` as a dependency in the umbrella app.

## D03 — Config at `~/.phylofetch/config.json`
**Decision**: JSON config in user home; no environment variables for tool paths.  
**Rationale**: Reproducible, self-contained, inspectable. Tool paths are set once via the UI and persisted.

## D04 — RunManager for all external tool calls
**Decision**: Every subprocess call goes through `RunManager.run()`, which stores command, stdout, stderr, tool versions, and environment in a per-run folder.  
**Rationale**: Full reproducibility. Every result can be traced back to the exact command and software versions that produced it.

## D05 — NCBI bracket-style FASTA headers
**Decision**: Rich provenance embedded in FASTA description field: `[key=value]` format after the sequence ID.  
**Rationale**: Parseable by downstream tools (ID before whitespace is unchanged), self-documenting, survives file transfers. NCBI convention understood by the field.

## D06 — Per-locus extraction log alongside FASTAs
**Decision**: `LOCUS_extraction.log` written into the same directory as the output FASTAs.  
**Rationale**: The output directory is self-contained. Anyone can open a locus folder and know exactly what produced it without access to the central run database.

## D07 — Genus-agnostic design
**Decision**: No Alternaria-specific defaults, taxon IDs, or locus names hardcoded.  
**Rationale**: phylofetch is for fungal genomics broadly. Users specify their own reference sequences and organism names. The NCBI default is "Fungi".

## D08 — BLAST pairwise alignment viewer (no extra deps)
**Decision**: Loci Extraction shows extracted vs reference alignment using BLAST `-outfmt 0` output in `st.code()`.  
**Rationale**: Requires no additional Python packages (BLAST is already required). Familiar format for biologists. Interactive MSA viewer (pymsaviz) deferred to Alignment Prep page.

## D09 — MACSE optional with graceful fallback
**Decision**: MACSE is optional; users configure a JAR path. Pages skip MACSE sections if JAR is not configured.  
**Rationale**: MACSE cannot be installed via conda and requires manual JAR download. Codon-aware alignment is a power-user feature; standard MAFFT alignment works for most use cases.

## D10 — IQ-TREE2 in-app (page 5, stub then full)
**Decision**: Wrap IQ-TREE2 in the Tree Visualization page with RunManager integration.  
**Rationale**: Keeps the workflow in one tool. Users can run IQ-TREE2 without leaving phylofetch, and the command is logged. Bayesian posteriors imported post-hoc (BEAST standalone is fine).

## D11 — LXD-001 fix: ITSx command hygiene
**Decision**: Remove `--multi_out T` flag (invalid in all ITSx versions); clear stale output files before rerun; always return last 40 log lines.  
**Rationale**: `--multi_out` caused silent failures. Stale files caused false-positive success detection. Log lines surfaced on failure enable UI-level diagnostics without requiring users to find log files.

## D12 — LXD-002 fix: BLAST HSP grouping by reference
**Decision**: `select_best_locus_group()` groups HSPs by `qseqid` (reference accession) first, then picks the single best reference by total bitscore before selecting the best `(contig, strand)` group.  
**Rationale**: Without this fix, HSPs from different reference sequences on the same contig could be merged into one false CDS. The fix ensures only one reference contributes to each extraction. Regression test added.

## D13 — Plug-in architecture via pip
**Decision**: phylofetch is designed to be `pip install`-able as a dependency of the umbrella Fungal Toolkit.  
**Rationale**: Clean separation of concerns. The umbrella app imports from `phylofetch.*` without copying code. Config namespaces don't collide (`~/.phylofetch/` vs `~/.fungal_toolkit/`).

## D14 — BUSCO v4/v5 + Compleasm auto-detection
**Decision**: `scan_busco_run()` auto-detects BUSCO format by looking for `short_summary*.txt` (BUSCO) vs `summary.txt` (Compleasm).  
**Rationale**: Both tools are in active use; supporting both without requiring the user to specify the format reduces friction.

## D15 — Supermatrix missing-taxon padding
**Decision**: Taxa absent from a locus get gap characters (`-` by default) in the concatenated supermatrix.  
**Rationale**: Standard practice for multi-locus supermatrix phylogenetics. IQ-TREE2 handles missing data; the partition file is still correct. User can choose `N` or `?` as alternatives.
