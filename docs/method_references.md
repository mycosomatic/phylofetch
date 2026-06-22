# phylofetch — methods references

Literature backing the workflow (genome assembly → marker extraction → codon-aware alignment
→ multilocus / phylogenomic tree), for defending methodological choices in write-ups.

**Verification key:**
- ✅ **verified** = citation details fetched from the source during research (2026-06-21).
- ⚠️ **confirm** = correct paper, but exact volume/pages/authors/DOI not re-fetched — verify
  before formal citation.

---

## A. The workflow precedent (genome-mined markers → codon-aware alignment → phylogeny)

### ✅ UnFATE — closest published analogue to this whole workflow
Ametrano CG, Jensen J, Lumbsch HT, Grewe F (2025). *UnFATE: A Comprehensive Probe Set and
Bioinformatics Pipeline for Phylogeny Reconstruction and Multilocus Barcoding of Filamentous
Ascomycetes (Ascomycota, Pezizomycotina).* Systematic Biology 74(5):740–757.
- DOI: https://doi.org/10.1093/sysbio/syaf011
- Open copy: https://pmc.ncbi.nlm.nih.gov/articles/PMC12699997/
- Supports: extracting single-copy ortholog markers from **genome assemblies / WGS**;
  **MACSE2 codon-aware alignment**; IQ-TREE + ASTRAL (concatenation + coalescent);
  multilocus **barcoding vs a genome database, outperforming ITS** — for filamentous
  Ascomycetes (includes *Alternaria* / Dothideomycetes).

### ✅ Genome-mining markers for non-model fungi
Feau N, Decourcelle T, Husson C, Desprez-Loustau M-L, Dutech C (2011). *Finding single copy
genes out of sequenced genomes for multilocus phylogenetics in non-model fungi.* PLoS ONE
6(4):e18803.
- DOI: https://doi.org/10.1371/journal.pone.0018803
- Open copy: https://pmc.ncbi.nlm.nih.gov/articles/PMC3076447/
- Supports: mining single-copy protein-coding markers from fungal genomes for phylogenetics;
  genome-derived markers match/beat ITS resolution.

### ✅ BUSCO orthologs → genome-scale species trees (the phylogenomics path)
Tsai C-H, Stajich JE (2025). *Phyling: phylogenetic inference from annotated genomes.*
bioRxiv 2025.07.30.666921 (preprint).
- DOI: https://doi.org/10.1101/2025.07.30.666921
- Open copy: https://pmc.ncbi.nlm.nih.gov/articles/PMC12324385/
- Supports: BUSCO single-copy orthologs from genomes → concatenation / ASTRAL species trees.

---

## B. Tool / method references (for a Methods section)

### ✅ MACSE (codon-aware alignment of coding sequences)
Ranwez V, Harispe S, Delsuc F, Douzery EJP (2011). *MACSE: Multiple Alignment of Coding
SEquences Accounting for Frameshifts and Stop Codons.* PLoS ONE 6(9):e22594.
- DOI: https://doi.org/10.1371/journal.pone.0022594
- https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0022594

Ranwez V, Douzery EJP, Cambon C, Chantret N, Delsuc F (2018). *MACSE v2: Toolkit for the
Alignment of Coding Sequences Accounting for Frameshifts and Stop Codons.* Molecular Biology
and Evolution 35(10):2582–2584.
- DOI: https://doi.org/10.1093/molbev/msy159
- Open copy: https://pmc.ncbi.nlm.nih.gov/articles/PMC6188553/

### ⚠️ Codon-structure-preserving alignment & by-hand curation (Codon Tip Prep companions, D-022)
phylofetch's Codon Tip Prep page (D-022) frames comparison tips into intron-free, frame-pinned
CDS so they line up codon-by-codon with the extracted isolate loci. It runs **no aligner** and is
agnostic to how you align/curate afterwards. Codon structure must be preserved either by a
**codon-aware aligner** or **by hand** — a plain nucleotide aligner (e.g. MAFFT on CDS) does not
guarantee it. Options (all optional companions, none a phylofetch dependency):
- **MACSE v2** — codon-aware, models frameshifts/stops (cited above); best for messy ORFs.
- **Translate → align AA → back-translate** (gap pattern from the protein alignment stamped onto
  the *original* nucleotides, gaps forced onto codon boundaries; no nucleotides invented):
  - ⚠️ **PAL2NAL** — Suyama M, Torrents D, Bork P (2006). *PAL2NAL: robust conversion of protein
    sequence alignments into the corresponding codon alignments.* Nucleic Acids Res 34:W609–W612.
    DOI: https://doi.org/10.1093/nar/gkl315
  - ⚠️ **TranslatorX** — Abascal F, Zardoya R, Telford MJ (2010). *TranslatorX: multiple alignment
    of nucleotide sequences guided by amino acid translations.* Nucleic Acids Res 38:W7–W13.
    DOI: https://doi.org/10.1093/nar/gkq291
- **Interactive frame-aware editors** (Mesquite alternatives that show the live AA translation so
  reading frame is visible while editing):
  - ⚠️ **AliView** — Larsson A (2014). *AliView: a fast and lightweight alignment viewer and editor
    for large datasets.* Bioinformatics 30(22):3276–3278. DOI:
    https://doi.org/10.1093/bioinformatics/btu531
  - ⚠️ **SeaView** — Gouy M, Guindon S, Gascuel O (2010). *SeaView version 4.* Molecular Biology
    and Evolution 27(2):221–224. DOI: https://doi.org/10.1093/molbev/msp259
  - ⚠️ **Jalview** — Waterhouse AM, Procter JB, Martin DMA, Clamp M, Barton GJ (2009). *Jalview
    Version 2.* Bioinformatics 25(9):1189–1191. DOI: https://doi.org/10.1093/bioinformatics/btp033
- Boundary aids carried into the alignment with **no extra tooling**: the full-gene FASTA is
  written exons-UPPER / introns-lower (case is inert to the aligners/tree tools but visible in any
  editor and tracks the sequence as gaps move); the per-sequence GFF3 from extraction imports into
  Geneious as an exon/intron feature track. (NB: intron positions are not always homologous across
  taxa — intron gain/loss — so per-sequence case-masking is the honest representation, not a single
  reference ruler.)

### ✅ Exonerate (protein-to-genome / spliced alignment)
Slater GSC, Birney E (2005). *Automated generation of heuristics for biological sequence
comparison.* BMC Bioinformatics 6:31.
- DOI: https://doi.org/10.1186/1471-2105-6-31
- Open copy: https://pmc.ncbi.nlm.nih.gov/articles/PMC553969/

### ✅ ITSx (rDNA region extraction)
Bengtsson-Palme J, Ryberg M, Hartmann M, et al. (2013). *Improved software detection and
extraction of ITS1 and ITS2 from ribosomal ITS sequences of fungi and other eukaryotes for
analysis of environmental sequencing data.* Methods in Ecology and Evolution 4(10):914–919.
- DOI: https://doi.org/10.1111/2041-210X.12073
- https://besjournals.onlinelibrary.wiley.com/doi/10.1111/2041-210X.12073

### ✅ BUSCO (single-copy ortholog assessment / sets)
Simão FA, Waterhouse RM, Ioannidis P, Kriventseva EV, Zdobnov EM (2015). *BUSCO: assessing
genome assembly and annotation completeness with single-copy orthologs.* Bioinformatics
31(19):3210–3212.
- DOI: https://doi.org/10.1093/bioinformatics/btv351
- https://academic.oup.com/bioinformatics/article/31/19/3210/211866
- ⚠️ confirm — BUSCO v5: Manni M, Berkeley MR, Seppey M, Simão FA, Zdobnov EM (2021).
  Molecular Biology and Evolution 38(10):4647–4654. DOI: https://doi.org/10.1093/molbev/msab199

---

## C. Component precedent — protein-guided marker extraction (Exonerate) in phylogenetics

### ⚠️ HybPiper (canonical target-capture pipeline that uses Exonerate to extract markers)
Johnson MG, Gardner EM, Liu Y, Medina R, Goffinet B, Shaw AJ, Zerega NJC, Wickett NJ (2016).
*HybPiper: Extracting coding sequence and introns for phylogenetics from high-throughput
sequencing reads using target enrichment.* Applications in Plant Sciences 4(7):1600016.
- DOI: https://doi.org/10.3732/apps.1600016
- ⚠️ confirm — cited from general knowledge; details not re-fetched this session.

---

## D. *Alternaria* / Dothideomycetes — genus-specific genome & multi-locus phylogeny

### ⚠️ confirm exact details (titles/journals from literature search; DOIs not all re-fetched)
- Woudenberg JHC, Groenewald JZ, Binder M, Crous PW (2013). *Alternaria redefined.* Studies
  in Mycology 75:171–212. DOI: https://doi.org/10.3114/sim0015
  (the ITS + GAPDH + RPB2 + TEF1 + Alt a1 + endoPG + OPA10-2 multilocus framework, integrated
  with genome/transcriptome data).
- *Phylogenomic analyses of* Alternaria *section* Alternaria*: A high-resolution, genome-wide
  study of lineage sorting and gene tree discordance.* Mycologia (2021).
  https://www.tandfonline.com/doi/full/10.1080/00275514.2021.1950456
- *Alternaria: update on species limits, evolution, multi-locus phylogeny, and classification.*
  Studies in Fungi (2023). https://doi.org/10.48130/SIF-2023-0001

---

## One-line defense

The approach is not novel/risky: **UnFATE (Syst. Biol. 2025)** is a peer-reviewed pipeline doing
genome-mined markers → MACSE codon alignment → multilocus phylogeny + barcoding for this exact
clade, and **Feau et al. (2011)** established genome-mining of phylogenetic markers in non-model
fungi over a decade ago. phylofetch adds packaging, per-locus provenance, and type-material
integration — not a new method.
