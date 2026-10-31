# Computational Biology & Microbial Functional Genomics Toolkit

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python](https://img.shields.io/badge/Python-3.10-3776AB.svg?logo=python&logoColor=white)](https://www.python.org/)
[![R](https://img.shields.io/badge/R-%3E%3D4.2.0-276DC3.svg?logo=r&logoColor=white)](https://www.r-project.org/)
[![Bioconductor](https://img.shields.io/badge/Bioconductor-3.18-green.svg)](https://bioconductor.org/)

An integrated computational biology repository containing reproducible analytical pipelines for **transcriptomic microarray normalization and paired differential expression analysis**, alongside **high-throughput BLAST pairwise alignment parsing for microbial variant discovery**
---

## Table of Contents
- [Project Architecture](#project-architecture)
- [Module 1: Transcriptomics & Differential Expression](#module-1-transcriptomics--differential-expression-01_transcriptomics_microarray)
  - [Overview & Scientific Background](#overview--scientific-background)
  - [Analytical Workflow](#analytical-workflow)
  - [Candidate Biomarker Discovery](#candidate-biomarker-discovery)
  - [Visualizations](#visualizations)
  - [Usage](#usage-transcriptomics)
- [Module 2: BLAST Sequence & Variant Parsing](#module-2-blast-sequence--variant-parsing-02_sequence_variant_parsing)
  - [Overview & Scientific Background](#overview--scientific-background-1)
  - [Key Features](#key-features)
  - [Usage & Real Benchmark](#usage--real-benchmark)
- [Environment Setup & Installation](#environment-setup--installation)
- [Conference Citations](#conference-citations)
- [License & Contact](#license--contact)

---

## Project Architecture

```text
bioinformatics-genomics-toolkit/
├── .gitignore
├── LICENSE
├── README.md
├── environment.yml
│
├── 01_transcriptomics_microarray/
│   ├── README.md
│   ├── normalize_affy_geo.R              # RMA background correction & GEO annotation
│   ├── diff_expr_limma.R                 # Paired limma empirical Bayes modeling
│   ├── data/
│   │   └── .gitkeep                      # Raw .CEL / .tar files directory
│   └── results/
│       ├── .gitkeep
│       ├── volcano_plot_GSE79973.png     # Publication-grade volcano plot
│       ├── heatmap_top_DEGs_GSE79973.png # Clustered heatmap across 10 matched pairs
│       └── GSE79973_DEG_limma_results.tsv# Full eBayes statistical summary table
│
└── 02_sequence_variant_parsing/
    ├── README.md
    ├── blast_variant_extractor.py        # Automated NCBI BLAST report parser
    ├── example_data/
    │   └── HPC2_AlpA.txt                 # Authentic H. pylori clinical contig alignment
    └── results/
        ├── .gitkeep
        └── HPC2_AlpA_extracted.fasta     # Standardized FASTA sequence output
