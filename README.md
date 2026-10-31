# Computational Biology & Microbial Functional Genomics Toolkit

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python](https://img.shields.io/badge/Python-3.10-3776AB.svg?logo=python&logoColor=white)](https://www.python.org/)
[![R](https://img.shields.io/badge/R-%3E%3D4.2.0-276DC3.svg?logo=r&logoColor=white)](https://www.r-project.org/)
[![Bioconductor](https://img.shields.io/badge/Bioconductor-3.18-green.svg)](https://bioconductor.org/)

An integrated computational biology repository containing reproducible analytical pipelines for **transcriptomic microarray normalization and paired differential expression analysis**, alongside **high-throughput BLAST pairwise alignment parsing for microbial variant discovery**
---

## Table of Contents
## <a id="project-architecture"></a>Project Architecture
...

## <a id="module-1-transcriptomics--differential-expression"></a>Module 1: Transcriptomics & Differential Expression (`01_transcriptomics_microarray`)

### <a id="module-1-overview--scientific-background"></a>Overview & Scientific Background
...

### <a id="analytical-workflow"></a>Analytical Workflow
...

### <a id="candidate-biomarker-discovery"></a>Candidate Biomarker Discovery
...

### <a id="visualizations"></a>Visualizations
...

### <a id="usage-transcriptomics"></a>Usage: Transcriptomics


## <a id="module-2-blast-sequence--variant-parsing"></a>Module 2: BLAST Sequence & Variant Parsing (`02_sequence_variant_parsing`)

### <a id="module-2-overview--scientific-background"></a>Overview & Scientific Background
...

### <a id="key-features"></a>Key Features
...

### <a id="usage--real-benchmark"></a>Usage & Real Benchmark

<img src="01_transcriptomics_microarray/results/volcano_plot_GSE79973.png" width="850">

<img src="01_transcriptomics_microarray/results/heatmap_top_DEGs_GSE79973.png" width="850">
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
