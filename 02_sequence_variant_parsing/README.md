# BLAST Pairwise Alignment Parser & Microbial Variant Assembler

[![Python](https://img.shields.io/badge/Python-3.10-3776AB.svg?logo=python&logoColor=white)](https://www.python.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

A command-line Python utility engineered to parse NCBI pairwise BLAST reports (Format 0), automate multi-block alignment assembly, extract genomic coordinate metadata, manage indels, and export standardized FASTA sequences for downstream variant profiling and phylogenetics[cite: 1, 11].

---

## Scientific Background

*Helicobacter pylori* outer membrane proteins (OMPs) mediate host epithelial cell adhesion, colonization persistence, and virulence modulation[cite: 8]. During genomic research at Shahid Beheshti University, 31 whole-genome sequenced (WGS) clinical *H. pylori* isolates from Iranian patients were screened to discover novel allelic variants across five major adhesin loci: *alpA*, *alpB*, *babA*, *hopZ*, and *sabA*[cite: 8].

Reference gene sequences were obtained from *H. pylori* strain MT5135 (*alpA*, *alpB*, *babA*, *sabA*) and strain 26695 (*hopZ*) from NCBI databases[cite: 8]. All assembly contigs were searched using BLAST with stringent thresholds of $\ge 95\%$ coverage and $\ge 85\%$ sequence identity[cite: 8].

### Genomic Discovery Impact
This automated sequence parsing pipeline facilitated the screening, verification, and assembly of **76 novel non-redundant alleles** deposited into the NCBI GenBank repository[cite: 9]:
* **`alpA`**: 30 novel alleles identified[cite: 9]
* **`alpB`**: 27 novel alleles identified[cite: 9]
* **`hopZ`**: 28 novel alleles identified[cite: 9]
* **`sabA`**: 21 novel alleles identified[cite: 9]
* **`babA`**: 5 novel alleles identified[cite: 9]

This work was presented as a poster at the *6th International & 18th Iranian Genetics Congress (GC2024)* (Abstract G-1201)[cite: 5, 8].

---

## Directory Architecture

```text
02_sequence_variant_parsing/
├── README.md
├── blast_variant_extractor.py       # Modular CLI parsing engine
├── example_data/
│   └── HPC2_AlpA.txt                # Clinical isolate WGS assembly BLAST alignment
└── results/
    ├── .gitkeep
    └── HPC2_AlpA_extracted.fasta    # Standardized FASTA sequence output
