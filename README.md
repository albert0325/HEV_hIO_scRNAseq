# Hepatitis E virus replication is maintained in proliferative cells within the intestinal crypt

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![GEO](https://img.shields.io/badge/GEO-GSE303209-blue.svg)](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE303209)
[![Journal](https://img.shields.io/badge/Journal-Science%20Advances-red.svg)](https://www.science.org/journal/sciadv)
[![R](https://img.shields.io/badge/R-4.5.2-276DC3.svg)](https://www.r-project.org/)

---

## Overview

This repository contains all bioinformatic analysis code accompanying the manuscript:

> **Hepatitis E virus replication is maintained in proliferative cells within the intestinal crypt**  
> *Science Advances* (2026)  

The study demonstrates that HEV predominantly infects proliferative transit-amplifying (TA) cells and intestinal stem cells (ISCs) within the crypts of human pluripotent stem cell-derived intestinal organoids (hIOs), using single-cell RNA sequencing (scRNA-seq), stemness scoring, and pseudotime trajectory analysis.

---


## Data Availability

Raw and processed scRNA-seq data are publicly available on NCBI Gene Expression Omnibus (GEO):

**Accession: [GSE303209](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE303209)**

The `data/` folder in this repository does **not** track large data files. Please download and organising the data locally to reproduce all analyses.

---

## Requirements

### R version

Analyses were performed under **R 4.5.2** (2025-10-31) on macOS Sequoia 15.7.4
(x86\_64-apple-darwin20).

### Core packages

| Package | Tested version | Source | Purpose |
|---|---|---|---|
| [Seurat](https://satijalab.org/seurat/) | 5.4.0 | CRAN | scRNA-seq processing, SCTransform integration, UMAP |
| [SeuratObject](https://satijalab.org/seurat/) | 5.3.0 | CRAN | Seurat S4 class infrastructure |
| [SingleR](https://bioconductor.org/packages/SingleR/) | 2.12.0 | Bioconductor | Reference-based cell type annotation |
| [SummarizedExperiment](https://bioconductor.org/packages/SummarizedExperiment/) | 1.40.0 | Bioconductor | Reference matrix container for SingleR |
| [CytoTRACE2](https://cytotrace.stanford.edu/) | 1.1.0 | GitHub | Developmental potency scoring |
| [monocle3](https://cole-trapnell-lab.github.io/monocle3/) | 1.4.26 | GitHub | Pseudotime trajectory inference |
| [ComplexHeatmap](https://bioconductor.org/packages/ComplexHeatmap/) | 2.26.1 | Bioconductor | SingleR score heatmap |
| [circlize](https://jokergoo.github.io/circlize_book/) | 0.4.17 | CRAN | Colour scale for ComplexHeatmap |
| [ggpubr](https://rpkgs.datanovia.com/ggpubr/) | 0.6.2 | CRAN | Statistical annotation on plots |
| [rstatix](https://rpkgs.datanovia.com/rstatix/) | 0.7.3 | CRAN | Wilcoxon tests within dplyr pipelines |
| [tidyverse](https://www.tidyverse.org/) | 2.0.0 | CRAN | Data wrangling and ggplot2 |
| [ggplot2](https://ggplot2.tidyverse.org/) | 4.0.2 | CRAN | Core plotting |
| [readxl](https://readxl.tidyverse.org/) | 1.4.5 | CRAN | Reading marker gene Excel file |
| [patchwork](https://patchwork.data-imaginist.com/) | 1.3.2 | CRAN | Figure panel composition |
| [RColorBrewer](https://CRAN.R-project.org/package=RColorBrewer) | 1.1-3 | CRAN | Colour palettes |
| [future](https://future.futureverse.org/) | 1.69.0 | CRAN | Parallelisation for SCTransform |
| [viridisLite](https://sjmgarnier.github.io/viridisLite/) | 0.4.3 | CRAN | Viridis colour scales |

### Installation

```r
# CRAN packages
install.packages(c(
  "Seurat", "tidyverse", "ggplot2", "patchwork",
  "ggpubr", "rstatix", "readxl", "circlize",
  "RColorBrewer", "future", "viridisLite", "remotes"
))

# Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c(
  "SingleR", "SummarizedExperiment",
  "ComplexHeatmap", "BiocParallel"
))

# GitHub packages
remotes::install_github("cole-trapnell-lab/monocle3")
remotes::install_github("digitalcytometry/cytotrace2", subdir = "cytotrace2")
```


---

## Reproducing the Analysis

1. **Clone this repository**
   ```bash
   git clone https://github.com/albert0325/HEV_hIO_scRNAseq.git
   cd HEV_hIO_scRNAseq
   ```

2. **Download the data**  
   Download the count matrices from GEO accession [GSE303209](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE303209) and place them in the `data/` directory.

3. **Run scripts in order from an R session opened at the repository root**
   ```r
   source("R/01_quality_control.R")    # produces: hev object, outputs/S9A.tiff
   source("R/02_annotation.R")         # produces: fig3A–D, S9B–C
   source("R/03_stemness_pseudotime.R")# produces: S11A–C
   source("R/04_ISG_response.R")       # produces: boxplot_sig_schoggins.pdf
   ```
   Each script sources `R/utils.R` and writes all figures to `outputs/`
   (created automatically if absent).

---

## Citation

If you use this code or data, please cite:

```
[Citation to be added upon publication — DOI will be inserted here]
```

---

## License

This project is licensed under the MIT License — see the [LICENSE](LICENSE) file for details.

---

## Contact

**Albert Li (李律)**  
Institute for Pharmacy and Molecular Biotechnology (IPMB)  
BioQuant, Heidelberg University, Germany  

📧  albert0325162@gmail.com  
🌐 https://albert0325.github.io/
---

### Code & Data

For technical questions about the code or analyses, please open a GitHub Issue:  
https://github.com/albert0325/HEV_hIO_scRNAseq/issues  

---

### Journal Correspondence

For journal-related matters, please contact the corresponding author via the journal.

---

### Laboratory Contacts

**Dao Thi Lab**  
Dr. Viet Loan Dao Thi  
📧 VietLoan.DaoThi@med.uni-heidelberg.de  
🌐 https://daothilab.com/

**Computational Regulatory Omics Lab**  
Prof. Dr. Carl Herrmann  
📧 carl.herrmann@bioquant.uni-heidelberg.de  
🌐 https://www.hdsu.org/


