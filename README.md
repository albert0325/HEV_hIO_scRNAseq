# Hepatitis E virus replication is maintained in proliferative cells within the intestinal crypt

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![GEO](https://img.shields.io/badge/GEO-GSE303209-blue.svg)](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE303209)
[![Journal](https://img.shields.io/badge/Journal-Science%20Advances-red.svg)](https://www.science.org/doi/full/10.1126/sciadv.aeb2333)
[![R](https://img.shields.io/badge/R-4.5.2-276DC3.svg)](https://www.r-project.org/)

---

## Overview

This repository contains scripts accompanying the manuscript:

> **Hepatitis E virus replication is maintained in proliferative cells within the intestinal crypt**  
> Prallet et al., Sci. Adv. 12, eaeb2333 (2026)   

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

---
## Citation

If you use this code or data, please cite:

```
Hepatitis E virus replication is maintained inproliferative cells within the intestinal crypt. 
Prallet et al., Sci. Adv. 12, eaeb2333 (2026)
(https://doi.org/10.1126/sciadv.aeb2333)

```

---

## License

This project is licensed under the MIT License — see the [LICENSE](LICENSE) file for details.




