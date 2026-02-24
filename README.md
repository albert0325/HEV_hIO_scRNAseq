# Hepatitis E virus replication is maintained in proliferative cells within the intestinal crypt

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![GEO](https://img.shields.io/badge/GEO-GSE303209-blue.svg)](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE303209)
[![Journal](https://img.shields.io/badge/Journal-Science%20Advances-red.svg)](https://www.science.org/journal/sciadv)
[![R](https://img.shields.io/badge/R-%3E%3D4.3-276DC3.svg)](https://www.r-project.org/)

---

## Overview

This repository contains all bioinformatic analysis code accompanying the manuscript:

> **Hepatitis E virus replication is maintained in proliferative cells within the intestinal crypt**  
> *Science Advances* (2026)  
> Computational Regulatory Omics Lab (CROmLab) — Institute for Pharmacy and Molecular Biotechnology, Heidelberg University & BioQuant, Heidelberg, Germany

The study demonstrates that HEV predominantly infects proliferative transit-amplifying (TA) cells and intestinal stem cells (ISCs) within the crypts of human pluripotent stem cell-derived intestinal organoids (hIOs), using single-cell RNA sequencing (scRNA-seq), stemness scoring, and pseudotime trajectory analysis.

---


## Data Availability

Raw and processed scRNA-seq data are publicly available on NCBI Gene Expression Omnibus (GEO):

**Accession: [GSE303209](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE303209)**

The `data/` folder in this repository does **not** track large data files. Please see [`data/README_data.md`](data/README_data.md) for detailed instructions on downloading and organising the data locally to reproduce all analyses.

---

## Requirements

### R version
R ≥ 4.3.0 is recommended.

### Core packages

| Package | Version tested | Purpose |
|---|---|---|
| [Seurat](https://satijalab.org/seurat/) | ≥ 5.0 | scRNA-seq processing, clustering, UMAP |
| [Monocle3](https://cole-trapnell-lab.github.io/monocle3/) | ≥ 1.3 | Pseudotime trajectory analysis |
| [SingleR](https://bioconductor.org/packages/SingleR/) | ≥ 2.4 | Automated cell type annotation |
| [UCell](https://bioconductor.org/packages/UCell/) | ≥ 2.6 | Stemness scoring |
| [tidyverse](https://www.tidyverse.org/) | ≥ 2.0 | Data wrangling and ggplot2 visualisation |
| [patchwork](https://patchwork.data-imaginist.com/) | ≥ 1.2 | Figure panel composition |
| [SeuratDisk](https://github.com/mojaveazure/seurat-disk) | ≥ 0.0.0.9021 | h5Seurat / h5ad I/O |

### Installation

```r
# CRAN packages
install.packages(c("Seurat", "tidyverse", "patchwork", "remotes"))

# Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c("SingleR", "UCell", "celldex"))

# Monocle3 (GitHub)
remotes::install_github("cole-trapnell-lab/monocle3")

# SeuratDisk (GitHub)
remotes::install_github("mojaveazure/seurat-disk")
```

---

## Reproducing the Analysis

1. **Clone this repository**
   ```bash
   git clone https://github.com/CROmLab/HEV-intestinal-organoids-scRNAseq.git
   cd HEV-intestinal-organoids-scRNAseq
   ```

2. **Download the data**  
   Follow the instructions in [`data/README_data.md`](data/README_data.md) to download the count matrices from GEO accession [GSE303209](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE303209) and place them in the `data/` directory.

3. **Run the scripts in order**  
   Scripts are numbered and should be run sequentially:
   ```
   Fig3_01_quality_control.R
   Fig3_02_annotation.R
   Fig3_03_stemness_pseudotime.R
   Fig3_04_HEV_tropism.R
   Fig3_05_ISG_response.R
   ```
   Each script saves intermediate `.rds` objects that are loaded by the next script, so the order matters.

4. **Outputs**  
   Figures and tables are written to an `output/` directory that is created automatically on the first run.

---

## Key Findings Recapitulated by This Code

- **17,902 single cells** passing QC from HEVcc-infected hIOs at 14 days post-infection
- Cell type annotation into ISCs, TA cells, enterocytes, goblet, enteroendocrine, and Paneth cells
- **500 HEV RNA-positive cells** (2.78% infection rate), enriched in TA cells and ISCs
- Stemness scoring confirming higher stemness in TA/ISC populations
- Pseudotime trajectory recapitulating ISC → TA → differentiated lineage hierarchy
- ISG upregulation restricted to HEV-infected enterocytes

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

**Computational Regulatory Omics Lab (CROmLab)**  
Institute for Pharmacy and Molecular Biotechnology, Heidelberg University & BioQuant  
Heidelberg, Germany  

For questions related to the bioinformatic analysis, please open a [GitHub Issue](../../issues) or contact the corresponding author via the journal.
