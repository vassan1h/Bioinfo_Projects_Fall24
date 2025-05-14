# HD-BrainSeq Shiny App

An interactive Shiny application for exploratory analysis of mRNA-Seq expression data from human post-mortem BA9 brain tissue. This app enables bioinformatic assessment of differential gene expression in Huntington’s Disease (HD) vs. control samples.

## 🔬 Features

- **Sample Metadata Visualization**
  - Summary tables, raw metadata view, histograms of continuous variables
- **RNA-Seq Count QC & Filtering**
  - Variance and sparsity-based gene filtering
  - Diagnostic scatter plots
  - Log10-transformed heatmap
  - Interactive PCA with axis selection
- **Differential Expression Analysis (DE)**
  - Filterable DESeq2 results
  - Histograms for p-values and log2 fold changes
  - Interactive volcano plot with customizable thresholds and color
- **Pathway Enrichment (ReactomePA)**
  - Ensembl-to-Entrez ID conversion via `org.Hs.eg.db`
  - Functional enrichment via `ReactomePA`
  - Barplots, dotplots, and network plots of enriched pathways

## 📁 Input File Requirements

| Tab | File Type | Format |
|-----|-----------|--------|
| Samples | Sample metadata | `.csv` with clinical or experimental metadata |
| Counts | Raw or normalized matrix | `.csv` or `.tsv`, genes × samples |
| DE | DESeq2 results | `.csv` with `log2FoldChange`, `pvalue`, `padj` columns |
| Reactome | Same DE file as above | Must use Ensembl IDs for gene identifiers |

## Getting Started

### 1. Install Required Packages

```r
install.packages(c("shiny", "bslib", "ggplot2", "tidyverse", "forcats", "dplyr",
                   "colourpicker", "DT", "patchwork", "gplots", "purrr"))

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("ReactomePA", "clusterProfiler", "biomaRt",
                       "org.Hs.eg.db", "rtracklayer", "fgsea"))
