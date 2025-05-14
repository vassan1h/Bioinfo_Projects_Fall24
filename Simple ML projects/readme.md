# Project 1: Genome Utility Functions

A collection of basic bioinformatics utility functions for working with genome sequences. Developed as part of an introductory project for understanding algorithmic operations in computational genomics.

## Functionality - Project - 1 :

The project includes standalone implementations of:

- **Hamming Distance**  
  Compute the number of differences between two equal-length strings.

- **Subsequence Detection**  
  Determine if one sequence is a subsequence of another.

- **Skew Calculation**  
  Track the GC skew across a genome to identify replication origins.

- **Minimum Skew Position**  
  Find the index where GC skew is minimized.

- **Pattern Matching**  
  Search for all exact matches of a pattern within a genome.

- **Frequent Words**  
  Identify the most frequent k-mers in a sequence.

Each function is implemented from scratch with an emphasis on clarity and algorithmic logic.

---

# Project 2: Replication of RNA-Seq Figures on Transcriptional Memory

This project replicates key figures from the study:

**Nickel-induced transcriptional memory in lung epithelial cells promotes interferon signaling upon nicotine exposure**  
📄 DOI: [10.1016/j.taap.2023.116753](https://doi.org/10.1016/j.taap.2023.116753)

## 🎯 Objective

To reproduce and validate RNA-seq–based visualizations (PCA, fold-change bar plots) from the study, which investigates how prior nickel exposure primes lung epithelial cells to mount a heightened interferon signaling response upon nicotine exposure.

## 🧬 Biological Context

- **System**: BEAS-2B lung epithelial cells
- **Treatments**: 
  - UT (Untreated)
  - UT-NTE (Nicotine-exposed)
  - NiW (Nickel pretreated and washed out)
  - NiW-NTE (Nickel-pretreated then nicotine-exposed)
- **Focus**: Transcriptional memory and IFN signaling genes

## 📈 Replicated Figures

- **Figure 2A (PCA Plot)**:  
  PCA of normalized RNA-seq data to visualize condition-specific clustering. This project uses a centered log-ratio (CLR) transformation (Python) vs. VST (DESeq2) in the original paper.

- **Figure 3F (Bar Plots of IFN Pathway Genes)**:  
  Bar plots for selected interferon pathway genes showing fold changes and significance across treatment conditions.

## Methods

- **Data Source**: Public RNA-seq counts and metadata (link below)
- **Preprocessing**: 
  - Filtering low-expression genes
  - Metadata mapping
- **Analysis Tools**: 
  - `PyDESeq2` for differential expression
  - `scikit-learn` for PCA
  - `matplotlib`/`seaborn` for plotting

## Some Insights

- Nickel pretreatment causes persistent transcriptional shifts (visible in PCA).
- Interferon signaling is more strongly induced by nicotine in previously nickel-exposed cells.
- Genes like **HLA-F**, **MX1**, **IFIM1**, and **IRF6** show condition-specific regulation, reinforcing the hypothesis of toxicant-induced memory.

## Usage

You know it.

## 👨‍🔬 Author

Vassanth M., MS Bioinformatics @BU 
