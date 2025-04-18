# Single-cell RNA-seq Analysis Portfolio

## 📊 PBMC Dataset Analysis with Seurat

This repository contains a comprehensive single-cell RNA sequencing (scRNA-seq) analysis workflow of Peripheral Blood Mononuclear Cells (PBMCs) using the Seurat package in R. The analysis demonstrates proficiency in computational biology techniques and provides a reproducible workflow for scRNA-seq data processing.

![UMAP_Visualization](../output/images/pbmc3k_umap.jpg)
*UMAP visualization of PBMC cell types identified through unsupervised clustering and differential expression analysis.*

## 🧬 Analysis Overview

This project demonstrates the complete workflow for analyzing single-cell RNA sequencing data, from quality control to cell type annotation:

1. **Data Import & Quality Control**
   - Loading 10X Genomics data
   - QC metrics calculation and filtering
   - Removal of low-quality cells and potential doublets

2. **Normalization & Feature Selection**
   - Log-normalization of count data
   - Identification of highly variable genes
   - Data scaling and mitochondrial percentage regression

3. **Dimensionality Reduction**
   - Principal Component Analysis (PCA)
   - Elbow plot assessment for optimal dimensions
   - UMAP embedding for visualization

4. **Clustering & Cell Type Identification**
   - Graph-based clustering using Louvain algorithm
   - Marker gene identification through differential expression
   - Cell type annotation based on canonical markers

5. **Visualization & Result Interpretation**
   - Gene expression visualization (violin plots, feature plots)
   - Marker gene heatmaps
   - Annotated UMAP visualizations

## 🛠️ Technologies & Skills Demonstrated

- **Programming**: R
- **Packages**: Seurat, dplyr, patchwork, presto, ggplot2
- **Analysis Techniques**: 
  - Single-cell data processing
  - Dimensionality reduction
  - Unsupervised clustering
  - Differential gene expression analysis
  - Data visualization
  - Cell type annotation

## 📋 Repository Structure

```
├── scripts/
│   ├── 01_qc_normalization.R       # Quality control and data normalization
│   └── 02_clustering_annotation.R  # Clustering and cell type annotation
├── output/
│   ├── images/                     # Visualization outputs
│   └── data/                       # Processed data objects
└── README.md                       # Project documentation
```

## 🚀 How to Run

1. Clone this repository:
   ```bash
   git clone https://github.com/aditya19991010/single_cell_transcriptomics.git
   cd single_cell_transcriptomics
   ```

2. Install required R packages:
   ```R
   if (!requireNamespace("BiocManager", quietly = TRUE))
     install.packages("BiocManager")
   
   pkg_list <- c("Seurat", "dplyr", "patchwork")
   for (pkg in pkg_list) {
     if (!requireNamespace(pkg, quietly = TRUE))
       BiocManager::install(pkg)
   }
   ```

3. Run the analysis scripts in order:
   ```R
   source("scripts/01_qc_normalization.R")
   source("scripts/02_clustering_annotation.R")
   ```

## 📈 Results & Findings

The analysis successfully identified 9 distinct cell populations in the PBMC dataset:
- Naive CD4+ T cells
- CD14+ Monocytes
- Memory CD4+ T cells
- B cells
- CD8+ T cells
- FCGR3A+ Monocytes
- Natural Killer (NK) cells
- Dendritic Cells (DCs)
- Platelets

Each population was identified based on the expression of canonical marker genes and validated through differential expression analysis.

## 📚 References

- Stuart, T., Butler, A., Hoffman, P. et al. Comprehensive Integration of Single-Cell Data. Cell (2019).
- 10X Genomics. PBMC3k Dataset. https://support.10xgenomics.com/single-cell-gene-expression/datasets

## 🤝 Contact Information

Feel free to reach out if you have any questions or would like to discuss potential collaborations!

- **Email**: aditya19991010@gmail.com


---

*This project was created as part of my bioinformatics portfolio to demonstrate proficiency in single-cell transcriptomics analysis.*
