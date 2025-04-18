#################################################
# Single-cell RNA-seq Analysis using Seurat
# This script analyzes a PBMC dataset with quality control,
# normalization, clustering, and visualization
#################################################

# Install required packages
# Note: It's good practice to check if packages are already installed before installing
pkg_list <- c("Seurat", "dplyr", "patchwork")

# Fixed the typo in BiocManager and added proper error checking
for (pkg in pkg_list) {
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  if (!requireNamespace(pkg, quietly = TRUE))
    BiocManager::install(pkg)
}

# Optional package installations (commented out)
# install.packages('ragg')  # For high-quality graphics
# 
# install.packages('devtools')
# devtools::install_github('immunogenomics/presto')  # For fast Wilcoxon tests

# Load the required packages
for (item in pkg_list) {
  library(item, character.only = TRUE) 
}

# Load presto for fast differential expression analysis
library('presto')

#################################################
# DATA IMPORT
#################################################
# Load the PBMC dataset (Peripheral Blood Mononuclear Cells)
# You may need to adjust the path to match your file location
pbmc.data <- Read10X(data.dir = "/home/group_nithya01/NR_Aditya/Lab_work/project/scRNA_demo/pbmc3k_filtered_gene_bc_matrices/filtered_gene_bc_matrices/hg19/")

# Check dimensions (genes x cells)
pbmc.data@Dim 

# Example: Looking at expression of specific genes in the first 30 cells
# CD3D (T cells), TCL1A (B cells), MS4A1 (B cells)
pbmc.data[c("CD3D", "TCL1A", "MS4A1"), 1:30] 

# Comparing memory usage: dense vs sparse matrix
# Sparse matrices are much more memory efficient for scRNA-seq data
dense.size <- object.size(as.matrix(pbmc.data))
print(dense.size, units = 'Mb') # Size of the dense matrix

sparse.size <- object.size(pbmc.data)
print(sparse.size, units = "Mb") # Size of the sparse matrix

# How many times larger is the dense representation?
print(dense.size/sparse.size)

#################################################
# QUALITY CONTROL
#################################################
# Create a Seurat object
# min.cells = 3: Keep genes expressed in at least 3 cells
# min.features = 200: Keep cells with at least 200 detected genes
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)
print(pbmc)

# Show QC metrics for the first 5 cells
head(pbmc@meta.data, 5)

# Calculate mitochondrial gene percentage
# High mitochondrial gene content often indicates dying/stressed cells
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
head(pbmc)
class(pbmc)

# Visualize QC metrics
# nFeature_RNA: Number of genes detected per cell
# nCount_RNA: Total number of molecules detected per cell
# percent.mt: Percentage of mitochondrial genes
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# Examine relationships between QC metrics
# Cell complexity (genes vs UMIs) and mitochondrial percentage
plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# Filter cells based on QC metrics:
# - Keep cells with 200-2500 features (genes)
# - Keep cells with <5% mitochondrial genes
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

#################################################
# NORMALIZATION
#################################################
# Log-normalization (standard for scRNA-seq)
pbmc_ln <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)

# Alternative: CPM normalization (counts per million)
pbmc_cpm <- NormalizeData(pbmc, normalization.method = "RC", scale.factor = 1e6)

# Proceed with log-normalized data
pbmc = pbmc_ln

#################################################
# FEATURE SELECTION
#################################################
# Identify highly variable genes (features)
# These genes will be used for downstream analysis
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(pbmc), 10)

# Plot variable features with and without labels
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

# Reset any existing graphics device
dev.off()

#################################################
# DATA SCALING
#################################################
# Scale data for all genes (can be memory intensive for large datasets)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)

# Regress out the effect of mitochondrial percentage 
# This helps to remove technical variation
pbmc <- ScaleData(pbmc, vars.to.regress = "percent.mt")

#################################################
# DIMENSIONALITY REDUCTION
#################################################
# Run PCA on variable features
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))

# Examine PCA results
print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)

# Visualize the top features for PC1 and PC2
VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")

# PCA plot
DimPlot(pbmc, reduction = "pca") + NoLegend()

# Heatmap for PC1
DimHeatmap(pbmc, dims = 1, cells = 500, balanced = TRUE)

# Heatmaps for PCs 1-15
DimHeatmap(pbmc, dims = 1:15, cells = 500, balanced = TRUE)

# Determine dimensionality for downstream analysis
ElbowPlot(pbmc)

#################################################
# CLUSTERING
#################################################
# Construct a KNN graph based on the euclidean distance in PCA space
# Using the first 10 PCs based on elbow plot
pbmc <- FindNeighbors(pbmc, dims = 1:10)

# Apply Louvain algorithm for community detection (clustering)
# Resolution parameter controls clustering granularity
pbmc <- FindClusters(pbmc, resolution = 0.5)

#################################################
# VISUALIZATION
#################################################
# Run UMAP for visualization
pbmc <- RunUMAP(pbmc, dims = 1:10)

# Plot UMAP with clusters
# Each cluster represents a putative cell type
DimPlot(pbmc, reduction = "umap")

# Note: The next steps would typically include:
# 1. Identifying marker genes for each cluster
# 2. Annotating cell types based on known markers
# 3. Differential expression analysis between conditions
# 4. Trajectory analysis or additional visualizations
