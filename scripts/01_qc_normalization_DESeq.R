#################################################
# Single-cell RNA-seq Analysis using Seurat
# This script analyzes a PBMC dataset with quality control,
# normalization, clustering, and visualization
#################################################

# Create necessary directories
dir.create("output/images", recursive = TRUE, showWarnings = FALSE)
dir.create("output/data", recursive = TRUE, showWarnings = FALSE)

#install packages
pkg_list <- c("Seurat", "dplyr", "patchwork")

for (pkg in pkg_list) {
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  if (!requireNamespace(pkg, quietly = TRUE))
    BiocManager::install(pkg)
}

#importing packages
for (item in pkg_list) {
  library(item, character.only = TRUE) 
}

# Set the data directory - adjust as needed
data_dir <- "/home/group_nithya01/NR_Aditya/Lab_work/project/scRNA_demo/pbmc3k_filtered_gene_bc_matrices/filtered_gene_bc_matrices/hg19/"

#import dataset
# Load the PBMC dataset
pbmc.data <- Read10X(data.dir = data_dir)
pbmc.data@Dim #genes bs

# Check some marker genes in first 30 cells
pbmc.data[c("CD3D", "TCL1A", "MS4A1"), 1:30]

# Compare dense vs sparse matrix size
dense.size <- object.size(as.matrix(pbmc.data))
sparse.size <- object.size(pbmc.data)
print(paste("Dense matrix size:", format(dense.size, units = 'Mb'), "MB"))
print(paste("Sparse matrix size:", format(sparse.size, units = "Mb"), "MB"))
print(paste("Ratio:", round(dense.size/sparse.size, 2), "times larger"))

##QC
# Create Seurat object
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)
print(pbmc)

# Calculate mitochondrial gene percentage
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

# Save QC metrics plot
png("output/images/qc_violin_plots.png", width = 900, height = 600, res = 100)
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()

# Save feature correlation plots
png("output/images/feature_correlations.png", width = 1000, height = 500, res = 100)
plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
print(plot1 + plot2)
dev.off()

# Filter cells based on QC metrics
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

# Normalize data
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)

# Find variable features
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

# Save variable feature plot
png("output/images/variable_features.png", width = 1000, height = 500, res = 100)
top10 <- head(VariableFeatures(pbmc), 10)
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
print(plot1 + plot2)
dev.off()

# Scale data
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- ScaleData(pbmc, vars.to.regress = "percent.mt")

# Run PCA
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))

# Save PCA plots
png("output/images/pca_loadings.png", width = 800, height = 600, res = 100)
VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")
dev.off()

png("output/images/pca_plot.png", width = 800, height = 600, res = 100)
DimPlot(pbmc, reduction = "pca") + NoLegend()
dev.off()

# Save elbow plot
png("output/images/elbow_plot.png", width = 800, height = 600, res = 100)
ElbowPlot(pbmc)
dev.off()

# Save PC heatmaps
png("output/images/pc1_heatmap.png", width = 800, height = 1000, res = 100)
DimHeatmap(pbmc, dims = 1, cells = 500, balanced = TRUE)
dev.off()

png("output/images/pcs_heatmap.png", width = 1200, height = 1800, res = 100)
DimHeatmap(pbmc, dims = 1:15, cells = 500, balanced = TRUE)
dev.off()

# Save preprocessed data
saveRDS(pbmc, file = "output/data/pbmc_preprocessed.rds")

print("QC and normalization complete. Preprocessed data saved to output/data/pbmc_preprocessed.rds")


# Note: The next steps would typically include:
# 1. Identifying marker genes for each cluster
# 2. Annotating cell types based on known markers
# 3. Differential expression analysis between conditions
# 4. Trajectory analysis or additional visualizations
