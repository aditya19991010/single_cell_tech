
#install packages
pkg_list <- c("Seurat", "dplyr", "patchwork")


for (pkg in pkg_list) {
  if (!require('BiocMangaer', quietly=True))
    installed.packages('BiocManager')
  BiocManager::install(pkg)
}
# install.packages('ragg')
# 
# install.packages('devtools')
# devtools::install_github('immunogenomics/presto')

#importing packages
for (item in pkg_list) {
  library(item, character.only = TRUE) 
}

library('presto')

#import dataset
# Load the PBMC dataset
pbmc.data <- Read10X(data.dir = "/home/group_nithya01/NR_Aditya/Lab_work/project/scRNA_demo/pbmc3k_filtered_gene_bc_matrices/filtered_gene_bc_matrices/hg19/")
pbmc.data@Dim #genes bs
pbmc.data[c("CD3D", "TCL1A", "MS4A1"), 1:30] #### each number represents the number of molecules detected

dense.size <- object.size(as.matrix(pbmc.data))
print(dense.size, units = 'Mb') #size of the matrix

sparse.size <- object.size(pbmc.data)
print(sparse.size, units = "Mb")

print(dense.size/sparse.size)

##QC
### Initialize the Seurat object with the raw (non-normalized data).
##################################################################################
#The number of unique genes detected in each cell.
#   Low-quality cells or empty droplets will often have very few genes
#   Cell doublets or multiplets may exhibit an aberrantly high gene count
#Similarly, the total number of molecules detected within a cell (correlates strongly with unique genes)
#The percentage of reads that map to the mitochondrial genome
#   Low-quality / dying cells often exhibit extensive mitochondrial contamination
#   We calculate mitochondrial QC metrics with the PercentageFeatureSet() function, which calculates the percentage of counts originating from a set of features
#   We use the set of all genes starting with MT- as a set of mitochondrial genes
##################################################################################

pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)
pbmc

# Show QC metrics for the first 5 cells
head(pbmc@meta.data, 5)

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
head(pbmc)
class(pbmc)

# Visualize QC metrics as a violin plot
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

#### FeatureScatter is typically used to visualize feature-feature relationships, but can be used
#### for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
# Shows Pearson correlation

plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

pbmc_ln <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc_cpm <- NormalizeData(pbmc, normalization.method = "RC", scale.factor = 1e6)

pbmc = pbmc_ln

pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(pbmc), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

dev.off()
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)

pbmc <- ScaleData(pbmc, vars.to.regress = "percent.mt")


pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))

# Examine and visualize PCA results a few different ways
print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)

VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")

DimPlot(pbmc, reduction = "pca") + NoLegend()


DimHeatmap(pbmc, dims = 1, cells = 500, balanced = TRUE)

DimHeatmap(pbmc, dims = 1:15, cells = 500, balanced = TRUE)

ElbowPlot(pbmc)

pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)


pbmc <- RunUMAP(pbmc, dims = 1:10)

# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(pbmc, reduction = "umap")

