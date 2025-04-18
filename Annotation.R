#################################################
# MARKER GENE IDENTIFICATION AND CELL TYPE ANNOTATION
#################################################

# Find all marker genes for cluster 2
# This identifies genes differentially expressed in cluster 2 compared to all other clusters
cluster2.markers <- FindMarkers(pbmc, ident.1 = 2)
head(cluster2.markers, n = 5)  # Display top 5 markers

# Find markers distinguishing cluster 5 from clusters 0 and 3 specifically
# Useful for understanding differences between specific cell populations
cluster5.markers <- FindMarkers(pbmc, ident.1 = 5, ident.2 = c(0, 3))
head(cluster5.markers, n = 5)  # Display top 5 markers

# Find markers for every cluster compared to all remaining cells
# only.pos = TRUE means we only report genes upregulated in each cluster
pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE)
head(pbmc.markers)  # Display top markers across all clusters

# Filter for highly expressed markers (log2 fold-change > 1)
# Group by cluster to organize markers by cell population
pbmc.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1)

# Find markers for cluster 0 using ROC test
# The ROC test is an alternative to the default Wilcoxon rank sum test
# logfc.threshold = 0.25 means we only consider genes with log2 fold-change > 0.25
cluster0.markers <- FindMarkers(pbmc, ident.1 = 0, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)

#################################################
# VISUALIZING MARKER GENES
#################################################

# Violin plots show expression distribution across clusters
# MS4A1 and CD79A are B cell markers
VlnPlot(pbmc, features = c("MS4A1", "CD79A"))

# You can plot raw counts (instead of normalized) on a log scale
# NKG7 is an NK cell marker, PF4 is a platelet marker
VlnPlot(pbmc, features = c("NKG7", "PF4"), slot = "counts", log = TRUE)

# Feature plots show expression on UMAP visualization
# This plots key markers for major immune cell types:
# - MS4A1: B cells
# - GNLY: NK cells
# - CD3E: T cells
# - CD14: Monocytes
# - FCER1A: Dendritic cells
# - FCGR3A: Non-classical monocytes
# - LYZ: Monocytes
# - PPBP: Platelets
# - CD8A: CD8+ T cells
FeaturePlot(pbmc, features = c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A", "LYZ", "PPBP", "CD8A"))

# Create a heatmap of top 10 markers per cluster
# slice_head(n = 10) selects the top 10 markers for each cluster
pbmc.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10

# Create a heatmap showing expression of top markers across clusters
DoHeatmap(pbmc, features = top10$gene) + NoLegend()

#################################################
# CELL TYPE ANNOTATION
#################################################

# Assign cell type labels to clusters based on marker gene expression
# This is based on canonical markers for immune cell populations:
new.cluster.ids <- c("Naive CD4 T", "CD14+ Mono", "Memory CD4 T", "B", "CD8 T", "FCGR3A+ Mono",
                     "NK", "DC", "Platelet")
names(new.cluster.ids) <- levels(pbmc)

# Apply the new cell type labels to the Seurat object
pbmc <- RenameIdents(pbmc, new.cluster.ids)

# Generate UMAP plot with cell type labels
DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()

#################################################
# FINAL VISUALIZATION AND EXPORT
#################################################

# Load ggplot2 for enhanced visualization options
library(ggplot2)

# Create a publication-quality UMAP plot with customized formatting
plot <- DimPlot(pbmc, reduction = "umap", label = TRUE, label.size = 4.5) + 
  xlab("UMAP 1") + 
  ylab("UMAP 2") +
  theme(axis.title = element_text(size = 18), 
        legend.text = element_text(size = 18)) + 
  guides(colour = guide_legend(override.aes = list(size = 10)))

# Save the UMAP visualization as a JPEG file
# Make sure the directory "../output/images/" exists
ggsave(filename = "../output/images/pbmc3k_umap.jpg", height = 7, width = 12, plot = plot, quality = 50)

# Note: After this analysis, you might want to:
# 1. Export the annotated Seurat object for future use
# 2. Perform additional analyses like trajectory inference or cell-cell interaction analysis
# 3. Compare cell type proportions between different conditions if you have multiple samples
# 4. Investigate specific gene modules or pathways of interest
