#################################################
# MARKER GENE IDENTIFICATION AND CELL TYPE ANNOTATION
#################################################


# Load required libraries
library(Seurat)
library(dplyr)
library(patchwork)
library(ggplot2)
library(presto)

# Create necessary directories
dir.create("output/images", recursive = TRUE, showWarnings = FALSE)
dir.create("output/data", recursive = TRUE, showWarnings = FALSE)

# Load preprocessed data
pbmc <- readRDS("output/data/pbmc_preprocessed.rds")

# Cluster cells
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)

# Run UMAP
pbmc <- RunUMAP(pbmc, dims = 1:10)

# Save initial UMAP clustering
png("output/images/umap_clusters.png", width = 900, height = 700, res = 100)
DimPlot(pbmc, reduction = "umap")
dev.off()

# Find markers for cluster 2
cluster2.markers <- FindMarkers(pbmc, ident.1 = 2)
write.csv(cluster2.markers, "output/data/cluster2_markers.csv")

# Find markers distinguishing cluster 5 from clusters 0 and 3
cluster5.markers <- FindMarkers(pbmc, ident.1 = 5, ident.2 = c(0, 3))
write.csv(cluster5.markers, "output/data/cluster5_vs_0_3_markers.csv")

# Find markers for all clusters
pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE)
write.csv(pbmc.markers, "output/data/all_cluster_markers.csv")

# Save filtered markers with high fold change
high_fc_markers <- pbmc.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1)
write.csv(high_fc_markers, "output/data/high_fc_markers.csv")

# Find markers for cluster 0 using ROC test
cluster0.markers <- FindMarkers(pbmc, ident.1 = 0, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)
write.csv(cluster0.markers, "output/data/cluster0_markers_roc.csv")

# Save violin plots for key markers
png("output/images/bcell_markers_vlnplot.png", width = 900, height = 600, res = 100)
VlnPlot(pbmc, features = c("MS4A1", "CD79A"))
dev.off()

png("output/images/innate_markers_vlnplot.png", width = 900, height = 600, res = 100)
VlnPlot(pbmc, features = c("NKG7", "PF4"), slot = "counts", log = TRUE)
dev.off()

# Save feature plots for key cell type markers
png("output/images/immune_markers_featureplot.png", width = 1200, height = 1200, res = 100)
FeaturePlot(pbmc, features = c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A", "LYZ", "PPBP", "CD8A"))
dev.off()

# Generate top 10 markers per cluster for heatmap
top10 <- pbmc.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup()
write.csv(top10, "output/data/top10_markers_per_cluster.csv")

# Save heatmap of top markers
png("output/images/top_markers_heatmap.png", width = 1200, height = 1800, res = 100)
DoHeatmap(pbmc, features = top10$gene) + NoLegend()
dev.off()

# Annotate cell types
new.cluster.ids <- c("Naive CD4 T", "CD14+ Mono", "Memory CD4 T", "B", "CD8 T", "FCGR3A+ Mono",
                     "NK", "DC", "Platelet")
names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.cluster.ids)

# Save basic UMAP with cell type labels
png("output/images/umap_celltypes.png", width = 900, height = 700, res = 100)
DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
dev.off()

# Create publication-quality UMAP plot
plot <- DimPlot(pbmc, reduction = "umap", label = TRUE, label.size = 4.5) + 
  xlab("UMAP 1") + 
  ylab("UMAP 2") +
  theme(axis.title = element_text(size = 18), 
        legend.text = element_text(size = 18)) + 
  guides(colour = guide_legend(override.aes = list(size = 10)))

# Save high-quality UMAP
ggsave(filename = "output/images/pbmc3k_umap.jpg", height = 7, width = 12, plot = plot, quality = 50)
ggsave(filename = "output/images/pbmc3k_umap.png", height = 7, width = 12, plot = plot)

# Save final annotated Seurat object
saveRDS(pbmc, file = "output/data/pbmc_final.rds")

print("Clustering and annotation complete. Final data saved to output/data/pbmc_final.rds")
# Note: After this analysis, you might want to:
# 1. Export the annotated Seurat object for future use
# 2. Perform additional analyses like trajectory inference or cell-cell interaction analysis
# 3. Compare cell type proportions between different conditions if you have multiple samples
# 4. Investigate specific gene modules or pathways of interest
