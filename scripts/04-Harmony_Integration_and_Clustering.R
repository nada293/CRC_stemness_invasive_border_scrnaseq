# ==============================================================================
# Stage 4: Harmony_Integration_and_Clustering
# Project: GSE144735 - Single-Cell RNA Sequencing Analysis
# Author: Nada Qenawy
# ==============================================================================

# 1. Load Required Libraries
library(Seurat)
library(harmony)
library(dplyr)
library(stringr)
library(ggplot2)

# 2. Extract Metadata Variables (patient_id & tissue_site)
# Format of orig.ident: "KUL01-B", "KUL01-T", "KUL01-N", etc.
seurat_obj$patient_id  <- sub("-.*", "", seurat_obj$orig.ident)
seurat_obj$tissue_site <- sub(".*-", "", seurat_obj$orig.ident)

# Verify new metadata columns
head(seurat_obj@meta.data[, c("orig.ident", "patient_id", "tissue_site")])

# 3. Create Before-Integration UMAP Plot for Comparison
p1 <- DimPlot(seurat_obj, reduction = "umap", group.by = "orig.ident") + 
  ggtitle("Before Integration")

# 4. Perform Harmony Integration Across Patients
seurat_obj <- RunHarmony(
  seurat_obj,
  group.by.vars    = "patient_id",  # Preserves biological variance (Border vs Core)
  reduction.use    = "pca",
  plot_convergence = TRUE,
  nclust           = 50,
  max_iter         = 10,
  early_stop       = TRUE
)

# 5. Run UMAP on Harmony Embeddings
seurat_obj <- RunUMAP(
  seurat_obj, 
  reduction      = "harmony", 
  dims           = 1:15, 
  reduction.name = "umap_after"
)

# 6. Perform Clustering on Integrated Data
seurat_obj <- FindNeighbors(seurat_obj, reduction = "harmony", dims = 1:15)
seurat_obj <- FindClusters(seurat_obj, resolution = 0.5)

# 7. Generate Post-Integration Visualizations
p2 <- DimPlot(seurat_obj, reduction = "umap_after", group.by = "patient_id") + 
  ggtitle("After Integration (by Patient)")

p3 <- DimPlot(seurat_obj, reduction = "umap_after", group.by = "orig.ident") + 
  ggtitle("After Integration (by Sample)")

p_clusters <- DimPlot(
  seurat_obj, 
  reduction = "umap_after", 
  group.by  = "seurat_clusters", 
  label     = TRUE, 
  repel     = TRUE
) + ggtitle("Harmony Integrated Clusters")

# Display plots
print(p1 + p2)
print(p_clusters)

# 8. Export & Save Figures
if (!dir.exists("figures")) dir.create("figures")

ggsave("figures/Harmony_Integration_Comparison.png", plot = p1 + p2, width = 10, height = 5, dpi = 300)
ggsave("figures/Harmony_Clusters.png", plot = p_clusters, width = 7, height = 6, dpi = 300)

# 9. Save Integrated Seurat Object
saveRDS(seurat_obj, file = "seurat_obj_harmony.rds")

cat("\n[SUCCESS] Harmony Integration & Clustering completed and saved successfully!\n")