# ==============================================================================
# Stage 5 : Canonical and DEG Annotation
# Project: GSE144735 - Single-Cell RNA Sequencing Analysis
# Author: Nada Qenawy
# ==============================================================================

library(Seurat)
library(dplyr)
library(ggplot2)

# ---- 0. Paths & Directories Setup ----
project_dir <- "D:/github/CRC_stemness_invasive_border_scrnaseq"
results_dir <- file.path(project_dir, "results")
figures_dir <- file.path(project_dir, "figures")

dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(figures_dir, showWarnings = FALSE, recursive = TRUE)

# ---- 1. Load Integrated Seurat Object ----
seurat_obj <- readRDS(file.path(results_dir, "07_harmony_integrated_unannotated.rds"))

DefaultAssay(seurat_obj) <- "RNA"
Idents(seurat_obj) <- "seurat_clusters"


# ==============================================================================
# FIGURE 1: Canonical Marker Genes Across Clusters (DotPlot)
# ==============================================================================

# Define canonical cell type marker gene list
canonical_markers <- list(
  "T cells"            = c("CD3D", "CD3E", "TRAC", "CD247"),
  "NK cells"           = c("NKG7", "KLRD1", "GNLY", "XCL1", "XCL2"),
  "B cells"            = c("MS4A1", "CD79A", "CD74", "CD37", "CD19"),
  "Plasma cells"       = c("JCHAIN", "MZB1", "IGHG1", "IGKC"),
  "Myeloid cells"      = c("LST1", "TYROBP", "AIF1", "LILRB1"),
  "Monocytes"          = c("LYZ", "S100A8", "S100A9", "FCN1", "CD14"),
  "Macrophages"        = c("C1QA", "C1QB", "C1QC", "APOE", "FOLR2"),
  "Dendritic cells"    = c("FCER1A", "CD1C", "CLEC10A", "HLA-DRA"),
  "Mast cells"         = c("TPSAB1", "TPSB2", "KIT", "MS4A2"),
  "Platelets"          = c("PPBP", "PF4", "NRGN"),
  "Epithelial cells"   = c("EPCAM", "KRT8", "KRT18", "KRT19"),
  "Fibroblasts / CAFs" = c("COL1A1", "COL1A2", "COL3A1", "DCN", "LUM"),
  "Endothelial cells"  = c("PECAM1", "VWF", "KDR", "EMCN"),
  "Pericytes"          = c("RGS5", "CSPG4", "MCAM", "NOTCH3")
)

# Filter markers available in dataset
canonical_markers_available <- lapply(canonical_markers, function(g) {
  g[g %in% rownames(seurat_obj)]
})
canonical_markers_available <- canonical_markers_available[sapply(canonical_markers_available, length) > 0]

# Generate and format Figure 1 DotPlot
p_fig1_dotplot <- DotPlot(
  seurat_obj,
  features = canonical_markers_available,
  assay = "RNA",
  dot.scale = 6,
  col.min = -2.5,
  col.max = 2.5
) +
  RotatedAxis() +
  scale_color_gradient2(low = "lightgrey", mid = "blue", high = "red") +
  ggtitle("Canonical Marker Genes Across Clusters") +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 8),
    strip.text.x = element_text(size = 9, face = "bold")
  )

print(p_fig1_dotplot)

ggsave(
  filename = file.path(figures_dir, "Canonical_Marker_annotation.png"),
  plot = p_fig1_dotplot, width = 18, height = 8, dpi = 300
)


# ==============================================================================
# FIGURE 2: Cell-Type Annotation & UMAP
# ==============================================================================

# Mapping cluster IDs (0 - 17) to cell type annotations
corrected_canonical_labels <- c(
  "0"  = "T cells",
  "1"  = "Epithelial cells",
  "2"  = "Plasma cells",
  "3"  = "Myeloid / Macrophages",
  "4"  = "NK / T cells",
  "5"  = "B cells",
  "6"  = "Fibroblasts / CAFs",
  "7"  = "Epithelial cells",
  "8"  = "Endothelial cells",
  "9"  = "Epithelial (Cycling)",
  "10" = "Fibroblasts / CAFs",
  "11" = "Pericytes",
  "12" = "Fibroblasts / CAFs",
  "13" = "Fibroblasts / CAFs",
  "14" = "Stromal / Fibroblasts",
  "15" = "Mast cells",
  "16" = "Endothelial cells",
  "17" = "B cells / Plasma cells"
)

# Assign canonical annotations to metadata
cluster_ids <- as.character(seurat_obj$seurat_clusters)
seurat_obj$canonical_annotation <- unname(corrected_canonical_labels[cluster_ids])
Idents(seurat_obj) <- "canonical_annotation"

reduction_name <- ifelse("umap_after" %in% names(seurat_obj@reductions), "umap_after", "umapafter")

p_fig2_umap <- DimPlot(
  seurat_obj,
  reduction = reduction_name,
  group.by = "canonical_annotation",
  label = TRUE,
  repel = TRUE,
  pt.size = 0.3
) +
  ggtitle("Annotated Cell Types") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    legend.title = element_blank()
  )

print(p_fig2_umap)

ggsave(
  filename = file.path(figures_dir, "Annotated_Cell_types.png"),
  plot = p_fig2_umap, width = 12, height = 8, dpi = 300
)


# ==============================================================================
# FIGURE 3: DEG-Based Annotation UMAP
# ==============================================================================

p_fig3_umap <- DimPlot(
  seurat_obj,
  reduction = reduction_name,
  group.by = "canonical_annotation",
  label = TRUE,
  label.size = 4.5,
  label.box = FALSE,
  repel = TRUE,
  pt.size = 0.4
) +
  labs(
    title = "2. DEG-Based Annotation",
    x = "umapafter_1",
    y = "umapafter_2"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 11),
    legend.title = element_blank()
  )

print(p_fig3_umap)

ggsave(
  filename = file.path(figures_dir, "DEG_Based_Annotation.png"),
  plot = p_fig3_umap, width = 10, height = 8, dpi = 300
)


# ==============================================================================
# FIGURE 4: Top Differentially Expressed Genes Across Clusters (DEG DotPlot)
# ==============================================================================

DefaultAssay(seurat_obj) <- "RNA"
Idents(seurat_obj) <- "seurat_clusters"

# Compute differentially expressed genes across clusters
deg_clusters <- FindAllMarkers(
  seurat_obj,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25,
  test.use = "wilcox"
)

# Export DEG results to CSV
write.csv(
  deg_clusters,
  file.path(results_dir, "08_cluster_DEGs.csv"),
  row.names = FALSE
)

# Target gene list for Top DEG plot
top_deg_genes_exact <- c(
  "MS4A1", "BANK1", "VPREB3", "C10orf99", "FABP1", "MT1G", "ASPM", 
  "GTSE1", "NEK2", "CD40LG", "TNFRSF4", "IL7R", "COL10A1", 
  "RP11-400N13.3", "COL11A1", "RERGL", "COX4I2", "HIGD1B", "TM4SF18", 
  "SOX18", "MYCT1", "KRT24", "PCOLCE2", "FIGF", "NRXN1", "CDH19", 
  "MYOT", "FCN1", "FPR3", "CLEC10A", "IGHA2", "IGHA1", "JCHAIN", 
  "TPSAB1", "TPSB2", "CPA3", "GZMH", "XCL2", "GNLY"
)

available_deg_genes <- top_deg_genes_exact[top_deg_genes_exact %in% rownames(seurat_obj)]

p_fig4_deg_dotplot <- DotPlot(
  seurat_obj,
  features = available_deg_genes,
  group.by = "seurat_clusters",
  dot.scale = 5,
  col.min = -1,
  col.max = 2
) +
  RotatedAxis() +
  scale_color_gradient2(low = "royalblue", mid = "purple", high = "firebrick3") +
  ggtitle("Top Differentially Expressed Genes Across Clusters") +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 8, face = "italic"),
    axis.text.y = element_text(size = 9, face = "bold"),
    panel.grid.major = element_line(color = "grey92", linewidth = 0.3)
  )

print(p_fig4_deg_dotplot)

ggsave(
  filename = file.path(figures_dir, "Top_DEG.png"),
  plot = p_fig4_deg_dotplot, width = 18, height = 8, dpi = 300
)

# ---- Save Final Processed Object ----
saveRDS(seurat_obj, file.path(results_dir, "09_final_annotated_seurat_obj.rds"))