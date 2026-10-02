# ==============================================================================
# Stage 6 : Downstream Analysis
# Project: GSE144735 - Single-Cell RNA Sequencing Analysis
# Author: Nada Qenawy
# ==============================================================================

# ---- 0. Load Required Libraries ----
library(Seurat)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(ggrepel)

# ---- 1. Environment & Paths Setup ----
project_dir <- "D:/github/CRC_stemness_invasive_border_scrnaseq"
results_dir <- file.path(project_dir, "results")
figures_dir <- file.path(project_dir, "figures")

dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(figures_dir, showWarnings = FALSE, recursive = TRUE)

# ---- 2. Load Final Annotated Seurat Object ----
seurat_obj <- readRDS(file.path(results_dir, "09_final_annotated_seurat_obj.rds"))

# ---- 3. Subset Malignant Epithelial Cells ----
epithelial_cells <- subset(
  seurat_obj, 
  subset = canonical_annotation %in% c("Epithelial cells", "Epithelial (Cycling)")
)

DefaultAssay(epithelial_cells) <- "RNA"

# ==============================================================================
# PART A: ISC Stemness Module Scoring & Spatial Comparison (Border vs. Core)
# ==============================================================================

# ---- A1. Calculate Intestinal Stem Cell (ISC) Signature Score ----
isc_signature_genes <- c(
  "LGR5", "ASCL2", "OLFM4", "SMOC2", "RGMB", 
  "EPHB2", "SOX9", "BMI1", "MSI1", "PROM1"
)

isc_genes_available <- isc_signature_genes[isc_signature_genes %in% rownames(epithelial_cells)]

epithelial_cells <- AddModuleScore(
  object = epithelial_cells,
  features = list(isc_genes_available),
  name = "ISC_Stemness_Score",
  assay = "RNA"
)

# Fix module score column name
colnames(epithelial_cells@meta.data)[colnames(epithelial_cells@meta.data) == "ISC_Stemness_Score1"] <- "ISC_Stemness_Score"

# ---- A2. Subset Border (B) vs. Core (T) Epithelial Cells ----
border_vs_core <- subset(epithelial_cells, subset = tissue_site %in% c("B", "T"))

border_vs_core$tissue_site <- factor(
  border_vs_core$tissue_site, 
  levels = c("B", "T"),
  labels = c("Invasive Border (B)", "Tumor Core (T)")
)

# ---- A3. Statistical Violin Plot (Stemness: Border vs. Core) ----
p_stemness_comp <- ggplot(
  border_vs_core@meta.data, 
  aes(x = tissue_site, y = ISC_Stemness_Score, fill = tissue_site)
) +
  geom_violin(trim = FALSE, alpha = 0.6) +
  geom_boxplot(width = 0.15, outlier.shape = NA, color = "black") +
  scale_fill_manual(values = c("Invasive Border (B)" = "#E41A1C", "Tumor Core (T)" = "#377EB8")) +
  stat_compare_means(method = "wilcox.test", label = "p.signif", label.x = 1.5) +
  labs(
    title = "Intestinal Stem Cell (ISC) Signature: Border vs. Core",
    x = "Tissue Site",
    y = "ISC Module Score"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 13),
    legend.position = "none",
    axis.text = element_text(size = 11, color = "black"),
    axis.title = element_text(size = 12, face = "bold")
  )

print(p_stemness_comp)
ggsave(
  filename = file.path(figures_dir, "ISC_Stemness_Border_vs_Core.png"),
  plot = p_stemness_comp, width = 6, height = 6, dpi = 300
)

# ---- A4. UMAP Feature Plot for Stemness Score ----
reduction_name <- ifelse("umap_after" %in% names(epithelial_cells@reductions), "umap_after", "umapafter")

p_stemness_umap <- FeaturePlot(
  epithelial_cells, 
  features = "ISC_Stemness_Score",
  reduction = reduction_name,
  cols = c("lightgrey", "yellow", "red")
) +
  ggtitle("ISC Stemness Score Expression (Epithelial Cells)") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 13))

print(p_stemness_umap)
ggsave(
  filename = file.path(figures_dir, "ISC_Stemness_UMAP.png"),
  plot = p_stemness_umap, width = 8, height = 7, dpi = 300
)

# ==============================================================================
# PART B: Differential Expression Analysis (Border vs. Core) & Volcano Plot
# ==============================================================================

# ---- B1. Perform Differential Expression Analysis ----
epithelial_bt <- subset(
  seurat_obj, 
  subset = canonical_annotation %in% c("Epithelial cells", "Epithelial (Cycling)") & tissue_site %in% c("B", "T")
)

DefaultAssay(epithelial_bt) <- "RNA"
Idents(epithelial_bt) <- "tissue_site"

deg_border_vs_core <- FindMarkers(
  epithelial_bt,
  ident.1 = "B", # Invasive Border
  ident.2 = "T", # Tumor Core
  logfc.threshold = 0.25,
  min.pct = 0.1,
  test.use = "wilcox"
)

# Classify DEG status
deg_border_vs_core$gene <- rownames(deg_border_vs_core)
deg_border_vs_core <- deg_border_vs_core %>%
  mutate(
    diffexpressed = case_when(
      avg_log2FC > 0.5 & p_val_adj < 0.05 ~ "Up in Border",
      avg_log2FC < -0.5 & p_val_adj < 0.05 ~ "Up in Core",
      TRUE ~ "Not Significant"
    )
  )

# Export complete DEG results
write.csv(
  deg_border_vs_core,
  file.path(results_dir, "10_Epithelial_DEGs_Border_vs_Core.csv"),
  row.names = FALSE
)

# ---- B2. Generate Publication-Ready Volcano Plot ----
top_border_genes <- deg_border_vs_core %>%
  filter(diffexpressed == "Up in Border") %>%
  arrange(p_val_adj) %>%
  head(15)

p_volcano <- ggplot(deg_border_vs_core, aes(x = avg_log2FC, y = -log10(p_val_adj), color = diffexpressed)) +
  geom_point(alpha = 0.7, size = 1.8) +
  scale_color_manual(
    values = c("Up in Border" = "#E41A1C", "Up in Core" = "#377EB8", "Not Significant" = "grey70")
  ) +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", color = "black", linewidth = 0.4) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black", linewidth = 0.4) +
  geom_text(
    data = top_border_genes,
    aes(label = gene),
    size = 3.2,
    vjust = -0.5,
    fontface = "bold.italic",
    show.legend = FALSE
  ) +
  labs(
    title = "Epithelial DEGs: Invasive Border (B) vs. Tumor Core (T)",
    x = "log2 Fold Change (Border / Core)",
    y = "-log10 Adjusted P-value",
    color = "Status"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 13),
    axis.text = element_text(size = 10, color = "black"),
    axis.title = element_text(size = 11, face = "bold"),
    legend.position = "top"
  )

print(p_volcano)
ggsave(
  filename = file.path(figures_dir, "Volcano_Epithelial_Border_vs_Core.png"),
  plot = p_volcano, width = 8, height = 7, dpi = 300
)

# ==============================================================================
# End of Downstream Script
# ==============================================================================