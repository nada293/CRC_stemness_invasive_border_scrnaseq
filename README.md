# Stem-like Malignant Cell States at the Invasive Border in Colorectal Cancer (scRNA-seq Re-analysis)

## The Question
Do malignant epithelial cells at the tumor invasive border express a stronger intestinal stem cell signature than those at the tumor core in colorectal cancer?

## The Data
- **Source:** Lee HO, Hong Y, Etlioglu HE, Cho YB, et al. "Lineage-dependent gene expression programs influence the immune landscape of colorectal cancer." *Nature Genetics*, 2020.
- **Accession:** [GSE144735](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE144735) (KUL3 Belgian cohort)
- **Size:** 27,414 cells, 6 patients, 3 tissue types (Tumor Core / Tumor Border / Normal Mucosa)
- **Download:** Supplementary files on the GEO page above (raw UMI count matrix + annotation file). Not included in this repository.

## What We Did
1. **Preprocessing & Quality Control:** Standard Seurat workflow including cell filtering, normalization, high-variance gene selection, and dimensional reduction (PCA & UMAP).
2. **Cell Type Annotation:** Annotated major tumor microenvironment populations (Epithelial, Stromal, Immune) using canonical biological markers and differential expression.
3. **Stemness Module Scoring:** Scored malignant epithelial cells for an Intestinal Stem Cell (ISC) marker gene signature (`LGR5`, `OLFM4`, `ASCL2`, `SMOC2`, `EPHB2`, `SOX9`, `BMI1`, `MSI1`, `PROM1`, `RGMB`).
4. **Spatial Comparative Analysis:** Evaluated stemness signature differences between the Invasive Border (**B**) and Tumor Core (**T**) using Wilcoxon rank-sum testing ($p < 0.0001$).
5. **Differential Expression Analysis:** Identified top driver genes enriched at the invasive margin (highlighting **`OLFM4`** with $log2FC > 2.5$) and visualized results via volcano plotting.

## How to Run It
1. Clone this repository:
   ```bash
   git clone [https://github.com/YOUR_USERNAME/CRC_stemness_invasive_border_scrnaseq.git](https://github.com/YOUR_USERNAME/CRC_stemness_invasive_border_scrnaseq.git)
   cd CRC_stemness_invasive_border_scrnaseq
   
   Run the main processing and figures script:

source("scripts/01_visualization_figures1_4.R")
Run the downstream stemness analysis and volcano plot pipeline:

source("scripts/02_downstream_analysis.R")

Requirements
($\ge 4.0$)R Packages: Seurat, dplyr, ggplot2, ggpubr, ggrepel

Results
Elevated Invasive Border Stemness: Malignant epithelial cells at the tumor invasive border show a statistically significant upregulation of the ISC stemness program compared to the tumor core ($p < 0.0001$).Marker Upregulation: OLFM4 emerged as the predominant driver gene at the border, co-upregulated with microenvironmental markers such as PLCB4, PIGR, and SPINK1

Structure
CRC_stemness_invasive_border_scrnaseq/
├── scripts/
│   ├── 01_visualization_figures1_4.R
│   └── 02_downstream_analysis.R
├── figures/
│   ├── Canonical_Marker_annotation.png
│   ├── Annotated_Cell_types.png
│   ├── DEG_Based_Annotation.png
│   ├── Top_DEG.png
│   ├── ISC_Stemness_Border_vs_Core.png
│   ├── ISC_Stemness_UMAP.png
│   └── Volcano_Epithelial_Border_vs_Core.png
├── results/
│   └── 10_Epithelial_DEGs_Border_vs_Core.csv
└── README.md 


Team & Contributions
— Nada Qenawy - Lead Bioinformatician / Primary ContributorExecuted
core scRNA-seq workflow: Data loading, Normalization, Feature Selection, Scaling, PCA, and UMAP Clustering.
Led cell-type annotation, DEG validation, and marker analysis across all clusters.
Developed and executed the entire downstream analytics pipeline (02_downstream_analysis.R), including Intestinal Stem Cell (ISC) module scoring and statistical testing ($p < 0.0001$).
Generated all publication-ready visualizations (Figures 1–4, ISC Violin/UMAP plots, and Epithelial Border vs. Core Volcano Plot).
Managed repository structure, code organization, and final technical documentation.

Marwan Waleed — Repository Setup & QC
Initialized the GitHub repository structure.
Assisted with initial Quality Control (QC) parameters.

Mohamed Jassem — Annotation Support
Assisted in validating cell-type annotations and marker expression analysis.