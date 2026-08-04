

# Install once if needed
install.packages("Seurat")
install.packages("patchwork")
install.packages("Matrix")
BiocManager::install("GEOquery")
library(Seurat)

#Download and Read Raw scRNA-seq Data for pRCC Samples

library(GEOquery)

# Download GEO metadata
gse152938 <- getGEO("GSE152938", GSEMatrix = TRUE)
gse_metadata <- pData(gse152938[[1]])

# View sample titles to find pRCC
View(gse_metadata[, c("title", "geo_accession")])




library(Seurat)
library(Matrix)
library(patchwork)


data <- Read10X(data.dir = "E:/LOY_KIRP/PRCC/")

seurat_obj <- CreateSeuratObject(counts = data, project = "KIRP_pRCC", min.cells = 3, min.features = 200)

#Quality Control & Normalization

# Mitochondrial % (human genes start with MT-)
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")

# Filter poor-quality cells
seurat_obj <- subset(seurat_obj, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 10)

# Normalize
seurat_obj <- NormalizeData(seurat_obj)
seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)
seurat_obj <- ScaleData(seurat_obj)


#Dimensionality Reduction & Clustering
seurat_obj <- RunPCA(seurat_obj)
ElbowPlot(seurat_obj)

seurat_obj <- FindNeighbors(seurat_obj, dims = 1:20)
seurat_obj <- FindClusters(seurat_obj, resolution = 0.5)
seurat_obj <- RunUMAP(seurat_obj, dims = 1:20)

DimPlot(seurat_obj, reduction = "umap", label = TRUE)









##Reference 

library(Seurat)
library(dplyr)
# Replace with your actual normal sample directory
normal_path <- "E:/LOY_KIRP/Normal Reference/"

# Read the 10X normal sample
normal_seurat <- Read10X(data.dir = normal_path) %>%
  CreateSeuratObject(project = "Normal_Ref")

# Normalize and scale
normal_seurat <- NormalizeData(normal_seurat)
normal_seurat <- ScaleData(normal_seurat)

# Define Y-chromosome genes used for LOY scoring
loy_genes <- c("RPS4Y1","ZFY","CDHL3","TBL1Y","USP9Y",
               "DDX3Y","UTY","TMSB4Y","NLGN4Y","HSFY2","KDM5D",
               "EIF1AY","RBMY1A1","PRY2")
loy_genes <- loy_genes[loy_genes %in% rownames(normal_seurat)]

# Calculate LOY module score
normal_seurat <- AddModuleScore(normal_seurat, features = list(loy_genes), name = "LOY")


LOY_cutoff <- median(normal_seurat$LOY1, na.rm = TRUE)
print(LOY_cutoff)






##prcc
# Define Y-chromosome genes used for LOY scoring
loy_genes <- c("RPS4Y1","ZFY","CDHL3","TBL1Y","USP9Y",
               "DDX3Y","UTY","TMSB4Y","NLGN4Y","HSFY2","KDM5D",
               "EIF1AY","RBMY1A1","PRY2")
loy_genes <- loy_genes[loy_genes %in% rownames(seurat_obj)]

# Add LOY module score (Y-linked gene expression)
seurat_obj <- AddModuleScore(seurat_obj, features = list(loy_genes), name = "LOY")

# Classify LOY status based on LOY cutoff from reference

seurat_obj$LOY_status <- ifelse(seurat_obj$LOY1 < LOY_cutoff, "LOY", "Y_present")
table(seurat_obj$LOY_status)




##Angle from all cell types
library(SingleR)
library(celldex)

# Reference covering all major human cell types
ref <- celldex::HumanPrimaryCellAtlasData()

# Use log-normalized expression data
data_matrix <- GetAssayData(seurat_obj, slot = "data")

# Run SingleR annotation
singleR_results <- SingleR(test = data_matrix,
                           ref = ref,
                           labels = ref$label.main)

# Add to metadata
seurat_obj$SingleR_annot <- singleR_results$labels


# Already done
seurat_obj$LOY_status <- ifelse(seurat_obj$LOY1 < LOY_cutoff, "LOY", "Y_present")

library(ggplot2)
library(dplyr)
LOY_plot_data <- seurat_obj@meta.data %>%
  group_by(SingleR_annot, LOY_status) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(SingleR_annot) %>%
  mutate(freq = n / sum(n))
ggplot(LOY_plot_data, aes(x = SingleR_annot, y = freq, fill = LOY_status)) +
  geom_bar(stat = "identity", position = "fill") +
  labs(
    title = "Proportion of LOY Cells by Broad Cell Type",
    x = "Cell Type",
    y = "Fraction"
  ) +
  scale_fill_manual(values = c("LOY" = "firebrick", "Y_present" = "steelblue")) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

##Run hoina
LOY_plot_data <- seurat_obj@meta.data %>%
  group_by(immune_simple, LOY_status) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(immune_simple) %>%
  mutate(freq = n / sum(n))

ggplot(LOY_plot_data, aes(x = immune_simple, y = freq, fill = LOY_status)) +
  geom_bar(stat = "identity", position = "fill") +
  labs(
    title = "Proportion of LOY Cells by Broad Cell Type",
    x = "Cell Type",
    y = "Fraction"
  ) +
  scale_fill_manual(values = c("LOY" = "firebrick", "Y_present" = "steelblue")) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))



seurat_obj$celltype_broad <- dplyr::case_when(
  grepl("T cell", seurat_obj$SingleR_annot, ignore.case = TRUE) ~ "T",
  grepl("B cell", seurat_obj$SingleR_annot, ignore.case = TRUE) ~ "B",
  grepl("NK", seurat_obj$SingleR_annot, ignore.case = TRUE) ~ "NK",
  grepl("Monocyte|Macrophage", seurat_obj$SingleR_annot, ignore.case = TRUE) ~ "Monocyte",
  grepl("Dendritic", seurat_obj$SingleR_annot, ignore.case = TRUE) ~ "DC",
  grepl("Epithelial", seurat_obj$SingleR_annot, ignore.case = TRUE) ~ "Epithelial",
  grepl("Fibroblast", seurat_obj$SingleR_annot, ignore.case = TRUE) ~ "Fibroblast",
  TRUE ~ "Other"
)
library(dplyr)
library(ggplot2)

LOY_plot_data <- seurat_obj@meta.data %>%
  group_by(celltype_broad, LOY_status) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(celltype_broad) %>%
  mutate(freq = n / sum(n))

ggplot(LOY_plot_data, aes(x = celltype_broad, y = freq, fill = LOY_status)) +
  geom_bar(stat = "identity", position = "fill") +
  labs(
    title = "Proportion of LOY Cells per Cell Type",
    y = "Fraction",
    x = "Cell Type"
  ) +
  scale_fill_manual(values = c("LOY" = "firebrick", "Y_present" = "steelblue")) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))






seurat_obj$meta_celltype <- case_when(
  grepl("T cell|B cell|Monocyte|NK|Dendritic", seurat_obj$SingleR_annot, ignore.case = TRUE) ~ "Immune",
  grepl("Epithelial|Proximal tubule|Renal", seurat_obj$SingleR_annot, ignore.case = TRUE) ~ "Epithelial",
  grepl("Fibroblast|Endothelial|Stroma", seurat_obj$SingleR_annot, ignore.case = TRUE) ~ "Stromal",
  TRUE ~ "Other"
)
library(dplyr)
library(ggplot2)

LOY_plot_data <- seurat_obj@meta.data %>%
  group_by(meta_celltype, LOY_status) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(meta_celltype) %>%
  mutate(freq = n / sum(n))

ggplot(LOY_plot_data, aes(x = meta_celltype, y = freq, fill = LOY_status)) +
  geom_bar(stat = "identity", position = "fill") +
  labs(
    title = "Proportion of LOY Cells by Major Cell Type",
    x = "Broad Cell Type",
    y = "Fraction"
  ) +
  scale_fill_manual(values = c("LOY" = "firebrick", "Y_present" = "steelblue")) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))




# ANGLE from immune cells

# Install once
BiocManager::install("SingleR")
BiocManager::install("celldex")
install.packages("Seurat")
install.packages("SeuratDisk")  # Optional, if you want to save as h5Seurat

library(SingleR)
library(celldex)
library(Seurat)
library(SummarizedExperiment)


ref <- celldex::MonacoImmuneData()

# Use your normalized data
data_matrix <- GetAssayData(seurat_obj, slot = "data")  # log-normalized expression


# Annotate with SingleR
singleR_results <- SingleR(test = data_matrix,
                           ref = ref,
                           labels = ref$label.main)


# Add annotation to metadata
seurat_obj$SingleR_annot<- singleR_results$labels

# Visualize on UMAP

DimPlot(
  seurat_obj,
  group.by = "SingleR_annot",
  label = TRUE,
  repel = TRUE
) +
  ggtitle("Immune Cells Within the TME") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))



# Recode SingleR labels into immune_simple categories
seurat_obj$immune_simple <- dplyr::recode(seurat_obj$SingleR_annot,
                                          "CD4+ T cells" = "CD4+ T",
                                          "CD8+ T cells" = "CD8+ T",
                                          "Naive B cells" = "B",
                                          "Memory B cells" = "B",
                                          "NK cells" = "NK",
                                          "Monocytes" = "Monocyte",
                                          "Dendritic cells" = "DC",
                                          .default = "Other"
)



# Install if you haven't already
install.packages("dplyr")  # Or install.packages("tidyverse")

# Load it
library(dplyr)

LOY_plot_data <- seurat_obj@meta.data %>%
  group_by(immune_simple, LOY_status) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(immune_simple) %>%
  mutate(freq = n / sum(n))

library(ggplot2)

ggplot(LOY_plot_data, aes(x = immune_simple, y = freq, fill = LOY_status)) +
  geom_bar(stat = "identity", position = "fill") +
  labs(
    title = "Proportion of LOY Cells per Immune Subtype",
    y = "Fraction",
    x = "Immune Cell Type"
  ) +
  scale_fill_manual(values = c("LOY" = "firebrick", "Y_present" = "steelblue")) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))








## angle from immune, epithelial and stromal cells

library(SingleR)
library(celldex)

ref <- HumanPrimaryCellAtlasData()
# Use normalized data (logcounts)
data_matrix <- GetAssayData(seurat_obj, slot = "data")

# Run SingleR again
singleR_results_hpca <- SingleR(test = data_matrix,
                                ref = ref,
                                labels = ref$label.main)

# Add new annotations to your Seurat object
seurat_obj$SingleR_annot <- singleR_results_hpca$labels
table(seurat_obj$SingleR_annot)


# Tabulate counts
annot_counts <- table(seurat_obj$SingleR_annot)

# Keep labels with >50 cells
valid_annots <- names(annot_counts[annot_counts > 50])

# Subset Seurat object
seurat_filtered <- subset(seurat_obj, subset = SingleR_annot %in% valid_annots)
# UMAP: Cell Types




DimPlot(
  seurat_filtered,
  group.by = "SingleR_annot",
  label = TRUE,
  repel = TRUE,
  pt.size = 1.2
) +
  ggtitle("UMAP of Annotated Cell Types (n > 50)") +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
    legend.text = element_text(face = "bold", size = 13),
    legend.title = element_text(face = "bold", size = 14),
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    strip.text = element_text(face = "bold")
  )

# UMAP: LOY status
DimPlot(seurat_filtered, group.by = "LOY_status", label = FALSE, pt.size = 1.2) +
  ggtitle("UMAP: LOY Status Distribution") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    legend.text = element_text(size = 12),
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank()
  )





# Subset only LOY-positive cells
seurat_loy <- subset(seurat_filtered, subset = LOY_status == "LOY")

DimPlot(
  seurat_loy,
  group.by = "SingleR_annot",
  label = TRUE,
  repel = TRUE,
  pt.size = 1.2
) +
  ggtitle("UMAP: Annotated Cell Types in LOY Cells") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    legend.text = element_text(size = 12),
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank()
  )






LOY_plot_data <- seurat_filtered@meta.data %>%
  group_by(SingleR_annot, LOY_status) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(SingleR_annot) %>%
  mutate(freq = n / sum(n))

ggplot(LOY_plot_data, aes(x = SingleR_annot, y = freq, fill = LOY_status)) +
  geom_bar(stat = "identity", position = "fill") +
  ylab("Proportion") +
  xlab("Cell Type") +
  ggtitle("LOY Proportion Across Cell Types") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))





FeaturePlot(
  seurat_filtered,
  features = "LOY1",
  pt.size = 1.5,
  label = FALSE
) +
  ggtitle("LOY Score Across Cells") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank()
  )













#Visualize Your 6 Prognostic Genes
# Your 6 prognostic genes
prog_genes <- c("WASF1", "CHMP1A", "MTHFD2", "PRC1", "MXD3", "CDKN2C")

# Check if they are present
prog_genes <- prog_genes[prog_genes %in% rownames(seurat_obj)]

# UMAP Feature Plots
FeaturePlot(seurat_obj, features = "WASF1", pt.size = 1)
# UMAP Feature Plots
FeaturePlot(seurat_obj, features = "CHMP1A", pt.size = 1)
FeaturePlot(seurat_obj, features = "MTHFD2", pt.size = 1)
FeaturePlot(seurat_obj, features = "PRC1", pt.size = 1)
FeaturePlot(seurat_obj, features = "MXD3", pt.size = 1)
FeaturePlot(seurat_obj, features = "CDKN2C", pt.size = 1)

# Violin plots by LOY status
VlnPlot(seurat_obj, features = prog_genes, group.by = "LOY_status", pt.size = 0.001)



# Use the recoded immune_simple column from SingleR
VlnPlot(seurat_obj,
        features = prog_genes,
        group.by = "immune_simple",
        split.by = "LOY_status",
        pt.size = 0.1) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# Optional: DotPlot for gene expression summary
DotPlot(seurat_obj, features = prog_genes, group.by = "immune_simple", split.by = "LOY_status") +
  RotatedAxis()








library(ggplot2)
library(cowplot)

# Dummy data for legend
legend_data <- data.frame(
  x = 1:2,
  y = 1,
  group = c("Y_present", "LOY")
)

# Match colors (gray = Y_present, blue = LOY)
loy_colors <- c("Y_present" = "gray", "LOY" = "blue")

# Dummy plot to generate legend
legend_plot <- ggplot(legend_data, aes(x, y, color = group)) +
  geom_point(size = 5) +
  scale_color_manual(values = loy_colors, name = "LOY Status") +
  theme_void() +
  theme(legend.position = "right")

# Extract legend
legend <- cowplot::get_legend(legend_plot)

# Main DotPlot
dot <- DotPlot(seurat_obj,
               features = prog_genes,
               group.by = "immune_simple",
               split.by = "LOY_status") +
  RotatedAxis()

# Combine plot and legend
combined_plot <- cowplot::plot_grid(dot, legend, rel_widths = c(1, 0.2))

# Add title to the whole plot
final_plot <- cowplot::plot_grid(
  ggdraw() + draw_label("Expression of Prognostic Genes across Immune Cell Types by LOY Status", fontface = "bold", size = 14),
  combined_plot,
  ncol = 1,
  rel_heights = c(0.1, 1)
)

# Print final plot
print(final_plot)

















DotPlot(
  seurat_obj,
  features = prog_genes,
  group.by = "immune_simple",
  split.by = "LOY_status"
) +
  RotatedAxis() +
  ggtitle("Prognostic Gene Expression Across Immune Cell Types by LOY Status") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y = element_text(size = 10),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10),
    strip.text = element_text(size = 12, face = "bold")
  )



library(Seurat)
library(ggplot2)
library(cowplot)

# 1. Create dummy data to generate custom LOY legend
legend_data <- data.frame(
  x = 1:2,
  y = 1,
  group = c("Y_present", "LOY")
)

# 2. Set the colors matching your DotPlot (gray = Y_present, blue = LOY)
loy_colors <- c("Y_present" = "gray", "LOY" = "blue")

# 3. Dummy plot just to get the legend
legend_plot <- ggplot(legend_data, aes(x, y, color = group)) +
  geom_point(size = 5) +
  scale_color_manual(values = loy_colors, name = "LOY Status") +
  theme_void() +
  theme(legend.position = "right")

# 4. Extract legend
legend <- cowplot::get_legend(legend_plot)

# 5. Main DotPlot with your theme and title
dot <- DotPlot(
  seurat_obj,
  features = prog_genes,
  group.by = "immune_simple",
  split.by = "LOY_status"
) +
  RotatedAxis() +
  ggtitle("Prognostic Gene Expression Across Immune Cell Types by LOY Status") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y = element_text(size = 10),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10),
    strip.text = element_text(size = 12, face = "bold")
  )

# 6. Combine plot and legend
combined_plot <- cowplot::plot_grid(dot, legend, rel_widths = c(1, 0.15))

# 7. Display final plot
print(combined_plot)
