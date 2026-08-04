#--------------IMMUNE AGGRESSIVENESS--------------
library(dplyr)
library(tidyr)
library(ggpubr)
library(tibble)
# Define immune & EMT genes with annotations
immune_emt_genes <- tribble(
  ~SYMBOL,      ~Category,
  "CD274",      "Immune Checkpoint (PD-L1)",
  "CTLA4",      "Immune Checkpoint",
  "PDCD1",      "Immune Checkpoint (PD-1)",
  "LAG3",       "Immune Checkpoint",
  "TIGIT",      "Immune Checkpoint",
  "MMP9",       "EMT/Aggressiveness",
  "TWIST1",     "EMT/Aggressiveness",
  "SNAI1",      "EMT/Aggressiveness",
  "SNAI2",      "EMT/Aggressiveness",
  "CD44",       "EMT/Aggressiveness"
)

# Use expr_final (already clean, gene symbols as rownames)
gene_subset <- intersect(immune_emt_genes$SYMBOL, rownames(expr_final))

# Extract relevant expression
expr_df <- expr_final[gene_subset, ] %>%
  as.data.frame() %>%
  rownames_to_column("SYMBOL") %>%
  pivot_longer(cols = -SYMBOL, names_to = "sample", values_to = "expression")

# Add LOY_group and gene category
expr_df <- expr_df %>%
  left_join(male_combined_filtered[, c("sample_dot", "LOY_group_cutoff")],
            by = c("sample" = "sample_dot")) %>%
  left_join(immune_emt_genes, by = "SYMBOL")

# Check NA presence
table(is.na(expr_df$LOY_group_cutoff))

# Plot boxplots
ggboxplot(expr_df, x = "LOY_group_cutoff", y = "expression", fill = "LOY_group_cutoff",
          facet.by = "Category", color = "LOY_group_cutoff",
          palette = c("Y_low" = "#1f77b4", "Y_high" = "#ff7f0e"),
          add = "jitter", legend = "none") +
  stat_compare_means(method = "wilcox.test", label = "p.format") +
  theme_minimal(base_size = 13) +
  labs(title = "Expression of Immune & EMT Genes by LOY Group",
       x = "LOY Group", y = "Expression (TPM)", fill = "LOY Group")

##Boxplot of selected immune/EMT gene expression across LOY groups

library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)

# Define genes of interest (adjust names if needed)
genes_of_interest <- c("MMP9", "MMP2", "TWIST1", "SNAI1", "CD44", "CD274", "CTLA4", "PDCD1","PD-L1")

# Ensure genes are present in your expr_final (which has HGNC symbols as rownames)
valid_genes <- genes_of_interest[genes_of_interest %in% rownames(expr_final)]

# Subset and reshape expression matrix
exp_subset <- expr_final[valid_genes, , drop = FALSE]

exp_long <- as.data.frame(t(exp_subset)) %>%
  rownames_to_column("sample") %>%
  pivot_longer(cols = -sample, names_to = "gene", values_to = "expression")

# Merge LOY group info from `male_combined_filtered`
exp_long <- left_join(exp_long, male_combined_filtered[, c("sample_dot", "LOY_group_cutoff")],
                      by = c("sample" = "sample_dot"))

# Plot expression boxplot
ggplot(exp_long, aes(x = LOY_group_cutoff, y = expression, fill = LOY_group_cutoff)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.8) +
  geom_jitter(width = 0.2, alpha = 0.5) +
  facet_wrap(~ gene, scales = "free_y") +
  stat_compare_means(method = "wilcox.test", label = "p.format", label.y.npc = "top") +
  scale_fill_manual(values = c("Y_low" = "#1f77b4", "Y_high" = "#ff7f0e")) +
  theme_minimal(base_size = 12) +
  labs(
    title = "Expression of EMT/Immune Genes by LOY Group",
    x = "LOY Group",
    y = "Expression Level",
    fill = "LOY Group"
  )

## BARPLOT FOR IMMUNE EMT GENE
library(ggplot2)
library(dplyr)

# Define genes of interest for DEG check
immune_aggressive_genes <- c("PDCD1", "CD274", "CTLA4", "LAG3", "TIGIT", "HAVCR2",
                             "MMP9", "MMP2", "TWIST1", "SNAI1", "ZEB1", "CD44")

# Add gene symbols as a column if missing
deg_results$GeneID <- rownames(deg_results)

# Subset DEGs of interest
deg_interest <- deg_results %>%
  filter(GeneID %in% immune_aggressive_genes)

# Stop if none found (optional)
if (nrow(deg_interest) == 0) {
  stop("No immune/EMT genes found in DEG results.")
}

# Plot barplot: version 1 (logFC > 0 = up in Y_low)
ggplot(deg_interest, aes(x = reorder(GeneID, logFC), y = logFC, fill = logFC > 0)) +
  geom_col(width = 0.6) +
  geom_text(aes(label = paste0("p=", signif(adj.P.Val, 2))),
            hjust = ifelse(deg_interest$logFC > 0, -0.1, 1.1),
            size = 3.2, color = "black") +
  coord_flip() +
  scale_fill_manual(
    values = c("TRUE" = "#1b9e77", "FALSE" = "#d95f02"),
    labels = c("Down in Y_low", "Up in Y_low")
  ) +
  theme_minimal(base_size = 13) +
  labs(
    title = "Expression of Immune & EMT Genes (Y_low vs Y_high)",
    x = "Gene Symbol",
    y = "log2 Fold Change",
    fill = "Expression in Y_low"
  )




#-----------immune infiltration by xCELL method------------
# Install devtools if not already installed
if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools")
}
devtools::install_github('dviraran/xCell')
library(xCell)
xcell_scores <- xCellAnalysis(expr_final)
xcell_df <- as.data.frame(t(xcell_scores))
xcell_df$sample <- rownames(xcell_df)
xcell_df$sample_dash <- gsub("\\.", "-", xcell_df$sample)

# Merge with your LOY classification
xcell_merged <- left_join(xcell_df, male_combined, by = "sample_dash")
library(tidyr)
library(dplyr)
library(tidyr)

# Keep only numeric columns + LOY_group
long_xcell <- xcell_merged %>%
  select(where(is.numeric), LOY_group) %>%  # keeps only numeric columns and LOY_group
  pivot_longer(cols = -LOY_group, names_to = "cell_type", values_to = "score")

library(dplyr)
library(ggplot2)
library(ggpubr)

# Step 1: Compute per-cell Wilcoxon test
cell_pvals <- long_xcell %>%
  group_by(cell_type) %>%
  summarise(p_value = wilcox.test(score ~ LOY_group)$p.value) %>%
  mutate(adj_p = p.adjust(p_value, method = "BH")) %>%
  arrange(adj_p)
# Step 2 (modified): Choose top N excluding LOY_score
top_cells <- cell_pvals %>%
  filter(cell_type != "LOY_score") %>%
  slice_head(n = 12) %>%
  pull(cell_type)
# Step 3: Filter
top_long_xcell <- filter(long_xcell, cell_type %in% top_cells)

# Step 4: Plot
ggplot(top_long_xcell, aes(x = LOY_group, y = score, fill = LOY_group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.8, width = 0.5) +
  geom_jitter(width = 0.2, size = 1.2, alpha = 0.5) +
  facet_wrap(~ cell_type, scales = "free_y", ncol = 4) +
  stat_compare_means(method = "wilcox.test", label = "p.format", size = 4) +
  scale_fill_manual(values = c("Y_low" = "#1f77b4", "Y_high" = "#ff7f0e")) +
  theme_minimal(base_size = 13) +
  theme(
    strip.text = element_text(face = "bold", size = 12),
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  labs(
    title = "Top Differential Immune Cells by LOY Group",
    x = "LOY Group",
    y = "xCell Score"
  )

install.packages("fmsb")
library(fmsb)
library(fmsb)
library(dplyr)
library(scales)

# Step 1: Compute group medians for top 10 immune cells
radar_cells <- cell_pvals %>%
  filter(cell_type != "LOY_score") %>%
  slice_head(n = 10) %>%
  pull(cell_type)

radar_data <- long_xcell %>%
  filter(cell_type %in% radar_cells) %>%
  group_by(LOY_group, cell_type) %>%
  summarise(median_score = median(score), .groups = "drop") %>%
  pivot_wider(names_from = cell_type, values_from = median_score)

# Step 2: Format for fmsb (min/max + data)
radar_plot_data <- rbind(
  rep(1, length(radar_cells)),  # Max values
  rep(0, length(radar_cells)),  # Min values
  radar_data[match(c("Y_high", "Y_low"), radar_data$LOY_group), radar_cells]
)
rownames(radar_plot_data) <- c("max", "min", "Y_high", "Y_low")

# Step 3: Radar plot
colors_border <- c("#FF7F0E", "#1F77B4")
radarchart(
  radar_plot_data,
  axistype = 1,
  pcol = colors_border, plwd = 3, plty = 1,
  cglcol = "grey", cglty = 1, axislabcol = "black", cglwd = 0.8,
  vlcex = 0.8,
  title = "Median xCell Scores by LOY Group"
)
legend("topright", legend = c("Y_high", "Y_low"), bty = "n", pch = 20, col = colors_border, text.col = "black", cex = 1.2)
dev.off()  # Reset current plotting device
radarchart(
  radar_plot_data,
  axistype = 1,
  pcol = colors_border, plwd = 3, plty = 1,
  cglcol = "grey", cglty = 1, axislabcol = "black", cglwd = 0.8,
  vlcex = 0.8,
  title = "Median xCell Scores by LOY Group"
)
legend("topright", legend = c("Y_high", "Y_low"), bty = "n", pch = 20, col = colors_border, text.col = "black", cex = 1.2)
# Get all cell types from xCell output
all_cells <- colnames(xcell_df)
# Manually define immune cell types based on xCell documentation or common knowledge
immune_cells <- c(
  "B-cells", "CD4+ T-cells", "CD8+ T-cells", "DC", "Macrophages",
  "Monocytes", "NK cells", "NKT", "Neutrophils", "Tregs",
  "Mast cells", "Th1 cells", "Th2 cells", "Memory B-cells", "Memory CD4+ T-cells",
  "Memory CD8+ T-cells", "Plasma cells"
)

# Find which immune cells are present in your data
immune_cells_present <- intersect(immune_cells, all_cells)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)

# Subset xcell_df for immune cells + sample info
immune_xcell_df <- xcell_df[, c("sample", immune_cells_present)]

# Pivot longer
long_immune <- immune_xcell_df %>%
  pivot_longer(cols = -sample, names_to = "cell_type", values_to = "score")

# Add LOY group info
long_immune <- left_join(long_immune, male_combined_filtered[, c("sample_dot", "LOY_group_cutoff")], 
                         by = c("sample" = "sample_dot"))

# Filter out any NA LOY groups if exist
long_immune <- filter(long_immune, !is.na(LOY_group_cutoff))

# Now plot boxplots for immune cells
ggplot(long_immune, aes(x = LOY_group_cutoff, y = score, fill = LOY_group_cutoff)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.8, width = 0.6) +
  geom_jitter(width = 0.2, size = 1.2, alpha = 0.6) +
  facet_wrap(~ cell_type, scales = "free_y", ncol = 4) +
  stat_compare_means(method = "wilcox.test", label = "p.format", size = 4) +
  scale_fill_manual(values = c("Y_low" = "#1f77b4", "Y_high" = "#ff7f0e")) +
  theme_minimal(base_size = 13) +
  theme(strip.text = element_text(face = "bold", size = 12),
        legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Immune Cell xCell Scores by LOY Group",
       x = "LOY Group",
       y = "xCell Score")
library(pheatmap)
library(dplyr)

# 1. Define immune cells (adjust as needed)
immune_cells <- c(
  "B-cells", "CD4+ T-cells", "CD8+ T-cells", "DC", "Macrophages",
  "Monocytes", "NK cells", "NKT", "Neutrophils", "Tregs",
  "Mast cells", "Th1 cells", "Th2 cells", "Memory B-cells", "Memory CD4+ T-cells",
  "Memory CD8+ T-cells", "Plasma cells"
)
immune_cells_present <- intersect(immune_cells, colnames(xcell_df))

# 2. Subset xCell immune cell scores (samples x immune cells)
expr_sub <- xcell_df[, immune_cells_present]
rownames(expr_sub) <- xcell_df$sample  # samples as rows

# 3. Transpose so rows = immune cells, columns = samples
expr_sub_t <- t(expr_sub)

# 4. Prepare LOY group annotation for columns
annotation_col <- male_combined_filtered %>%
  filter(sample_dot %in% colnames(expr_sub_t)) %>%
  select(sample_dot, LOY_group_cutoff) %>%
  distinct()

# Order samples so Y_low first, then Y_high
annotation_col <- annotation_col %>%
  arrange(LOY_group_cutoff)

# Reorder columns in expr_sub_t to match annotation
expr_sub_t <- expr_sub_t[, annotation_col$sample_dot]

rownames(annotation_col) <- annotation_col$sample_dot
annotation_col <- annotation_col[, "LOY_group_cutoff", drop=FALSE]

# 5. Color palette for LOY group annotation
ann_colors <- list(
  LOY_group_cutoff = c("Y_low" = "#1f77b4", "Y_high" = "#ff7f0e")
)

# 6. Scale rows (immune cells) for heatmap clarity
expr_scaled <- t(scale(t(expr_sub_t)))

# 7. Draw heatmap
pheatmap(expr_scaled,
         cluster_rows = TRUE,       # cluster immune cells
         cluster_cols = FALSE,      # no clustering columns to keep groups separate
         annotation_col = annotation_col,
         annotation_colors = ann_colors,
         show_colnames = FALSE,
         fontsize_row = 10,
         main = "Heatmap of Immune Cell xCell Scores by LOY Group")

library(ggridges)
library(ggplot2)

ggplot(long_immune, aes(x = score, y = cell_type, fill = LOY_group_cutoff)) +
  geom_density_ridges(alpha = 0.7, scale = 1.2) +
  scale_fill_manual(values = c("Y_low" = "#1f77b4", "Y_high" = "#ff7f0e")) +
  theme_minimal(base_size = 14) +
  labs(title = "Ridgeline Plot of Immune Cell Scores by LOY Group",
       x = "xCell Score", y = "Immune Cell Type", fill = "LOY Group")
immune_summary <- long_immune %>%
  group_by(cell_type, LOY_group_cutoff) %>%
  summarise(median_score = median(score), .groups = "drop")

immune_pvals <- long_immune %>%
  group_by(cell_type) %>%
  summarise(p_value = wilcox.test(score ~ LOY_group_cutoff)$p.value) %>%
  mutate(signif = p.adjust(p_value, method = "BH") < 0.05)

bubble_data <- left_join(immune_summary, immune_pvals, by = "cell_type")

ggplot(bubble_data, aes(x = cell_type, y = LOY_group_cutoff, size = median_score, color = -log10(p_value))) +
  geom_point(alpha = 0.8) +
  scale_color_viridis_c() +
  coord_flip() +
  labs(title = "Bubble Plot of Immune Cell Scores and Significance",
       x = "Immune Cell Type", y = "LOY Group", size = "Median Score", color = "-log10(p-value)") +
  theme_minimal()

install.packages("GGally")
library(GGally)

# Prepare data: samples + immune cell scores + LOY group
immune_data <- xcell_df[, c("sample", immune_cells_present)]
immune_data <- left_join(immune_data, male_combined_filtered[, c("sample_dot", "LOY_group_cutoff")],
                         by = c("sample" = "sample_dot"))
immune_data <- na.omit(immune_data)

ggparcoord(data = immune_data, columns = 2:(1 + length(immune_cells_present)),
           groupColumn = "LOY_group_cutoff",
           scale = "uniminmax",
           showPoints = TRUE,
           alphaLines = 0.5) +
  scale_color_manual(values = c("Y_low" = "#1f77b4", "Y_high" = "#ff7f0e")) +
  theme_minimal() +
  labs(title = "Parallel Coordinate Plot of Immune Cell Profiles by LOY Group",
       x = "Immune Cell Types", y = "Scaled xCell Score")


#-------------TME-----------
immune_cells <- c("CD4+ T-cells", "CD8+ T-cells", "Tregs", "Th1 cells", "Th2 cells", "B-cells", "Plasma cells", 
                  "NK cells", "Macrophages M1", "Macrophages M2", "DC", "aDC", "cDC", "Neutrophils")
stromal_cells <- c("Fibroblasts", "Endothelial cells", "Pericytes", "MSC")
tme_scores <- c("ImmuneScore", "StromaScore", "MicroenvironmentScore")

# Immune subset
immune_data <- long_xcell %>% filter(cell_type %in% immune_cells)
stromal_data <- long_xcell %>% filter(cell_type %in% stromal_cells)
tme_score_data <- long_xcell %>% filter(cell_type %in% tme_scores)
# Immune
ggplot(immune_data, aes(x = LOY_group, y = score, fill = LOY_group)) +
  geom_boxplot(outlier.size = 0.4) +
  stat_compare_means(method = "wilcox.test", label = "p.format", size = 3) +
  facet_wrap(~ cell_type, scales = "free_y", ncol = 4) +
  scale_fill_manual(values = c("Y_low" = "#1f77b4", "Y_high" = "#ff7f0e")) +
  labs(title = "Immune Cell Enrichment by LOY Group") +
  theme_bw(base_size = 11)

# Stromal
ggplot(stromal_data, aes(x = LOY_group, y = score, fill = LOY_group)) +
  geom_boxplot(outlier.size = 0.4) +
  stat_compare_means(method = "wilcox.test", label = "p.format", size = 3) +
  facet_wrap(~ cell_type, scales = "free_y", ncol = 2) +
  scale_fill_manual(values = c("Y_low" = "#1f77b4", "Y_high" = "#ff7f0e")) +
  labs(title = "Stromal Cell Enrichment by LOY Group") +
  theme_bw(base_size = 11)

# TME scores
ggplot(tme_score_data, aes(x = LOY_group, y = score, fill = LOY_group)) +
  geom_boxplot(outlier.size = 0.4) +
  stat_compare_means(method = "wilcox.test", label = "p.format", size = 4) +
  facet_wrap(~ cell_type, scales = "free_y") +
  scale_fill_manual(values = c("Y_low" = "#1f77b4", "Y_high" = "#ff7f0e")) +
  labs(title = "Tumor Microenvironment Scores by LOY Group") +
  theme_bw(base_size = 11)
library(dplyr)
library(tibble)

# Step 0: Define top TME cell columns
top_tme_cells <- c("CD4+ T-cells", "CD8+ T-cells", "Tregs", "Th1 cells", "Th2 cells", "B-cells", "Plasma cells", 
                   "NK cells", "Macrophages M1", "Macrophages M2", "DC", "aDC", "cDC", "Neutrophils",
                   "Fibroblasts", "Endothelial cells", "Pericytes", "MSC")

# Step 1: Subset matrix
tme_matrix <- xcell_merged %>%
  select(sample_dash, all_of(top_tme_cells), LOY_group) %>%
  column_to_rownames("sample_dash")


# 2. Extract annotation for rows (samples)
annotation_row <- data.frame(LOY_group = tme_matrix$LOY_group)
rownames(annotation_row) <- rownames(tme_matrix)

# 3. Drop LOY_group column from matrix
tme_matrix <- tme_matrix %>% select(-LOY_group)

# 4. Scale and plot heatmap
pheatmap(t(scale(t(tme_matrix))),
         annotation_row = annotation_row,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
         show_rownames = FALSE,
         clustering_method = "ward.D2")


# Load libraries
library(dplyr)
library(tibble)
library(pheatmap)

# Define TME cell types to include
top_tme_cells <- c("CD4+ T-cells", "CD8+ T-cells", "Tregs", "Th1 cells", "Th2 cells", "B-cells", "Plasma cells", 
                   "NK cells", "Macrophages M1", "Macrophages M2", "DC", "aDC", "cDC", "Neutrophils",
                   "Fibroblasts", "Endothelial cells", "Pericytes", "MSC")

# Prepare matrix
tme_matrix <- xcell_merged %>%
  select(sample_dash, all_of(top_tme_cells), LOY_group) %>%
  column_to_rownames("sample_dash")

# Extract annotation
annotation_row <- data.frame(LOY_group = tme_matrix$LOY_group)
rownames(annotation_row) <- rownames(tme_matrix)

# Remove LOY_group from matrix for heatmap
tme_matrix <- tme_matrix %>% select(-LOY_group)

# Scale
tme_matrix_scaled <- t(scale(t(tme_matrix)))

# Order samples by LOY group
ordered_samples <- annotation_row %>%
  arrange(LOY_group) %>%
  rownames()

# Reorder matrix and annotation
tme_matrix_scaled <- tme_matrix_scaled[ordered_samples, ]
annotation_row <- annotation_row[ordered_samples, , drop = FALSE]

# Add gaps between Y_high and Y_low
group_changes <- which(diff(as.numeric(factor(annotation_row$LOY_group))) != 0)

# Plot heatmap
pheatmap(tme_matrix_scaled,
         annotation_row = annotation_row,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
         show_rownames = FALSE,
         cluster_rows = TRUE,
         cluster_cols = FALSE,  # disable column clustering
         gaps_row = group_changes,  # visually split Y_high and Y_low
         fontsize_row = 7)

library(dplyr)
library(pheatmap)
library(tibble)

# 1. Subset the matrix (samples as rows, TME scores as columns)
tme_matrix <- xcell_merged %>%
  select(sample_dash, all_of(top_tme_cells), LOY_group) %>%
  column_to_rownames("sample_dash")

# 2. Extract annotation for rows (samples)
annotation_row <- data.frame(LOY_group = tme_matrix$LOY_group)
rownames(annotation_row) <- rownames(tme_matrix)

# 3. Order samples by LOY_group
ordered_samples <- annotation_row %>%
  arrange(LOY_group) %>%
  rownames()

# 4. Reorder matrix and annotation accordingly
tme_matrix <- tme_matrix[ordered_samples, ]
annotation_row <- annotation_row[ordered_samples, , drop = FALSE]

# 5. Drop LOY_group column from the matrix
tme_matrix <- tme_matrix %>% select(-LOY_group)

# 6. Scale and plot heatmap
pheatmap(
  t(scale(t(tme_matrix))),
  annotation_row = annotation_row,
  color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
  show_rownames = FALSE,
  clustering_method = "ward.D2",
  cluster_rows = FALSE  # Don't cluster rows, keep Y_high and Y_low grouped
)


library(tidyr)
immune_summary <- xcell_merged %>%
  select(all_of(immune_cells), LOY_group) %>%
  pivot_longer(cols = -LOY_group, names_to = "Cell_Type", values_to = "Score") %>%
  group_by(Cell_Type, LOY_group) %>%
  summarise(
    n = n(),
    mean_score = mean(Score, na.rm = TRUE),
    median_score = median(Score, na.rm = TRUE),
    sd = sd(Score, na.rm = TRUE)
  ) %>%
  arrange(Cell_Type, LOY_group)
library(ggplot2)

ggplot(immune_summary, aes(x = Cell_Type, y = mean_score, fill = LOY_group)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_minimal(base_size = 12) +
  coord_flip() +
  labs(title = "Mean Immune Cell Infiltration by LOY Group",
       x = "Immune Cell Type", y = "Mean Infiltration Score") +
  scale_fill_manual(values = c("Y_low" = "#1f77b4", "Y_high" = "#ff7f0e"))
write.csv(immune_summary, "Immune_Infiltration_LOY_Summary.csv", row.names = FALSE)

# Number of immune cells

library(dplyr)
library(tidyr)

immune_cells <- c("CD4+ T-cells", "CD8+ T-cells", "Tregs", "Th1 cells", "Th2 cells", 
                  "B-cells", "Plasma cells", "NK cells", "Macrophages M1", "Macrophages M2", 
                  "DC", "aDC", "cDC", "Neutrophils")

immune_summary <- xcell_merged %>%
  select(all_of(immune_cells), LOY_group) %>%
  pivot_longer(cols = -LOY_group, names_to = "Cell_Type", values_to = "Score") %>%
  group_by(Cell_Type, LOY_group) %>%
  summarise(
    n = n(),                           # Number of samples per group per cell type
    mean_score = mean(Score, na.rm=TRUE),
    median_score = median(Score, na.rm=TRUE),
    sd_score = sd(Score, na.rm=TRUE)
  ) %>%
  arrange(Cell_Type, LOY_group)

print(immune_summary)
write.csv(immune_summary, file = "immune_cell_abundance_by_LOY_group.csv", row.names = FALSE)
immune_counts <- xcell_merged %>%
  select(all_of(immune_cells), LOY_group) %>%
  group_by(LOY_group) %>%
  summarise(across(all_of(immune_cells), ~ sum(.x, na.rm = TRUE))) %>%
  pivot_longer(-LOY_group, names_to = "Cell_Type", values_to = "Estimated_Cell_Abundance")

print(immune_counts)
immune_positive_counts <- xcell_merged %>%
  select(all_of(immune_cells), LOY_group) %>%
  pivot_longer(cols = -LOY_group, names_to = "Cell_Type", values_to = "Score") %>%
  group_by(LOY_group, Cell_Type) %>%
  summarise(
    n_samples_positive = sum(Score > 0, na.rm = TRUE),
    total_samples = n()
  )

print(immune_positive_counts)

write.csv(immune_counts,"immune counts.csv")
write.csv(immune_positive_counts,"immune positive counts.csv")


xcell_df <- as.data.frame(t(xcell_scores))
xcell_df <- rownames_to_column(xcell_df, var = "sample")

subtype_markers <- xcell_df %>%
  select(sample, `Macrophages`, `CD8+ T-cells`, `NK cells`, `Tregs`, `aDC`, `B-cells`)
library(dplyr)

xcell_df <- xcell_df %>%
  mutate(
    immune_subtype = case_when(
      Macrophages > 0.4 & Tregs > 0.3 ~ "C1_Wound_Healing",
      `CD8+ T-cells` > 0.4 & `NK cells` > 0.3 ~ "C2_IFNG_Dominant",
      `B-cells` > 0.3 & `Th2 cells` > 0.3 ~ "C3_Inflammatory",
      Macrophages > 0.4 & `CD8+ T-cells` < 0.2 ~ "C4_Lymph_Depleted",
      ImmuneScore < 0.1 ~ "C5_Immuno_Quiet",
      Tregs > 0.4 ~ "C6_TGFb_Dominant",
      TRUE ~ "Unclassified"
    )
  )
table(xcell_df$immune_subtype)
library(ggplot2)

ggplot(xcell_df, aes(x = immune_subtype, fill = immune_subtype)) +
  geom_bar() +
  theme_minimal() +
  labs(title = "Immune Subtype Distribution (xCell-based)",
       x = "Subtype", y = "Number of Samples")


#--------------Cybersort------------
# Load immune infiltration results
cibersortx_results <- read.csv("CIBERSORTx_Job5_Results.csv", row.names = 1, check.names = FALSE)
# Convert to dash format to match metadata
cibersortx_results$sample_dash <- gsub("\\.", "-", rownames(cibersortx_results))
# Merge with LOY metadata
immune_df <- merge(cibersortx_results, male_combined, by = "sample_dash")

# Check dimensions
dim(immune_df)
library(tidyverse)
library(ggpubr)

# Step 1: Identify immune cell columns
immune_cols <- colnames(cibersortx_results)[1:22]  # adjust if needed

# Step 2: Merge CIBERSORTx results with LOY metadata
immune_df <- merge(cibersortx_results, male_combined, by = "sample_dash")

# Step 3: Reshape to long format
immune_long <- immune_df %>%
  select(all_of(c("sample_dash", immune_cols, "LOY_group_cutoff"))) %>%
  pivot_longer(
    cols = all_of(immune_cols),
    names_to = "Cell_Type",
    values_to = "Fraction"
  )
# Wilcoxon test by cell type
sig_cells <- immune_long %>%
  group_by(Cell_Type) %>%
  summarise(p_value = wilcox.test(Fraction ~ LOY_group_cutoff)$p.value) %>%
  mutate(p_adj = p.adjust(p_value, method = "fdr")) %>%
  filter(p_adj < 0.05)

# Filter long data to only significant cell types
sig_immune_long <- immune_long %>%
  filter(Cell_Type %in% sig_cells$Cell_Type)

# Plot
ggplot(sig_immune_long, aes(x = LOY_group_cutoff, y = Fraction, fill = LOY_group_cutoff)) +
  geom_boxplot(outlier.size = 0.5) +
  stat_compare_means(method = "wilcox.test", label = "p.format") +
  facet_wrap(~ Cell_Type, scales = "free", ncol = 3) +
  scale_fill_manual(values = c("Y_low" = "#1f77b4", "Y_high" = "#ff7f0e")) +
  theme_minimal() +
  labs(title = "Significantly Different Immune Cell Types by LOY Group",
       x = "", y = "Fraction") +
  theme(strip.text = element_text(size = 10))

library(ggplot2)
library(ggpubr)

# Plot all immune cells (faceted by cell type)
ggplot(immune_long, aes(x = LOY_group_cutoff, y = Fraction, fill = LOY_group_cutoff)) +
  geom_boxplot(outlier.size = 0.5) +
  stat_compare_means(method = "wilcox.test", label = "p.format", size = 3) +
  facet_wrap(~ Cell_Type, scales = "free_y", ncol = 4) +
  scale_fill_manual(values = c("Y_low" = "#1f77b4", "Y_high" = "#ff7f0e")) +
  theme_minimal(base_size = 12) +
  theme(
    strip.text = element_text(size = 10, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  ) +
  labs(
    title = "CIBERSORTx Immune Cell Fractions by LOY Group (All 22 Cell Types)",
    x = "LOY Group (Optimal Cutoff)",
    y = "Estimated Cell Fraction"
  )

ggplot(immune_long, aes(x = LOY_group_cutoff, y = Fraction, fill = LOY_group_cutoff)) +
  geom_boxplot(outlier.size = 0.5) +
  stat_compare_means(
    method = "wilcox.test",
    label = "p.format",
    size = 4,          # Make p-value text larger
    label.y.npc = "top", # Position at top of each facet
    fontface = "bold"  # Make it bold
  ) +
  facet_wrap(~ Cell_Type, scales = "free_y", ncol = 4) +
  scale_fill_manual(values = c("Y_low" = "#1f77b4", "Y_high" = "#ff7f0e")) +
  theme_minimal(base_size = 12) +
  theme(
    strip.text = element_text(size = 10, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  ) +
  labs(
    title = "CIBERSORTx Immune Cell Fractions by LOY Group (All 22 Cell Types)",
    x = "LOY Group (Optimal Cutoff)",
    y = "Estimated Cell Fraction"
  )
library(dplyr)
library(ggpubr)

# Compute p-values per cell type
pval_df <- immune_long %>%
  group_by(Cell_Type) %>%
  summarise(p.value = wilcox.test(Fraction ~ LOY_group_cutoff)$p.value) %>%
  mutate(label = paste0("p = ", signif(p.value, 2)))
# Join p-values to plotting data
immune_long <- immune_long %>%
  left_join(pval_df, by = "Cell_Type")
library(ggplot2)

# Define a fixed Y-position per cell type for p-value text
y_positions <- immune_long %>%
  group_by(Cell_Type) %>%
  summarise(y_pos = max(Fraction, na.rm = TRUE) + 0.02)

immune_long <- immune_long %>%
  left_join(y_positions, by = "Cell_Type")

# Plot
ggplot(immune_long, aes(x = LOY_group_cutoff, y = Fraction, fill = LOY_group_cutoff)) +
  geom_boxplot(outlier.size = 0.5) +
  geom_text(aes(y = y_pos, label = label), size = 3.5, fontface = "bold") +  # Clear p-values
  facet_wrap(~ Cell_Type, scales = "free_y", ncol = 4) +
  scale_fill_manual(values = c("Y_low" = "#1f77b4", "Y_high" = "#ff7f0e")) +
  theme_minimal(base_size = 12) +
  theme(
    strip.text = element_text(size = 10, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  ) +
  labs(
    title = "CIBERSORTx Immune Cell Fractions by LOY Group (All 22 Cell Types)",
    x = "LOY Group (Optimal Cutoff)",
    y = "Estimated Cell Fraction"
  )

ggplot(immune_long, aes(x = LOY_group_cutoff, y = Fraction, fill = LOY_group_cutoff)) +
  geom_boxplot(outlier.size = 0.5) +
  geom_text(aes(y = y_pos, label = label), size = 3.5, fontface = "bold") +  # Only this shows p-values
  facet_wrap(~ Cell_Type, scales = "free_y", ncol = 4) +
  scale_fill_manual(values = c("Y_low" = "#1f77b4", "Y_high" = "#ff7f0e")) +
  theme_minimal(base_size = 12) +
  theme(
    strip.text = element_text(size = 10, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  ) +
  labs(
    title = "CIBERSORTx Immune Cell Fractions by LOY Group (All 22 Cell Types)",
    x = "LOY Group (Optimal Cutoff)",
    y = "Estimated Cell Fraction"
  )
#--------------TME estimate-------------
if (!requireNamespace("estimate", quietly = TRUE)) {
  install.packages("estimate")
}

library(estimate)
library(tidyverse)

# expr_final: your expression matrix (genes × samples)
# Ensure gene symbols are unique
expr_estimate <- expr_final[!duplicated(rownames(expr_final)), ]
# Save in GCT format
gct_file <- "expression_for_estimate.gct"
write.table(cbind(NAME=rownames(expr_estimate), Description="NA", expr_estimate),
            file = gct_file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
# Add GCT file headers manually (required for ESTIMATE)
# Load file
lines <- readLines(gct_file)
n_genes <- nrow(expr_estimate)
n_samples <- ncol(expr_estimate)

# Insert required GCT header lines
header <- c("#1.2", paste(n_genes, n_samples, sep = "\t"))
writeLines(c(header, lines), gct_file)
library(estimate)

# Run ESTIMATE scoring
estimate_scores <- estimateScore(gct_file, output.ds = "estimate_scores.gct", platform = "illumina")
# Load results
estimate_result <- read.table("estimate_scores.gct", skip = 2, header = TRUE, sep = "\t", check.names = FALSE)
rownames(estimate_result) <- estimate_result$NAME
estimate_result <- estimate_result[, -c(1, 2)]  # remove NAME and Description columns
estimate_result <- t(estimate_result)  # transpose: now samples as rows
# Create dataframe
estimate_df <- as.data.frame(estimate_result)
estimate_df$sample <- rownames(estimate_df)

# Merge using sample ID
estimate_merged <- estimate_df %>%
  left_join(male_combined_filtered[, c("sample_dot", "LOY_group_cutoff")],
            by = c("sample" = "sample_dot"))
# Optional: rename LOY group column for easier plotting
estimate_merged <- estimate_merged %>%
  rename(LOY_group = LOY_group_cutoff)
library(ggpubr)
library(tidyr)

# Convert to long format
estimate_long <- estimate_merged %>%
  pivot_longer(cols = c("StromalScore", "ImmuneScore", "ESTIMATEScore"),
               names_to = "ScoreType", values_to = "Score")

# Plot violin + boxplot + p-values
ggplot(estimate_long, aes(x = LOY_group, y = Score, fill = LOY_group)) +
  geom_violin(trim = FALSE, alpha = 0.8) +
  geom_boxplot(width = 0.1, outlier.shape = NA, alpha = 0.6) +
  geom_jitter(width = 0.15, size = 0.5, alpha = 0.4) +
  stat_compare_means(method = "wilcox.test", label = "p.format", size = 4) +
  facet_wrap(~ ScoreType, scales = "free_y") +
  scale_fill_manual(values = c("Y_high" = "#D73027", "Y_low" = "#4575B4")) +
  theme_minimal(base_size = 14) +
  labs(title = "ESTIMATE Scores by LOY Group", x = "LOY Group", y = "Score") +
  theme(legend.position = "none")
