# risk immune infiltration 
# Install if needed
if (!require("xCell")) BiocManager::install("xCell")
library(xCell)

# Use expression matrix with genes as rownames and samples as columns
# Make sure gene symbols match xCell expectations
xcell_results <- xCellAnalysis(expr_final)

# Transpose to sample × cell type format
xcell_df <- as.data.frame(t(xcell_results))
xcell_df$sample <- rownames(xcell_df)

# Merge with risk group
xcell_df <- left_join(xcell_df, surv_subset[, c("sample", "risk_group")], by = "sample")
library(dplyr)

immune_cells <- setdiff(colnames(xcell_df), c("sample", "risk_group"))

wilcox_results <- lapply(immune_cells, function(cell) {
  test <- wilcox.test(xcell_df[[cell]] ~ xcell_df$risk_group)
  data.frame(CellType = cell, 
             p.value = test$p.value, 
             median_High = median(xcell_df[xcell_df$risk_group == "High", cell], na.rm = TRUE),
             median_Low = median(xcell_df[xcell_df$risk_group == "Low", cell], na.rm = TRUE))
})

wilcox_df <- do.call(rbind, wilcox_results)
wilcox_df$adj.p <- p.adjust(wilcox_df$p.value, method = "fdr")
significant_cells <- wilcox_df %>% filter(adj.p < 0.05)

# View top results
head(significant_cells[order(significant_cells$adj.p), ])
library(pheatmap)

# Subset significant immune cell types
sig_cells <- significant_cells$CellType
heatmap_matrix <- xcell_df[, sig_cells]
rownames(heatmap_matrix) <- xcell_df$sample

# Optional: Z-score scale
heatmap_matrix_scaled <- t(scale(t(heatmap_matrix)))


# Create annotation for risk group
annotation_col <- data.frame(RiskGroup = xcell_df$risk_group)
rownames(annotation_col) <- xcell_df$sample

# Plot heatmap
pheatmap(heatmap_matrix_scaled,
         annotation_row = annotation_col,
         show_rownames = FALSE,
         clustering_method = "ward.D2",
         color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
         main = "Immune Cell Infiltration (Significant by Risk Group)")
# 1. Ensure you have the same row order between heatmap and annotation
annotation_row <- data.frame(RiskGroup = xcell_df$risk_group)
rownames(annotation_row) <- xcell_df$sample

# 2. Subset the matrix to only significant immune cells
sig_cells <- significant_cells$CellType
heatmap_matrix <- xcell_df[, sig_cells]
rownames(heatmap_matrix) <- xcell_df$sample

# 3. Order samples by risk group
ordered_samples <- xcell_df %>%
  arrange(risk_group) %>%
  pull(sample)


heatmap_matrix_ordered <- heatmap_matrix[ordered_samples, ]
annotation_ordered <- annotation_row[ordered_samples, , drop = FALSE]

# 4. Scale matrix (optional but preferred for heatmaps)
heatmap_matrix_scaled <- t(scale(t(heatmap_matrix_ordered)))

# 5. Plot heatmap with reordered samples
library(pheatmap)

pheatmap(heatmap_matrix_scaled,
         annotation_row = annotation_ordered,
         cluster_cols = TRUE,
         cluster_rows = FALSE,  # Don't cluster rows since we order by risk
         color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
         main = "Immune Cell Infiltration (Separated by Risk Group)",
         show_rownames = FALSE)
# Ensure both dataframes have consistent sample names
common_samples <- intersect(xcell_df$sample, surv_subset$sample)

# Subset and sort both dataframes by sample
xcell_matched <- xcell_df %>% filter(sample %in% common_samples) %>% arrange(sample)
surv_matched  <- surv_subset %>% filter(sample %in% common_samples) %>% arrange(sample)

# Double check alignment
stopifnot(all(xcell_matched$sample == surv_matched$sample))

# Now run correlation safely
risk_corr <- sapply(immune_cells, function(cell) {
  cor.test(xcell_matched[[cell]], surv_matched$risk_score)$estimate
})

# Convert to data frame
risk_corr_df <- data.frame(CellType = immune_cells, Correlation = risk_corr)
risk_corr_df <- risk_corr_df[order(abs(risk_corr_df$Correlation), decreasing = TRUE), ]

# View top results
head(risk_corr_df)
write.csv(risk_corr_df, "Immune_vs_RiskScore_Correlation.csv", row.names = FALSE)
# Store both estimate and p-value
risk_corr_list <- lapply(immune_cells, function(cell) {
  test <- cor.test(xcell_matched[[cell]], surv_matched$risk_score)
  data.frame(CellType = cell,
             Correlation = test$estimate,
             P_value = test$p.value)
})

# Combine into a single data frame
risk_corr_df <- do.call(rbind, risk_corr_list)

# Adjust p-values if desired
risk_corr_df$adj_p <- p.adjust(risk_corr_df$P_value, method = "fdr")

# Sort by absolute correlation
risk_corr_df <- risk_corr_df[order(abs(risk_corr_df$Correlation), decreasing = TRUE), ]

# View top hits
head(risk_corr_df)
library(ggplot2)

# Top 3 correlated cell types
top_cells <- head(risk_corr_df$CellType, 3)

# Plot
for (cell in top_cells) {
  p <- ggplot(data = xcell_matched, aes(x = .data[[cell]], y = surv_matched$risk_score)) +
    geom_point(alpha = 0.6) +
    geom_smooth(method = "lm", se = TRUE, color = "red") +
    labs(title = paste0("Correlation: ", cell),
         x = paste(cell, "Score"),
         y = "Risk Score") +
    theme_minimal()
  print(p)
}
library(ggplot2)

# Top 3 immune cell types most correlated with risk score
top_cells <- head(risk_corr_df$CellType, 3)

# Plot with correlation and p-value in the title
for (cell in top_cells) {
  x <- xcell_matched[[cell]]
  y <- surv_matched$risk_score
  
  # Perform correlation test
  test <- cor.test(x, y)
  cor_val <- round(test$estimate, 3)
  p_val <- signif(test$p.value, 3)
  
  # Create scatterplot with correlation + p-value
  p <- ggplot(data = xcell_matched, aes(x = .data[[cell]], y = y)) +
    geom_point(alpha = 0.6) +
    geom_smooth(method = "lm", se = TRUE, color = "red") +
    labs(title = paste0(cell, "\nCorrelation = ", cor_val, ", p = ", p_val),
         x = paste(cell, "Score"),
         y = "Risk Score") +
    theme_minimal()
  
  print(p)
}
# Significantly correlated immune cells (adjusted p < 0.05)
sig_cells <- risk_corr_df %>%
  filter(adj_p < 0.05) %>%
  arrange(adj_p) %>%
  pull(CellType)
library(ggpubr)

for (cell in sig_cells) {
  p <- ggplot(xcell_matched, aes(x = surv_matched$risk_group, y = .data[[cell]], fill = surv_matched$risk_group)) +
    geom_boxplot(alpha = 0.8) +
    geom_jitter(width = 0.2, size = 0.7, alpha = 0.5) +
    labs(title = paste("Cell Type:", cell),
         x = "Risk Group",
         y = paste(cell, "Score")) +
    theme_minimal() +
    theme(legend.position = "none") +
    scale_fill_manual(values = c("low" = "#1f77b4", "high" = "#d62728")) +
    stat_compare_means(method = "wilcox.test", label = "p.format")
  
  print(p)
}

library(tidyverse)
install.packages("effsize")
library(effsize)     # for cohen.d
library(ggpubr)      # for stat_compare_means
library(dplyr)
library(tidyr)

# Filter for significant cells
sig_cells <- risk_corr_df %>% filter(adj_p < 0.05) %>% pull(CellType)

# Create long-form data for ggplot
boxplot_data <- xcell_matched %>%
  select(all_of(sig_cells)) %>%
  mutate(risk_group = surv_matched$risk_group) %>%
  pivot_longer(cols = all_of(sig_cells), names_to = "CellType", values_to = "Score")
# Compute stats
stats_df <- boxplot_data %>%
  group_by(CellType) %>%
  summarise(
    p_value = wilcox.test(Score ~ risk_group)$p.value,
    cohens_d = cohen.d(Score ~ risk_group)$estimate
  ) %>%
  mutate(p_label = paste0("p = ", signif(p_value, 2),
                          "\nCohen's d = ", round(cohens_d, 2)))
# Join stats to plotting data
boxplot_data <- left_join(boxplot_data, stats_df, by = "CellType")

# Faceted boxplot with stats
ggplot(boxplot_data, aes(x = risk_group, y = Score, fill = risk_group)) +
  geom_boxplot(alpha = 0.8, outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 0.5, alpha = 0.5) +
  facet_wrap(~ CellType, scales = "free_y") +
  labs(title = "Significantly Correlated Immune Cell Types by Risk Group",
       x = "Risk Group", y = "Cell Score") +
  scale_fill_manual(values = c("low" = "#1f77b4", "high" = "#d62728")) +
  geom_text(data = stats_df, aes(x = 1.5, y = Inf, label = p_label),
            inherit.aes = FALSE, vjust = 1.5, size = 3) +
  theme_minimal() +
  theme(legend.position = "none",
        strip.text = element_text(size = 10))


# Match expression data with risk group samples
expr_risk <- expr_final[, surv_subset$sample]
install.packages("dplyr")
library(dplyr)

# Install and load xCell
if (!requireNamespace("xCell", quietly = TRUE)) {
  remotes::install_github("dviraran/xCell")
}
library(xCell)

# Run xCell (use log2 TPM + 1 values if not already normalized)
xcell_results <- xCellAnalysis(expr_risk)

# Transpose to get sample-wise data
xcell_df <- as.data.frame(t(xcell_results))
xcell_df$sample <- rownames(xcell_df)

# Merge with risk group info
immune_merged <- left_join(xcell_df, surv_subset[, c("sample", "risk_group")], by = "sample")

# Melt for visualization
library(reshape2)
long_df <- melt(immune_merged, id.vars = c("sample", "risk_group"), variable.name = "Cell_Type", value.name = "Score")

# Boxplot
library(ggplot2)
ggplot(long_df, aes(x = risk_group, y = Score, fill = risk_group)) +
  geom_boxplot() +
  facet_wrap(~ Cell_Type, scales = "free_y") +
  theme_bw() +
  labs(title = "Immune Cell Infiltration by Risk Group", y = "xCell Score")

# Calculate p-values first
pval_df <- long_df %>%
  group_by(Cell_Type) %>%
  summarise(p = wilcox.test(Score ~ risk_group)$p.value) %>%
  mutate(adj_p = p.adjust(p, method = "BH")) %>%
  filter(adj_p < 0.05)  # significant only

# Filter long_df to show only significant ones
long_df_sig <- long_df %>% filter(Cell_Type %in% pval_df$Cell_Type)
install.packages("ggplot2")   # Run this if not already installed
library(ggplot2)
library(ggpubr)  # Required for stat_compare_means()

ggplot(long_df_sig, aes(x = risk_group, y = Score, fill = risk_group)) +
  geom_boxplot(outlier.shape = NA) +
  stat_compare_means(method = "wilcox.test", label = "p.format") +
  facet_wrap(~ Cell_Type, scales = "free_y") +
  theme_bw() +
  labs(
    title = "Significant Immune Cell Differences by Risk Group",
    y = "xCell Score", x = "Risk Group"
  ) +
  scale_fill_manual(values = c("Low" = "#1f77b4", "High" = "#ff7f0e"))


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


#--------------RISK------------
# Add risk group to immune_long
# Convert rownames to a column for merging
cibersortx_results <- read.csv("CIBERSORTx_Job5_Results.csv", row.names = 1, check.names = FALSE)
cibersortx_df <- cibersortx_results %>% 
  tibble::rownames_to_column("sample_dot")

# Ensure sample IDs in risk data use dot format (if needed)
surv_subset$sample_dot <- gsub("-", ".", surv_subset$sample_dash)

# Merge immune data with risk group info
immune_risk_df <- merge(cibersortx_df, 
                        surv_subset[, c("sample_dot", "risk_group")], 
                        by = "sample_dot")

# Check merged data
head(immune_risk_df)
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggpubr)
immune_cells <- immune_risk_df %>%
  select(where(is.numeric)) %>%
  colnames()

immune_long <- immune_risk_df %>%
  pivot_longer(
    cols = all_of(immune_cells),
    names_to = "immune_cell",
    values_to = "fraction"
  )
sapply(immune_risk_df, class)
library(ggplot2)
library(ggpubr)
library(dplyr)

# Remove unwanted columns before pivoting (if not already done)
immune_risk_df_filtered <- immune_risk_df %>%
  select(-c(`P-value`, Correlation, RMSE, sample_dash))  # exclude non-immune columns

immune_cells <- immune_risk_df_filtered %>%
  select(-sample_dot, -risk_group) %>%
  colnames()

immune_long <- immune_risk_df_filtered %>%
  pivot_longer(
    cols = all_of(immune_cells),
    names_to = "immune_cell",
    values_to = "fraction"
  )

# Reorder immune_cell factor to control facet order if desired (optional)
immune_long$immune_cell <- factor(immune_long$immune_cell,
                                  levels = unique(immune_long$immune_cell))

# Plot with improved aesthetics
ggplot(immune_long, aes(x = risk_group, y = fraction, fill = risk_group)) +
  geom_boxplot(outlier.shape = NA, width = 0.7) +
  geom_jitter(width = 0.15, alpha = 0.6, size = 1.5) +
  facet_wrap(~ immune_cell, scales = "free_y", ncol = 3) +  # fewer columns, more rows
  stat_compare_means(
    method = "wilcox.test",
    label = "p.format",
    label.size = 4,
    vjust = 1.5
  ) +
  scale_fill_manual(values = c("Low" = "#4DBBD5", "High" = "#E64B35")) +
  theme_bw(base_size = 15) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
    strip.text = element_text(size = 14, face = "bold"),
    legend.position = "right",
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold")
  ) +
  labs(
    title = "Cibersort Immune Cell infiltrationby Risk Group",
    x = "Risk Group",
    y = "Fraction"
  )

#------------------TME------------
# Install if not already installed
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

# Merge with your clinical/risk data (surv_subset)
estimate_merged <- left_join(estimate_df, surv_subset[, c("sample", "risk_group")], by = "sample")
library(ggpubr)

# Melt the data
estimate_long <- estimate_merged %>%
  pivot_longer(cols = c("StromalScore", "ImmuneScore", "ESTIMATEScore"), 
               names_to = "ScoreType", values_to = "Score")

# Plot
ggplot(estimate_long, aes(x = risk_group, y = Score, fill = risk_group)) +
  geom_boxplot(alpha = 0.8, outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.5, size = 1) +
  stat_compare_means(method = "wilcox.test", label = "p.format") +
  facet_wrap(~ ScoreType, scales = "free_y") +
  scale_fill_manual(values = c("Low" = "#1f77b4", "High" = "#ff7f0e")) +
  theme_minimal() +
  labs(title = "TME Scores by Risk Group (ESTIMATE)", x = "Risk Group", y = "Score")
library(ggplot2)
library(ggpubr)
library(dplyr)

# Filter out NA values in risk_group
estimate_long_filtered <- estimate_long %>%
  filter(!is.na(risk_group))

# Create violin plot with p-values
ggplot(estimate_long_filtered, aes(x = risk_group, y = Score, fill = risk_group)) +
  geom_violin(trim = FALSE, alpha = 0.8) +
  geom_boxplot(width = 0.1, outlier.shape = NA, alpha = 0.6) +
  geom_jitter(width = 0.2, size = 1, alpha = 0.5) +
  stat_compare_means(method = "wilcox.test", label = "p.format", size = 4) +
  facet_wrap(~ ScoreType, scales = "free_y") +
  scale_fill_manual(values = c("Low" = "#1f77b4", "High" = "#ff7f0e")) +
  theme_minimal(base_size = 14) +
  labs(
    title = "TME Scores by Risk Group (ESTIMATE, Violin Plot)",
    x = "Risk Group", y = "Score"
  ) +
  theme(
    strip.text = element_text(size = 12, face = "bold"),
    axis.text.x = element_text(size = 12),
    legend.position = "none"
  )
ggplot(estimate_long_filtered, aes(x = risk_group, y = Score, fill = risk_group)) +
  geom_violin(trim = FALSE, alpha = 0.9) +
  geom_boxplot(width = 0.1, outlier.shape = NA, alpha = 0.6, color = "black") +
  geom_jitter(width = 0.2, size = 1, alpha = 0.5) +
  stat_compare_means(method = "wilcox.test", label = "p.format", size = 4) +
  facet_wrap(~ ScoreType, scales = "free_y") +
  scale_fill_manual(values = c("High" = "#D7263D", "Low" = "#1B9AAA")) +  # Red and Blue
  theme_minimal(base_size = 14) +
  labs(
    title = "TME Scores by Risk Group (ESTIMATE, Violin Plot)",
    x = "Risk Group", y = "Score"
  ) +
  theme(
    strip.text = element_text(size = 12, face = "bold"),
    axis.text.x = element_text(size = 12),
    legend.position = "none"
  )
ggplot(estimate_long_filtered, aes(x = risk_group, y = Score, fill = risk_group)) +
  geom_violin(trim = FALSE, alpha = 0.9) +
  geom_boxplot(width = 0.1, outlier.shape = NA, alpha = 0.6, color = "black", size = 1.2) +  # <— bold lines
  geom_jitter(width = 0.2, size = 1, alpha = 0.5) +
  stat_compare_means(method = "wilcox.test", label = "p.format", size = 4) +
  facet_wrap(~ ScoreType, scales = "free_y") +
  scale_fill_manual(values = c("High" = "#D7263D", "Low" = "#1B9AAA")) +  # Custom fill
  theme_minimal(base_size = 14) +
  labs(
    title = "TME Scores by Risk Group (ESTIMATE, Violin Plot)",
    x = "Risk Group", y = "Score"
  ) +
  theme(
    strip.text = element_text(size = 12, face = "bold"),
    axis.text.x = element_text(size = 12),
    legend.position = "none"
  )
library(ggplot2)
library(ggpubr)
library(tidyr)
library(dplyr)
# Merge ESTIMATE results with metadata by sample name

# Now you can reshape
estimate_long <- estimate_corr %>%
  pivot_longer(cols = c("StromalScore", "ImmuneScore", "ESTIMATEScore"),
               names_to = "ScoreType", values_to = "Score")

# Reshape to long format for ggplot
estimate_long <- estimate_corr %>%
  pivot_longer(cols = c("StromalScore", "ImmuneScore", "ESTIMATEScore"),
               names_to = "ScoreType", values_to = "Score")

# Scatter plot with regression and correlation
ggplot(estimate_long, aes(x = risk_score, y = Score)) +
  geom_point(alpha = 0.6, color = "#1B9AAA") +
  geom_smooth(method = "lm", se = TRUE, color = "#D7263D", size = 1) +
  stat_cor(method = "pearson", label.x.npc = "left", label.y.npc = "top", size = 4) +
  facet_wrap(~ ScoreType, scales = "free_y") +
  theme_minimal(base_size = 14) +
  labs(title = "Correlation Between Risk Score and ESTIMATE Scores",
       x = "Risk Score", y = "TME Score") +
  theme(strip.text = element_text(size = 12, face = "bold"))
# Calculate correlation coefficients
cor_res <- sapply(c("StromalScore", "ImmuneScore", "ESTIMATEScore"), function(x) {
  cor.test(estimate_corr[[x]], estimate_corr$risk_score)$estimate
})

cor_df <- data.frame(
  ScoreType = names(cor_res),
  Correlation = as.numeric(cor_res)
)

ggplot(cor_df, aes(x = ScoreType, y = Correlation, fill = ScoreType)) +
  geom_bar(stat = "identity", width = 0.6, color = "black") +
  geom_text(aes(label = round(Correlation, 2)), vjust = -0.5, size = 5) +
  scale_fill_manual(values = c("StromalScore" = "#FFA07A", 
                               "ImmuneScore" = "#87CEEB", 
                               "ESTIMATEScore" = "#90EE90")) +
  ylim(min(cor_df$Correlation) - 0.1, 1) +
  theme_minimal(base_size = 14) +
  labs(title = "Pearson Correlation Between Risk Score and TME Components",
       x = "TME Score Type", y = "Pearson Correlation Coefficient") +
  theme(legend.position = "none",
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 13, face = "bold"))



#-------------IMVIGOR210 -------------
# Install Bioconductor (if not already installed)
