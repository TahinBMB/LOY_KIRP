######Training set



setwd("E:/LOY_KIRP")






library(dplyr)
library(biomaRt)
library(tibble)  # for column_to_rownames()

#---------Step 1: Load data-------

metadat <- read.delim("TCGA-KIRP.clinical.tsv")
exp_data <- read.delim("TCGA-KIRP.star_tpm.tsv")
surv_data <- read.delim("TCGA-KIRP.survival.tsv")


#------Preprocessing------

metadat_filtered <- metadat %>%
  filter(sample_type.samples != "Solid Tissue Normal") %>%
  filter(gender.demographic == "male")

# For matching, metadata sample IDs use dash "-", expression data column names use dots "."
# So convert metadata sample IDs from dash to dot format for matching expression columns
male_sample_ids <- gsub("-", ".", metadat_filtered$sample)

# --- Step 3: Prepare expression matrix ---

# Set Ensembl_ID column as rownames
exp_data <- exp_data %>% column_to_rownames(var = "Ensembl_ID")

# Make sure expression sample column names use dot format to match metadata
colnames(exp_data) <- gsub("-", ".", colnames(exp_data))

# Subset expression data to keep only male tumor samples
exp_data_male <- exp_data[, colnames(exp_data) %in% male_sample_ids]

# --- Step 4: Clean Ensembl IDs by removing version suffixes ---

ensembl_ids_clean <- sub("\\..*", "", rownames(exp_data_male))

# Add clean Ensembl IDs as a new column for deduplication
exp_data_male$ensembl_id_clean <- ensembl_ids_clean

# --- Step 5: Handle duplicate genes ---

# Calculate mean expression per gene (row) across all samples
expr_cols <- setdiff(colnames(exp_data_male), "ensembl_id_clean")
exp_data_male$mean_expr <- rowMeans(exp_data_male[, expr_cols], na.rm = TRUE)

# Keep the row with highest mean expression for each clean Ensembl ID
exp_data_unique <- exp_data_male %>%
  arrange(desc(mean_expr)) %>%
  distinct(ensembl_id_clean, .keep_all = TRUE)

# Set clean Ensembl IDs as rownames
rownames(exp_data_unique) <- exp_data_unique$ensembl_id_clean

# Remove helper columns
exp_data_unique <- exp_data_unique[, expr_cols]

# --- Step 6: Map Ensembl IDs to HGNC gene symbols ---

# Connect to Ensembl
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# Get mapping for clean Ensembl IDs
gene_map <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol"),
  filters = "ensembl_gene_id",
  values = rownames(exp_data_unique),
  mart = mart
)

# Add Ensembl IDs as a column for merge
exp_data_unique$ensembl_id <- rownames(exp_data_unique)

# Merge expression data with gene symbols
expr_with_symbols <- merge(
  exp_data_unique,
  gene_map,
  by.x = "ensembl_id",
  by.y = "ensembl_gene_id",
  all.x = TRUE
)

# --- Step 7: Remove rows without gene symbols and duplicates ---

# Remove rows missing gene symbols or empty
expr_with_symbols <- expr_with_symbols[!is.na(expr_with_symbols$hgnc_symbol) & expr_with_symbols$hgnc_symbol != "", ]

# Remove duplicated gene symbols, keep first occurrence
expr_with_symbols <- expr_with_symbols[!duplicated(expr_with_symbols$hgnc_symbol), ]

# Set gene symbols as rownames
rownames(expr_with_symbols) <- expr_with_symbols$hgnc_symbol

# Remove helper columns
expr_final <- expr_with_symbols[, !(colnames(expr_with_symbols) %in% c("ensembl_id", "hgnc_symbol"))]

# --- Final: expr_final is your clean expression matrix ready for analysis ---

dim(expr_final)
head(expr_final[, 1:5])


expr_final
metadat_filtered

#-------loy ssgsea--------
#3. Compute ssGSEA LOY scores (for male samples only)
LOY_set <- c("RPS4Y1","ZFY","CDHL3","TBL1Y","USP9Y",
             "DDX3Y","UTY","TMSB4Y","NLGN4Y","HSFY2","KDM5D",
             "EIF1AY","RBMY1A1","PRY2")

# Custom ssGSEA function
ssgsea <- function(X, gene_sets, alpha = 0.25, scale = T, norm = F, single = T) {
  row_names <- rownames(X)
  num_genes <- nrow(X)
  gene_sets <- lapply(gene_sets, function(genes) {which(row_names %in% genes)})
  R <- matrixStats::colRanks(X, preserveShape = T, ties.method = 'average')
  es <- apply(R, 2, function(R_col) {
    gene_ranks <- order(R_col, decreasing = TRUE)
    es_sample <- sapply(gene_sets, function(gene_set_idx) {
      indicator_pos <- gene_ranks %in% gene_set_idx
      indicator_neg <- !indicator_pos
      rank_alpha <- (R_col[gene_ranks] * indicator_pos) ^ alpha
      step_cdf_pos <- cumsum(rank_alpha) / sum(rank_alpha)
      step_cdf_neg <- cumsum(indicator_neg) / sum(indicator_neg)
      step_cdf_diff <- step_cdf_pos - step_cdf_neg
      if (scale) step_cdf_diff <- step_cdf_diff / num_genes
      if (single) sum(step_cdf_diff) else step_cdf_diff[which.max(abs(step_cdf_diff))]
    })
    unlist(es_sample)
  })
  if (length(gene_sets) == 1) es <- matrix(es, nrow = 1)
  if (norm) es <- es / diff(range(es))
  rownames(es) <- names(gene_sets)
  colnames(es) <- colnames(X)
  return(es)
}
library(tibble)

# Run ssGSEA
LOY_ssgsea_male <- ssgsea(as.matrix(expr_final), LOY_set)
LOY_ssgsea_male <- na.omit(LOY_ssgsea_male)
LOY_means_male <- colMeans(LOY_ssgsea_male)

LOY_df_male <- data.frame(
  sample = names(LOY_means_male),
  LOY_score = as.numeric(LOY_means_male)
)
LOY_df_male$sample_dash <- gsub("\\.", "-", LOY_df_male$sample)



#------survival by optimal cutoff-------

# Make sure survminer and survival are loaded
library(survival)
library(survminer)

# Step 1: Prepare the survival dataframe (if not already)
male_surv <- left_join(LOY_df_male, metadat_filtered, by = c("sample_dash" = "sample"))
male_surv <- left_join(male_surv, surv_data, by = c("sample_dash" = "sample"))

# Remove missing survival data
male_surv_clean <- male_surv %>%
  filter(!is.na(OS.time), !is.na(OS), !is.na(LOY_score))

# Step 2: Determine optimal cutoff for LOY_score
cutpoint <- surv_cutpoint(
  male_surv_clean,
  time = "OS.time",
  event = "OS",
  variables = "LOY_score",
  minprop = 0.1  # ensures each group has at least 10% of samples
)

summary(cutpoint)  # Check the optimal cutoff value
# Correct extraction of cutoff value
cutoff_value <- cutpoint$cutpoint["LOY_score", "cutpoint"]

#  Classify groups based on LOY_score
male_surv_clean <- male_surv_clean %>%
  mutate(LOY_group = ifelse(LOY_score > cutoff_value, "Y_high", "Y_low"))

# Set factor levels to ensure Y_high appears first (basal group)
male_surv_clean$LOY_group <- factor(male_surv_clean$LOY_group, levels = c("Y_high", "Y_low"))

# Create boxplot
# Load required libraries
library(ggplot2)
library(ggpubr)  # for stat_compare_means

# Create boxplot with p-value
ggplot(male_surv_clean, aes(x = LOY_group, y = LOY_score, fill = LOY_group)) +
  geom_boxplot(outlier.shape = NA, width = 0.6) +
  geom_jitter(width = 0.2, alpha = 0.6, size = 2) +
  scale_fill_manual(values = c("Y_high" = "steelblue", "Y_low" = "firebrick")) +
  labs(title = "LOY Score by Group",
       x = "LOY Group",
       y = "LOY Score") +
  stat_compare_means(method = "wilcox.test", 
                     label.y = max(male_surv_clean$LOY_score) * 1.05,
                     size = 5) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none")

##enhanced boxplot
library(ggplot2)
library(ggpubr)

ggplot(male_surv_clean, aes(x = LOY_group, y = LOY_score, fill = LOY_group)) +
  geom_boxplot(
    outlier.shape = NA, 
    width = 0.6,
    alpha = 0.8,
    lwd = 0.5,           # Thinner boxplot lines
    fatten = 1.2         # Slightly thicker median line
  ) +
  geom_jitter(
    width = 0.2, 
    alpha = 0.5, 
    size = 2,
    shape = 21,          # Filled circles with outline
    color = "black"      # Outline color
  ) +
  scale_fill_manual(
    values = c("Y_high" = "#3B7EA1", "Y_low" = "#D55E00"),  # Colorblind-friendly palette
    labels = c("Y_high" = "Y_high", "Y_low" = "Y_low")  # More descriptive labels
  ) +
  scale_x_discrete(
    labels = c("Y_high" = "Y_high", "Y_low" = "Y_low")  # Consistent labeling
  ) +
  labs(
    title = "Association Between LOY Groups and LOY Scores",
    x = "LOY Group",
    y = "LOY Score (normalized)",
    caption = "Wilcoxon rank-sum test; p-value shown"
  ) +
  stat_compare_means(
    method = "wilcox.test", 
    label = "p.format",                # Show just the p-value
    label.x = 1.5,                     # Center the p-value
    label.y = max(male_surv_clean$LOY_score) * 1.08,
    size = 4.5,
    bracket.size = 0.5,                # Thinner bracket
    tip.length = 0.02                  # Shorter bracket tips
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "none",
    panel.grid.major = element_line(linewidth = 0.25),  # Thinner grid lines
    panel.grid.minor = element_blank(),
    axis.line = element_line(linewidth = 0.5),         # Add axis lines
    axis.ticks = element_line(linewidth = 0.5),        # Add ticks
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    axis.title = element_text(face = "bold"),
    plot.caption = element_text(size = 10, hjust = 0)
  ) +
  coord_cartesian(ylim = c(min(male_surv_clean$LOY_score), 
                           max(male_surv_clean$LOY_score) * 1.1))  # Adjust y-axis limits

#------survival by optimal cutoff-------

# Make sure survminer and survival are loaded
library(survival)
library(survminer)

# Step 1: Prepare the survival dataframe (if not already)
male_surv <- left_join(LOY_df_male, metadat_filtered, by = c("sample_dash" = "sample"))
male_surv <- left_join(male_surv, surv_data, by = c("sample_dash" = "sample"))

# Remove missing survival data
male_surv_clean <- male_surv %>%
  filter(!is.na(OS.time), !is.na(OS), !is.na(LOY_score))

# Step 2: Determine optimal cutoff for LOY_score
cutpoint <- surv_cutpoint(
  male_surv_clean,
  time = "OS.time",
  event = "OS",
  variables = "LOY_score",
  minprop = 0.1  # ensures each group has at least 10% of samples
)

summary(cutpoint)  # Check the optimal cutoff value

# Step 3: Categorize samples based on optimal cutoff
male_surv_cut <- surv_categorize(cutpoint)




# Step 4: Fit survival model using optimal cutoff
surv_fit_optimal <- survfit(Surv(OS.time, OS) ~ LOY_score, data = male_surv_cut)

# Step 5: Plot KM curve
ggsurvplot(
  surv_fit_optimal,
  data = male_surv_cut,
  pval = TRUE,
  risk.table = FALSE,
  legend.title = "LOY Group (Optimal Cut)",
  legend.labs = c("Y_low", "Y_high"),
  title = "Overall Survival by Optimal LOY Cutoff",
  xlab = "Time (Days)",
  ylab = "Survival Probability",
  surv.scale = "percent",  # Shows y-axis as %
  ggtheme = theme_minimal(base_size = 14),  # Clean base theme
  legend = "top",  # Position the legend
  font.main = c(16, "bold"),
  font.x = c(14),
  font.y = c(14),
  font.tickslab = c(12),
  pval.size = 5,
  size = 1.2,  # Line thickness
  linetype = "solid"
)



#----------------------deg--------------

library(limma)


# Match samples using dot-format sample names
male_surv_clean$sample_dot <- gsub("-", ".", male_surv_clean$sample_dash)

# Ensure overlap
common_samples <- intersect(male_surv_clean$sample_dot, colnames(expr_final))

# Subset expression data and survival metadata
expr_final_subset <- expr_final[, common_samples]
surv_subset <- male_surv_clean %>% filter(sample_dot %in% common_samples)

# Match order
surv_subset <- surv_subset[match(colnames(expr_final_subset), surv_subset$sample_dot), ]

# Sanity check
stopifnot(all(colnames(expr_final_subset) == surv_subset$sample_dot))
# Apply cutoff and define group
LOY_cutoff <- 0.02317572
surv_subset <- surv_subset %>%
  mutate(
    LOY_group = ifelse(LOY_score > LOY_cutoff, "Y_high", "Y_low"),
    LOY_group = factor(LOY_group, levels = c("Y_high", "Y_low"))
  )
group <- factor(surv_subset$LOY_group, levels = c("Y_high", "Y_low"))  # Y_high = baseline

design <- model.matrix(~ group)
colnames(design)
fit <- lmFit(expr_final_subset, design)
fit <- eBayes(fit)
deg_results <- topTable(fit, coef = "groupY_low", number = Inf, adjust.method = "BH")

head(deg_results)

write.csv(deg_results, "DEG_Ylow_vs_Yhigh_optimalCutoff.csv", row.names = TRUE)



library(ggplot2)
library(ggrepel)  # For better text label placement
library(dplyr)

# Mark significance and direction
deg_results <- deg_results %>%
  mutate(
    significant = adj.P.Val < 0.05 & abs(logFC) > 0.5,
    regulation = case_when(
      adj.P.Val < 0.05 & logFC > 0.5  ~ "Upregulated in Y_low",
      adj.P.Val < 0.05 & logFC < -0.5 ~ "Downregulated in Y_low",
      TRUE ~ "Not Significant"
    )
  )

# Select top 10 up and top 10 down genes by adj.P.Val and abs(logFC) for labeling
top_up <- deg_results %>%
  filter(regulation == "Upregulated in Y_low") %>%
  arrange(adj.P.Val, desc(logFC)) %>%
  head(10)

top_down <- deg_results %>%
  filter(regulation == "Downregulated in Y_low") %>%
  arrange(adj.P.Val, logFC) %>%
  head(10)

top_genes <- bind_rows(top_up, top_down)

# Plot


library(ggplot2)
library(ggrepel)

ggplot(deg_results, aes(x = logFC, y = -log10(adj.P.Val), color = regulation)) +
  geom_point(alpha = 0.7, size = 1.8) +
  scale_color_manual(
    values = c(
      "Upregulated in Y_low" = "#1f78b4",  # Blue
      "Downregulated in Y_low" = "#e31a1c",  # Red
      "Not Significant" = "grey70"
    )
  ) +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", color = "black", size = 0.5) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black", size = 0.5) +
  geom_text_repel(
    data = top_genes,
    aes(label = rownames(top_genes)),
    size = 4,
    box.padding = 0.4,
    point.padding = 0.3,
    max.overlaps = Inf,
    segment.color = "grey50",
    fontface = "bold"  # Bold font here
  ) +
  labs(
    title = "Differential Expression Y_low vs Y_high",
    x = expression(Log[2]~Fold~Change~"(Y_low vs Y_high)"),
    y = expression(-log[10]~"Adjusted P-value"),
    color = "Gene Regulation"
  ) +
  theme_minimal(base_size = 15) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 17),
    legend.position = "right",
    legend.title = element_text(size = 13, face = "bold"),
    legend.text = element_text(size = 12),
    axis.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 12)
  )

#------------heatmap--------------


library(ComplexHeatmap)
library(circlize)
library(dplyr)
library(tibble)
library(RColorBrewer)


# Step 1: Select top 20 DEGs (10 up + 10 down)
top_degs_20 <- bind_rows(top_up, top_down)
top_genes <- rownames(top_degs_20)

# Step 2: Subset expression matrix to top genes and samples (matching order)
heatmap_mat <- expr_final[top_genes, surv_subset$sample_dot]

# Step 3: Scale rows (genes) for visualization (z-score)
heatmap_mat_scaled <- t(scale(t(as.matrix(heatmap_mat))))

# Step 4: Prepare group factor for column splitting
group_factor <- surv_subset$LOY_group[match(colnames(heatmap_mat_scaled), surv_subset$sample_dot)]

# Step 5: Create top annotation (colored bar for LOY groups)
col_ha <- HeatmapAnnotation(
  Group = group_factor,
  col = list(Group = c("Y_high" = "steelblue", "Y_low" = "firebrick")),
  show_annotation_name = TRUE
)

# Step 6: Define color function for heatmap
col_fun <- colorRamp2(c(-2, 0, 2), c("navy", "white", "firebrick3"))

# Step 7: Plot heatmap with columns split by LOY group
Heatmap(
  heatmap_mat_scaled,
  name = "Z-score",
  col = col_fun,
  top_annotation = col_ha,
  show_row_names = TRUE,
  show_column_names = FALSE,
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  column_split = group_factor,
  column_title = "Top 20 DEGs Expression by LOY Group",
  row_names_gp = gpar(fontsize = 11, fontface = "bold"),
  heatmap_legend_param = list(title = "Expression (Z-score)")
)





sig_degs <- deg_results %>% filter(adj.P.Val < 0.05 & abs(logFC) > 0.5)
up_genes <- rownames(sig_degs)[sig_degs$logFC > 0.5]
down_genes <- rownames(sig_degs)[sig_degs$logFC < -0.5]




library(clusterProfiler)
library(org.Hs.eg.db)

# Convert gene symbols to Entrez IDs
entrez_up <- bitr(up_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
entrez_down <- bitr(down_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

# GO BP enrichment for upregulated genes
ego_up <- enrichGO(
  gene = entrez_up$ENTREZID,
  OrgDb = org.Hs.eg.db,
  ont = "BP",
  pAdjustMethod = "BH",
  qvalueCutoff = 0.05,
  readable = TRUE
)

# Plot top 10 enriched GO BP terms


library(clusterProfiler)
library(ggplot2)

# Generate the barplot and customize fonts
barplot(ego_up, showCategory = 10, title = "GO Biological Process Enriched in Y_low") +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    axis.title.x = element_text(face = "bold", size = 12),
    axis.title.y = element_text(face = "bold", size = 12),
    axis.text.x = element_text(face = "bold", size = 10),
    axis.text.y = element_text(face = "bold", size = 10),
    legend.title = element_text(face = "bold", size = 12),
    legend.text = element_text(face = "bold", size = 10)
  )



# GO BP enrichment for downregulated genes
ego_down <- enrichGO(
  gene = entrez_down$ENTREZID,
  OrgDb = org.Hs.eg.db,
  ont = "BP",
  pAdjustMethod = "BH",
  qvalueCutoff = 0.05,
  readable = TRUE
)

# Plot top 10 enriched GO BP terms
barplot(
  ego_down,
  showCategory = 10,
  title = "GO Biological Process Enriched in Y_high"
) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    axis.title.x = element_text(face = "bold", size = 12),
    axis.title.y = element_text(face = "bold", size = 12),
    axis.text.x = element_text(face = "bold", size = 10),
    axis.text.y = element_text(face = "bold", size = 10),
    legend.title = element_text(face = "bold", size = 12),
    legend.text = element_text(face = "bold", size = 10)
  )






##dorkar nai

dotplot(
  ego_down,
  showCategory = 10,
  title = "GO Biological Process Enriched in Downregulated Genes",
  font.size = 10
) +
  theme(
    axis.text.y = element_text(size = 10, hjust = 1),
    plot.title = element_text(size = 14, face = "bold"),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9)
  )




dotplot(
  ego_up,
  showCategory = 10,
  title = "GO Biological Process Enriched in Upregulated Genes",
  font.size = 10
) +
  theme(
    axis.text.y = element_text(size = 10, hjust = 1),
    plot.title = element_text(size = 14, face = "bold"),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9)
  )





###net er jonno hoina

kk_up <- enrichKEGG(
  gene = entrez_up$ENTREZID,
  organism = "hsa",
  pAdjustMethod = "BH",
  qvalueCutoff = 0.05
)
dotplot(kk_up, showCategory = 10, title = "KEGG Pathways Enriched in Y_low")+
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    axis.title.x = element_text(face = "bold", size = 12),
    axis.title.y = element_text(face = "bold", size = 12),
    axis.text.x = element_text(face = "bold", size = 10),
    axis.text.y = element_text(face = "bold", size = 10),
    legend.title = element_text(face = "bold", size = 12),
    legend.text = element_text(face = "bold", size = 10)
  )




kk_down <- enrichKEGG(
  gene = entrez_down$ENTREZID,
  organism = "hsa",
  pAdjustMethod = "BH",
  qvalueCutoff = 0.05
)
dotplot(kk_down, showCategory = 10, title = "KEGG Pathways Enriched in Downregulated Genes")+
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    axis.title.x = element_text(face = "bold", size = 12),
    axis.title.y = element_text(face = "bold", size = 12),
    axis.text.x = element_text(face = "bold", size = 10),
    axis.text.y = element_text(face = "bold", size = 10),
    legend.title = element_text(face = "bold", size = 12),
    legend.text = element_text(face = "bold", size = 10)
  )





#BiocManager::install("ReactomePA")
library(ReactomePA)




reactome_up <- enrichPathway(
  gene = entrez_up$ENTREZID,
  organism = "human",
  pAdjustMethod = "BH",
  qvalueCutoff = 0.05,
  readable = TRUE
)

dotplot(reactome_up, showCategory = 10, title = "Reactome Pathways Enriched in Upregulated Genes")+
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    axis.title.x = element_text(face = "bold", size = 12),
    axis.title.y = element_text(face = "bold", size = 12),
    axis.text.x = element_text(face = "bold", size = 10),
    axis.text.y = element_text(face = "bold", size = 10),
    legend.title = element_text(face = "bold", size = 12),
    legend.text = element_text(face = "bold", size = 10)
  )

# Load required library
library(enrichplot)
# Compute term similarity
reactome_up_sim <- pairwise_termsim(reactome_up)

# Now plot
emapplot(
  reactome_up_sim,
  showCategory = 30,
  color = "p.adjust",
  layout = "kk"
) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    legend.title = element_text(face = "bold", size = 12),
    legend.text = element_text(face = "bold", size = 10)
  )

# Compute term similarity
reactome_up_sim <- pairwise_termsim(reactome_up)

# Enhanced emapplot with bold edges
emapplot(
  reactome_up_sim,
  showCategory = 30,
  color = "p.adjust",
  layout = "kk",
  cex_line = 1.5,          # BOLD EDGE width
  line_color = "black"     # EDGE COLOR
) +
  ggtitle("Enriched Reactome Pathways in Y_low Group") +
  xlab("Semantic Space X") +
  ylab("Semantic Space Y") +
  scale_color_gradient(low = "red", high = "blue", name = "Adjusted p-value") +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    axis.title.x = element_text(face = "bold", size = 12),
    axis.title.y = element_text(face = "bold", size = 12),
    axis.text.x = element_text(face = "bold", size = 10),
    axis.text.y = element_text(face = "bold", size = 10),
    legend.title = element_text(face = "bold", size = 12),
    legend.text = element_text(face = "bold", size = 10)
  )



library(enrichplot)
library(ggplot2)

emapplot(
  reactome_up_sim,
  showCategory = 30,
  color = "p.adjust",
  layout = "kk"
) +
  ggtitle("Enriched Reactome Pathways in Y_low Group") +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    legend.title = element_text(face = "bold", size = 12),
    legend.text = element_text(face = "bold", size = 10)
  )


reactome_down <- enrichPathway(
  gene = entrez_down$ENTREZID,
  organism = "human",
  pAdjustMethod = "BH",
  qvalueCutoff = 0.05,
  readable = TRUE
)

barplot(reactome_down, showCategory = 10, title = "Reactome Pathways Enriched in Downregulated Genes")+
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    axis.title.x = element_text(face = "bold", size = 12),
    axis.title.y = element_text(face = "bold", size = 12),
    axis.text.x = element_text(face = "bold", size = 10),
    axis.text.y = element_text(face = "bold", size = 10),
    legend.title = element_text(face = "bold", size = 12),
    legend.text = element_text(face = "bold", size = 10)
  )


# Compute term similarity for downregulated gene enrichment
reactome_down_sim <- pairwise_termsim(reactome_down)

emapplot(
  reactome_down_sim,
  showCategory = 30,
  color = "p.adjust",
  layout = "kk"
) +
  ggtitle("Enriched Reactome Pathways in Y_high Group") +
  xlab("Semantic Space X") +
  ylab("Semantic Space Y") +
  scale_color_gradient(low = "red", high = "blue", name = "Adjusted p-value") +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    axis.title.x = element_text(face = "bold", size = 12),
    axis.title.y = element_text(face = "bold", size = 12),
    axis.text.x = element_text(face = "bold", size = 10),
    axis.text.y = element_text(face = "bold", size = 10),
    legend.title = element_text(face = "bold", size = 12),
    legend.text = element_text(face = "bold", size = 10)
  )



#-----------hallmark---------------
library(clusterProfiler)

library(msigdbr)

# Prepare ranked gene list for GSEA: all genes ordered by logFC
gene_list <- deg_results$logFC
names(gene_list) <- rownames(deg_results)
gene_list <- sort(gene_list, decreasing = TRUE)

# Load Hallmark gene sets from MSigDB
hallmark <- msigdbr(species = "Homo sapiens", category = "H")

gsea_res <- GSEA(
  geneList = gene_list,
  TERM2GENE = hallmark[, c("gs_name", "gene_symbol")],
  pAdjustMethod = "BH",
  minGSSize = 15,
  maxGSSize = 500,
  verbose = FALSE,
  eps = 0  # For better precision of small p-values
)

# Load the required library
library(enrichplot)


dotplot(gsea_res,
        showCategory = 25,
        color = "p.adjust",
        x = "NES",
        title = "Top 25 Hallmark Pathways (GSEA)") +
  theme_classic(base_size = 14) +
  scale_color_gradient(low = "red", high = "blue") +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(face = "bold")
  )





#----------selected hallmark analysis in Ylow---------------

# 1. Filter GSEA result for NES > 0 and adjusted p < 0.05
gsea_pos <- gsea_res@result
gsea_pos <- gsea_pos[gsea_pos$NES > 0 & gsea_pos$p.adjust < 0.05, ]
library(dplyr)

top3_pos <- gsea_pos %>%
  arrange(desc(NES)) %>%
  slice_head(n = 3)

##gsea plot
library(enrichplot)
library(dplyr)

# Filter positive NES and adj p < 0.05
gsea_pos <- gsea_res@result %>%
  filter(NES > 0, p.adjust < 0.05) %>%
  arrange(desc(NES))

# Select top 3 pathways
top3_pos <- head(gsea_pos$ID, 3)

# Plot combined GSEA enrichment curves for top 3 positive pathways
gseaplot2(gsea_res, geneSetID = top3_pos, 
          title = "GSEA Hallmark Pathways Enriched in Y_low", 
          base_size = 14)





# 2. Select hallmark pathways of interest
target_hallmarks <- c("HALLMARK_MITOTIC_SPINDLE", "HALLMARK_E2F_TARGETS","HALLMARK_G2M_CHECKPOINT")

# 3. Get gene symbols for those Hallmark pathways
hallmark_selected <- msigdbr(species = "Homo sapiens", category = "H") %>%
  dplyr::filter(gs_name %in% target_hallmarks)

# 4. Save the hallmark genes to CSV
write.csv(hallmark_selected, "selected_hallmark_genes.csv", row.names = FALSE)

# 5. Get unique list of hallmark gene symbols
hallmark_genes <- unique(hallmark_selected$gene_symbol)


# 6. Identify DEGs from your results (e.g., using FDR < 0.05 as threshold)
deg_sig <- deg_results[deg_results$adj.P.Val < 0.05, ]
deg_gene_names <- rownames(deg_sig)

# 7. Find intersection between DEGs and selected hallmark genes
common_genes <- intersect(deg_gene_names, hallmark_genes)

# 8. Save overlapping genes
common_df <- data.frame(Gene = common_genes)
write.csv(common_df, "DEG_overlap_with_hallmark.csv", row.names = FALSE)






#----------Downstream--------

gene_list <- read.csv("DEG_overlap_with_hallmark.csv")$Gene
common_genes <- intersect(gene_list, rownames(expr_final))
expr_subset <- expr_final[common_genes, ]

# Match samples between expression and survival
common_samples <- intersect(colnames(expr_subset), male_surv_clean$sample)
expr_subset <- expr_subset[, common_samples]
surv_subset <- male_surv_clean %>% filter(sample %in% common_samples)
surv_subset <- surv_subset[match(common_samples, surv_subset$sample), ]

# Check alignment
stopifnot(all(colnames(expr_subset) == surv_subset$sample))



library(survival)

expr_t <- t(expr_subset)
univ_results <- lapply(colnames(expr_t), function(gene) {
  fit <- coxph(Surv(OS.time, OS) ~ expr_t[, gene], data = surv_subset)
  summary(fit)$coefficients
})

univ_df <- do.call(rbind, univ_results)
univ_df <- data.frame(Gene = colnames(expr_t), univ_df)
univ_df <- univ_df[, c("Gene", "coef", "exp.coef.", "se.coef.", "Pr...z..")]
colnames(univ_df) <- c("Gene", "coef", "HR", "SE", "p.value")

# Filter significant genes
sig_genes <- univ_df %>% filter(p.value < 0.05)
print(sig_genes)


# Save significant univariate Cox results as CSV
write.csv(sig_genes, file = "Significant_Univariate_Cox_Genes.csv", row.names = FALSE)




library(glmnet)

# Build input matrix and survival object
x <- as.matrix(expr_t[, sig_genes$Gene])
y <- Surv(surv_subset$OS.time, surv_subset$OS)

# LASSO
cvfit <- cv.glmnet(x, y, family = "cox", alpha = 1)

plot(cvfit, xlab = "log(Lambda)", ylab = "Mean CV Error")



# Extract selected genes
lasso_coef <- coef(cvfit, s = "lambda.min")
selected_genes <- rownames(lasso_coef)[lasso_coef[,1] != 0]
selected_genes


##bootstap validation
library(glmnet)
library(survival)
library(ggplot2)
library(dplyr)

set.seed(123)

x <- as.matrix(expr_t[, sig_genes$Gene])
y <- Surv(surv_subset$OS.time, surv_subset$OS)

nboot <- 1000  # number of bootstrap samples
selected_matrix <- matrix(0, nrow = length(sig_genes$Gene), ncol = nboot,
                          dimnames = list(sig_genes$Gene, NULL))

for (i in seq_len(nboot)) {
  # Bootstrap sample indices
  boot_idx <- sample(seq_len(nrow(x)), replace = TRUE)
  
  x_boot <- x[boot_idx, ]
  y_boot <- y[boot_idx]
  
  # Fit LASSO Cox on bootstrap sample
  cv_boot <- cv.glmnet(x_boot, y_boot, family = "cox", alpha = 1, nfolds = 5)
  coef_boot <- coef(cv_boot, s = "lambda.min")
  
  # Mark selected genes (nonzero coefficients)
  selected_matrix[, i] <- as.numeric(coef_boot != 0)
}

# Calculate selection frequency for each gene
selection_freq <- rowMeans(selected_matrix)

selection_df <- data.frame(
  Gene = names(selection_freq),
  Frequency = selection_freq
) %>% arrange(desc(Frequency))

# View top stable genes
head(selection_df, 20)




library(glmnet)

# Fit LASSO Cox model
cvfit <- cv.glmnet(x, y, family = "cox", alpha = 1)

# Extract coefficients at lambda.min
coef_min <- coef(cvfit, s = "lambda.min")
selected_genes <- rownames(coef_min)[which(coef_min != 0)]
selected_genes <- selected_genes[selected_genes != "(Intercept)"]

# Plot coefficient paths with thicker lines
plot(cvfit$glmnet.fit, xvar = "lambda", label = FALSE, lwd = 2)


# Extract coefficient values at lambda.min
coefs <- coef(cvfit$glmnet.fit)
lambda_min_idx <- which.min(abs(cvfit$lambda - cvfit$lambda.min))
nz <- which(coefs[, lambda_min_idx] != 0)

# Coordinates for labels
x_pos <- log(cvfit$lambda.min)
y_pos <- coefs[nz, lambda_min_idx]
labels <- rownames(coefs)[nz]

# To reduce label overlap, jitter the y position slightly
set.seed(123)  # for reproducibility
y_jitter <- jitter(y_pos, amount = 0.02)

# Add bold black text labels, slightly to the right of the points
# Add bold black text labels, slightly to the right of the points
text(
  x = x_pos ,   # slight shift right of the vertical line
  y = y_jitter,
  labels = labels,
  pos = 4,           # labels to the right
  cex = 0.7,        # slightly smaller font size
  font = 2,          # bold font
  col = "black"
)


# Optional: add points at coefficients at lambda.min for visual cues
points(
  rep(x_pos, length(y_pos)),
  y_pos,
  pch = 19,
  col = "red",
  cex = 1.2
)




#Selection frequency bar plot from bootstrap
library(ggplot2)

ggplot(selection_df, aes(x = reorder(Gene, Frequency), y = Frequency)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +
  labs(title = "Bootstrap Selection Frequency of Genes in LASSO Model",
       x = "Gene", y = "Selection Frequency") +
  theme_minimal(base_size = 14)


##Build risk score and plot survival curves by median cutoff

# Final model coefficients at lambda.min
final_coef <- coef(cvfit, s = "lambda.min")
selected_genes <- rownames(final_coef)[final_coef[,1] != 0]

# Calculate risk score for each patient
risk_score <- as.vector(x[, selected_genes] %*% final_coef[selected_genes, 1])

# Add risk score to survival dataframe
surv_subset$risk_score <- risk_score

# Dichotomize risk groups (median cutoff)
surv_subset$risk_group <- ifelse(surv_subset$risk_score > median(risk_score), "High", "Low")

#saving as csv
# Build data frame with sample ID, risk score, and risk group
risk_output <- data.frame(
  sample = rownames(surv_subset),
  risk_score = surv_subset$risk_score,
  risk_group = surv_subset$risk_group
)

# Add expression values of selected genes (optional)
risk_output <- cbind(risk_output, x[, selected_genes])

# Save as CSV
write.csv(risk_output, "selected_genes_risk_score_groups.csv", row.names = FALSE)



##save as csv for selected genes only
# Extract coefficients for the 6 selected genes
final_coef <- coef(cvfit, s = "lambda.min")
selected_genes <- rownames(final_coef)[final_coef[,1] != 0]
selected_genes <- selected_genes[selected_genes != "(Intercept)"]

# Get coefficients as numeric vector
selected_coef <- as.numeric(final_coef[selected_genes, 1])
names(selected_coef) <- selected_genes

# Compute risk score using the selected genes only
risk_score <- as.vector(x[, selected_genes] %*% selected_coef)

# Create a dataframe with coefficients and the risk score
output_df <- data.frame(
  gene = selected_genes,
  coefficient = selected_coef
)

# Append risk score (same for all genes) as a new column, if you want to include it this way
# This shows the risk score per sample in a separate CSV
write.csv(output_df, "selected_genes_coefficients.csv", row.names = FALSE)



library(survminer)

# KM plot
fit_km <- survfit(Surv(OS.time, OS) ~ risk_group, data = surv_subset)

ggsurvplot(fit_km,
           data = surv_subset,
           risk.table = FALSE,
           pval = TRUE,
           conf.int = FALSE,
           palette = c("#E7B800", "#2E9FDF"),
           title = "Survival curves by LASSO risk group")

table(surv_subset$LOY_group, surv_subset$risk_group)
ylow_counts <- table(surv_subset$LOY_group, surv_subset$risk_group)["Y_low", ]
print(ylow_counts)



# Make sure this runs successfully first
counts <- as.data.frame(table(surv_subset$LOY_group, surv_subset$risk_group))
colnames(counts) <- c("LOY_group", "Risk_group", "Count")
print(counts)
str(counts)

ggplot(counts, aes(x = LOY_group, y = Count, fill = Risk_group)) +
  geom_bar(stat = "identity", position = "fill") +  # position = "fill" makes it proportional
  scale_fill_manual(values = c("Low" = "#1B9AAA", "High" = "#D7263D")) +
  labs(
    title = "Proportion of Risk Groups within LOY Groups",
    x = "LOY Group",
    y = "Proportion",
    fill = "Risk Group"
  ) +
  theme_minimal(base_size = 16) +
  theme(
    plot.title = element_text(face = "bold", size = 18, hjust = 0.5),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(face = "bold"),
    legend.title = element_text(face = "bold")
  )



library(survival)
library(survminer)

# Final model coefficients at lambda.min
final_coef <- coef(cvfit, s = "lambda.min")
selected_genes <- rownames(final_coef)[final_coef[, 1] != 0]

# Calculate risk score
risk_score <- as.vector(x[, selected_genes] %*% final_coef[selected_genes, 1])
surv_subset$risk_score <- risk_score

# Median-based risk group
surv_subset$risk_group <- ifelse(surv_subset$risk_score > median(risk_score), "High", "Low")

# Kaplan-Meier fit
fit_km <- survfit(Surv(OS.time, OS) ~ risk_group, data = surv_subset)

# Publication-quality plot
ggsurvplot(
  fit_km,
  data = surv_subset,
  risk.table = FALSE,
  pval = TRUE,
  conf.int = FALSE,
  palette = c("#D7263D", "#1B9AAA"),  # Better contrast than yellow/blue
  size = 2,  # Thicker lines
  xlab = "Time (days)",
  ylab = "Overall Survival Probability",
  title = "Kaplan-Meier Survival Curves by LASSO Risk Group",
  
  ggtheme = theme_minimal(base_size = 16),
  font.main = c(20, "bold", "black"),
  font.x = c(18, "bold"),
  font.y = c(18, "bold"),
  font.tickslab = c(16, "bold"),
  font.legend = c(16, "bold"),
  legend.title = "Risk Group",
  legend.labs = c("High", "Low")
)


#--------------BOXPLOT OF SELECTED GENES----------------------



library(dplyr)
library(reshape2)
library(ggpubr)

# Your LASSO genes
lasso_genes <- c("WASF1", "MTHFD2", "PRC1", "MXD3", "CDKN2C", "CHMP1A")

# Prepare group dataframe
group_df <- male_surv %>%
  select(sample, LOY_score) %>%
  mutate(
    LOY_group_opt = ifelse(LOY_score > median(LOY_score, na.rm = TRUE), "Y_high", "Y_low"),
    sample_dot = sample
  )

# Ensure the genes are present in expression matrix
valid_genes <- intersect(lasso_genes, rownames(expr_final))

# Subset expression data for valid genes and samples
expr_subset <- expr_final[valid_genes, group_df$sample_dot]

# Convert to long format
expr_long <- melt(as.matrix(expr_subset))
colnames(expr_long) <- c("Gene", "Sample", "Expression")
expr_long$Sample <- as.character(expr_long$Sample)

# Merge with LOY group data
expr_long <- left_join(expr_long, group_df, by = c("Sample" = "sample_dot"))

# Plot all genes together with facets
p_all <- ggboxplot(expr_long,
                   x = "LOY_group_opt", y = "Expression", fill = "LOY_group_opt",
                   facet.by = "Gene", scales = "free_y",   # free y-scale per gene
                   palette = c("#1f77b4", "#ff7f0e")) +
  stat_compare_means(method = "wilcox.test", label = "p.format") +
  labs(
    x = "LOY Group",
    y = "Expression (TPM)",
    title = "Expression of LASSO Genes (Y_high vs Y_low)"
  ) +
  theme_minimal(base_size = 14) +  # Increase base font size
  theme(
    legend.position = "bottom",
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
    axis.title = element_text(face = "bold", size = 14),
    axis.text = element_text(face = "bold", size = 12),
    legend.title = element_text(face = "bold", size = 14),
    legend.text = element_text(face = "bold", size = 12),
    strip.text = element_text(face = "bold", size = 14),  # facet label bold
    plot.margin = margin(20, 20, 30, 20)  # top, right, bottom, left margins
  )

print(p_all)



#------------------OS BY EXPRESSION OF SELECTED GENES----------------






library(survival)
library(survminer)
library(dplyr)

for (gene in selected_genes) {
  gene_expr <- expr_final[gene, , drop = TRUE]
  
  gene_surv_df <- data.frame(
    Sample = names(gene_expr),
    Expression = as.numeric(gene_expr)
  ) %>%
    mutate(sample_dash = gsub("\\.", "-", Sample)) %>%
    left_join(male_surv_clean, by = "sample_dash") %>%
    filter(!is.na(OS.time), !is.na(OS), !is.na(Expression))
  
  if (nrow(gene_surv_df) < 10) next  # Skip small groups
  
  gene_surv_df <- gene_surv_df %>%
    mutate(ExprGroup = ifelse(Expression > median(Expression, na.rm = TRUE), "High", "Low"),
           ExprGroup = factor(ExprGroup, levels = c("Low", "High")))
  
  concord <- survConcordance(Surv(OS.time, OS) ~ Expression, data = gene_surv_df)
  c_index <- round(concord$concordance, 3)
  
  fit <- survfit(Surv(OS.time, OS) ~ ExprGroup, data = gene_surv_df)
  
  # Enhanced plot
  plot <- ggsurvplot(
    fit,
    data = gene_surv_df,
    pval = TRUE,
    conf.int = FALSE,
    risk.table = FALSE,
    palette = c("#1f77b4", "#ff7f0e"),  # Blue & Orange
    title = paste("Overall Survival by", gene, "Expression"),
    subtitle = paste("Concordance Index:", c_index),
    xlab = "Time (days)",
    ylab = "Survival Probability",
    legend.title = gene,
    legend.labs = c("Low Expression", "High Expression"),
    risk.table.title = "Number at Risk",
    risk.table.height = 0.25,
    surv.median.line = "hv",
    ggtheme = theme_minimal(base_size = 16) +
      theme(
        plot.title = element_text(face = "bold", hjust = 0.5, size = 18),
        plot.subtitle = element_text(face = "bold", hjust = 0.5, size = 14),
        axis.title = element_text(face = "bold", size = 16),
        axis.text = element_text(face = "bold", size = 14),
        legend.title = element_text(face = "bold", size = 14),
        legend.text = element_text(face = "bold", size = 13),
        strip.text = element_text(face = "bold", size = 13)
      )
  )
  
  print(plot)
  readline(prompt = "Press [Enter] to view next KM plot...")
}




#---------lasso summery------
# Safely extract coefficients
lasso_coef <- coef(cvfit, s = "lambda.min")

# Check if it's not NULL
if (!is.null(lasso_coef)) {
  lasso_df <- as.data.frame(as.matrix(lasso_coef))
  lasso_df$Gene <- rownames(lasso_df)
  colnames(lasso_df)[1] <- "Coefficient"
  
  # Filter only non-zero coefficients
  lasso_df <- lasso_df[lasso_df$Coefficient != 0, ]
  
  # List of selected genes
  selected_genes <- lasso_df$Gene
  print(selected_genes)
} else {
  cat("No coefficients found at lambda.min. Check your cvfit or input data.\n")
}


library(survival)
library(dplyr)

# Step 1: Extract LASSO-selected genes
lasso_coef <- coef(cvfit, s = "lambda.min")
lasso_df <- as.data.frame(as.matrix(lasso_coef))
lasso_df$Gene <- rownames(lasso_df)
colnames(lasso_df)[1] <- "Coefficient"
selected_genes <- lasso_df %>% filter(Coefficient != 0) %>% pull(Gene)

# Step 2: Prepare Cox model data
cox_data <- cbind(surv_subset[, c("OS.time", "OS")], expr_t[, selected_genes])

# Step 3: Fit multivariate Cox model
cox_formula <- as.formula(paste("Surv(OS.time, OS) ~", paste(selected_genes, collapse = " + ")))
cox_fit <- coxph(cox_formula, data = cox_data)

# Step 4: Get summary
cox_summary <- summary(cox_fit)

# Step 5: Extract individual gene significance
cox_table <- as.data.frame(cox_summary$coefficients)
cox_table$HR <- exp(cox_table$coef)
cox_table$Gene <- rownames(cox_table)
cox_table <- cox_table[, c("Gene", "HR", "coef", "se(coef)", "z", "Pr(>|z|)")]
colnames(cox_table) <- c("Gene", "HR", "Coefficient", "SE", "Z", "p_value")

# Step 6: Save to CSV
write.csv(cox_table, "LASSO_Cox_Summary.csv", row.names = FALSE)

# Step 7: Overall model statistics
overall_stats <- data.frame(
  Concordance_Index = round(cox_summary$concordance[1], 3),
  Concordance_SE = round(cox_summary$concordance[2], 3),
  Likelihood_Ratio_Test = cox_summary$logtest["test"],
  LRT_p_value = cox_summary$logtest["pvalue"],
  Wald_Test = cox_summary$waldtest["test"],
  Wald_p_value = cox_summary$waldtest["pvalue"],
  Score_Test = cox_summary$sctest["test"],
  Score_p_value = cox_summary$sctest["pvalue"],
  AIC = AIC(cox_fit)
)

# Save overall model summary
write.csv(overall_stats, "LASSO_Cox_Overall_Model_Stats.csv", row.names = FALSE)

# Print to console
print(cox_table)
print(overall_stats)

# Save cox_table as CSV
write.csv(cox_table, "cox_lasso_genes_results.csv", row.names = FALSE)

# Convert overall_stats to data frame (if it’s a single-row matrix)
overall_stats_df <- as.data.frame(t(overall_stats))

# Add row names as a column (optional, for clarity)
overall_stats_df$Metric <- rownames(overall_stats_df)

overall_stats_df <- overall_stats_df[, c("Metric", "C")]

# Rename the value column for clarity
colnames(overall_stats_df) <- c("Metric", "Value")

# Save overall_stats as CSV
write.csv(overall_stats_df, "cox_model_overall_stats.csv", row.names = FALSE)



#-----------------------FOREST PLOT UNIVARIATE------------------------

install.packages(c("survminer", "forestplot", "rms", "timeROC"))


library(forestplot)

# Prepare table text
hr_ci <- paste0(round(cox_table$HR, 2), " (", 
                round(cox_table$HR / exp(cox_table$SE), 2), "-", 
                round(cox_table$HR * exp(cox_table$SE), 2), ")")

tabletext <- cbind(
  c("Gene", cox_table$Gene),
  c("HR (95% CI)", hr_ci),
  c("p-value", signif(cox_table$p_value, 3))
)

# Forest plot





library(forestplot)
library(grid)  # For unit()

# Assuming you have columns: Gene, HR, Lower CI, Upper CI, and optionally p-value
# You can add a p-value column to `tabletext` if you want

tabletext <- cbind(
  c("Gene", as.character(cox_table$Gene)),
  c("HR (95% CI)", paste0(
    formatC(cox_table$HR, format = "f", digits = 2), " (",
    formatC(cox_table$HR / exp(cox_table$SE), format = "f", digits = 2), "-",
    formatC(cox_table$HR * exp(cox_table$SE), format = "f", digits = 2), ")"
  ))
)

forestplot(
  labeltext = tabletext,
  mean = c(NA, cox_table$HR),
  lower = c(NA, cox_table$HR / exp(cox_table$SE)),
  upper = c(NA, cox_table$HR * exp(cox_table$SE)),
  zero = 1,
  boxsize = 0.3,
  lineheight = unit(10, "mm"),
  graph.pos = 2,
  col = fpColors(
    box = "#1f77b4",       # Blue boxes
    line = "#1f77b4",      # Blue lines
    summary = "black"
  ),
  xlab = expression("Hazard Ratio (log scale)"),
  title = "Hazard Ratios of LASSO Genes",
  txt_gp = fpTxtGp(
    label = gpar(fontface = "bold", cex = 1.2),     # Left labels
    ticks = gpar(cex = 1.0),                        # X-axis ticks
    xlab = gpar(fontface = "bold", cex = 1.2),      # X-axis label
    title = gpar(fontface = "bold", cex = 1.4)      # Title
  ),
  lwd.ci = 2,
  ci.vertices = TRUE,
  ci.vertices.height = 0.1,
  line.margin = 0.1,
  clip = c(0.1, 10),         # HR axis limits (adjust if needed)
  is.summary = c(TRUE, rep(FALSE, nrow(cox_table)))
)







#--------------------ROC UNIVARIATE----------------


# Get linear predictor (risk score) from Cox model
# Assuming expr_t is transposed expression matrix with samples as rows and genes as columns
# Subset expr_t for the LASSO genes and samples in surv_subset

genes <- c("WASF1", "MTHFD2", "PRC1", "MXD3", "CDKN2C", "CHMP1A")

# Subset expression matrix for these genes and matching samples
expr_for_pred <- expr_t[surv_subset$sample, genes]

# Add these expression columns to surv_subset
surv_pred_data <- cbind(surv_subset, expr_for_pred)

# Now predict linear predictor (risk score)
lp <- predict(cox_fit, newdata = surv_pred_data, type = "lp")


library(timeROC)

# Time-dependent ROC (1, 3, 5 year)
time_roc <- timeROC(T = surv_subset$OS.time,
                    delta = surv_subset$OS,
                    marker = lp,
                    cause = 1,
                    times = c(365, 1095, 1825),
                    iid = TRUE)

plot(time_roc, time = 365, col = "red", title = FALSE, lwd = 2.5, cex.axis = 1.2, cex.lab = 1.3)
plot(time_roc, time = 1095, add = TRUE, col = "blue", lwd = 2.5)
plot(time_roc, time = 1825, add = TRUE, col = "darkgreen", lwd = 2.5)

title(main = "Time-dependent ROC Curves for LASSO-Based Prognostic Model", cex.main = 1.5, font.main = 2)

auc_vals <- round(time_roc$AUC, 3)
legend("bottomright",
       legend = c(paste0("1 Year AUC = ", auc_vals[1]),
                  paste0("3 Year AUC = ", auc_vals[2]),
                  paste0("5 Year AUC = ", auc_vals[3])),
       col = c("red", "blue", "darkgreen"),
       lty = 1,
       lwd = 2.5,
       bty = "n",
       cex = 1.2,
       text.font = 2)











#----------------HEATMAP OF LASSO GENES BASED ON RISK SCORES---------------

library(pheatmap)

# Extract expression for selected genes (genes x samples)
expr_selected <- expr_t[, selected_genes]

# Order samples by risk group (Low then High)
ord <- order(surv_subset$risk_group)
expr_ordered <- expr_selected[ord, ]

# Annotation for samples (risk group)
annotation_col <- data.frame(RiskGroup = factor(surv_subset$risk_group[ord]))
rownames(annotation_col) <- rownames(expr_ordered)

pheatmap(t(expr_ordered), # transpose so genes in rows, samples in columns
         annotation_col = annotation_col,
         show_rownames = TRUE,
         show_colnames = FALSE,
         cluster_rows = TRUE,
         cluster_cols = FALSE,
         main = "Heatmap of LASSO-selected Gene Expression")



#-----------------PCA------------------------




library(ggplot2)
library(ggfortify)

# PCA on expression matrix (samples x genes)
pca_res <- prcomp(expr_selected, scale. = TRUE)

# Create dataframe for plotting
pca_df <- data.frame(pca_res$x[, 1:2],
                     RiskGroup = surv_subset$risk_group)

# Enhanced, publication-ready PCA plot
ggplot(pca_df, aes(x = PC1, y = PC2, color = RiskGroup, fill = RiskGroup)) +
  geom_point(size = 3, alpha = 0.8) +
  stat_ellipse(type = "norm", level = 0.95, size = 1, alpha = 0.3) +   # 95% confidence ellipse
  scale_color_manual(values = c("High" = "#E64B35", "Low" = "#4DBBD5")) +
  scale_fill_manual(values = c("High" = "#E64B3555", "Low" = "#4DBBD555")) +  # semi-transparent fills
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 20, hjust = 0.5),
    plot.subtitle = element_text(face = "bold", size = 16, hjust = 0.5),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 13),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title = element_text(face = "bold", size = 14),
    axis.text = element_text(size = 12)
  ) +
  labs(
    title = "PCA of LASSO-selected Genes Expression",
    subtitle = "Samples clustered by Risk Group",
    x = paste0("PC1 (", round(summary(pca_res)$importance[2,1] * 100, 1), "% variance)"),
    y = paste0("PC2 (", round(summary(pca_res)$importance[2,2] * 100, 1), "% variance)"),
    color = "Risk Group",
    fill = "Risk Group"
  )


#vital status by risk grp

library(ggplot2)
library(dplyr)

# Assume you have `risk_group` and `vital_status` in `surv_subset`
# Example: risk_group = "High" / "Low", vital_status = "Alive" / "Dead"

# Count and calculate proportion
plot_data <- surv_subset %>%
  count(risk_group, vital_status.demographic) %>%
  group_by(risk_group) %>%
  mutate(prop = n / sum(n))

# Plot
ggplot(plot_data, aes(x = risk_group, y = prop, fill = vital_status.demographic)) +
  geom_bar(stat = "identity", position = "fill") +
  labs(x = "Risk Group", y = "Proportion", fill = "Vital Status") +
  scale_fill_manual(values = c("Alive" = "#4CAF50", "Dead" = "#F44336")) +
  theme_minimal() +
  ggtitle("Vital Status by Risk Group")



ggplot(surv_subset, aes(x = vital_status.demographic, y = risk_score, fill = vital_status.demographic)) +
  geom_violin(trim = FALSE) +
  geom_boxplot(width = 0.2, outlier.shape = NA, color = "black") +
  labs(x = "Vital Status", y = "Risk Score") +
  scale_fill_manual(values = c("Alive" = "#2196F3", "Dead" = "#E91E63")) +
  theme_minimal() +
  ggtitle("Risk Score by Vital Status")






# Calculate risk score (linear predictor) for each patient
risk_score <- predict(cox_fit, newdata = surv_pred_data, type = "lp")

# Combine with survival info
risk_df <- data.frame(
  Sample = rownames(surv_pred_data),
  Risk_Score = risk_score,
  Vital_Status = ifelse(surv_pred_data$OS == 1, "Dead", "Alive"),
  OS_Time = surv_pred_data$OS.time
)

# Sort by risk score
risk_df <- risk_df[order(risk_df$Risk_Score), ]
risk_df$Patient_ID <- 1:nrow(risk_df)

# Plot
library(ggplot2)
ggplot(risk_df, aes(x = Patient_ID, y = Risk_Score, color = Vital_Status)) +
  geom_point(size = 2.5) +
  scale_color_manual(values = c("Alive" = "blue", "Dead" = "red")) +
  labs(title = "Risk Score Distribution by Vital Status",
       x = "Patient (Ranked by Risk Score)",
       y = "Risk Score") +
  theme_minimal(base_size = 14)



ggplot(risk_df, aes(x = Patient_ID, y = OS_Time, color = Vital_Status)) +
  geom_point(size = 2.5) +
  scale_color_manual(values = c("Alive" = "blue", "Dead" = "red")) +
  labs(
       x = "Patient (Ranked by Risk Score)",
       y = "Survival Time (days)") +
  theme_minimal(base_size = 14)



# Add risk group
risk_df$Risk_Group <- ifelse(risk_df$Risk_Score >= median(risk_df$Risk_Score), "High Risk", "Low Risk")

# Example facet by Risk Group
ggplot(risk_df, aes(x = Patient_ID, y = Risk_Score, color = Risk_Group)) +
  geom_point(size = 2.5) +
  facet_wrap(~ Risk_Group, scales = "free_x") +
  scale_color_manual(values = c("Low Risk" = "darkgreen", "High Risk" = "firebrick")) +
  labs(title = "Risk Score by Group",
       x = "Patient (Ranked)",
       y = "Risk Score") +
  theme_minimal(base_size = 14)
