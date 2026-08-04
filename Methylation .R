
######Training set



setwd("E:/LOY_KIRP")






library(dplyr)
library(biomaRt)
library(tibble)  # for column_to_rownames()

#---------Step 1: Load data-------

metadat <- read.delim("GDC KIRP Clinical.tsv")
exp_data <- read.delim("GDC KIRP TPM.tsv")
surv_data <- read.delim("GDC KIRP Surv.tsv")
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
    labels = c("Y_high" = "High LOY", "Y_low" = "Low LOY")  # More descriptive labels
  ) +
  scale_x_discrete(
    labels = c("Y_high" = "High LOY", "Y_low" = "Low LOY")  # Consistent labeling
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
  palette = c("#1f77b4", "#ff7f0e"),
  title = "Overall Survival by Optimal LOY Cutoff"
)






#----------------------deg--------------

library(limma)

####deg

cut_os <- surv_cutpoint(
  male_surv_clean,
  time = "OS.time",
  event = "OS",
  variables = "LOY_score",
  minprop = 0.1
)
male_combined <- male_surv_clean  # or `male_surv_cut` if you want the categorized version directly


# Step 1: Assign LOY group using optimal cutoff
optimal_cutoff <- cut_os$cutpoint["LOY_score", "cutpoint"]
male_combined$LOY_group_cutoff <- ifelse(
  male_combined$LOY_score >= optimal_cutoff,
  "Y_high",
  "Y_low"
)
 
# Step 2: Harmonize sample IDs
male_combined$sample_dot <- gsub("-", ".", male_combined$sample_dash)
common_samples <- intersect(male_combined$sample_dot, colnames(expr_final))

# Step 3: Filter expression and metadata
expr_mat_filtered <- expr_final[, common_samples]
male_combined_filtered <- male_combined[male_combined$sample_dot %in% common_samples, ]
male_combined_filtered <- male_combined_filtered[match(common_samples, male_combined_filtered$sample_dot), ]

stopifnot(all(colnames(expr_mat_filtered) == male_combined_filtered$sample_dot))

# Step 4: DEG analysis (Y_high is baseline)
group_factor <- factor(male_combined_filtered$LOY_group_cutoff, levels = c("Y_high", "Y_low"))

design <- model.matrix(~ 0 + group_factor)
colnames(design) <- levels(group_factor)

fit <- lmFit(expr_mat_filtered, design)
contrast.matrix <- makeContrasts(Y_low_vs_Y_high = Y_low - Y_high, levels = design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

deg_results <- topTable(fit2, number = Inf, adjust.method = "BH")
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
ggplot(deg_results, aes(x = logFC, y = -log10(adj.P.Val), color = regulation)) +
  geom_point(alpha = 0.7, size = 1.5) +
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
    size = 3.5,
    box.padding = 0.3,
    point.padding = 0.2,
    max.overlaps = Inf,
    segment.color = "grey50"
  ) +
  labs(
    title = "Volcano Plot: Differential Expression Y_low vs Y_high",
    x = expression(Log[2]~Fold~Change~"(Y_low vs Y_high)"),
    y = expression(-log[10]~"Adjusted P-value"),
    color = "Gene Regulation"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.position = "right",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 11),
    axis.title = element_text(size = 13),
    axis.text = element_text(size = 11)
  )



#------------heatmap--------------


# Required libraries
library(ComplexHeatmap)
library(circlize)
library(dplyr)
library(RColorBrewer)
library(tibble)

#--- Step 1: Select Top 20 DEGs (10 up + 10 down) ---

top_degs_20 <- bind_rows(top_up, top_down)
top_deg_names <- rownames(top_degs_20)

#--- Step 2: Subset expression data for these genes ---

# Subset expr_final matrix to the 20 DEGs
heatmap_matrix <- expr_final[top_deg_names, common_samples]

# Optional: z-score normalize genes (rows) for better heatmap visualization
heatmap_matrix_scaled <- t(scale(t(as.matrix(heatmap_matrix))))

#--- Step 3: Create annotation for LOY groups ---

# Create LOY group annotation from matched metadata
group_annotation <- male_combined_filtered %>%
  dplyr::select(sample_dot, LOY_group) %>%
  column_to_rownames("sample_dot")

# Define colors
group_colors <- list(LOY_group = c("Y_low" = "#1f78b4", "Y_high" = "#ff7f0e"))

# Create annotation object
col_ha <- HeatmapAnnotation(
  df = group_annotation,
  col = group_colors,
  annotation_legend_param = list(title = "LOY Group")
)

#--- Step 4: Plot Heatmap ---

Heatmap(
  heatmap_matrix_scaled,
  name = "Z-score",
  top_annotation = col_ha,
  show_row_names = TRUE,
  show_column_names = FALSE,
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  column_split = group_annotation$LOY_group,
  row_names_gp = gpar(fontsize = 10),
  column_title = "Top 20 DEGs by LOY Group",
  heatmap_legend_param = list(title = "Expression Z-score")
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
barplot(ego_up, showCategory = 10, title = "GO Biological Process Enriched in Upregulated Genes")





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
  title = "GO Biological Process Enriched in Downregulated Genes"
)




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
barplot(kk_up, showCategory = 10, title = "KEGG Pathways Enriched in Upregulated Genes")



kk_down <- enrichKEGG(
  gene = entrez_down$ENTREZID,
  organism = "hsa",
  pAdjustMethod = "BH",
  qvalueCutoff = 0.05
)
barplot(kk_down, showCategory = 10, title = "KEGG Pathways Enriched in Downregulated Genes")




#BiocManager::install("ReactomePA")
library(ReactomePA)




reactome_up <- enrichPathway(
  gene = entrez_up$ENTREZID,
  organism = "human",
  pAdjustMethod = "BH",
  qvalueCutoff = 0.05,
  readable = TRUE
)

dotplot(reactome_up, showCategory = 10, title = "Reactome Pathways Enriched in Upregulated Genes")





reactome_down <- enrichPathway(
  gene = entrez_down$ENTREZID,
  organism = "human",
  pAdjustMethod = "BH",
  qvalueCutoff = 0.05,
  readable = TRUE
)

dotplot(reactome_down, showCategory = 10, title = "Reactome Pathways Enriched in Downregulated Genes")





#-----------hallmark---------------

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
        color = "p.adjust",           # color by adjusted p-value
        x = "NES",                 # use Normalized Enrichment Score on x-axis
        title = "Top 25 Hallmark Pathways (GSEA)") +
  ggplot2::theme_classic() +
  ggplot2::scale_color_gradient(low = "red", high = "blue")



#----------selected hallmark analysis in Ylow---------------

# 1. Filter GSEA result for NES > 0 and adjusted p < 0.05
gsea_pos <- gsea_res@result
gsea_pos <- gsea_pos[gsea_pos$NES > 0 & gsea_pos$p.adjust < 0.05, ]

gsea_neg <- gsea_res@result %>%
  filter(NES < 0 & p.adjust < 0.05)



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

plot(cvfit, main = "Cross-Validation Curve (LASSO)", xlab = "log(Lambda)", ylab = "Mean CV Error")
abline(v = log(cvfit$lambda.min), col = "red", lty = 2)
abline(v = log(cvfit$lambda.1se), col = "blue", lty = 2)
legend("topright", legend = c("lambda.min", "lambda.1se"),
       col = c("red", "blue"), lty = 2, bty = "n")


# Extract selected genes
lasso_coef <- coef(cvfit, s = "lambda.min")
selected_genes <- rownames(lasso_coef)[lasso_coef[,1] != 0]
selected_genes


# Save selected LASSO genes to CSV
write.csv(data.frame(Gene = selected_genes),
          file = "LASSO_Selected_Genes.csv",
          row.names = FALSE)


#--------------BOXPLOT OF SELECTED GENES----------------------
# Load required packages
library(dplyr)
library(reshape2)
library(ggpubr)

# Define the LASSO gene list
lasso_genes <- c("WASF1", "MTHFD2", "PRC1", "MXD3", "CDKN2C", "CHMP1A")

# Extract LOY group info from metadata
group_df <- male_combined_filtered %>%
  dplyr::select(sample, LOY_group_cutoff) %>%
  dplyr::rename(
    sample_dot = sample,
    LOY_group_opt = LOY_group_cutoff
  )

# Ensure the genes are present in expression matrix
valid_genes <- intersect(lasso_genes, rownames(expr_final))

# Subset expression data for valid LASSO genes and samples
expr_subset <- expr_final[valid_genes, group_df$sample_dot]

# Convert to long format
expr_long <- melt(as.matrix(expr_subset))
colnames(expr_long) <- c("Gene", "Sample", "Expression")

# Convert Sample ID to character for joining
expr_long$Sample <- as.character(expr_long$Sample)

# Merge with LOY group data
expr_long <- left_join(expr_long, group_df, by = c("Sample" = "sample_dot"))

# Create boxplot
ggboxplot(expr_long,
          x = "LOY_group_opt", y = "Expression", fill = "LOY_group_opt",
          facet.by = "Gene", scales = "free_y",
          palette = c("#1f77b4", "#ff7f0e")) +
  stat_compare_means(method = "wilcox.test", label = "p.format") +
  labs(
    x = "LOY Group",
    y = "Expression (TPM)",
    title = "Expression of LASSO Genes (Y_high vs Y_low)"
  ) +
  theme_minimal()






# Load required packages
library(dplyr)
library(reshape2)
library(ggpubr)

# Define the LASSO gene list
lasso_genes <- c("WASF1", "MTHFD2", "PRC1", "MXD3", "CDKN2C", "CHMP1A")

# Extract LOY group info from metadata
group_df <- male_combined_filtered %>%
  dplyr::select(sample, LOY_group_cutoff) %>%
  dplyr::rename(
    sample_dot = sample,
    LOY_group_opt = LOY_group_cutoff
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

# ---- Plot one boxplot per gene ----
for (g in valid_genes) {
  gene_data <- expr_long %>% filter(Gene == g)
  
  p <- ggboxplot(gene_data,
                 x = "LOY_group_opt", y = "Expression", fill = "LOY_group_opt",
                 palette = c("#1f77b4", "#ff7f0e")) +
    stat_compare_means(method = "wilcox.test", label = "p.format") +
    labs(
      x = "LOY Group",
      y = "Expression (TPM)",
      title = paste("Expression of", g, "(Y_high vs Y_low)")
    ) +
    theme_minimal()
  
  print(p)  # Show the plot
  
  
}




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
  
  plot <- ggsurvplot(
    fit,
    data = gene_surv_df,
    pval = TRUE,
    conf.int = TRUE,
    risk.table = TRUE,
    palette = c("#1f77b4", "#ff7f0e"),
    title = paste("Overall Survival by", gene, "Expression"),
    subtitle = paste("Concordance Index:", c_index),
    xlab = "Time (days)",
    ylab = "Survival Probability",
    legend.title = gene,
    legend.labs = c("Low Expression", "High Expression"),
    risk.table.title = "Number at risk",
    risk.table.height = 0.25,
    surv.median.line = "hv",
    ggtheme = theme_minimal(base_size = 14)
  )
  
  print(plot)  # Display in RStudio Viewer
  # Optional: pause between plots
  readline(prompt = "Press [enter] to view next KM plot...")
}














library(dplyr)
library(survival)
library(survminer)



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
write.csv(cox_table, "LASSO_Multivariate_Cox_Summary.csv", row.names = FALSE)

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

#------km-------
# Use the same selected genes and their coefficients
lasso_df <- as.data.frame(as.matrix(coef(cvfit, s = "lambda.min")))
lasso_df$Gene <- rownames(lasso_df)
colnames(lasso_df)[1] <- "Coefficient"
lasso_df <- lasso_df[lasso_df$Coefficient != 0, ]

# Get expression data for selected genes
expr_selected <- expr_t[, lasso_df$Gene]

# Compute risk score = sum(expression * coefficient) for each patient
risk_score <- as.numeric(as.matrix(expr_selected) %*% lasso_df$Coefficient)
names(risk_score) <- rownames(expr_selected)

# Combine with survival data
risk_df <- data.frame(
  sample = names(risk_score),
  risk_score = risk_score
)

# Add LOY-compatible sample ID if needed
risk_df$sample_dash <- gsub("\\.", "-", risk_df$sample)

# Merge with survival info
merged_risk <- left_join(risk_df, male_surv_clean, by = "sample_dash")

# Remove missing OS/time
merged_risk <- merged_risk %>% filter(!is.na(OS.time) & !is.na(OS))

# Create risk group (High if above median)
median_cutoff <- median(merged_risk$risk_score, na.rm = TRUE)
merged_risk <- merged_risk %>%
  mutate(RiskGroup = ifelse(risk_score > median_cutoff, "High", "Low")) %>%
  mutate(RiskGroup = factor(RiskGroup, levels = c("Low", "High")))


# Fit KM survival model
fit_risk <- survfit(Surv(OS.time, OS) ~ RiskGroup, data = merged_risk)

# Plot
library(survminer)
ggsurvplot(
  fit_risk,
  data = merged_risk,
  pval = TRUE,
  conf.int = TRUE,
  risk.table = TRUE,
  palette = c("#1f77b4", "#ff7f0e"),
  title = "Kaplan-Meier Survival by LASSO Risk Score",
  xlab = "Time (days)",
  ylab = "Survival Probability",
  legend.title = "Risk Group",
  legend.labs = c("Low Risk", "High Risk"),
  surv.median.line = "hv",
  ggtheme = theme_minimal(base_size = 14)
)

median_cutoff <- median(risk_df$risk_score, na.rm = TRUE)
risk_df$risk_group <- ifelse(risk_df$risk_score > median_cutoff, "High Risk", "Low Risk")
summary(risk_df$risk_group)
table(risk_df$risk_group)
# Merge to bring in LOY group info
merged_df <- left_join(risk_df, male_combined[, c("sample_dash", "LOY_group_cutoff")], by = "sample_dash")
write.csv(merged_df, "merged_df.csv", row.names = FALSE)

# Check for NAs just in case
table(is.na(merged_df$LOY_group_cutoff))

# Cross-tabulate
table(merged_df$risk_group, merged_df$LOY_group_cutoff)


library(ggplot2)

ggplot(merged_df, aes(x = LOY_group_cutoff, fill = risk_group)) +
  geom_bar(position = "fill") +
  ylab("Proportion") +
  scale_fill_manual(values = c("High Risk" = "red", "Low Risk" = "steelblue")) +
  theme_minimal() +
  ggtitle("Risk Group Proportion by LOY Group")





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
forestplot(labeltext = tabletext,
           mean = c(NA, cox_table$HR),
           lower = c(NA, cox_table$HR / exp(cox_table$SE)),
           upper = c(NA, cox_table$HR * exp(cox_table$SE)),
           zero = 1,
           boxsize = 0.2,
           lineheight = unit(8, "mm"),
           col = fpColors(box = "darkblue", line = "blue", summary = "black"),
           title = "Hazard Ratios of LASSO Genes")


#------------------KM------------------------


library(survival)
library(survminer)

top_gene <- cox_table$Gene[which.min(cox_table$p_value)]
expr_vec <- expr_t[, top_gene]

# Median cutoff
km_group <- ifelse(expr_vec > median(expr_vec), "High", "Low")
km_group <- factor(km_group, levels = c("Low", "High"))

surv_obj <- Surv(surv_subset$OS.time, surv_subset$OS)

fit <- survfit(surv_obj ~ km_group)

# KM Plot
ggsurvplot(fit,
           data = data.frame(km_group),
           pval = TRUE,
           risk.table = TRUE,
           title = paste("Kaplan-Meier for", top_gene),
           palette = c("red", "blue"))



#--------------------ROC UNIVARIATE----------------

library(timeROC)

# Construct linear predictor from model
cox_fit <- coxph(Surv(OS.time, OS) ~ ., data = data_nomogram)
lp <- predict(cox_fit, type = "lp")

# Time-dependent ROC (1, 3, 5 year)
time_roc <- timeROC(T = surv_subset$OS.time,
                    delta = surv_subset$OS,
                    marker = lp,
                    cause = 1,
                    times = c(365, 1095, 1825),
                    iid = TRUE)

# Plot
plot(time_roc, time = 365, col = "red", title = TRUE)
plot(time_roc, time = 1095, add = TRUE, col = "blue")
plot(time_roc, time = 1825, add = TRUE, col = "darkgreen")
legend("bottomright", legend = c("1 Year", "3 Year", "5 Year"),
       col = c("red", "blue", "darkgreen"), lty = 1)





#-----------------PCA------------------------


library(ggplot2)
library(ggfortify)

# PCA on expression matrix (samples x genes)
pca_res <- prcomp(expr_selected, scale. = TRUE)

# Create dataframe for plotting
pca_df <- data.frame(pca_res$x[, 1:2],
                     RiskGroup = surv_subset$risk_group)

# Publication-ready PCA plot with ellipses
ggplot(pca_df, aes(x = PC1, y = PC2, color = RiskGroup, fill = RiskGroup)) +
  geom_point(size = 3, alpha = 0.8) +
  stat_ellipse(type = "norm", level = 0.95, size = 1, alpha = 0.3) +   # 95% confidence ellipse
  scale_color_manual(values = c("High" = "#E64B35", "Low" = "#4DBBD5")) +
  scale_fill_manual(values = c("High" = "#E64B3555", "Low" = "#4DBBD555")) +  # semi-transparent fills
  theme_minimal(base_size = 14) +
  theme(
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 11),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title = element_text(face = "bold")
  ) +
  labs(
    title = "PCA of LASSO-selected Genes Expression",
    subtitle = "Samples clustered by Risk Group",
    x = paste0("PC1 (", round(summary(pca_res)$importance[2,1] * 100, 1), "% variance)"),
    y = paste0("PC2 (", round(summary(pca_res)$importance[2,2] * 100, 1), "% variance)"),
    color = "Risk Group",
    fill = "Risk Group"
  )




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







###-----multivariate cox-------


library(dplyr)
library(survival)

# Recreate relevant variables
surv_subset <- surv_subset %>%
  mutate(
    stage = as.factor(ajcc_pathologic_stage.diagnoses),
    smoking = as.numeric(cigarettes_per_day.exposures)
  )

# Remove NAs to avoid model fitting issues
surv_model_data <- surv_subset %>%
  filter(!is.na(stage), !is.na(smoking), !is.na(risk_score))

# Fit multivariate Cox model with risk_score
multi_cox <- coxph(Surv(OS.time, OS) ~ risk_score + stage + smoking,
                   data = surv_model_data)

# Summarize Cox results
summary_cox <- summary(multi_cox)

# Build summary dataframe
cox_df <- data.frame(
  Variable = rownames(summary_cox$coefficients),
  HR = round(summary_cox$coefficients[, "exp(coef)"], 3),
  CI.lower = round(summary_cox$conf.int[, "lower .95"], 3),
  CI.upper = round(summary_cox$conf.int[, "upper .95"], 3),
  p.value = signif(summary_cox$coefficients[, "Pr(>|z|)"], 3)
)

# Clean variable names
cox_df$Variable <- gsub("risk_score", "Risk Score", cox_df$Variable)
cox_df$Variable <- gsub("stage", "Stage: ", cox_df$Variable)
cox_df$Variable <- gsub("smoking", "Smoking", cox_df$Variable)

# Save results
write.csv(cox_df, "Multivariate_Cox_Results.csv", row.names = FALSE)


library(ggplot2)
library(readr)

cox_df <- read_csv("Multivariate_Cox_Results.csv")

# Order for plotting
cox_df <- cox_df %>%
  arrange(desc(HR)) %>%
  mutate(Variable = factor(Variable, levels = rev(Variable)))

# Add label column with HR and p-value
cox_df$label <- paste0("HR=", cox_df$HR, ", p=", cox_df$p.value)

# Draw forest plot
ggplot(cox_df, aes(x = HR, y = Variable)) +
  geom_point(color = "blue", size = 3) +
  geom_errorbarh(aes(xmin = CI.lower, xmax = CI.upper), height = 0.2, color = "gray30") +
  geom_vline(xintercept = 1, linetype = "dashed", color = "red") +
  xlab("Hazard Ratio (HR)") + ylab("") +
  theme_minimal(base_size = 14) +
  ggtitle("Multivariate Cox Regression - Forest Plot") +
  scale_x_log10() +
  geom_text(aes(label = label), hjust = -0.1, size = 3.5) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = 0.5))



library(timeROC)

# Define time points in days
time_points <- c(365, 1095, 1825)

# Compute time-dependent ROC
roc_obj <- timeROC(
  T = surv_model_data$OS.time,
  delta = surv_model_data$OS,
  marker = surv_model_data$risk_score,
  cause = 1,
  times = time_points,
  iid = TRUE
)

# Plot the 1-year ROC curve
plot(roc_obj, time = 365, col = "blue", title = FALSE, lwd = 2)

# Add 3-year and 5-year curves
plot(roc_obj, time = 1095, add = TRUE, col = "green", lwd = 2)
plot(roc_obj, time = 1825, add = TRUE, col = "red", lwd = 2)

# Add AUC values as text
text(0.65, 0.20, paste0("AUC at 1-year: ", round(roc_obj$AUC[1], 3)), col = "blue", cex = 1)
text(0.65, 0.15, paste0("AUC at 3-year: ", round(roc_obj$AUC[2], 3)), col = "green", cex = 1)
text(0.65, 0.10, paste0("AUC at 5-year: ", round(roc_obj$AUC[3], 3)), col = "red", cex = 1)

# Add legend
legend("bottomright",
       legend = c("1-year", "3-year", "5-year"),
       col = c("blue", "green", "red"),
       lty = 1, lwd = 2, cex = 0.9)

# Add title
title("Time-dependent ROC Curves with AUC")



#-------DNA Methylation--------
# Load methylation data (rows = probes, columns = samples)
library(data.table)
meth_data <- fread("TCGA-KIRP.methylation450.tsv", data.table = FALSE)
rownames(meth_data) <- meth_data[,1]
meth_data <- meth_data[,-1]

#subset only tumors
tumor_cols <- grep("-01A", colnames(meth_data))  # TCGA tumor barcode
meth_data <- meth_data[, tumor_cols]

# Install if not already
if (!requireNamespace("IlluminaHumanMethylation450kanno.ilmn12.hg19", quietly = TRUE)) {
  BiocManager::install("IlluminaHumanMethylation450kanno.ilmn12.hg19")
}
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)

# Extract annotation
ann <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)

# Load libraries
# Load dplyr just in case
library(data.table)
library(dplyr)
library(reshape2)
library(ggplot2)
library(ggpubr)

# 1. Load methylation data (make sure path and filename are correct)
meth_data <- fread("TCGA-KIRP.methylation450.tsv", data.table = FALSE)
rownames(meth_data) <- meth_data[,1]
meth_data <- meth_data[,-1]

# Subset tumor samples only (TCGA barcode with -01A)
tumor_cols <- grep("-01A", colnames(meth_data))
meth_data <- meth_data[, tumor_cols]

# 2. Load Illumina450k annotation (install if needed)
if (!requireNamespace("IlluminaHumanMethylation450kanno.ilmn12.hg19", quietly = TRUE)) {
  BiocManager::install("IlluminaHumanMethylation450kanno.ilmn12.hg19")
}
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
ann <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)

# Convert annotation to tibble for dplyr compatibility
ann_df <- as_tibble(as.data.frame(ann))

# 3. Define genes of interest
genes_of_interest <- c("MTHFD2", "WASF1", "CHMP1A", "MDX3", "PRC1", "CDKN2C")

# 4. Extract promoter probes for these genes (TSS200 or TSS1500)
promoter_probes <- ann_df %>%
  filter(grepl("TSS200|TSS1500", UCSC_RefGene_Group)) %>%
  filter(sapply(strsplit(UCSC_RefGene_Name, ";"), function(x) any(x %in% genes_of_interest))) %>%
  dplyr::select(Name, UCSC_RefGene_Name)


probe_ids <- promoter_probes$Name

# 5. Subset methylation data to promoter probes
meth_subset <- meth_data[probe_ids, , drop = FALSE]

# 6. Match methylation sample names with your LOY group sample names
# (Assuming male_combined is your metadata dataframe with LOY_group and sample_dash columns)
# Make sure sample names match in format, e.g. both TCGA-XX-XXXX-01A

matched_samples <- intersect(colnames(meth_subset), male_combined$sample_dash)

# Filter methylation and metadata for matched samples only
meth_subset <- meth_subset[, matched_samples, drop = FALSE]
male_combined_matched <- male_combined %>% filter(sample_dash %in% matched_samples)

# 7. Melt methylation matrix to long format for plotting
melted <- melt(as.matrix(meth_subset))
# Explicitly call reshape2’s melt on your matrix
library(reshape2)  # make sure reshape2 is loaded

melted <- reshape2::melt(as.matrix(meth_subset),
                         varnames = c("Probe", "Sample"),
                         value.name = "BetaValue")

colnames(melted) <- c("Probe", "Sample", "Beta")

# 8. Join gene annotation to melted data
melted$Probe <- as.character(melted$Probe)
melted <- left_join(melted, promoter_probes, by = c("Probe" = "Name"))

# 9. Add LOY group info from matched metadata
melted$LOY_group <- male_combined_matched$LOY_group[match(melted$Sample, male_combined_matched$sample_dash)]

# 10. Clean NA Beta values if any
melted_clean <- melted %>% filter(!is.na(Beta))

# 11. Simplify gene names (take first gene if multiple)
melted_clean <- melted_clean %>%
  mutate(Gene = sapply(strsplit(UCSC_RefGene_Name, ";"), `[`, 1))

# 12. Plot promoter methylation grouped by LOY and faceted by gene
print(
  ggplot(melted_clean, aes(x = LOY_group, y = Beta, fill = LOY_group)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.2, alpha = 0.5, size = 1) +
    facet_wrap(~ Gene, scales = "free_y") +
    stat_compare_means(method = "wilcox.test", label = "p.signif") +
    labs(title = "Promoter Methylation in Y_high vs Y_low Groups",
         y = "Beta Value", x = "LOY Group") +
    theme_bw()
)
##Compare gene body methylation between LOY-high vs LOY-low.

# Assume ann_df is your annotation data.frame
gene_body_probes <- ann_df %>%
  filter(grepl("Body", UCSC_RefGene_Group)) %>%
  filter(sapply(strsplit(as.character(UCSC_RefGene_Name), ";"), function(genes) any(genes %in% genes_of_interest))) %>%
  dplyr::select(Name, UCSC_RefGene_Name)

beta_gene_body <- meth_data[rownames(meth_data) %in% gene_body_probes$Name, ]

# Melt and annotate
library(dplyr)
library(reshape2)

# Add probe and gene name
beta_df <- cbind(Probe = rownames(beta_gene_body), as.data.frame(beta_gene_body)) %>%
  left_join(gene_body_probes, by = c("Probe" = "Name"))

# Reshape
# Explicitly call reshape2::melt on your data.frame
melted_body <- reshape2::melt(
  beta_df,
  id.vars      = c("Probe", "UCSC_RefGene_Name"),
  variable.name = "Sample",
  value.name    = "Beta"
)

# Compute average gene body methylation per gene per sample
gene_body_meth <- melted_body %>%
  group_by(UCSC_RefGene_Name, Sample) %>%
  summarise(Mean_Beta = mean(Beta, na.rm = TRUE), .groups = "drop")

# Add LOY group info
gene_body_meth$LOY_group <- male_combined$LOY_group[
  match(gene_body_meth$Sample, male_combined$sample_dash)
]

# Clean
gene_body_meth <- gene_body_meth %>% filter(!is.na(Mean_Beta), !is.na(LOY_group))

# Plot
library(ggplot2)
library(ggpubr)

ggplot(gene_body_meth, aes(x = LOY_group, y = Mean_Beta, fill = LOY_group)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.5, size = 1) +
  facet_wrap(~ UCSC_RefGene_Name, scales = "free_y") +
  stat_compare_means(method = "wilcox.test", label = "p.signif") +
  labs(title = "Gene Body Methylation by LOY Group", y = "Mean Beta", x = "LOY Group") +
  theme_bw()

library(dplyr)

gene_body_meth %>%
  group_by(UCSC_RefGene_Name, LOY_group) %>%
  summarise(
    mean_beta = mean(Mean_Beta, na.rm = TRUE),
    median_beta = median(Mean_Beta, na.rm = TRUE),
    sd_beta = sd(Mean_Beta, na.rm = TRUE),
    n = n()
  ) %>%
  print(n = Inf)
####takde only sig probs
run_probe_test <- function(df) {
  # Skip if Probe is NA or empty
  if(is.na(df$Probe[1]) || df$Probe[1] == "") {
    return(data.frame(
      Probe = NA_character_,
      UCSC_RefGene_Name = ifelse(length(df$UCSC_RefGene_Name) > 0, df$UCSC_RefGene_Name[1], NA_character_),
      p_value = NA_real_,
      mean_beta_high = NA_real_,
      mean_beta_low = NA_real_,
      delta_beta = NA_real_
    ))
  }
  
  # Remove NAs in LOY_group and Beta
  valid_idx <- !is.na(df$LOY_group) & !is.na(df$Beta)
  df_sub <- df[valid_idx, ]
  
  # Make factor of LOY_group in subset
  lo <- factor(df_sub$LOY_group)
  
  # Check exactly two groups present
  if(length(levels(lo)) != 2 || length(unique(lo)) != 2) {
    return(data.frame(
      Probe = df$Probe[1],
      UCSC_RefGene_Name = df$UCSC_RefGene_Name[1],
      p_value = NA_real_,
      mean_beta_high = NA_real_,
      mean_beta_low = NA_real_,
      delta_beta = NA_real_
    ))
  }
  
  # Run test safely
  test <- wilcox.test(df_sub$Beta ~ lo)
  
  mean_high <- mean(df_sub$Beta[lo == "Y_high"], na.rm = TRUE)
  mean_low <- mean(df_sub$Beta[lo == "Y_low"], na.rm = TRUE)
  
  data.frame(
    Probe = df$Probe[1],
    UCSC_RefGene_Name = df$UCSC_RefGene_Name[1],
    p_value = test$p.value,
    mean_beta_high = mean_high,
    mean_beta_low = mean_low,
    delta_beta = mean_high - mean_low
  )
}

# Assuming you have a long‐format data frame named `melted_clean` (or `melted_body`)
# with columns: Probe, UCSC_RefGene_Name, LOY_group, Beta

# 1. Split into a list by Probe
probe_groups <- split(melted_clean, melted_clean$Probe)

# 2. Run your Wilcoxon‐based tester on each probe
probe_pvals <- do.call(rbind, lapply(probe_groups, run_probe_test))

# 3. Adjust p‐values for multiple testing
probe_pvals$adj_p_value <- p.adjust(probe_pvals$p_value, method = "BH")

# 4. Look at only the significant probes
sig_probes <- probe_pvals %>% filter(!is.na(adj_p_value) & adj_p_value < 0.05)

print(sig_probes)

nrow(sig_probes)
sig_probes
##location of the probes
# 1. Extract the significant probe IDs
sig_ids <- sig_probes$Probe

# 2. Pull their annotation from ann_df
# Make sure dplyr is loaded
library(dplyr)

probe_annotation <- ann_df %>%
  filter(Name %in% sig_ids) %>%
  dplyr::select(
    Probe               = Name,
    Gene                = UCSC_RefGene_Name,
    Region              = UCSC_RefGene_Group,
    CpG_Island_Relation = Relation_to_Island,
    Chromosome          = chr,
    Coordinate          = pos
  )

print(probe_annotation)
write.csv(probe_annotation, "Significant_Probe_Annotation.csv", row.names = FALSE)


####
library(ggplot2)
library(ggpubr)
library(dplyr)

# Define the CHMP1A probe
chmp1a_probe <- "cg05441039"

# Merge LOY group from metadata (male_combined) into methylation data (melted_body)
chmp1a_df <- melted_body %>%
  filter(Probe == chmp1a_probe) %>%
  mutate(LOY_group = male_combined$LOY_group[match(Sample, male_combined$sample_dash)]) %>%
  filter(!is.na(LOY_group)) %>%
  mutate(LOY_group = factor(LOY_group, levels = c("Y_low", "Y_high")))

# Plot
ggplot(chmp1a_df, aes(x = LOY_group, y = Beta, fill = LOY_group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.8) +
  geom_jitter(width = 0.1, alpha = 0.5) +
  scale_fill_manual(values = c("#E69F00", "#56B4E9")) +
  stat_compare_means(method = "wilcox.test", label = "p.format") +
  labs(
    title = paste("Methylation of CHMP1A Probe:", chmp1a_probe),
    x = "LOY Group", y = "Beta"
  ) +
  theme_minimal(base_size = 14)
####
# Define the CHMP1A probe
mthfd2_probe <- "cg22704057"

mthfd2_df <- melted_clean %>% 
  filter(Probe == mthfd2_probe) %>%
  mutate(LOY_group = factor(LOY_group, levels = c("Y_low", "Y_high")))

ggplot(mthfd2_df, aes(x = LOY_group, y = Beta, fill = LOY_group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.8) +
  geom_jitter(width = 0.1, alpha = 0.5) +
  scale_fill_manual(values = c("#E69F00", "#56B4E9")) +
  stat_compare_means(method = "wilcox.test", label = "p.format") +
  labs(
    title = paste("Methylation of MTHFD2 Probe:", mthfd2_probe),
    x = "LOY Group", y = "Beta"
  ) +
  theme_minimal(base_size = 14)

##
# Define PRC1 probe
prc1_probe <- "cg01407062"

# Filter promoter methylation data for PRC1 probe
prc1_df <- melted_clean %>%
  filter(Probe == prc1_probe) %>%
  mutate(LOY_group = factor(LOY_group, levels = c("Y_low", "Y_high")))

# Plot
ggplot(prc1_df, aes(x = LOY_group, y = Beta, fill = LOY_group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.8) +
  geom_jitter(width = 0.1, alpha = 0.5) +
  scale_fill_manual(values = c("#E69F00", "#56B4E9")) +
  stat_compare_means(method = "wilcox.test", label = "p.format") +
  labs(
    title = paste("Methylation of PRC1 Probe:", prc1_probe),
    x = "LOY Group", y = "Beta"
  ) +
  theme_minimal(base_size = 14)

library(dplyr)

# Define your genes of interest:
genes_interest <- c("WASF1", "MXD3", "CDKN2C")

# 1. Filter the probe results dataframe for these genes
filtered_probes <- probe_results_df %>%
  filter(grepl(paste(genes_interest, collapse = "|"), UCSC_RefGene_Name))

# 2. For each gene, keep the probe with the smallest adjusted p-value (most significant)
top_probes <- filtered_probes %>%
  group_by(Gene = sapply(strsplit(UCSC_RefGene_Name, ";"), `[`, 1)) %>%
  slice_min(order_by = adj_p_value, n = 1, with_ties = FALSE) %>%
  ungroup()

# 3. Select and print relevant columns
top_probes %>%
  dplyr::select(Probe, UCSC_RefGene_Name, mean_beta_high, mean_beta_low, delta_beta, p_value, adj_p_value) -> df_to_print

print(df_to_print, n = Inf)

# Filter all probes annotated as MXD3
mxd3_probes <- probe_pvals %>% 
  filter(grepl("MXD3", UCSC_RefGene_Name))

# Pick the probe with the smallest p-value
top_mxd3_probe <- mxd3_probes %>% 
  slice_min(order_by = p_value, n = 1, with_ties = FALSE)

# Now combine with your previous top_probes for WASF1 and CDKN2C
top_probes_all <- bind_rows(
  top_probes,  # previous WASF1 and CDKN2C
  top_mxd3_probe
)

# Select relevant columns and print
top_probes_all_df <- as_tibble(top_probes_all)

top_probes_all_df %>%
  dplyr::select(Probe, UCSC_RefGene_Name, mean_beta_high, mean_beta_low, delta_beta, p_value, adj_p_value) %>%
  base::print(n = Inf)

library(ggplot2)
library(dplyr)
library(ggpubr)

# Define probes for WASF1 and CDKN2C
probes <- c("cg05715287", "cg01145389", "cg09487733")

# Function to plot methylation for a given probe
plot_methylation <- function(probe_id) {
  df <- melted_clean %>%
    filter(Probe == probe_id) %>%
    mutate(LOY_group = factor(LOY_group, levels = c("Y_low", "Y_high")))
  
  gene_name <- unique(df$UCSC_RefGene_Name)
  
  ggplot(df, aes(x = LOY_group, y = Beta, fill = LOY_group)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.8) +
    geom_jitter(width = 0.1, alpha = 0.5) +
    scale_fill_manual(values = c("#E69F00", "#56B4E9")) +
    stat_compare_means(method = "wilcox.test", label = "p.format") +
    labs(
      title = paste("Methylation of Probe:", probe_id, "Gene(s):", gene_name),
      x = "LOY Group", y = "Beta"
    ) +
    theme_minimal(base_size = 14)
}

# Plot for each probe
for (p in probes) {
  print(plot_methylation(p))
}

#
# Load dplyr
library(dplyr)

# Define your 3 probes
your_probes <- c("cg05715287", "cg01145389", "cg09487733")

# Check if probes exist in ann_df
valid_probes <- your_probes[your_probes %in% ann_df$Name]
missing_probes <- setdiff(your_probes, valid_probes)

# Report missing probes
if (length(missing_probes) > 0) {
  message("Missing probe(s) in annotation: ", paste(missing_probes, collapse = ", "))
}

probe_annotation <- ann_df %>%
  filter(Name %in% valid_probes) %>%
  dplyr::select(
    Probe               = Name,
    Gene                = UCSC_RefGene_Name,
    Region              = UCSC_RefGene_Group,
    CpG_Island_Relation = Relation_to_Island,
    Chromosome          = chr,
    Coordinate          = pos
  )

# Print and write
print(probe_annotation)
write.csv(probe_annotation, "Three_Probes_Annotation.csv", row.names = FALSE)


##
library(dplyr)
library(ggplot2)
library(ggpubr)

# 1. Use your significant probe IDs
sig_ids <- sig_probes$Probe  # cg01407062, cg05441039, cg22704057

# 2. Filter your melted_long (or melted_body) for just those probes
#    and bring in LOY_group from male_combined
plot_df <- melted_body %>%
  filter(Probe %in% sig_ids) %>%
  mutate(
    LOY_group = male_combined$LOY_group_cutoff[match(Sample, male_combined$sample_dash)]
  ) %>%
  filter(!is.na(LOY_group)) %>%
  mutate(
    Probe = factor(Probe, levels = sig_ids),                     # keep original order
    LOY_group = factor(LOY_group, levels = c("Y_low","Y_high"))  # ensure consistent ordering
  )

# 3. Make one boxplot per probe, comparing Y_low vs Y_high
ggboxplot(
  plot_df,
  x = "LOY_group",
  y = "Beta",
  fill = "LOY_group",
  facet.by = "Probe",
  scales = "free_y",
  palette = c("Y_low"  = "#1f77b4", 
              "Y_high" = "#ff7f0e")
) +
  stat_compare_means(
    aes(group = LOY_group),
    method = "wilcox.test", 
    label = "p.signif"
  ) +
  labs(
    title = "Methylation at Significant CpGs by LOY Status",
    x = "LOY Group",
    y = "Beta Value"
  ) +
  theme_minimal()

