# 15_pathway_effect_heatmap_partial.R
#
# Regenerates the cross-analysis pathway-effect heatmap (formerly
# 11_pathway_effect_heatmap.R) using the EXACT partial-tAge decomposition
# (14a/14b/14d) instead of the pathway-restricted-rerun approach -- same
# unbiased selection principle (recurrence-filtered, hierarchically
# clustered), now on numbers verified exact to floating-point precision.

source("R/config.R")
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(pheatmap)
})

df <- read.csv(file.path(RERUN_DIR, "partial_tage_ALL.csv"))
df <- df[df$model == "yugene" & df$analysis %in% c("meta_analysis", "temporal_pooled"), ]
df$label <- factor(df$label, levels = c("CICQ", "SSCQ", "RS", "SIPS", "OIS",
                                         "Fibroblast", "Keratinocyte", "Melanocyte"))
df$pathway_short <- gsub("^HALLMARK_", "", df$pathway)

mat_d <- df %>% select(pathway_short, label, cohens_d) %>%
  pivot_wider(names_from = label, values_from = cohens_d) %>% as.data.frame()
rownames(mat_d) <- mat_d$pathway_short; mat_d$pathway_short <- NULL
mat_d <- as.matrix(mat_d)

mat_p <- df %>% select(pathway_short, label, p_adj) %>%
  pivot_wider(names_from = label, values_from = p_adj) %>% as.data.frame()
rownames(mat_p) <- mat_p$pathway_short; mat_p$pathway_short <- NULL
mat_p <- as.matrix(mat_p)[rownames(mat_d), colnames(mat_d)]

n_sig <- rowSums(mat_p < 0.05, na.rm = TRUE)
cat("Distribution of n groups significant (of 8), all 50 pathways:\n"); print(table(n_sig))
keep <- names(n_sig)[n_sig >= 6]
cat(sprintf("Kept %d / %d Hallmark pathways (significant in >=6 of 8 groups)\n", length(keep), nrow(mat_d)))

mat_d_keep <- mat_d[keep, , drop = FALSE]
mat_p_keep <- mat_p[keep, , drop = FALSE]

sig_stars <- matrix("", nrow = nrow(mat_p_keep), ncol = ncol(mat_p_keep), dimnames = dimnames(mat_p_keep))
sig_stars[mat_p_keep < 0.001] <- "***"
sig_stars[mat_p_keep >= 0.001 & mat_p_keep < 0.01] <- "**"
sig_stars[mat_p_keep >= 0.01 & mat_p_keep < 0.05] <- "*"

col_annotation <- data.frame(
  Analysis = c(rep("Meta-analysis (vs Proliferating)", 5), rep("Temporal (vs own baseline)", 3))
)
rownames(col_annotation) <- colnames(mat_d_keep)

cap <- min(10, ceiling(quantile(abs(mat_d_keep), 0.97, na.rm = TRUE)))
breaks <- seq(-cap, cap, length.out = 101)

png(file.path(RERUN_DIR, "pathway_effect_heatmap.png"), width = 10, height = 12, units = "in", res = 300)
pheatmap(
  mat_d_keep,
  color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
  breaks = breaks,
  display_numbers = sig_stars,
  number_color = "black",
  fontsize_number = 10,
  cluster_cols = FALSE,
  cluster_rows = TRUE,
  clustering_method = "average",
  gaps_col = 5,
  annotation_col = col_annotation,
  na_col = "grey85",
  main = "Hallmark pathways moving tAge -- EXACT partial decomposition (Cohen's d, yugene_diff model)\nRows: significant (padj<0.05) in >=6/8 groups. */**/*** = padj<.05/.01/.001",
  fontsize_row = 9,
  fontsize_col = 10,
  angle_col = 45,
  border_color = "grey60"
)
dev.off()
cat(sprintf("Saved -> %s\n", file.path(RERUN_DIR, "pathway_effect_heatmap.png")))

cat("\n== Cross-analysis correlation (Cohen's d, yugene model, all 50 pathways) ==\n")
cor_mat <- cor(mat_d, use = "pairwise.complete.obs", method = "spearman")
print(round(cor_mat[c("Fibroblast", "Keratinocyte", "Melanocyte"), c("CICQ", "SSCQ", "RS", "SIPS", "OIS")], 2))
