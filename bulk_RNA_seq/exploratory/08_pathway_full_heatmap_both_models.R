# 08_pathway_full_heatmap_both_models.R
#
# Figure for section 2.1.5.2 (Pathway-level Transcriptomic Age Differences):
# comprehensive partial-tAge heatmap, meta-analysis only, ALL 50 pathways
# (every pathway reaches significance in at least one of the 5 conditions x
# 2 models = 10 tests, so no pathway is dropped), both EN models shown
# side-by-side per condition so model-(dis)agreement is visible directly in
# the figure rather than only in partial_tage_ALL.csv. Complements
# 06_pathway_effect_heatmap.R (yugene-only, recurrence-filtered cross-
# analysis view) and 07_pathway_divergence_meta_conditions.R (yugene-only,
# meta-analysis divergence ranking) rather than replacing them.

source("R/config.R")
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(pheatmap)
})

df <- read.csv(file.path(RERUN_DIR, "partial_tage_ALL.csv"), check.names = FALSE)
meta <- df[df$analysis == "meta_analysis", ]
meta$pathway_short <- gsub("^HALLMARK ", "", meta$pathway)
meta$col <- paste(meta$label, meta$model, sep = " | ")
col_order <- as.vector(t(outer(c("CICQ", "SSCQ", "RS", "SIPS", "OIS"), c("scaled", "yugene"),
                                FUN = function(a, b) paste(a, b, sep = " | "))))

mat_d <- meta %>% select(pathway_short, col, cohens_d) %>%
  pivot_wider(names_from = col, values_from = cohens_d) %>% as.data.frame()
rownames(mat_d) <- mat_d$pathway_short; mat_d$pathway_short <- NULL
mat_d <- as.matrix(mat_d)[, col_order]

mat_p <- meta %>% select(pathway_short, col, p_adj) %>%
  pivot_wider(names_from = col, values_from = p_adj) %>% as.data.frame()
rownames(mat_p) <- mat_p$pathway_short; mat_p$pathway_short <- NULL
mat_p <- as.matrix(mat_p)[rownames(mat_d), col_order]

sig_stars <- matrix("", nrow = nrow(mat_p), ncol = ncol(mat_p), dimnames = dimnames(mat_p))
sig_stars[mat_p < 0.001] <- "***"
sig_stars[mat_p >= 0.001 & mat_p < 0.01] <- "**"
sig_stars[mat_p >= 0.01 & mat_p < 0.05] <- "*"

col_annotation <- data.frame(
  Condition = rep(c("CICQ", "SSCQ", "RS", "SIPS", "OIS"), each = 2),
  Model = rep(c("scaled_diff", "yugene_diff"), times = 5)
)
rownames(col_annotation) <- col_order
display_col_labels <- rep(c("scaled", "yugene"), times = 5)

cap <- min(6, ceiling(quantile(abs(mat_d), 0.97, na.rm = TRUE)))
breaks <- seq(-cap, cap, length.out = 101)

png(file.path(RERUN_DIR, "figure_pathway_heatmap_all_both_models.png"),
    width = 11, height = 15, units = "in", res = 300)
pheatmap(
  mat_d,
  color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
  breaks = breaks,
  display_numbers = sig_stars,
  number_color = "black",
  fontsize_number = 8,
  cluster_cols = FALSE,
  cluster_rows = TRUE,
  clustering_method = "average",
  gaps_col = seq(2, 8, by = 2),
  annotation_col = col_annotation,
  labels_col = display_col_labels,
  na_col = "grey85",
  main = paste0("All 50 pathways, partial-tAge decomposition (Cohen's d), meta-analysis vs Proliferating\n",
                "Both EN models shown per condition (scaled/yugene). */**/*** = padj<.05/.01/.001, BH-adjusted within analysis x model"),
  fontsize_row = 7,
  fontsize_col = 9,
  angle_col = 0,
  border_color = "grey70"
)
dev.off()
cat(sprintf("Saved -> %s\n", file.path(RERUN_DIR, "figure_pathway_heatmap_all_both_models.png")))

n_sig <- rowSums(mat_p < 0.05, na.rm = TRUE)
cat("\nDistribution of n significant tests (of 10: 5 conditions x 2 models), all 50 pathways:\n")
print(table(n_sig))
