# 07_pathway_divergence_meta_conditions.R
#
# Which Hallmark pathways differ MOST between the 5 meta-analysis conditions
# themselves (CICQ/SSCQ/RS/SIPS/OIS vs each other, not vs Proliferating in
# general), from the exact partial-tAge decomposition (01-05). Divergence =
# range of Cohen's d across the 5 conditions (max-min), a topic-blind ranking.

source("R/config.R")
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(pheatmap)
})

df <- read.csv(file.path(RERUN_DIR, "partial_tage_ALL.csv"))
meta <- df[df$analysis == "meta_analysis" & df$model == "yugene", ]
meta$label <- factor(meta$label, levels = c("CICQ", "SSCQ", "RS", "SIPS", "OIS"))
meta$pathway_short <- gsub("^HALLMARK ", "", meta$pathway)

mat_d <- meta %>% select(pathway_short, label, cohens_d) %>%
  pivot_wider(names_from = label, values_from = cohens_d) %>% as.data.frame()
rownames(mat_d) <- mat_d$pathway_short; mat_d$pathway_short <- NULL
mat_d <- as.matrix(mat_d)

mat_p <- meta %>% select(pathway_short, label, p_adj) %>%
  pivot_wider(names_from = label, values_from = p_adj) %>% as.data.frame()
rownames(mat_p) <- mat_p$pathway_short; mat_p$pathway_short <- NULL
mat_p <- as.matrix(mat_p)[rownames(mat_d), colnames(mat_d)]

divergence <- apply(mat_d, 1, function(x) max(x) - min(x))
keep <- names(sort(divergence, decreasing = TRUE))[1:20]
cat("Top 20 pathways by divergence (max-min Cohen's d) across CICQ/SSCQ/RS/SIPS/OIS:\n")
print(round(divergence[keep], 2))

mat_d_keep <- mat_d[keep, , drop = FALSE]
mat_p_keep <- mat_p[keep, , drop = FALSE]

sig_stars <- matrix("", nrow = nrow(mat_p_keep), ncol = ncol(mat_p_keep), dimnames = dimnames(mat_p_keep))
sig_stars[mat_p_keep < 0.001] <- "***"
sig_stars[mat_p_keep >= 0.001 & mat_p_keep < 0.01] <- "**"
sig_stars[mat_p_keep >= 0.01 & mat_p_keep < 0.05] <- "*"

cap <- min(8, ceiling(quantile(abs(mat_d_keep), 0.97, na.rm = TRUE)))
breaks <- seq(-cap, cap, length.out = 101)

png(file.path(RERUN_DIR, "pathway_divergence_meta_conditions.png"), width = 8, height = 10, units = "in", res = 300)
pheatmap(
  mat_d_keep,
  color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
  breaks = breaks,
  display_numbers = sig_stars,
  number_color = "black",
  fontsize_number = 11,
  cluster_cols = FALSE,
  cluster_rows = TRUE,
  clustering_method = "average",
  na_col = "grey85",
  main = "Top 20 pathways separating arrest conditions from each other\n(meta-analysis, partial-tAge decomposition, Cohen's d yugene_diff). */**/*** = padj<.05/.01/.001",
  fontsize_row = 10,
  fontsize_col = 12,
  angle_col = 0,
  border_color = "grey60"
)
dev.off()
cat(sprintf("\nSaved -> %s\n", file.path(RERUN_DIR, "pathway_divergence_meta_conditions.png")))
