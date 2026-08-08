# 11_pathway_effect_heatmap.R
#
# One clean, data-driven (not narrative-driven) figure showing which MSigDB
# Hallmark pathways actually move tAge, across BOTH analyses side by side:
#   - meta_analysis: CICQ / SSCQ / RS / SIPS / OIS vs pooled Proliferating
#     (recount3 cross-study fibroblast comparison)
#   - temporal: Fibroblast / Keratinocyte / Melanocyte, pooled irradiated
#     (4+10+20 days) vs each cell type's own 'none' baseline
#
# UNBIASED PATHWAY SELECTION (not the DNA-damage/mTOR/hypoxia/interferon
# subset picked earlier for narrative purposes): a pathway is kept if it's
# BH-significant (p_adj < 0.05, yugene_diff model) in at least 2 of the 8
# groups -- a purely data-driven recurrence filter, not a topic-based one.
# Rows are then hierarchically clustered so the figure's grouping comes from
# the data, not from a preconceived stress-pathway narrative.

source("R/config.R")
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(pheatmap)
})

df <- read.csv(file.path(RERUN_DIR, "pathway_tage_score_all.csv"))
df <- df[df$status == "OK", ]

# Consistent, readable group labels/order per analysis
df$group <- factor(df$group, levels = c("CICQ", "SSCQ", "RS", "SIPS", "OIS",
                                         "Fibroblast", "Keratinocyte", "Melanocyte"))
df$pathway_short <- gsub("^HALLMARK_", "", df$pathway)

# ── Effect-size matrix (Cohen's d, yugene_diff model) ───────────────────────
mat_d <- df %>%
  select(pathway_short, group, cohens_d_yugene) %>%
  pivot_wider(names_from = group, values_from = cohens_d_yugene) %>%
  as.data.frame()
rownames(mat_d) <- mat_d$pathway_short
mat_d$pathway_short <- NULL
mat_d <- as.matrix(mat_d)

# ── Significance matrix, for the recurrence filter and star overlay ────────
mat_p <- df %>%
  select(pathway_short, group, p_adj_yugene) %>%
  pivot_wider(names_from = group, values_from = p_adj_yugene) %>%
  as.data.frame()
rownames(mat_p) <- mat_p$pathway_short
mat_p$pathway_short <- NULL
mat_p <- as.matrix(mat_p)
mat_p <- mat_p[rownames(mat_d), colnames(mat_d)]

# Data-driven recurrence filter: significant in >=6 of the 8 groups (broad
# consensus, not a single lucky test) -- chosen by checking the actual
# distribution of n_sig across all 50 pathways first (0-7 sig groups out of
# 8; >=2 kept 47/50, useless for a figure; >=6 kept a readable 16/50 while
# still being a purely count-based, topic-blind cutoff).
n_sig <- rowSums(mat_p < 0.05, na.rm = TRUE)
keep <- names(n_sig)[n_sig >= 6]
cat(sprintf("Kept %d / %d Hallmark pathways (significant in >=6 of 8 groups)\n", length(keep), nrow(mat_d)))

mat_d_keep <- mat_d[keep, , drop = FALSE]
mat_p_keep <- mat_p[keep, , drop = FALSE]

sig_stars <- matrix("", nrow = nrow(mat_p_keep), ncol = ncol(mat_p_keep),
                     dimnames = dimnames(mat_p_keep))
sig_stars[mat_p_keep < 0.001] <- "***"
sig_stars[mat_p_keep >= 0.001 & mat_p_keep < 0.01] <- "**"
sig_stars[mat_p_keep >= 0.01 & mat_p_keep < 0.05] <- "*"

col_annotation <- data.frame(
  Analysis = c(rep("Meta-analysis (vs Proliferating)", 5), rep("Temporal (vs own baseline)", 3))
)
rownames(col_annotation) <- colnames(mat_d_keep)

# Cap color scale symmetrically, robust to a couple of extreme outlier cells
cap <- min(8, ceiling(quantile(abs(mat_d_keep), 0.97, na.rm = TRUE)))
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
  main = "Hallmark pathways moving tAge (Cohen's d, yugene_diff model)\nRows: significant (padj<0.05) in >=6/8 groups. */**/*** = padj<.05/.01/.001",
  fontsize_row = 9,
  fontsize_col = 10,
  angle_col = 45,
  border_color = "grey60",
  legend = TRUE,
  annotation_legend = TRUE,
  annotation_names_col = FALSE
)
dev.off()
cat(sprintf("Saved -> %s\n", file.path(RERUN_DIR, "pathway_effect_heatmap.png")))

# ── Quantify: does temporal Fibroblast resemble meta-analysis RS/SIPS/OIS? ──
# Both are fibroblast senescence signals -- a natural cross-analysis check.
cat("\n== Cross-analysis correlation (Cohen's d, yugene_diff, all 50 pathways -- not just the filtered subset) ==\n")
cor_mat <- cor(mat_d, use = "pairwise.complete.obs", method = "spearman")
print(round(cor_mat[c("Fibroblast", "Keratinocyte", "Melanocyte"), c("CICQ", "SSCQ", "RS", "SIPS", "OIS")], 2))
