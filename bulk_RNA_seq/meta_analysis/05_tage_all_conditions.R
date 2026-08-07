# 05_tage_all_conditions.R
#
# tAge (transcriptomic age, Gladyshev-Lab/tAge EN elastic-net clock) for every
# arrest/senescence condition in the regenerated recount3 meta-analysis
# object: Proliferating, Contact-inhibited CQ, Serum-starved CQ, Replicative
# CS (RS), Oncogene-induced CS (OIS), Stress-induced CS (SIPS). All samples
# are Fibroblast by construction (see 01_build_cs_cq_object.R).
#
# Same pipeline as ../../scrna_seq/GSE226225/final_analysis/tage_analysis.R
# and the earlier ad hoc bulk run (species='human', EN models, control group
# = Proliferating), just applied to the full 6-condition object instead of
# the Proliferating+CQ-only cq_samples.rds used before CS data was available.

source("R/config.R")

suppressPackageStartupMessages({
  library(SummarizedExperiment)
  library(tAge)
  library(Biobase)
  library(dplyr)
  library(ggplot2)
})

SCRNA_DIR <- "/home/ro/APFS_copy/root/Backup/Documents/modules/scrna_seq/GSE226225/final_analysis"
MODEL_DIR <- file.path(SCRNA_DIR, "tAge_models")
Sys.setenv(RETICULATE_PYTHON = file.path(SCRNA_DIR, ".venv/bin/python"))
model_paths <- list(
  scaled_diff = file.path(MODEL_DIR, "EN_Chronoage_Multispecies_Multitissue_scaleddiff.pkl"),
  yugene_diff = file.path(MODEL_DIR, "EN_Chronoage_Multispecies_Multitissue_yugenediff.pkl")
)

rse <- readRDS(file.path(RERUN_DIR, "cs_cq_all_study_processed.rds"))
cat(sprintf("Loaded regenerated object: %d genes x %d samples\n", nrow(rse), ncol(rse)))
print(table(colData(rse)$cell_substate))

eset <- ExpressionSet(
  assayData = as.matrix(assay(rse)),
  phenoData = AnnotatedDataFrame(as.data.frame(colData(rse)))
)

cat("\n[meta-analysis, all 6 conditions] preprocessing...\n")
tAge_eset <- tAge_preprocessing(
  eset, species = "human", gene_mapping_type = "Gene.Symbol",
  control_group_column = "cell_state", control_group_label = "Proliferating",
  verbose = TRUE, count_threshold = 10, percent_threshold = 20
)

cat("\n[meta-analysis, all 6 conditions] predicting tAge...\n")
tage_result <- predict_tAge(tAge_eset, model_paths, species = "human", mode = "EN")

# Tidy condition labels to match the paper's short names
tage_result$condition <- tage_result$cell_substate
tage_result$condition[tage_result$condition == "Replicative CS"] <- "RS"
tage_result$condition[tage_result$condition == "Oncogene-induced CS"] <- "OIS"
tage_result$condition[tage_result$condition == "Stress-induced CS"] <- "SIPS"
tage_result$condition[tage_result$condition == "Contact_inhibited CQ"] <- "CICQ"
tage_result$condition[tage_result$condition == "Serum_starved CQ"] <- "SSCQ"
tage_result$condition <- factor(tage_result$condition,
                                 levels = c("Proliferating", "CICQ", "SSCQ", "RS", "SIPS", "OIS"))

write.csv(tage_result, file.path(RERUN_DIR, "tage_all_conditions.csv"), row.names = FALSE)

cat("\n== tAge spread by condition (both EN models) ==\n")
spread_stats <- tage_result %>%
  group_by(condition) %>%
  summarise(
    n = n(),
    scaled_diff_median = median(scaled_diff_EN_tAge), scaled_diff_mean = mean(scaled_diff_EN_tAge),
    scaled_diff_sd = sd(scaled_diff_EN_tAge), scaled_diff_min = min(scaled_diff_EN_tAge), scaled_diff_max = max(scaled_diff_EN_tAge),
    yugene_diff_median = median(yugene_diff_EN_tAge), yugene_diff_mean = mean(yugene_diff_EN_tAge),
    yugene_diff_sd = sd(yugene_diff_EN_tAge), yugene_diff_min = min(yugene_diff_EN_tAge), yugene_diff_max = max(yugene_diff_EN_tAge)
  )
print(spread_stats)
write.csv(spread_stats, file.path(RERUN_DIR, "tage_spread_by_condition.csv"), row.names = FALSE)

# Wilcoxon vs Proliferating, both models, BH-corrected across the 10 tests (5 conditions x 2 models)
cat("\n== Wilcoxon rank-sum vs Proliferating (BH-corrected across 5 conditions x 2 models) ==\n")
prolif <- tage_result[tage_result$condition == "Proliferating", ]
other_conditions <- setdiff(levels(tage_result$condition), "Proliferating")
wilcox_rows <- do.call(rbind, lapply(other_conditions, function(cond) {
  grp <- tage_result[tage_result$condition == cond, ]
  p_scaled <- wilcox.test(grp$scaled_diff_EN_tAge, prolif$scaled_diff_EN_tAge)$p.value
  p_yugene <- wilcox.test(grp$yugene_diff_EN_tAge, prolif$yugene_diff_EN_tAge)$p.value
  data.frame(condition = cond, model = c("scaled_diff", "yugene_diff"), p = c(p_scaled, p_yugene))
}))
wilcox_rows$p.adj <- p.adjust(wilcox_rows$p, method = "BH")
print(wilcox_rows)
write.csv(wilcox_rows, file.path(RERUN_DIR, "tage_wilcoxon_vs_proliferating.csv"), row.names = FALSE)

# All-pairs comparison (every condition vs every other condition, not just vs
# Proliferating), both models, BH-corrected across all 15 pairs x 2 models = 30 tests.
cat("\n== Wilcoxon rank-sum, ALL PAIRS of conditions (BH-corrected across 15 pairs x 2 models) ==\n")
all_conditions <- levels(tage_result$condition)
pairs <- combn(all_conditions, 2, simplify = FALSE)
pairwise_rows <- do.call(rbind, lapply(pairs, function(pr) {
  g1 <- tage_result[tage_result$condition == pr[1], ]
  g2 <- tage_result[tage_result$condition == pr[2], ]
  p_scaled <- wilcox.test(g1$scaled_diff_EN_tAge, g2$scaled_diff_EN_tAge)$p.value
  p_yugene <- wilcox.test(g1$yugene_diff_EN_tAge, g2$yugene_diff_EN_tAge)$p.value
  diff_scaled <- median(g1$scaled_diff_EN_tAge) - median(g2$scaled_diff_EN_tAge)
  diff_yugene <- median(g1$yugene_diff_EN_tAge) - median(g2$yugene_diff_EN_tAge)
  data.frame(
    condition_1 = pr[1], condition_2 = pr[2],
    model = c("scaled_diff", "yugene_diff"),
    median_diff_1_minus_2 = c(diff_scaled, diff_yugene),
    p = c(p_scaled, p_yugene)
  )
}))
pairwise_rows$p.adj <- p.adjust(pairwise_rows$p, method = "BH")
pairwise_rows <- pairwise_rows[order(pairwise_rows$p.adj), ]
print(pairwise_rows)
write.csv(pairwise_rows, file.path(RERUN_DIR, "tage_pairwise_all_conditions.csv"), row.names = FALSE)

# Heatmap of pairwise median tAge differences (scaled_diff model), significance-masked
pairwise_wide <- pairwise_rows[pairwise_rows$model == "scaled_diff", ]
sig_symbol <- function(p) ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", "")))
pairwise_wide$sig <- sig_symbol(pairwise_wide$p.adj)
p_heat <- ggplot(pairwise_wide, aes(x = condition_2, y = condition_1, fill = median_diff_1_minus_2)) +
  geom_tile(colour = "black") +
  geom_text(aes(label = paste0(round(median_diff_1_minus_2, 1), sig)), size = 3.5) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red2", midpoint = 0,
                        name = "Median tAge\ndiff (1 - 2)") +
  theme_bw() + labs(x = NULL, y = NULL,
                     title = "Pairwise tAge differences between conditions (scaled_diff model)",
                     subtitle = "*/**/*** = BH-adjusted p < 0.05/0.01/0.001")
ggsave(file.path(RERUN_DIR, "tage_pairwise_heatmap.png"), p_heat, width = 7, height = 6, dpi = 300)

# Violin/boxplot of tAge spread by condition, both models
plot_df <- tage_result %>%
  dplyr::select(condition, scaled_diff_EN_tAge, yugene_diff_EN_tAge) %>%
  tidyr::pivot_longer(cols = c(scaled_diff_EN_tAge, yugene_diff_EN_tAge), names_to = "model", values_to = "tAge")
p <- ggplot(plot_df, aes(x = condition, y = tAge, fill = condition)) +
  geom_violin(alpha = 0.6, trim = FALSE) +
  geom_boxplot(width = 0.15, outlier.shape = NA, fill = "white") +
  facet_wrap(~model, ncol = 1, scales = "free_y") +
  theme_bw() + theme(legend.position = "none") +
  labs(x = NULL, y = "tAge", title = "Transcriptomic age by arrest/senescence condition (Fibroblast, recount3 meta-analysis)")
ggsave(file.path(RERUN_DIR, "tage_spread_by_condition.png"), p, width = 8, height = 8, dpi = 300)

cat(sprintf("\nDone. Outputs -> %s\n", RERUN_DIR))
