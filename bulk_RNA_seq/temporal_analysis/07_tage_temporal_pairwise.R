# 07_tage_temporal_pairwise.R
#
# All-pairs timepoint comparisons (not just each timepoint vs 'none'), per
# cell type, both EN models -- the temporal counterpart to
# meta_analysis/05_tage_all_conditions.R's tage_pairwise_all_conditions.csv.
# 04_tage_by_celltype.R only tested each timepoint against its own 'none'
# baseline; this adds 4_days vs 10_days, 4_days vs 20_days, and 10_days vs
# 20_days, so trajectory claims (e.g. a partial reversal by 20 days) rest on
# a direct test rather than being inferred from separate vs-baseline tests.
# Reads the already-computed tage_temporal_by_celltype.csv -- no re-run of
# tAge prediction.

source("R/config.R")

TIME_LEVELS <- c("none", "4_days", "10_days", "20_days")
CT_LEVELS <- c("Fibroblast", "Keratinocyte", "Melanocyte")

tage_temporal <- read.csv(file.path(RERUN_DIR, "tage_temporal_by_celltype.csv"))
tage_temporal$time_after_treatment <- factor(tage_temporal$time_after_treatment, levels = TIME_LEVELS)

pairs <- combn(TIME_LEVELS, 2, simplify = FALSE)

rows <- do.call(rbind, lapply(CT_LEVELS, function(ct) {
  df_ct <- tage_temporal[tage_temporal$cell_type == ct, ]
  do.call(rbind, lapply(pairs, function(pr) {
    g1 <- df_ct[df_ct$time_after_treatment == pr[1], ]
    g2 <- df_ct[df_ct$time_after_treatment == pr[2], ]
    p_scaled <- wilcox.test(g1$scaled_diff_EN_tAge, g2$scaled_diff_EN_tAge)$p.value
    p_yugene <- wilcox.test(g1$yugene_diff_EN_tAge, g2$yugene_diff_EN_tAge)$p.value
    diff_scaled <- median(g1$scaled_diff_EN_tAge) - median(g2$scaled_diff_EN_tAge)
    diff_yugene <- median(g1$yugene_diff_EN_tAge) - median(g2$yugene_diff_EN_tAge)
    data.frame(cell_type = ct, timepoint_1 = pr[1], timepoint_2 = pr[2],
               model = c("scaled_diff", "yugene_diff"),
               median_diff_1_minus_2 = c(diff_scaled, diff_yugene),
               p = c(p_scaled, p_yugene))
  }))
}))
# BH across the full family: 3 cell types x 6 pairs x 2 models = 36 tests.
rows$p.adj <- p.adjust(rows$p, method = "BH")
rows <- rows[order(rows$p.adj), ]
write.csv(rows, file.path(RERUN_DIR, "tage_temporal_pairwise_all_timepoints.csv"), row.names = FALSE)

cat(sprintf("Saved -> %s (%d rows, BH-corrected across %d tests)\n",
            file.path(RERUN_DIR, "tage_temporal_pairwise_all_timepoints.csv"), nrow(rows), nrow(rows)))
cat(sprintf("  %d of %d comparisons significant at padj<0.05\n", sum(rows$p.adj < 0.05), nrow(rows)))
cat("\n== All pairwise timepoint comparisons NOT involving 'none' (i.e. 4d/10d/20d vs each other) ==\n")
print(rows[rows$timepoint_1 != "none" & rows$timepoint_2 != "none",
           c("cell_type", "timepoint_1", "timepoint_2", "model", "median_diff_1_minus_2", "p.adj")],
      row.names = FALSE)

# Global Kruskal-Wallis across all 4 timepoints, one test per cell type x
# model (uncorrected; a single omnibus test per panel, not part of the
# pairwise family above). Carried over from the superseded
# 05_tage_spread_figure.R so this omnibus result isn't lost with it.
cat("\n== Kruskal-Wallis across all 4 timepoints, per cell type x model ==\n")
kw <- do.call(rbind, lapply(CT_LEVELS, function(ct) {
  d <- tage_temporal[tage_temporal$cell_type == ct, ]
  data.frame(cell_type = ct, model = c("scaled_diff", "yugene_diff"),
             kruskal_p = c(kruskal.test(d$scaled_diff_EN_tAge ~ d$time_after_treatment)$p.value,
                           kruskal.test(d$yugene_diff_EN_tAge ~ d$time_after_treatment)$p.value))
}))
print(kw, row.names = FALSE)
write.csv(kw, file.path(RERUN_DIR, "tage_temporal_kruskal.csv"), row.names = FALSE)
cat(sprintf("Saved -> %s\n", file.path(RERUN_DIR, "tage_temporal_kruskal.csv")))
