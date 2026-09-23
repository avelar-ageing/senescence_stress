# 07_tage_temporal_tests.R
#
# ALL temporal tAge significance testing, for ALL THREE clocks, in one place and
# under ONE BH family. Supersedes 07_tage_temporal_pairwise.R (chronological
# only, 36-test family) and the testing half of 08_mortality_temporal.py
# (mortality only, 18-test family), which produced two files and two different
# adjusted floors for the same 6v6 design - FDR 0.004 against 0.003 - purely
# because the families differed in size. Pooling gives one floor of 0.0034.
#
# 08_mortality_temporal.py still predicts the mortality values (it needs the
# model .pkl); it no longer tests them.
#
# Output: rerun_outputs/tage_temporal_tests.csv, with `model` and `test`
# columns, replacing tage_temporal_pairwise_all_timepoints.csv,
# tage_temporal_kruskal.csv and the pairwise/kruskal rows of
# mortality_temporal_pairwise.csv.
source("R/config.R")
source(file.path("temporal_analysis", "R_keratinocyte_batch.R"))

TIME_LEVELS <- c("none", "4_days", "10_days", "20_days")
CT_LEVELS   <- c("Fibroblast", "Keratinocyte", "Melanocyte")

# --- chronological clocks: batch-centre keratinocytes, as the DEG arm does ----
chron <- read.csv(file.path(RERUN_DIR, "tage_temporal_by_celltype.csv"))
chron <- batch_centre_keratinocytes(chron, c("scaled_diff_EN_tAge", "yugene_diff_EN_tAge"))
# --- mortality: already centred by 08_mortality_temporal.py -------------------
mort <- read.csv(file.path(RERUN_DIR, "mortality_temporal_by_celltype.csv"))
mort$time_after_treatment <- mort$timepoint

VALUES <- list(yugene_diff = list(d = chron, col = "yugene_diff_EN_tAge"),
               scaled_diff = list(d = chron, col = "scaled_diff_EN_tAge"),
               mortality   = list(d = mort,  col = "mortality_tAge"))

pairs <- combn(TIME_LEVELS, 2, simplify = FALSE)
floor_p <- function(a, b) 2 / choose(a + b, min(a, b))

rows <- do.call(rbind, lapply(names(VALUES), function(mdl) {
  d <- VALUES[[mdl]]$d; v <- VALUES[[mdl]]$col
  do.call(rbind, lapply(CT_LEVELS, function(ct) {
    dc <- d[d$cell_type == ct, ]
    do.call(rbind, lapply(pairs, function(pr) {
      x <- dc[[v]][dc$time_after_treatment == pr[1]]
      y <- dc[[v]][dc$time_after_treatment == pr[2]]
      data.frame(test = "pairwise_timepoints", model = mdl, cell_type = ct,
                 group_1 = pr[1], group_2 = pr[2], n_1 = length(x), n_2 = length(y),
                 median_1 = median(x), median_2 = median(y),
                 diff = median(y) - median(x),
                 p = suppressWarnings(wilcox.test(x, y)$p.value),
                 p_floor = floor_p(length(x), length(y)))
    }))
  }))
}))
# ONE BH family across all three clocks: 3 models x 3 cell types x 6 pairs = 54.
rows$p_adj <- p.adjust(rows$p, method = "BH")

omni <- do.call(rbind, lapply(names(VALUES), function(mdl) {
  d <- VALUES[[mdl]]$d; v <- VALUES[[mdl]]$col
  do.call(rbind, lapply(CT_LEVELS, function(ct) {
    dc <- d[d$cell_type == ct, ]
    k <- kruskal.test(dc[[v]] ~ factor(dc$time_after_treatment, levels = TIME_LEVELS))
    data.frame(test = "kruskal", model = mdl, cell_type = ct, group_1 = NA, group_2 = NA,
               n_1 = nrow(dc), n_2 = NA, median_1 = NA, median_2 = NA, diff = NA,
               p = k$p.value, p_floor = NA)
  }))
}))
omni$p_adj <- p.adjust(omni$p, method = "BH")   # its own family: 9 omnibus tests

out <- rbind(rows, omni)
write.csv(out, file.path(RERUN_DIR, "tage_temporal_tests.csv"), row.names = FALSE)
cat(sprintf("Saved -> tage_temporal_tests.csv (%d pairwise in one BH family + %d kruskal)\n",
            nrow(rows), nrow(omni)))
cat(sprintf("  pooled floor: raw %.5f -> BH %.5f (%d of %d tests tied at it)\n",
            min(rows$p), min(rows$p_adj), sum(rows$p_adj == min(rows$p_adj)), nrow(rows)))
cat(sprintf("  significant at p_adj < 0.05: %d of %d\n", sum(rows$p_adj < 0.05), nrow(rows)))
for (ct in CT_LEVELS) {
  cat(sprintf("\n-- %s --\n", ct))
  s <- rows[rows$cell_type == ct, c("model","group_1","group_2","diff","p_adj")]
  print(s[order(s$model, match(s$group_1,TIME_LEVELS), match(s$group_2,TIME_LEVELS)), ],
        row.names = FALSE, digits = 3)
}
