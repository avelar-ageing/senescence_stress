# 14c_analyze_partial_scores.R
#
# Analyzes the exact linear "partial tAge" decomposition (14a+14b) --
# computes Cohen's d / Wilcoxon per pathway per condition/cell-type from the
# EXACT per-sample partial contributions (no re-imputation, no re-
# normalization on a restricted gene set; verified exact reconstruction of
# the real tAge prediction to floating-point precision). Then compares
# against the earlier (flawed) pathway-restricted-rerun approach to show how
# much that approach's forced ~99% imputation distorted the numbers.

source("R/config.R")

PT <- file.path(RERUN_DIR, "partial_tage")
groups_meta <- read.csv(file.path(PT, "meta_groups.csv"))  # group values already match cell_substate spelling exactly

cohens_d <- function(x, y) {
  nx <- length(x); ny <- length(y)
  pooled_sd <- sqrt(((nx - 1) * var(x) + (ny - 1) * var(y)) / (nx + ny - 2))
  (mean(x) - mean(y)) / pooled_sd
}

score_pathways <- function(scores_csv, groups_df, test_label, control_label, label) {
  scores <- read.csv(scores_csv)
  m <- merge(scores, groups_df, by = "sample_id")
  test_rows <- m[m$group == test_label, ]
  ctrl_rows <- m[m$group == control_label, ]
  pathway_cols <- setdiff(colnames(scores), c("sample_id", "full_tAge_direct", "full_tAge_reconstructed"))

  do.call(rbind, lapply(pathway_cols, function(pw) {
    d <- cohens_d(test_rows[[pw]], ctrl_rows[[pw]])
    p <- wilcox.test(test_rows[[pw]], ctrl_rows[[pw]])$p.value
    data.frame(label = label, pathway = pw, cohens_d = d, wilcox_p = p,
               n_test = nrow(test_rows), n_control = nrow(ctrl_rows))
  }))
}

all_results <- list()

meta_conditions <- c("Contact_inhibited CQ" = "CICQ", "Serum_starved CQ" = "SSCQ",
                     "Replicative CS" = "RS", "Stress-induced CS" = "SIPS",
                     "Oncogene-induced CS" = "OIS")
for (model in c("scaled", "yugene")) {
  scores_csv <- file.path(PT, sprintf("meta_partial_scores_%s.csv", model))
  for (cond_value in names(meta_conditions)) {
    r <- score_pathways(scores_csv, groups_meta, cond_value, "Proliferating", meta_conditions[[cond_value]])
    r$model <- model
    all_results[[paste("meta", meta_conditions[[cond_value]], model)]] <- r
  }
}

for (ct in c("Fibroblast", "Keratinocyte", "Melanocyte")) {
  groups_ct <- read.csv(file.path(PT, sprintf("%s_groups.csv", ct)))
  for (model in c("scaled", "yugene")) {
    scores_csv <- file.path(PT, sprintf("%s_partial_scores_%s.csv", ct, model))
    r <- score_pathways(scores_csv, groups_ct, "irradiated", "none", ct)
    r$model <- model
    all_results[[paste("temporal", ct, model)]] <- r
  }
}

final <- do.call(rbind, all_results)
final$p_adj <- ave(final$wilcox_p, final$model, FUN = function(x) p.adjust(x, method = "BH"))
write.csv(final, file.path(RERUN_DIR, "partial_tage_scores_final.csv"), row.names = FALSE)

cat("== Top 10 partial-tAge findings overall (by |Cohen's d|, yugene model) ==\n")
yug <- final[final$model == "yugene", ]
yug <- yug[order(-abs(yug$cohens_d)), ]
print(head(yug[, c("label", "pathway", "cohens_d", "p_adj", "n_test", "n_control")], 10), row.names = FALSE)

cat("\n== Compare to the earlier PATHWAY-RESTRICTED-RERUN approach (08/09 scripts) for the same headline cells ==\n")
old <- read.csv(file.path(RERUN_DIR, "pathway_tage_score_all.csv"))
old <- old[old$status == "OK", ]
compare_targets <- list(
  c("Melanocyte", "HALLMARK_INTERFERON_ALPHA_RESPONSE"),
  c("Melanocyte", "HALLMARK_DNA_REPAIR"),
  c("OIS", "HALLMARK_NOTCH_SIGNALING"),
  c("RS", "HALLMARK_MITOTIC_SPINDLE"),
  c("SIPS", "HALLMARK_PI3K_AKT_MTOR_SIGNALING")
)
for (tgt in compare_targets) {
  lbl <- tgt[1]; pw <- tgt[2]
  old_row <- old[old$group == lbl & old$pathway == pw, ]
  new_row <- final[final$label == lbl & final$pathway == pw & final$model == "yugene", ]
  cat(sprintf("%-12s %-40s  OLD(pathway-restricted-rerun) d=%6.2f padj=%.1e  |  NEW(exact partial) d=%6.2f padj=%.1e\n",
              lbl, pw,
              if (nrow(old_row)) old_row$cohens_d_yugene[1] else NA, if (nrow(old_row)) old_row$p_adj_yugene[1] else NA,
              if (nrow(new_row)) new_row$cohens_d[1] else NA, if (nrow(new_row)) new_row$p_adj[1] else NA))
}
