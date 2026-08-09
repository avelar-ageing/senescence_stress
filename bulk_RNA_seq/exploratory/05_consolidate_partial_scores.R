# 05_consolidate_partial_scores.R
#
# Consolidates the exact partial-tAge decomposition (01/02/03/04) across
# EVERY group comparison computed: 5 meta-analysis conditions, 3 pooled
# temporal cell types (irradiated vs none), and 9 per-timepoint temporal
# comparisons (3 cell types x 3 timepoints, 6v6 vs none).

source("R/config.R")
PT <- file.path(RERUN_DIR, "partial_tage")

cohens_d <- function(x, y) {
  nx <- length(x); ny <- length(y)
  pooled_sd <- sqrt(((nx - 1) * var(x) + (ny - 1) * var(y)) / (nx + ny - 2))
  (mean(x) - mean(y)) / pooled_sd
}

score_pathways <- function(scores_csv, groups_df, test_label, control_label, label, analysis, cell_type = NA, timepoint = NA) {
  scores <- read.csv(scores_csv)
  m <- merge(scores, groups_df, by = "sample_id")
  test_rows <- m[m$group == test_label, ]
  ctrl_rows <- m[m$group == control_label, ]
  pathway_cols <- setdiff(colnames(scores), c("sample_id", "full_tAge_direct", "full_tAge_reconstructed"))

  do.call(rbind, lapply(pathway_cols, function(pw) {
    d <- cohens_d(test_rows[[pw]], ctrl_rows[[pw]])
    p <- tryCatch(wilcox.test(test_rows[[pw]], ctrl_rows[[pw]])$p.value, error = function(e) NA)
    data.frame(analysis = analysis, label = label, cell_type = cell_type, timepoint = timepoint,
               pathway = pw, cohens_d = d, wilcox_p = p, n_test = nrow(test_rows), n_control = nrow(ctrl_rows))
  }))
}

all_results <- list()

# ── Meta-analysis: 5 conditions ─────────────────────────────────────────────
groups_meta <- read.csv(file.path(PT, "meta_groups.csv"))
meta_conditions <- c("Contact_inhibited CQ" = "CICQ", "Serum_starved CQ" = "SSCQ",
                     "Replicative CS" = "RS", "Stress-induced CS" = "SIPS",
                     "Oncogene-induced CS" = "OIS")
for (model in c("scaled", "yugene")) {
  scores_csv <- file.path(PT, sprintf("meta_partial_scores_%s.csv", model))
  for (cond_value in names(meta_conditions)) {
    r <- score_pathways(scores_csv, groups_meta, cond_value, "Proliferating",
                         meta_conditions[[cond_value]], "meta_analysis")
    r$model <- model
    all_results[[paste("meta", meta_conditions[[cond_value]], model)]] <- r
  }
}

# ── Temporal, pooled (irradiated vs none) ───────────────────────────────────
for (ct in c("Fibroblast", "Keratinocyte", "Melanocyte")) {
  groups_ct <- read.csv(file.path(PT, sprintf("%s_groups.csv", ct)))
  for (model in c("scaled", "yugene")) {
    scores_csv <- file.path(PT, sprintf("%s_partial_scores_%s.csv", ct, model))
    r <- score_pathways(scores_csv, groups_ct, "irradiated", "none", ct, "temporal_pooled", cell_type = ct)
    r$model <- model
    all_results[[paste("temporal_pooled", ct, model)]] <- r
  }
}

# ── Temporal, per-timepoint (3 celltypes x 3 timepoints, 6v6) ───────────────
for (ct in c("Fibroblast", "Keratinocyte", "Melanocyte")) {
  for (tp in c("4_days", "10_days", "20_days")) {
    grp <- paste0(ct, "_", tp)
    groups_grp <- read.csv(file.path(PT, sprintf("%s_groups.csv", grp)))
    for (model in c("scaled", "yugene")) {
      scores_csv <- file.path(PT, sprintf("%s_partial_scores_%s.csv", grp, model))
      r <- score_pathways(scores_csv, groups_grp, tp, "none", grp, "temporal_bytimepoint",
                           cell_type = ct, timepoint = tp)
      r$model <- model
      all_results[[paste("temporal_bytimepoint", grp, model)]] <- r
    }
  }
}

final <- do.call(rbind, all_results)
final$p_adj <- ave(final$wilcox_p, paste(final$analysis, final$model), FUN = function(x) p.adjust(x, method = "BH"))
write.csv(final, file.path(RERUN_DIR, "partial_tage_ALL.csv"), row.names = FALSE)

cat(sprintf("Done. %d rows across %d group-comparisons x 50 pathways x 2 models.\n",
            nrow(final), length(unique(paste(final$analysis, final$label)))))
cat(sprintf("Saved -> %s\n", file.path(RERUN_DIR, "partial_tage_ALL.csv")))

cat("\n== Sanity: n rows per analysis ==\n")
print(table(final$analysis, final$model))
