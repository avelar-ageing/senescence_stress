# 18_partial_tage_within_study.R
#
# Re-estimates the META-ANALYSIS pathway contributions within study, so 2.1.5.2
# rests on the same de-confounded design as 2.1.5.1.
#
# WHY. meta_analysis/11 showed the untreated Proliferating controls are not
# exchangeable (study medians of whole-transcriptome tAge span 87.8 units on
# scaled_diff), and meta_analysis/13 showed that pooling across studies inflates
# or invents whole-transcriptome effects: CICQ falls from +47.6 to +5.9 and SSCQ
# from +21.7 to -3.4 once each sample is compared only to its own study's
# controls. The pathway contributions in partial_tage_ALL.csv were computed from
# exactly the same pooled contrasts, so they inherit exactly the same confound.
# A pathway contribution is a weighted sum over that set's genes of the same
# control-subtracted matrix, so nothing protects it.
#
# SCOPE: the meta-analysis only. The time course is a SINGLE study (ERP021140,
# all 72 samples), where every comparison is already within study and within cell
# type, so the temporal results need no correction on this account.
#
# METHOD: identical to meta_analysis/13, applied per pathway. Per study,
# d_s = mean(test) - mean(control); combined as sum(w_s d_s)/sum(w_s) with
# w_s = n_t n_c/(n_t + n_c); null by permuting the condition label WITHIN study.
# The permutation is shared across all 50 pathways for a given condition and
# model, which is both faster and correct - the pathways are re-scored under the
# same relabelling, preserving their mutual correlation.
#
# BH is applied across 50 pathways x 5 conditions within each model (250 tests),
# matching the existing convention in 05_consolidate_partial_scores.R.
#
# Output: rerun_outputs/partial_tage_within_study.csv

source("R/config.R")
suppressPackageStartupMessages(library(dplyr))

set.seed(1)
NPERM <- 10000
PT <- file.path(RERUN_DIR, "partial_tage")
CONDS <- c("Contact_inhibited CQ" = "CICQ", "Serum_starved CQ" = "SSCQ",
           "Replicative CS" = "RS", "Stress-induced CS" = "SIPS",
           "Oncogene-induced CS" = "OIS")

groups <- read.csv(file.path(PT, "meta_groups.csv"))
meta <- read.csv(file.path(RERUN_DIR, "sample_metadata_RERUN.csv"))
groups$study <- meta$study[match(groups$sample_id, meta$external_id)]
stopifnot(!any(is.na(groups$study)))

old <- read.csv(file.path(RERUN_DIR, "partial_tage_ALL.csv"))

res <- list()
for (model in c("scaled", "yugene")) {
  sc <- read.csv(file.path(PT, sprintf("meta_partial_scores_%s.csv", model)),
                 check.names = FALSE)
  pw_cols <- setdiff(colnames(sc), c("sample_id", "full_tAge_direct",
                                     "full_tAge_reconstructed"))
  rownames(sc) <- sc$sample_id
  S <- as.matrix(sc[, pw_cols, drop = FALSE])

  for (cond_value in names(CONDS)) {
    lab <- CONDS[[cond_value]]
    g <- groups[groups$group %in% c(cond_value, "Proliferating"), ]
    # keep only studies containing both arms
    keep <- names(which(tapply(g$group, g$study,
                               function(x) length(unique(x)) == 2)))
    g <- g[g$study %in% keep, ]
    blocks <- lapply(split(g, g$study), function(b) {
      A <- S[b$sample_id, , drop = FALSE]
      nt <- sum(b$group == cond_value); nc <- sum(b$group == "Proliferating")
      list(A = A, nt = nt, nc = nc, colsum = colSums(A),
           test_rows = which(b$group == cond_value),
           w = nt * nc / (nt + nc))
    })
    W <- sum(vapply(blocks, function(b) b$w, numeric(1)))

    # diff_s = (1/nt + 1/nc) * colSums(test rows) - (1/nc) * colSums(all rows)
    combine <- function(pick) {
      acc <- numeric(ncol(S))
      for (i in seq_along(blocks)) {
        b <- blocks[[i]]
        cs <- colSums(b$A[pick[[i]], , drop = FALSE])
        acc <- acc + b$w * ((1 / b$nt + 1 / b$nc) * cs - b$colsum / b$nc)
      }
      acc / W
    }
    obs <- combine(lapply(blocks, function(b) b$test_rows))
    cnt <- numeric(ncol(S))
    for (p in seq_len(NPERM)) {
      pick <- lapply(blocks, function(b) sample.int(b$nt + b$nc, b$nt))
      cnt <- cnt + (abs(combine(pick)) >= abs(obs))
    }
    p_perm <- (1 + cnt) / (NPERM + 1)

    # per-study sign consistency
    ds <- vapply(blocks, function(b) {
      (1 / b$nt + 1 / b$nc) * colSums(b$A[b$test_rows, , drop = FALSE]) -
        b$colsum / b$nc
    }, numeric(ncol(S)))
    npos <- rowSums(ds > 0)

    o <- old[old$analysis == "meta_analysis" & old$label == lab &
               old$model == model, ]
    res[[paste(model, lab)]] <- data.frame(
      analysis = "meta_analysis", label = lab, model = model,
      pathway = colnames(S), n_studies = length(blocks),
      contrib_diff_within_study = obs,
      contrib_diff_pooled = o$contrib_diff[match(colnames(S), o$pathway)],
      wilcox_p_pooled = o$wilcox_p[match(colnames(S), o$pathway)],
      p_adj_pooled = o$p_adj[match(colnames(S), o$pathway)],
      studies_positive = npos, p_perm = p_perm, row.names = NULL)
  }
  cat(sprintf("%s done\n", model))
}

final <- bind_rows(res)
final$p_perm_adj <- ave(final$p_perm, final$model,
                        FUN = function(x) p.adjust(x, "BH"))
write.csv(final, file.path(RERUN_DIR, "partial_tage_within_study.csv"),
          row.names = FALSE)

cat("\n== how many pathways are significant, pooled vs within study? ==\n")
print(final %>% group_by(model, label) %>%
        summarise(n_sig_pooled = sum(p_adj_pooled < 0.05, na.rm = TRUE),
                  n_sig_within = sum(p_perm_adj < 0.05),
                  n_lost = sum(p_adj_pooled < 0.05 & p_perm_adj >= 0.05, na.rm = TRUE),
                  n_gained = sum(p_adj_pooled >= 0.05 & p_perm_adj < 0.05, na.rm = TRUE),
                  n_sign_flip = sum(sign(contrib_diff_pooled) !=
                                    sign(contrib_diff_within_study)),
                  .groups = "drop"), n = 20)

cat("\n== pathways significant within study on BOTH models, per condition ==\n")
both <- final %>% filter(p_perm_adj < 0.05) %>%
  count(label, pathway) %>% filter(n == 2)
for (L in unique(both$label)) {
  cat(sprintf("\n-- %s (%d) --\n", L, sum(both$label == L)))
  cat(paste(" ", both$pathway[both$label == L]), sep = "\n")
}
cat(sprintf("\nSaved -> %s\n", file.path(RERUN_DIR, "partial_tage_within_study.csv")))
