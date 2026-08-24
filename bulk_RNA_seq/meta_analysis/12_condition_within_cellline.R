# 12_condition_within_cellline.R
#
# Do the condition effects survive when each condition is compared to
# Proliferating controls OF THE SAME CELL LINE?
#
# SUPERSEDED FOR THE HEADLINE ESTIMATE - see 13_condition_within_study.R. Cell
# line is the wrong stratum: a line is shared across studies (IMR90 appears in 15
# of them) while batch, passage, protocol and library prep are study-specific, so
# a within-line contrast can still be a cross-study contrast. Script 13 uses
# study, which holds line AND batch constant. This script is retained because it
# answers the narrower question of whether an effect reproduces on more than one
# genetic background, which study-level stratification does not show directly.
#
# WHY. 11_proliferating_baseline_heterogeneity.R showed the untreated controls are
# not exchangeable: among Proliferating samples alone, IMR90 sits at -15.3 tAge
# units and HCA2/BJ primary fibroblasts at +45.4 on scaled_diff, a 61-unit gap,
# study medians span 87.8 units, and cell-line medians span 70.7. Condition
# effects are 8-61 units, so the pooled comparison in 2.1.5.1 is only
# interpretable if condition and cell line are not confounded. They are: CICQ
# contains no IMR90 at all, OIS is 28/48 IMR90-derived, while Proliferating is
# 33/91 IMR90.
#
# CELL LINE HERE IS THE RESOLVED STRAIN, not the metadata's `cell_line`, which
# lumps 8 unrelated strains under "Primary" and hides IMR90-hTERT inside
# "IMR-90" (see script 10).
#
# LIMITS OF THIS TEST, stated up front. Stratifying costs power and not every
# condition x line cell exists: there is no CICQ IMR90 and no OIS HDF strain, so
# those conditions can never be compared on the same line. Where a stratum has
# n < 3 on either side it is skipped, which drops RS to a single line and drops
# every 2-vs-2 CICQ study. A surviving effect in >= 2 lines is evidence the
# effect is not a line artefact; a single-stratum effect is not.
#
# Output: rerun_outputs/condition_within_cellline.csv

source("R/config.R")
suppressPackageStartupMessages(library(dplyr))

d <- read.csv(file.path(RERUN_DIR, "tage_all_conditions.csv"))
d$condition <- factor(d$condition, levels = c("Proliferating", "CICQ", "SSCQ", "RS", "SIPS", "OIS"))
# Strain-level cell line from script 10. The metadata's own `cell_line` cannot be
# used here: "Primary" lumps 8 unrelated strains from 8 studies and "IMR-90"
# contains IMR90-hTERT, so stratifying on it does not hold cell line constant.
ann <- read.csv(file.path(RERUN_DIR, "immortalisation_annotation_corrected.csv"))
d$cell_line <- ann$cell_line_resolved[match(d$external_id, ann$external_id)]
MODELS <- c(scaled_diff = "scaled_diff_EN_tAge", yugene_diff = "yugene_diff_EN_tAge")
CONDS <- setdiff(levels(d$condition), "Proliferating")
MIN_N <- 3

cat("== condition x cell line availability (Proliferating = the control column) ==\n")
print(table(d$cell_line, d$condition))
cat("\n== condition x tissue ==\n")
print(table(d$tissue, d$condition))

strat_test <- function(strat_col, label) {
  rows <- do.call(rbind, lapply(sort(unique(d[[strat_col]])), function(L) {
    s <- d[d[[strat_col]] == L, ]
    ctrl <- s[s$condition == "Proliferating", ]
    if (nrow(ctrl) < MIN_N) return(NULL)
    do.call(rbind, lapply(CONDS, function(cond) {
      tst <- s[s$condition == cond, ]
      if (nrow(tst) < MIN_N) return(NULL)
      do.call(rbind, lapply(names(MODELS), function(mn) {
        v <- MODELS[[mn]]
        x <- tst[[v]]; y <- ctrl[[v]]
        data.frame(test = paste0("condition_within_", label), stratum = L,
                   condition = cond, model = mn,
                   n_test = length(x), n_control = length(y),
                   median_test = median(x), median_control = median(y),
                   diff_within = median(x) - median(y),
                   p = wilcox.test(x, y)$p.value)
      }))
    }))
  }))
  rows$p_adj <- ave(rows$p, rows$model, FUN = function(x) p.adjust(x, "BH"))
  rows
}

by_line <- strat_test("cell_line", "cellline")
by_tissue <- strat_test("tissue", "tissue")

# pooled effects for reference
pooled <- do.call(rbind, lapply(CONDS, function(cond) {
  do.call(rbind, lapply(names(MODELS), function(mn) {
    v <- MODELS[[mn]]
    data.frame(condition = cond, model = mn,
               diff_pooled = median(d[[v]][d$condition == cond]) -
                             median(d[[v]][d$condition == "Proliferating"]))
  }))
}))

for (mn in names(MODELS)) {
  cat(sprintf("\n\n===== %s : condition vs same-line Proliferating =====\n", mn))
  s <- by_line[by_line$model == mn, ]
  s <- merge(s, pooled[pooled$model == mn, c("condition", "diff_pooled")], by = "condition")
  print(s[order(s$condition, s$stratum),
          c("condition", "stratum", "n_test", "n_control", "diff_within",
            "diff_pooled", "p", "p_adj")], row.names = FALSE, digits = 3)

  cat(sprintf("\n-- %s : summary per condition --\n", mn))
  agg <- s %>% group_by(condition) %>%
    summarise(n_strata = n(),
              n_sig = sum(p_adj < 0.05),
              strata_sig = paste(stratum[p_adj < 0.05], collapse = ","),
              diff_pooled = first(diff_pooled),
              diff_within_median = median(diff_within),
              diff_within_range = sprintf("%.1f to %.1f", min(diff_within), max(diff_within)),
              .groups = "drop")
  print(as.data.frame(agg), row.names = FALSE, digits = 3)
}

cat("\n\n===== same, stratified by tissue of origin =====\n")
print(by_tissue[order(by_tissue$model, by_tissue$condition, by_tissue$stratum),
                c("stratum", "condition", "model", "n_test", "n_control",
                  "diff_within", "p", "p_adj")], row.names = FALSE, digits = 3)

out <- bind_rows(by_line, by_tissue)
write.csv(out, file.path(RERUN_DIR, "condition_within_cellline.csv"), row.names = FALSE)
cat(sprintf("\nSaved -> %s\n", file.path(RERUN_DIR, "condition_within_cellline.csv")))
