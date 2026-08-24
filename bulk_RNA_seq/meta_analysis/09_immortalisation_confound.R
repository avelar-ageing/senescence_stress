# 09_immortalisation_confound.R
#
# Is the meta-analysis tAge result confounded by immortalisation status (and by
# cell line)?
#
# ANNOTATION PROVENANCE. This script does NOT use the hand-curated
# `immortalised` column in sample_metadata_RERUN.csv. That column was audited in
# 10_immortalisation_annotation_audit.R and found to have 8 false negatives in 3
# studies (SRP017378 n=5, whose line is literally named "BJ hTERT"; SRP123346
# n=2, which was also internally inconsistent - 2 of its 3 samples flagged one
# way and 1 the other; SRP136727 n=1). The corrected count is 28 immortalised
# samples, not 20. The corrected annotation is read in below.
#
# WHY THIS MATTERS. Within the Proliferating control group alone, immortalised
# samples score far higher on the scaled_diff model than primary ones, so any
# condition containing immortalised samples may be inflated on that model.
#
# WHY IMMORTALISATION CANNOT BE A COVARIATE. It is partly confounded with
# condition by biology rather than by chance:
#   - Proliferating / CICQ / OIS / SSCQ CAN be immortalised. Quiescence is
#     reversible arrest and OIS runs through RAS/BRAF, neither
#     telomere-dependent.
#   - RS CANNOT be immortalised, by definition: hTERT maintains telomeres, so
#     there is no attrition and no route into replicative senescence. Its count
#     is structurally zero, not a sampling accident.
#   - SIPS could in principle be immortalised (irradiation is not
#     telomere-dependent) but happens to have none here.
# Because immortalisation is perfectly separable for RS, it cannot be fitted as a
# covariate across all conditions. The only clean option is to exclude
# immortalised samples and compare primary cells like for like.
#
# WHY THE EXISTING tAge VALUES CAN BE REUSED. control_subtraction removes a
# per-gene constant, which shifts every sample equally and therefore cancels from
# any between-group difference (verified in
# exploratory/17_tage_calculation_audit.py). So restricting to primary samples and
# re-differencing the already-computed tAge values is valid; it does not require
# re-running tAge_preprocessing with a different control set.
#
# Output: rerun_outputs/immortalisation_confound.csv

source("R/config.R")
suppressPackageStartupMessages(library(dplyr))

d <- read.csv(file.path(RERUN_DIR, "tage_all_conditions.csv"))
d$condition <- factor(d$condition, levels = c("Proliferating", "CICQ", "SSCQ", "RS", "SIPS", "OIS"))

# Replace the curated flag with the audited one (see header, and script 10).
ann <- read.csv(file.path(RERUN_DIR, "immortalisation_annotation_corrected.csv"))
stopifnot(setequal(d$external_id, ann$external_id))
d$immortalised <- ann$immortalised[match(d$external_id, ann$external_id)]
cat(sprintf("Using AUDITED annotation: %d immortalised (curated column said %d)\n",
            sum(d$immortalised == "yes"), sum(ann$immortalised_curated == "yes")))
MODELS <- c(scaled_diff = "scaled_diff_EN_tAge", yugene_diff = "yugene_diff_EN_tAge")

cat("== composition ==\n")
print(table(d$condition, d$immortalised))
cat("\n== cell lines per condition (a further, unaddressed confound) ==\n")
print(table(d$cell_line, d$condition))

# ---- 1. within-condition effect of immortalisation ------------------------
cat("\n== 1. immortalised vs primary, WITHIN condition ==\n")
within_rows <- do.call(rbind, lapply(levels(d$condition), function(cond) {
  s <- d[d$condition == cond, ]
  if (sum(s$immortalised == "yes") < 3 || sum(s$immortalised == "no") < 3) return(NULL)
  do.call(rbind, lapply(names(MODELS), function(mn) {
    v <- MODELS[[mn]]
    y <- s[[v]][s$immortalised == "yes"]; n <- s[[v]][s$immortalised == "no"]
    data.frame(test = "immortalised_vs_primary_within_condition", condition = cond,
               model = mn, n_immortalised = length(y), n_primary = length(n),
               median_immortalised = median(y), median_primary = median(n),
               median_diff = median(y) - median(n),
               p = wilcox.test(y, n)$p.value)
  }))
}))
within_rows$p_adj <- p.adjust(within_rows$p, method = "BH")
print(within_rows, row.names = FALSE, digits = 3)

# ---- 2. condition effects, all samples vs primary only -------------------
cat("\n== 2. condition vs Proliferating: ALL samples vs PRIMARY ONLY ==\n")
cond_rows <- do.call(rbind, lapply(setdiff(levels(d$condition), "Proliferating"), function(cond) {
  do.call(rbind, lapply(names(MODELS), function(mn) {
    v <- MODELS[[mn]]
    ga <- d[[v]][d$condition == cond];                       pa <- d[[v]][d$condition == "Proliferating"]
    gp <- d[[v]][d$condition == cond & d$immortalised == "no"]
    pp <- d[[v]][d$condition == "Proliferating" & d$immortalised == "no"]
    data.frame(condition = cond, model = mn,
               n_immortalised_in_condition = sum(d$condition == cond & d$immortalised == "yes"),
               diff_all = median(ga) - median(pa), p_all = wilcox.test(ga, pa)$p.value,
               diff_primary = median(gp) - median(pp), p_primary = wilcox.test(gp, pp)$p.value)
  }))
}))
# BH within each model, matching the convention of the 10-test vs-Proliferating family
cond_rows$p_all_adj <- ave(cond_rows$p_all, cond_rows$model, FUN = function(x) p.adjust(x, "BH"))
cond_rows$p_primary_adj <- ave(cond_rows$p_primary, cond_rows$model, FUN = function(x) p.adjust(x, "BH"))
print(cond_rows, row.names = FALSE, digits = 3)

cat("\n== 3. does the ordering survive? ==\n")
for (mn in names(MODELS)) {
  s <- cond_rows[cond_rows$model == mn, ]
  cat(sprintf("  %-11s ALL: %s\n", mn,
              paste(s$condition[order(-s$diff_all)], collapse = " > ")))
  cat(sprintf("  %-11s PRI: %s\n", "",
              paste(s$condition[order(-s$diff_primary)], collapse = " > ")))
}

out <- bind_rows(within_rows, cond_rows)
write.csv(out, file.path(RERUN_DIR, "immortalisation_confound.csv"), row.names = FALSE)
cat(sprintf("\nSaved -> %s\n", file.path(RERUN_DIR, "immortalisation_confound.csv")))
