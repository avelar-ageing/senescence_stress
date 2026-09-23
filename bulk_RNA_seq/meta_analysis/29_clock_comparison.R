# 29_clock_comparison.R
#
# The numbers behind "which clock is more reliable is not settled by dispersion"
# in section 2.1.5.
#
# WHY THIS EXISTS. That passage quoted six figures -- study-level control
# heterogeneity in control-SD units for two clocks, and the immortalisation
# artefact in control-SD units with a p-value for three -- and no script in the
# pipeline produced any of them. They were computed ad hoc in an earlier session.
# Two of the p-values were then left stale by the SRP089801 reclassification
# (meta_analysis/10), which changed the immortalised proliferating group from 21
# samples to 18. Reverse-engineering the old values confirms both the definitions
# below and that they predate that change: 31.14/27.35 = 1.14, 6.60/15.74 = 0.42
# and 0.54/0.71 = 0.76 are exactly the old gaps over the control SDs.
#
# DEFINITIONS, made explicit because "control standard deviations" is ambiguous.
# The denominator throughout is the SD of ALL 91 untreated proliferating samples
# on that clock -- not the within-study SD, and not the primary-only SD. Both
# numerators are then expressed in those units:
#   heterogeneity  = range of per-study median tAge among proliferating controls
#   immortalisation = median(immortalised controls) - median(primary controls)
# The p-value is a two-sided Mann-Whitney on the same two groups, and is the same
# quantity meta_analysis/16 section A reports unnormalised.
#
# The point of the comparison is that it does NOT select a clock: a scale-free
# model can look better or worse depending on which spread you divide by, so
# effect-over-control-SD is a circular selection rule. The one asymmetry that is
# testable is that the two scaled-difference models register the study-nested
# immortalisation artefact and YuGene does not.
#
# Output: rerun_outputs/clock_comparison.csv

source("R/config.R")
suppressPackageStartupMessages({ library(dplyr) })

d <- read.csv(file.path(RERUN_DIR, "tage_all_conditions.csv"))
ann <- read.csv(file.path(RERUN_DIR, "immortalisation_annotation_corrected.csv"))
mort <- read.csv(file.path(RERUN_DIR, "mortality_tage.csv"))
i <- match(d$external_id, ann$external_id); stopifnot(!any(is.na(i)))
d$immortalised <- ann$immortalised[i]
d$study        <- ann$study[i]
d$mortality_tAge <- mort$mortality_tAge[match(d$external_id, mort$external_id)]

# CLOCK ORDER (2026-09-01): YuGene, Scaled Difference, Mortality everywhere -
# figures, printed tables and the manuscript text all follow this order.
CLOCKS <- c(yugene_diff = "yugene_diff_EN_tAge",
            scaled_diff = "scaled_diff_EN_tAge",
            mortality   = "mortality_tAge")

p <- d[d$cell_state == "Proliferating", ]
cat(sprintf("untreated proliferating samples: %d (%d immortalised, %d primary)\n\n",
            nrow(p), sum(p$immortalised == "yes"), sum(p$immortalised == "no")))

rows <- lapply(names(CLOCKS), function(nm) {
  v <- CLOCKS[[nm]]
  x <- p[[v]]
  control_sd <- sd(x)                       # denominator: all 91 controls
  sm <- tapply(x, p$study, median)
  het_range <- diff(range(sm))
  y <- x[p$immortalised == "yes"]; n <- x[p$immortalised == "no"]
  gap <- median(y) - median(n)
  pv <- suppressWarnings(wilcox.test(y, n))$p.value
  data.frame(model = nm, n_controls = length(x), control_sd = control_sd,
             n_studies = length(sm),
             study_median_range = het_range,
             heterogeneity_in_control_sd = het_range / control_sd,
             n_immortalised = length(y), n_primary = length(n),
             immortalisation_gap = gap,
             immortalisation_in_control_sd = gap / control_sd,
             p = pv)
})
out <- bind_rows(rows)

cat("== study-level heterogeneity of untreated controls ==\n")
print(out[, c("model", "control_sd", "study_median_range", "heterogeneity_in_control_sd")],
      row.names = FALSE, digits = 3)
cat("\n== the study-nested immortalisation artefact, same denominator ==\n")
print(out[, c("model", "n_immortalised", "n_primary", "immortalisation_gap",
              "immortalisation_in_control_sd", "p")], row.names = FALSE, digits = 3)

cat("\n  Quote as: heterogeneity ")
cat(sprintf("%.1f (scaled) vs %.1f (YuGene) control SDs; artefact %+.2f, %+.2f and %+.2f\n",
            out$heterogeneity_in_control_sd[out$model == "scaled_diff"],
            out$heterogeneity_in_control_sd[out$model == "yugene_diff"],
            out$immortalisation_in_control_sd[out$model == "scaled_diff"],
            out$immortalisation_in_control_sd[out$model == "mortality"],
            out$immortalisation_in_control_sd[out$model == "yugene_diff"]))
cat(sprintf("  with p = %.3g (scaled), %.3g (mortality), %.2f (YuGene)\n",
            out$p[out$model == "scaled_diff"], out$p[out$model == "mortality"],
            out$p[out$model == "yugene_diff"]))

write.csv(out, file.path(RERUN_DIR, "clock_comparison.csv"), row.names = FALSE)
cat(sprintf("\nSaved -> %s\n", file.path(RERUN_DIR, "clock_comparison.csv")))
