# 11_proliferating_baseline_heterogeneity.R
#
# How much does tAge vary among the PROLIFERATING CONTROLS themselves?
#
# WHY. Every condition effect in 2.1.5.1 is a shift relative to the pooled
# Proliferating group. That is only interpretable against a baseline: if
# untreated proliferating fibroblasts already differ from one another by tens of
# tAge units depending on which study or cell line they came from, then a
# condition effect of comparable size is not distinguishable from batch
# structure. This script measures that noise floor and asks which sample
# attributes explain it.
#
# WHAT THE ZERO MEANS HERE. control_subtraction sets the pooled Proliferating
# per-gene median as the reference, so this group sits at ~0 BY CONSTRUCTION.
# That constrains only its CENTRE, not its SPREAD - the between-study and
# between-line dispersion measured here is real and is not an artefact of
# centring. What cannot be recovered is any absolute baseline difference that is
# common to all controls.
#
# STRATIFIERS AVAILABLE. All 230 samples are fibroblasts (cell_type is constant),
# so "cell type" cannot be tested. The available structure is cell_line (7
# levels), tissue of origin (Lung / Foreskin / Skin), immortalisation status
# (audited, see script 10), and study.
#
# Tests: Kruskal-Wallis per attribute, then pairwise Wilcoxon with BH within
# attribute x model. Study is reported descriptively - with 1-6 controls per
# study a per-study test is not meaningful, but the spread of study medians is
# the quantity of interest.
#
# Output: rerun_outputs/proliferating_baseline_heterogeneity.csv

source("R/config.R")
suppressPackageStartupMessages(library(dplyr))

d <- read.csv(file.path(RERUN_DIR, "tage_all_conditions.csv"))
ann <- read.csv(file.path(RERUN_DIR, "immortalisation_annotation_corrected.csv"))
d$immortalised <- ann$immortalised[match(d$external_id, ann$external_id)]

MODELS <- c(scaled_diff = "scaled_diff_EN_tAge", yugene_diff = "yugene_diff_EN_tAge")
p <- d[d$condition == "Proliferating", ]
cat(sprintf("Proliferating controls: n = %d\n\n", nrow(p)))

# ---- 0. the noise floor, for comparison against the condition effects --------
cat("== 0. overall spread of the control group ==\n")
floor_rows <- do.call(rbind, lapply(names(MODELS), function(mn) {
  v <- p[[MODELS[[mn]]]]
  sm <- tapply(v, p$study, median)
  cl <- tapply(v, p$cell_line, median)
  data.frame(model = mn, n = length(v),
             median = median(v), IQR = IQR(v), SD = sd(v),
             range = diff(range(v)),
             n_studies = length(sm),
             study_median_range = diff(range(sm)),
             study_median_IQR = IQR(sm),
             cell_line_median_range = diff(range(cl)))
}))
print(floor_rows, row.names = FALSE, digits = 3)
cat("\n  Compare these to the condition effects (8-61 tAge units): any effect not\n")
cat("  clearly exceeding the between-study spread of the controls is not\n")
cat("  separable from study-level batch structure.\n")

# ---- 1. Kruskal-Wallis per attribute ----------------------------------------
cat("\n== 1. does any attribute explain control tAge? (Kruskal-Wallis) ==\n")
ATTRS <- c("cell_line", "tissue", "immortalised", "study")
kw <- do.call(rbind, lapply(ATTRS, function(a) {
  do.call(rbind, lapply(names(MODELS), function(mn) {
    s <- p[!is.na(p[[a]]), ]
    keep <- names(which(table(s[[a]]) >= 3))
    s <- s[s[[a]] %in% keep, ]
    if (length(unique(s[[a]])) < 2) return(NULL)
    k <- kruskal.test(s[[MODELS[[mn]]]], factor(s[[a]]))
    data.frame(test = "kruskal_within_proliferating", attribute = a, model = mn,
               n_levels = length(unique(s[[a]])), n = nrow(s),
               statistic = unname(k$statistic), p = k$p.value)
  }))
}))
kw$p_adj <- ave(kw$p, kw$model, FUN = function(x) p.adjust(x, "BH"))
print(kw, row.names = FALSE, digits = 3)

# ---- 2. per-level medians ----------------------------------------------------
cat("\n== 2. control tAge by level (levels with n >= 3) ==\n")
lev <- do.call(rbind, lapply(ATTRS, function(a) {
  do.call(rbind, lapply(names(MODELS), function(mn) {
    s <- p[!is.na(p[[a]]), ]
    do.call(rbind, lapply(names(which(table(s[[a]]) >= 3)), function(L) {
      v <- s[[MODELS[[mn]]]][s[[a]] == L]
      data.frame(test = "level_median_within_proliferating", attribute = a,
                 model = mn, level = L, n = length(v),
                 median = median(v), IQR = IQR(v))
    }))
  }))
}))
for (mn in names(MODELS)) {
  cat(sprintf("\n-- %s --\n", mn))
  print(lev[lev$model == mn, c("attribute", "level", "n", "median", "IQR")],
        row.names = FALSE, digits = 3)
}

# ---- 3. pairwise cell_line and tissue contrasts ------------------------------
cat("\n== 3. pairwise contrasts among control cell lines / tissues ==\n")
pw <- do.call(rbind, lapply(c("cell_line", "tissue"), function(a) {
  do.call(rbind, lapply(names(MODELS), function(mn) {
    s <- p[!is.na(p[[a]]), ]
    L <- names(which(table(s[[a]]) >= 3))
    cmb <- if (length(L) >= 2) combn(L, 2) else matrix(nrow = 2, ncol = 0)
    if (!ncol(cmb)) return(NULL)
    do.call(rbind, lapply(seq_len(ncol(cmb)), function(i) {
      x <- s[[MODELS[[mn]]]][s[[a]] == cmb[1, i]]
      y <- s[[MODELS[[mn]]]][s[[a]] == cmb[2, i]]
      data.frame(test = "pairwise_within_proliferating", attribute = a, model = mn,
                 level_a = cmb[1, i], level_b = cmb[2, i],
                 n_a = length(x), n_b = length(y),
                 median_diff = median(x) - median(y),
                 p = wilcox.test(x, y)$p.value)
    }))
  }))
}))
pw$p_adj <- ave(pw$p, paste(pw$attribute, pw$model), FUN = function(x) p.adjust(x, "BH"))
print(pw[order(pw$attribute, pw$model, pw$p), ], row.names = FALSE, digits = 3)

out <- bind_rows(floor_rows %>% mutate(test = "control_spread"), kw, lev, pw)
write.csv(out, file.path(RERUN_DIR, "proliferating_baseline_heterogeneity.csv"),
          row.names = FALSE)
cat(sprintf("\nSaved -> %s\n",
            file.path(RERUN_DIR, "proliferating_baseline_heterogeneity.csv")))
