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
# levels), tissue of origin (Lung / Foreskin among controls; verified per strain,
# see meta_analysis/10), immortalisation status
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
# `cell_line` in the metadata lumps 8 strains under "Primary" and puts
# IMR90-hTERT inside "IMR-90"; use the strain-level column from script 10.
d$cell_line <- ann$cell_line_resolved[match(d$external_id, ann$external_id)]
# TISSUE FROM THE VERIFIED MAP, not from the metadata column (2026-08-31). The
# metadata called "Skin" the 6 samples now recorded as HCA2 and 4 of the HDF series;
# HDF series' own paper (Mitra 2018) both say foreskin. See the TISSUE_VERIFIED
# table in meta_analysis/10. After the correction, no Skin sample remains among the
# proliferating controls - HDF161, the only genuinely adult dermal strain, is in a
# study with no internal controls - so the tissue contrast here is Foreskin against
# Lung, and any three-way Lung/Foreskin/Skin split is an artefact of the old column.
d$tissue <- ann$tissue_verified[match(d$external_id, ann$external_id)]

MODELS <- c(yugene_diff = "yugene_diff_EN_tAge", scaled_diff = "scaled_diff_EN_tAge")
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

# ---------------------------------------------------------------------------
# 4. CELL LINE WITHIN TISSUE, PRIMARY CELLS ONLY (added 2026-08-31)
#
# Sections above show that tissue of origin separates the controls and that
# strains differ, but strain and tissue are nested, so neither result isolates
# line identity. This asks the question the nesting leaves open: among primary
# cells of the SAME tissue, do different strains still differ? Immortalised
# samples are excluded so that hTERT status cannot contribute.
#
# It matters because "background" has been used loosely for two things -- tissue
# of origin and the particular strain -- and only this comparison separates them.
# ---------------------------------------------------------------------------
cat("\n== 4. between-strain variation WITHIN tissue, primary controls only ==\n")
# all three clocks here, so mortality is merged in (script 11 otherwise runs the
# two chronological models only)
mt <- read.csv(file.path(RERUN_DIR, "mortality_tage.csv"))
p$mortality_tAge <- mt$mortality_tAge[match(p$external_id, mt$external_id)]
MODELS3 <- c(MODELS, mortality = "mortality_tAge")
pri <- p[p$immortalised == "no", ]
wt <- list()
for (ti in sort(unique(pri$tissue))) {
  q <- pri[pri$tissue == ti, ]
  keep <- names(which(table(q$cell_line) >= 3))
  if (length(keep) < 2) {
    cat(sprintf("  %-9s only %d strain(s) with n>=3 - not testable\n", ti, length(keep)))
    next
  }
  q <- q[q$cell_line %in% keep, ]
  for (mn in names(MODELS3)) {
    v <- MODELS3[[mn]]
    med <- tapply(q[[v]], q$cell_line, median)
    kw <- kruskal.test(q[[v]], factor(q$cell_line))
    wt[[length(wt) + 1]] <- data.frame(
      test = "within_tissue_between_strain_primary", tissue = ti, model = mn,
      n = nrow(q), n_strains = length(keep),
      strain_median_range = diff(range(med)),
      lowest = names(med)[which.min(med)],  lowest_median = min(med),
      highest = names(med)[which.max(med)], highest_median = max(med),
      statistic = unname(kw$statistic), p = kw$p.value)
  }
}
if (length(wt)) {
  wt <- bind_rows(wt); wt$p_adj <- p.adjust(wt$p, "BH")
  print(wt[, c("tissue", "model", "n", "n_strains", "strain_median_range",
               "lowest", "lowest_median", "highest", "highest_median", "p", "p_adj")],
        row.names = FALSE, digits = 3)
  cat("\n  For comparison, the BETWEEN-tissue difference in the same primary cells:\n")
  for (mn in names(MODELS3)) {
    v <- MODELS3[[mn]]
    f <- pri[[v]][pri$tissue == "Foreskin"]; l <- pri[[v]][pri$tissue == "Lung"]
    cat(sprintf("    %-12s foreskin %+.1f vs lung %+.1f -> %+.1f\n",
                mn, median(f), median(l), median(f) - median(l)))
  }
  cat("\n  => line identity is not absorbed by tissue: it is significant within lung on\n")
  cat("     all three clocks and within foreskin on the mortality clock. Tissue is the\n")
  cat("     larger effect on the scaled-difference and mortality clocks, but on YuGene\n")
  cat("     the within-lung spread exceeds the between-tissue difference.\n")
  write.csv(wt, file.path(RERUN_DIR, "within_tissue_between_strain.csv"), row.names = FALSE)
  cat(sprintf("\nSaved -> %s\n", file.path(RERUN_DIR, "within_tissue_between_strain.csv")))
}
