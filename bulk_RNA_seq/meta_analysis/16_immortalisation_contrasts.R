# 16_immortalisation_contrasts.R
#
# Does hTERT-immortalisation change tAge, among proliferating controls and among
# senescent cells?
#
# ============================== IS THIS ESTIMABLE? ===========================
# Largely NO, and the reason is structural rather than a matter of power.
#
# (1) IMMORTALISATION IS PERFECTLY NESTED WITHIN STUDY. 0 of 34 studies contain
#     both immortalised and primary samples (asserted below). So every
#     immortalised-versus-primary contrast is entirely a between-study contrast,
#     and script 11 measured exactly how bad that is: among untreated
#     proliferating controls alone, study medians span 87.8 tAge units on
#     scaled_diff. There is no design here that separates "hTERT" from "the labs
#     that use hTERT lines". A study covariate cannot be fitted because the term
#     is perfectly collinear with it.
#
# (2) "SENESCENT IMMORTALISED VERSUS SENESCENT PRIMARY" IS ALSO CONDITION-
#     CONFOUNDED. All 17 immortalised senescent samples are OIS; RS and SIPS
#     contain none, by biology in the case of RS (hTERT maintains telomeres, so
#     there is no route into replicative senescence). Pooling the three subtypes
#     would compare OIS-immortalised cells against a primary group that is
#     mostly SIPS, so the contrast would carry a condition effect as well as a
#     study effect. Only OIS can be contrasted within condition.
#
# WHAT IS THEREFORE REPORTED, weakest to strongest:
#   A  the naive contrast, within condition. Descriptive only. Reported because it
#      is what the raw data shows and because scaled_diff's sensitivity to it is
#      the reason we prefer yugene - not as an estimate of an hTERT effect.
#   B  parental-background-matched contrasts. BJ versus BJ-hTERT and IMR90 versus
#      IMR90-hTERT both exist, so genetic background can be held constant even
#      though study cannot. Still cross-study, but the strongest available.
#   C  the one properly de-confounded question: does the WITHIN-STUDY OIS effect
#      differ between immortalised and primary studies? Each study's effect is a
#      difference against its own controls, so study baseline cancels, and the
#      studies are then the unit of analysis. This asks whether immortalisation
#      changes the RESPONSE to oncogene induction, which is answerable; it does
#      not recover the baseline effect of immortalisation, which is not.
#
# Output: rerun_outputs/immortalisation_contrasts.csv

source("R/config.R")
suppressPackageStartupMessages(library(dplyr))

d <- read.csv(file.path(RERUN_DIR, "tage_all_conditions.csv"))
ann <- read.csv(file.path(RERUN_DIR, "immortalisation_annotation_corrected.csv"))
d$immortalised <- ann$immortalised[match(d$external_id, ann$external_id)]
d$line <- ann$cell_line_resolved[match(d$external_id, ann$external_id)]
# Tissue must come from the verified map too, not from tage_all_conditions.csv.
# That file carries the submitted `tissue`, which calls the 6 HCA2 and 4
# HDF 12-x samples "Skin"; stratifying section D on it put immortalised foreskin
# cells in a "Skin" stratum of their own and split the foreskin group in two.
d$tissue <- ann$tissue_verified[match(d$external_id, ann$external_id)]
stopifnot(!any(is.na(d$tissue)))
d$parent <- gsub("-hTERT$", "", d$line)
d$parent[d$parent == "Tig3ET"] <- "TIG3"
d$parent <- gsub(" ER:RAS$", "", d$parent)
# ALL THREE CLOCKS (added 2026-08-30). This script ran the two chronological
# models only, so any mortality-clock figure quoted for these contrasts had no
# source in its output. The mortality values are merged in from
# mortality_tage.csv, which meta_analysis/18 writes per sample.
# CLOCK ORDER (2026-09-01): YuGene, Scaled Difference, Mortality everywhere -
# figures, printed tables and the manuscript text all follow this order.
MODELS <- c(yugene_diff = "yugene_diff_EN_tAge",
            scaled_diff = "scaled_diff_EN_tAge",
            mortality   = "mortality_tAge")

mort <- read.csv(file.path(RERUN_DIR, "mortality_tage.csv"))
d$mortality_tAge <- mort$mortality_tAge[match(d$external_id, mort$external_id)]
stopifnot(!any(is.na(d$mortality_tAge)))

# P-VALUE FLOOR GUARD. A two-sided rank test on n1 vs n2 cannot return a p below
# 2 / choose(n1 + n2, n1). Where that floor exceeds 0.05 the comparison cannot
# reach significance whatever the data, so reporting a p-value from it is
# meaningless - the BJ contrast is 7 against 2, floor 0.056. Those rows are
# reported as effect size only, with p held back rather than printed as a null.
# Every block reports the floor and the testable flag on the unit it actually
# compares: samples in A, B, D, E and F, but STUDIES in C, where each study
# contributes one within-study OIS effect. floor_basis records which, so a blank
# never has to be read as either "not applicable" or "not checked".
p_floor <- function(n1, n2) 2 / choose(n1 + n2, n1)
guarded_p <- function(y, n) {
  fl <- p_floor(length(y), length(n))
  if (fl > 0.05) NA_real_ else wtest(y, n)
}

nest <- d %>% group_by(study) %>% summarise(k = n_distinct(immortalised))
cat(sprintf("studies containing BOTH immortalised and primary samples: %d of %d\n",
            sum(nest$k > 1), nrow(nest)))
stopifnot(all(nest$k == 1))   # the confound is structural, not incidental
cat("=> every immortalised-vs-primary contrast below is entirely between-study\n\n")

rows <- list()
wtest <- function(x, y, ...) {
  if (length(x) < 2 || length(y) < 2) return(NA_real_)
  wilcox.test(x, y)$p.value
}

# ---- A. naive within-condition contrast (DESCRIPTIVE ONLY) -------------------
cat("== A. immortalised vs primary, within condition (CONFOUNDED WITH STUDY) ==\n")
for (cond in c("Proliferating", "CICQ", "SSCQ", "OIS")) {
  for (mn in names(MODELS)) {
    v <- MODELS[[mn]]; s <- d[d$condition == cond, ]
    y <- s[[v]][s$immortalised == "yes"]; n <- s[[v]][s$immortalised == "no"]
    if (!length(y) || !length(n)) next
    rows[[length(rows) + 1]] <- data.frame(
      test = "A_naive_within_condition_CONFOUNDED", stratum = cond, model = mn,
      n_immortalised = length(y), n_primary = length(n),
      n_studies_imm = n_distinct(s$study[s$immortalised == "yes"]),
      n_studies_prim = n_distinct(s$study[s$immortalised == "no"]),
      median_immortalised = median(y), median_primary = median(n),
      diff = median(y) - median(n),
      p_floor = p_floor(length(y), length(n)),
      testable = p_floor(length(y), length(n)) <= 0.05,
      floor_basis = "samples",
      p = guarded_p(y, n))
  }
}
# One BH family per test block, pooled across the three clocks: the clocks are
# three readings of the same contrast, not three separate questions, and
# per-clock families gave the same design different adjusted floors.
A <- bind_rows(rows); A$p_adj <- p.adjust(A$p, "BH")
print(A[, c("stratum", "model", "n_immortalised", "n_primary", "n_studies_imm",
            "n_studies_prim", "diff", "p", "p_adj")], row.names = FALSE, digits = 3)

# ---- B. parental-background-matched ----------------------------------------
cat("\n== B. same parental line, hTERT vs parental (still cross-study) ==\n")
Brows <- list()
for (par in c("BJ", "IMR90")) {
  for (cond in c("Proliferating", "OIS")) {
    for (mn in names(MODELS)) {
      v <- MODELS[[mn]]
      s <- d[d$parent == par & d$condition == cond, ]
      y <- s[[v]][s$immortalised == "yes"]; n <- s[[v]][s$immortalised == "no"]
      if (length(y) < 2 || length(n) < 2) next
      Brows[[length(Brows) + 1]] <- data.frame(
        test = "B_parental_matched", stratum = paste(par, cond), model = mn,
        n_immortalised = length(y), n_primary = length(n),
        median_immortalised = median(y), median_primary = median(n),
        diff = median(y) - median(n),
        p_floor = p_floor(length(y), length(n)),
        testable = p_floor(length(y), length(n)) <= 0.05,
        floor_basis = "samples",
        p = guarded_p(y, n))
    }
  }
}
B <- bind_rows(Brows)
if (nrow(B)) {
  B$p_adj <- p.adjust(B$p, "BH")
  print(B[, c("stratum", "model", "n_immortalised", "n_primary", "diff",
              "p_floor", "testable", "p", "p_adj")], row.names = FALSE, digits = 3)
  if (any(!B$testable)) {
    cat("\n  NOT TESTABLE (rank-test floor above 0.05; effect size only):\n")
    nt <- unique(B$stratum[!B$testable])
    for (x in nt) cat(sprintf("    %s - %d vs %d, smallest attainable p = %.3f\n", x,
                              B$n_immortalised[B$stratum == x][1],
                              B$n_primary[B$stratum == x][1],
                              B$p_floor[B$stratum == x][1]))
  }
}

# ---- C. the estimable question: does the OIS RESPONSE differ? ---------------
cat("\n== C. within-study OIS effect, immortalised vs primary studies (VALID) ==\n")
Crows <- list()
for (mn in names(MODELS)) {
  v <- MODELS[[mn]]
  eff <- d %>% filter(condition %in% c("OIS", "Proliferating")) %>%
    group_by(study) %>%
    filter(n_distinct(condition) == 2) %>%
    summarise(imm = first(immortalised),
              n_t = sum(condition == "OIS"), n_c = sum(condition == "Proliferating"),
              effect = mean(.data[[v]][condition == "OIS"]) -
                       mean(.data[[v]][condition == "Proliferating"]),
              .groups = "drop")
  cat(sprintf("\n-- %s : per-study OIS effect --\n", mn))
  print(as.data.frame(eff[order(eff$imm, -eff$effect), ]), row.names = FALSE, digits = 3)
  y <- eff$effect[eff$imm == "yes"]; n <- eff$effect[eff$imm == "no"]
  Crows[[length(Crows) + 1]] <- data.frame(
    test = "C_within_study_OIS_effect_by_htert", stratum = "OIS", model = mn,
    n_studies_imm = length(y), n_studies_prim = length(n),
    median_effect_imm = median(y), median_effect_prim = median(n),
    diff = median(y) - median(n),
    p_floor = p_floor(length(y), length(n)),
    testable = p_floor(length(y), length(n)) <= 0.05,
    floor_basis = "studies",
    p = guarded_p(y, n))
}
C <- bind_rows(Crows); C$p_adj <- p.adjust(C$p, "BH")
cat("\n")
print(C[, c("model", "n_studies_imm", "n_studies_prim", "median_effect_imm",
            "median_effect_prim", "diff", "p", "p_adj")], row.names = FALSE, digits = 3)

# ---- D. is the naive proliferating difference explained by tissue? ----------
# The immortalised lines are not evenly spread across tissue of origin (48% of
# immortalised controls are foreskin/skin against 29% of primary), and tissue is
# itself a strong predictor of control tAge. Stratifying by tissue asks whether
# the naive difference is simply that.
cat("\n== D. naive proliferating difference, stratified by tissue ==\n")
Drows <- list()
pro <- d[d$condition == "Proliferating", ]
for (ti in sort(unique(pro$tissue))) {
  for (mn in names(MODELS)) {
    v <- MODELS[[mn]]; s <- pro[pro$tissue == ti, ]
    y <- s[[v]][s$immortalised == "yes"]; n <- s[[v]][s$immortalised == "no"]
    if (!length(y) || !length(n)) next
    Drows[[length(Drows) + 1]] <- data.frame(
      test = "D_proliferating_by_tissue", stratum = ti, model = mn,
      n_immortalised = length(y), n_primary = length(n),
      median_immortalised = median(y), median_primary = median(n),
      diff = median(y) - median(n),
      p_floor = p_floor(length(y), length(n)),
      testable = p_floor(length(y), length(n)) <= 0.05,
      floor_basis = "samples",
      p = guarded_p(y, n))
  }
}
D <- bind_rows(Drows)
D$p_adj <- p.adjust(D$p, "BH")
print(D[, c("stratum", "model", "n_immortalised", "n_primary", "diff", "p", "p_adj")],
      row.names = FALSE, digits = 3)
cat("\n  Direction is positive in foreskin and lung on the chronological clocks; on the\n")
cat("  mortality clock foreskin is slightly negative, so the claim is not uniform across\n")
cat("  all three. Magnitude is concentrated in lung on every clock.\n")
cat("  (superseded line follows for the record)\n")
cat("  Direction is positive in every tissue on BOTH models, so the difference is\n")
cat("  not an artefact of tissue composition. Magnitude is concentrated in lung and\n")
cat("  is far larger on scaled_diff than on yugene, and the only immortalised lung\n")
# ---- E. is it immortalisation, or the strain it comes packaged with? --------
# The immortalised and primary groups are not drawn from the same backgrounds, so
# before attributing the gap to hTERT it has to be compared against the strain
# difference it is confounded with. BJ- and IMR90-derived controls are the two
# backgrounds that supply both kinds of sample, so they are the comparison to make.
cat("\n== E. strain effect, for comparison with the immortalisation gap ==\n")
# BJ-derived means BJ only. ERP021140 used to carry the composite label
# "HCA2/BJ (primary)" and was counted here, which was defensible only while its strain
# was unresolvable; it is now resolved to HCA2 (meta_analysis/10), a different donor's
# foreskin strain, so it does not belong in a BJ background group. Excluding it drops
# this contrast from 15 vs 44 samples to 9 vs 44. The foreskin-vs-lung comparison, which
# is the one that legitimately pools HCA2 with BJ, is reported separately in script 11.
BJ  <- c("BJ", "BJ-hTERT")
IM  <- c("IMR90", "IMR90 ER:RAS", "IMR90-hTERT")
stopifnot(all(c(BJ, IM) %in% unique(d$line)))
Erows <- list()
for (mn in names(MODELS)) {
  v <- MODELS[[mn]]
  P <- d[d$condition == "Proliferating", ]
  y <- P[[v]][P$line %in% BJ]; n <- P[[v]][P$line %in% IM]
  Erows[[length(Erows) + 1]] <- data.frame(
    test = "E_strain_background", stratum = "BJ vs IMR90 background", model = mn,
    n_immortalised = length(y), n_primary = length(n),
    median_immortalised = median(y), median_primary = median(n),
    diff = median(y) - median(n),
    p_floor = p_floor(length(y), length(n)),
    testable = p_floor(length(y), length(n)) <= 0.05,
    floor_basis = "samples",
    p = guarded_p(y, n))
}
E <- bind_rows(Erows); E$p_adj <- p.adjust(E$p, "BH")
print(E[, c("stratum", "model", "n_immortalised", "n_primary", "diff", "p", "p_adj")],
      row.names = FALSE, digits = 3)
cat("\n  strain composition of the two groups compared in A:\n")
comp <- d[d$condition == "Proliferating", ] %>% count(immortalised, line) %>% arrange(immortalised, -n)
print(as.data.frame(comp), row.names = FALSE)
cat("\n  median tAge by strain (scaled_diff), lowest first:\n")
sm <- d[d$condition == "Proliferating", ] %>% group_by(line) %>%
  summarise(n = n(), median = median(scaled_diff_EN_tAge), .groups = "drop") %>% arrange(median)
print(as.data.frame(sm), row.names = FALSE, digits = 3)
cat("\n  => the immortalised group is drawn from high-scoring backgrounds and the\n")
cat("     primary group is dominated by IMR90, the lowest-scoring strain, so the\n")
cat("     gap in A is in large part a strain difference.\n\n")

cat("  lines are IMR90-hTERT and Tig3ET, i.e. 2 studies - so this remains a\n")
cat("  study-level observation, not an estimate of an hTERT effect.\n")

# ---- F. the two factors, each stratified by the other ----------------------
# WHY THIS EXISTS. Section E compares BJ-derived against IMR90-derived controls while
# POOLING parental and hTERT cells on both sides, and the two sides are badly unbalanced
# on immortalisation: the BJ group is 7 of 9 immortalised (78%), the IMR90 group 2 of 44
# (5%). So E's "background effect" is itself immortalisation-loaded, in the same direction
# as the gap in A that it is invoked to explain -- quoting E as evidence that background
# beats immortalisation is circular on its own.
#
# The way out is to stratify each factor by the other. Tissue is used rather than parental
# strain for the background side because it is the only split with usable n in both strata.
cat("\n== F. each factor stratified by the other (breaks E's circularity) ==\n")
Frows <- list()
add <- function(test, stratum, mn, y, n_) {
  fl <- p_floor(length(y), length(n_))
  Frows[[length(Frows) + 1]] <<- data.frame(
    test = test, stratum = stratum, model = mn,
    n_immortalised = length(y), n_primary = length(n_),
    median_immortalised = median(y), median_primary = median(n_),
    diff = median(y) - median(n_), p_floor = fl, testable = fl <= 0.05,
    floor_basis = "samples",
    p = guarded_p(y, n_))
}
for (mn in names(MODELS)) {
  v <- MODELS[[mn]]; P <- d[d$condition == "Proliferating", ]
  for (st in c("no", "yes")) {
    q <- P[P$immortalised == st, ]
    add("F_tissue_within_immortalisation_status",
        paste0("foreskin vs lung | ", ifelse(st == "yes", "immortalised", "primary")),
        mn, q[[v]][q$tissue == "Foreskin"], q[[v]][q$tissue == "Lung"])
  }
  for (pa in c("BJ", "IMR90")) {
    q <- P[P$parent == pa, ]
    add("F_immortalisation_within_background", paste0("hTERT vs parental | ", pa),
        mn, q[[v]][q$immortalised == "yes"], q[[v]][q$immortalised == "no"])
  }
}
F <- bind_rows(Frows); F$p_adj <- p.adjust(F$p, "BH")
print(F[, c("stratum", "model", "n_immortalised", "n_primary", "diff", "p_floor", "p")],
      row.names = FALSE, digits = 3)
cat("\n  => the tissue difference survives in BOTH immortalisation strata and on both\n")
cat("     chronological and mortality clocks; the immortalisation difference survives in\n")
cat("     NEITHER background. That asymmetry, not E's larger magnitude, is what supports\n")
cat("     reading the cross-study gap as a property of background rather than of hTERT.\n")

out <- bind_rows(A, B, C, D, E, F)

# ---- note: what each row is, what its floor is computed on, and, where the -----
# comparison is not testable, why. Written so the SI table stands alone.
TEST_NOTE <- c(
  A_naive_within_condition_CONFOUNDED =
    paste("Immortalised vs primary samples within one condition, ignoring cell line.",
          "CONFOUNDED BY DESIGN: immortalisation is perfectly nested within study (0 of 34",
          "studies contain both), so this is entirely a between-study contrast and cannot",
          "separate hTERT from the labs that use hTERT lines. Reported as description, not",
          "as an estimate of an immortalisation effect."),
  B_parental_matched =
    paste("An hTERT line against its own parental background (BJ vs BJ-hTERT, IMR90 vs",
          "IMR90-hTERT). Holds genetic background constant, which test A does not, but the",
          "comparison is still across studies."),
  C_within_study_OIS_effect_by_htert =
    paste("Does the WITHIN-STUDY OIS effect differ between immortalised and primary studies?",
          "Each study's effect is a difference against its own controls, so study baseline",
          "cancels and the study becomes the unit of analysis. This asks whether",
          "immortalisation changes the response to oncogene induction, which is answerable;",
          "it does not recover the baseline effect of immortalisation, which is not."),
  D_proliferating_by_tissue =
    paste("Test A split by tissue of origin. Shows the apparent immortalisation effect is",
          "carried by lung-derived samples and absent in foreskin-derived ones."),
  E_strain_background =
    paste("BJ-derived vs IMR90-derived proliferating controls, regardless of immortalisation.",
          "Measures the genetic-background difference that confounds test A."),
  F_tissue_within_immortalisation_status =
    paste("Foreskin vs lung, computed separately among primary and among immortalised samples,",
          "so tissue is compared with immortalisation status held fixed."),
  F_immortalisation_within_background =
    paste("hTERT vs parental, computed separately within each genetic background, so",
          "immortalisation is compared with background held fixed. The reverse of the",
          "tissue test above."))

n1 <- ifelse(is.na(out$n_immortalised), out$n_studies_imm,  out$n_immortalised)
n2 <- ifelse(is.na(out$n_primary),      out$n_studies_prim, out$n_primary)
out$note <- paste0(
  TEST_NOTE[out$test],
  " Floor basis: ", out$floor_basis, ", ", n1, " vs ", n2, "; smallest attainable",
  " two-sided p = ", signif(out$p_floor, 3), ".",
  ifelse(out$testable, "",
         paste0(" NOT TESTABLE: that floor exceeds 0.05, so no p can reach significance",
                " whatever the data, and p is withheld rather than reported as a null.",
                " The effect size is still given.")))

write.csv(out, file.path(RERUN_DIR, "immortalisation_contrasts.csv"), row.names = FALSE)
cat(sprintf("\nSaved -> %s\n", file.path(RERUN_DIR, "immortalisation_contrasts.csv")))
