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
d$parent <- gsub("-hTERT$", "", d$line)
d$parent[d$parent == "Tig3ET"] <- "TIG3"
d$parent <- gsub(" ER:RAS$", "", d$parent)
MODELS <- c(scaled_diff = "scaled_diff_EN_tAge", yugene_diff = "yugene_diff_EN_tAge")

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
      diff = median(y) - median(n), p = wtest(y, n))
  }
}
A <- bind_rows(rows); A$p_adj <- ave(A$p, A$model, FUN = function(x) p.adjust(x, "BH"))
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
        diff = median(y) - median(n), p = wtest(y, n))
    }
  }
}
B <- bind_rows(Brows)
if (nrow(B)) {
  B$p_adj <- ave(B$p, B$model, FUN = function(x) p.adjust(x, "BH"))
  print(B[, c("stratum", "model", "n_immortalised", "n_primary", "diff", "p", "p_adj")],
        row.names = FALSE, digits = 3)
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
    diff = median(y) - median(n), p = wtest(y, n))
}
C <- bind_rows(Crows); C$p_adj <- p.adjust(C$p, "BH")
cat("\n")
print(C[, c("model", "n_studies_imm", "n_studies_prim", "median_effect_imm",
            "median_effect_prim", "diff", "p", "p_adj")], row.names = FALSE, digits = 3)

out <- bind_rows(A, B, C)
write.csv(out, file.path(RERUN_DIR, "immortalisation_contrasts.csv"), row.names = FALSE)
cat(sprintf("\nSaved -> %s\n", file.path(RERUN_DIR, "immortalisation_contrasts.csv")))
