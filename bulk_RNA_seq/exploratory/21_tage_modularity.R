# 21_tage_modularity.R
#
# Is transcriptomic age MODULAR - i.e. does the composition of the tAge shift
# depend on cell type and timepoint rather than being one shared programme?
#
# WHY THIS SCRIPT EXISTS. The modularity argument was computed ad hoc in
# conversation and never scripted, so it was neither reproducible nor checkable.
# This reproduces it, and re-tests each claim.
#
# WHAT IS TESTED, and what is deliberately NOT.
#   1  Profile similarity: correlate each group's 50-set contribution profile with
#      every other group's. If tAge composition is cell-type-specific, pairs within
#      a cell type should agree more than pairs across cell types. Tested against a
#      permutation null that shuffles the cell-type labels of the groups.
#   2  Dimensionality: PCA on the profile matrix. One shared programme would put
#      nearly all variance on PC1.
#   3  Breadth: how many sets contribute in all three cell types.
#
#   DROPPED - the fourth claim, that groups with near-identical AGGREGATE tAge have
#   unrelated compositions (illustrated by keratinocyte 4 days +19.0 versus
#   melanocyte 20 days +19.5). That compares aggregate tAge ACROSS preprocessing
#   runs, which TAGE_CALCULATION_AUDIT.md section 5 forbids: gene retention,
#   measured coefficient weight and input scale all differ by run. It is not
#   recoverable and is not reported.
#
# WHY CORRELATIONS ARE SAFE ACROSS RUNS where magnitudes are not: Spearman
# correlation is invariant to the per-run scale factor that makes magnitudes
# incomparable. Everything here is rank-based for that reason.
#
# The temporal arm (3 cell types x 3 timepoints, one study) is primary. The
# meta-analysis conditions are reported alongside, using the WITHIN-STUDY
# contributions from exploratory/18, since its pooled contrasts are confounded.
#
# Output: rerun_outputs/tage_modularity.csv

source("R/config.R")
suppressPackageStartupMessages({ library(dplyr) })
set.seed(1)
NPERM <- 10000

df <- read.csv(file.path(RERUN_DIR, "partial_tage_ALL.csv"), check.names = FALSE)
tp <- df %>% filter(analysis == "temporal_bytimepoint")
res <- list()

spear <- function(M) suppressWarnings(cor(M, method = "spearman", use = "pairwise.complete.obs"))

for (mdl in c("scaled", "yugene")) {
  s <- tp %>% filter(model == mdl) %>%
    mutate(group = paste(cell_type, timepoint, sep = "_"))
  M <- reshape(s[, c("pathway", "group", "contrib_diff")], idvar = "pathway",
               timevar = "group", direction = "wide")
  rownames(M) <- M$pathway; M$pathway <- NULL
  colnames(M) <- gsub("^contrib_diff\\.", "", colnames(M))
  M <- as.matrix(M)
  ct <- sub("_.*$", "", colnames(M))

  # ---- 1. within- vs between-cell-type profile similarity ------------------
  R <- spear(M)
  ut <- which(upper.tri(R), arr.ind = TRUE)
  same <- ct[ut[, 1]] == ct[ut[, 2]]
  wi <- mean(R[upper.tri(R)][same]); bw <- mean(R[upper.tri(R)][!same])
  # permutation null: shuffle which cell type each group belongs to
  gap_null <- replicate(NPERM, {
    p <- sample(ct)
    sm <- p[ut[, 1]] == p[ut[, 2]]
    if (all(sm) || !any(sm)) return(NA_real_)
    mean(R[upper.tri(R)][sm]) - mean(R[upper.tri(R)][!sm])
  })
  gap_null <- gap_null[!is.na(gap_null)]
  p_gap <- (1 + sum(abs(gap_null) >= abs(wi - bw))) / (length(gap_null) + 1)
  res[[paste(mdl, "similarity")]] <- data.frame(
    test = "profile_similarity", model = mdl, n_groups = ncol(M),
    within_celltype_rho = wi, between_celltype_rho = bw, gap = wi - bw,
    p_perm = p_gap)

  # also within-cell-type pairs listed per cell type
  for (c1 in unique(ct)) {
    idx <- which(ct == c1)
    sub <- R[idx, idx]
    res[[paste(mdl, c1, "within")]] <- data.frame(
      test = "within_one_celltype", model = mdl, cell_type = c1,
      mean_rho = mean(sub[upper.tri(sub)]))
  }

  # ---- 2. dimensionality --------------------------------------------------
  pca <- prcomp(t(scale(M)), center = TRUE, scale. = FALSE)
  ve <- pca$sdev^2 / sum(pca$sdev^2)
  res[[paste(mdl, "pca")]] <- data.frame(
    test = "pca_variance", model = mdl,
    pc1 = ve[1], pc1_2 = sum(ve[1:2]), pc1_3 = sum(ve[1:3]),
    n_pc_for_80 = which(cumsum(ve) >= 0.8)[1])

  # ---- 2b. do the later PCs separate cell type? ---------------------------
  # PC1 carries 60-68%, which is the SHARED shift: every group moves the same
  # way on it. The modularity question is whether the residual structure is
  # cell-type organised, so report PC2/PC3 scores by cell type.
  sc <- as.data.frame(pca$x[, 1:3, drop = FALSE])
  sc$cell_type <- sub("_.*$", "", rownames(sc))
  agg <- sc %>% group_by(cell_type) %>%
    summarise(PC1 = mean(PC1), PC2 = mean(PC2), PC3 = mean(PC3), .groups = "drop")
  res[[paste(mdl, "pcscores")]] <- data.frame(
    test = "pc_scores_by_celltype", model = mdl, cell_type = agg$cell_type,
    PC1_mean = agg$PC1, PC2_mean = agg$PC2, PC3_mean = agg$PC3)
  # between-cell-type separation as a share of total spread on each PC
  for (k in 1:3) {
    v <- sc[[paste0("PC", k)]]
    res[[paste(mdl, "pcsep", k)]] <- data.frame(
      test = "pc_celltype_separation", model = mdl, PC = k,
      between_ss_frac = summary(aov(v ~ sc$cell_type))[[1]][["Sum Sq"]][1] /
                        sum(summary(aov(v ~ sc$cell_type))[[1]][["Sum Sq"]]))
  }

  # ---- 3. breadth ---------------------------------------------------------
  # Matches the claim in 2.1.5.4, which is about the POOLED temporal comparison
  # (irradiated vs none per cell type), not about per-timepoint tests. An earlier
  # version of this script counted a set as present in a cell type if it reached
  # significance at ANY timepoint, which is a far weaker criterion and gave 24-26
  # of 50 sets in all three cell types - not comparable to the published claim.
  pooled <- df %>% filter(analysis == "temporal_pooled", model == mdl, p_adj < 0.05) %>%
    distinct(pathway, cell_type) %>% count(pathway)
  res[[paste(mdl, "breadth")]] <- data.frame(
    test = "breadth_pooled", model = mdl,
    n_sets_in_3_celltypes = sum(pooled$n == 3),
    n_sets_in_2 = sum(pooled$n == 2), n_sets_in_1 = sum(pooled$n == 1),
    sets_in_3 = paste(pooled$pathway[pooled$n == 3], collapse = "; "))
}

# ---- meta-analysis conditions, within-study contributions -----------------
w <- read.csv(file.path(RERUN_DIR, "partial_tage_within_study.csv"), check.names = FALSE)
for (mdl in c("scaled", "yugene")) {
  s <- w %>% filter(model == mdl)
  M <- reshape(s[, c("pathway", "label", "contrib_diff_within_study")],
               idvar = "pathway", timevar = "label", direction = "wide")
  rownames(M) <- M$pathway; M$pathway <- NULL
  M <- as.matrix(M); colnames(M) <- gsub("^contrib_diff_within_study\\.", "", colnames(M))
  R <- spear(M)
  pca <- prcomp(t(scale(M)), center = TRUE, scale. = FALSE)
  ve <- pca$sdev^2 / sum(pca$sdev^2)
  res[[paste(mdl, "meta")]] <- data.frame(
    test = "meta_conditions", model = mdl, n_groups = ncol(M),
    mean_rho_between_conditions = mean(R[upper.tri(R)]),
    pc1 = ve[1], pc1_3 = sum(ve[1:3]), n_pc_for_80 = which(cumsum(ve) >= 0.8)[1])
  cat(sprintf("\n-- meta condition profile correlations (%s) --\n", mdl))
  print(round(R, 2))
}

# ---- 3b. breadth under the reporting rule actually used in the text --------
# The claim in 2.1.5.4 - "COMPLEMENT is the only set significant in all three cell
# types" - is about INTERPRETABLE sets significant on BOTH models, not about all 50
# on either model. Checked here explicitly, because the looser counts above (11-14
# of 50) would otherwise look like a contradiction.
rep <- read.csv(file.path(RERUN_DIR, "pathway_representation.csv"), check.names = FALSE)
interp <- rep$pathway[rep$tier == "INTERPRETABLE"]
pooled_all <- df %>% filter(analysis == "temporal_pooled", p_adj < 0.05) %>%
  distinct(pathway, cell_type, model) %>%
  count(pathway, model) %>% filter(n == 3) %>%
  count(pathway) %>% filter(n == 2)          # all 3 cell types, on both models
both3 <- sort(pooled_all$pathway)
res[["breadth_rule"]] <- data.frame(
  test = "breadth_reporting_rule", model = "both",
  n_all_sets = length(both3),
  n_interpretable = sum(both3 %in% interp),
  sets_in_3 = paste(intersect(both3, interp), collapse = "; "))
cat("\n=== 3b. all 3 cell types, BOTH models ===\n")
cat(sprintf("  any set: %d  -> %s\n", length(both3), paste(both3, collapse = "; ")))
cat(sprintf("  of these, INTERPRETABLE: %d -> %s\n", sum(both3 %in% interp),
            paste(intersect(both3, interp), collapse = "; ")))

out <- bind_rows(res)
write.csv(out, file.path(RERUN_DIR, "tage_modularity.csv"), row.names = FALSE)
cat("\n\n=== 1. profile similarity (temporal) ===\n")
print(out %>% filter(test == "profile_similarity") %>%
        select(model, n_groups, within_celltype_rho, between_celltype_rho, gap, p_perm),
      row.names = FALSE, digits = 3)
cat("\n=== within a single cell type ===\n")
print(out %>% filter(test == "within_one_celltype") %>%
        select(model, cell_type, mean_rho), row.names = FALSE, digits = 3)
cat("\n=== 2. dimensionality (temporal) ===\n")
print(out %>% filter(test == "pca_variance") %>%
        select(model, pc1, pc1_2, pc1_3, n_pc_for_80), row.names = FALSE, digits = 3)
cat("\n=== 2b. PC scores by cell type ===\n")
print(out %>% filter(test == "pc_scores_by_celltype") %>%
        select(model, cell_type, PC1_mean, PC2_mean, PC3_mean), row.names = FALSE, digits = 3)
cat("\n=== 2c. how much of each PC's spread is between cell types? ===\n")
print(out %>% filter(test == "pc_celltype_separation") %>%
        select(model, PC, between_ss_frac), row.names = FALSE, digits = 3)
cat("\n=== 3. breadth, pooled per cell type ===\n")
print(out %>% filter(test == "breadth_pooled") %>%
        select(model, n_sets_in_3_celltypes, n_sets_in_2, n_sets_in_1),
      row.names = FALSE, digits = 3)
cat("\nsets significant in all three cell types:\n")
for (i in which(out$test == "breadth_pooled")) {
  cat(sprintf("  %s: %s\n", out$model[i], out$sets_in_3[i]))
}
cat("\n=== meta-analysis conditions ===\n")
print(out %>% filter(test == "meta_conditions") %>%
        select(model, n_groups, mean_rho_between_conditions, pc1, pc1_3, n_pc_for_80),
      row.names = FALSE, digits = 3)
cat(sprintf("\nSaved -> %s\n", file.path(RERUN_DIR, "tage_modularity.csv")))
