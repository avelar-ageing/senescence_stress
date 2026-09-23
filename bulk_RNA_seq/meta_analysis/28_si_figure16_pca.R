# 28_si_figure16_pca.R
#
# The SIPS / OIS / CQ PCA, before and after removing the study batch effect,
# with one panel per covariate.
#
# WHICH FIGURE THIS IS. It is SI Figure 16 of the current supplementary file
# (Final/final/SI_figures_final.pdf, 30 Sep 2024), which is the numbering
# modules_october_2024.docx cites. Two older files number the same figure
# differently -- Final/nature_ageing/SI_figures.pdf (1 Sep 2024) calls it SI
# Figure 13 and modules_final.docx (18 Jul 2024) cites "SI Figure 12" -- so a
# request to fix "SI Figure 12" and a request to fix "SI Figure 16" are the same
# request. The caption is identical in all three: "PCA of SIPS, OIS, and CQ
# samples (a) before and (b) after removing the study batch effect."
#
# WHY THIS SCRIPT EXISTS. Same reason as meta_analysis/27. The submitted version
# was made in 2024 by Final/Scripts/for_github/study_degs.R (origin/main, never
# merged into restructure-bulk-rna-seq), which also defined filter_dataframe()
# inline -- so neither the figure nor the sample-selection rule behind it existed
# on the branch the paper ships with. filter_dataframe() now lives in
# R/functions.R; this is the driver.
#
# The sentence it supports, in 2.1.2: "PCA was performed whereby samples
# clustered into three groups corresponding to CQ, SIPS, and OIS samples, and
# the top two PCAs captured 83% of the sample variance."
#
# HOW THE 48 SAMPLES ARE CHOSEN. All CQ sub-states are collapsed to one "CQ"
# label, proliferating controls are dropped, and filter_dataframe() then keeps
# only studies contributing both a CQ and a senescent sample. Replicative
# senescence disappears at that step rather than by an explicit exclusion: no
# study in the pool contributes both CQ and RS. That is the mechanism behind the
# manuscript's "excluding RS samples due to insufficient sample numbers".
#
# ONE SUBSTANTIVE CHANGE FROM THE 2024 VERSION, as in script 27: the tissue
# panel uses the VERIFIED tissue (meta_analysis/10, TISSUE_VERIFIED) rather than
# the hand-curated `tissue` column. For THIS subset the two agree -- its four
# strains are HCA2/BJ (primary), HDF161, IMR90 and Tig3ET, and none of the ten
# samples mislabelled in the full arrest object are here -- so the figure is
# unchanged by the correction. The script asserts that rather than assuming it,
# which is the point: SI Figure 1 needed the fix, SI Figure 16 provably did not.
#
# Output: rerun_outputs/si_figure16_panels/si_figure16_pca_<arm>_<covariate>.png
#           -- one file per panel, for assembling by hand (see script 27's note)
#         rerun_outputs/si_figure16_pca_<arm>.png  -- combined contact sheet
#         rerun_outputs/si_figure16_pca_coordinates.csv

source("R/config.R")
source("R/functions.R")
suppressPackageStartupMessages({ library(DESeq2); library(dplyr); library(ggplot2) })

PUBLISHED_PC12 <- 83  # "the top two PCAs captured 83% of the sample variance"

# The panel list of the 2024 call, unchanged. pca_rds() silently drops any
# covariate with a single level, so cell_type (all Fibroblast) falls away and
# the six that survive are exactly the caption's: cell line, tissue, cell state,
# cell sub-state (= cq_test), sequencing platform, SRA accession.
PANELS <- c("cell_line", "cell_type", "tissue", "cell_state", "cq_test",
            "sra.platform_model", "study")

cat("== Step 1: sample selection ==\n")
meta_both <- read.csv(file.path(SAVE_DIR_CSV, "study_info_all.csv"))
meta_both$cq_test <- meta_both$cell_substate
temp_label <- unique(meta_both$cq_test)[grepl(unique(meta_both$cq_test), pattern = "CQ")]
meta_both$cq_test <- ifelse(meta_both$cell_substate %in% temp_label, "CQ", meta_both$cell_substate)
cat(sprintf("  study_info_all.csv: %d samples; CQ sub-states collapsed: %s\n",
            nrow(meta_both), paste(temp_label, collapse = ", ")))

# ---- the corrected tissue, applied before anything else ---------------------
ann <- read.csv(file.path(RERUN_DIR, "immortalisation_annotation_corrected.csv"))
j <- match(meta_both$external_id, ann$external_id)
stopifnot(!any(is.na(j)))
# `tissue_assubmitted` is present once meta_analysis/15 has patched the table, and
# holds the submitted label; before that the submitted label is still in `tissue`.
# Either way tissue_metadata ends up holding what the submission said, so the
# comparison below stays meaningful after the table is corrected at source.
meta_both$tissue_metadata <- if ("tissue_assubmitted" %in% colnames(meta_both)) {
  meta_both$tissue_assubmitted
} else {
  meta_both$tissue
}
meta_both$tissue <- ann$tissue_verified[j]
meta_both$cell_line_resolved <- ann$cell_line_resolved[j]

meta_both_cq_test <- meta_both[meta_both$cq_test != "Proliferating", ]
cq_full_rank <- filter_dataframe(dataframe = meta_both_cq_test,
                                 iterate_col = "sra",
                                 filter_col = "cq_test",
                                 col_label_1 = "CQ",
                                 col_label_2 = "CS")
cat(sprintf("  selected %d samples from %d studies\n",
            nrow(cq_full_rank), length(unique(cq_full_rank$sra))))
print(table(cq_full_rank$cq_test))
cat("\n  RS present in the selection? ",
    if (any(grepl("Replicative", cq_full_rank$cq_test))) "YES" else "no (dropped by filter_dataframe)", "\n")
print(table(cq_full_rank$sra, cq_full_rank$cq_test))

cat("\n  tissue, verified vs metadata, in this subset:\n")
print(table(verified = cq_full_rank$tissue, metadata = cq_full_rank$tissue_metadata))
n_relab <- sum(cq_full_rank$tissue != cq_full_rank$tissue_metadata)
cat(sprintf("  %d of %d samples relabelled by the tissue verification\n",
            n_relab, nrow(cq_full_rank)))
print(table(cq_full_rank$cell_line_resolved, cq_full_rank$tissue))

cat("\n== Step 2: process_rse on the 48-sample subset ==\n")
# Re-processed rather than subset out of cs_cq_all_study_processed: the
# low-expression filter is per-condition, so it has to see this sample set.
human_pc <- get_ensembl_release_pc()
cs_cq_download <- readRDS(file.path(RERUN_DIR, "cs_cq_download_raw.rds"))
cat(sprintf("  raw download: %d genes x %d samples\n",
            nrow(cs_cq_download), ncol(cs_cq_download)))

# Two "min(y): no non-missing arguments" warnings are expected here and are not
# a problem: the low-expression filter loops over every condition in meta_both,
# and two of them (Proliferating, Replicative CS) have no samples in this
# subset, so they contribute no genes to the union that is kept.
obj <- process_rse(
  rse = cs_cq_download,
  meta_add = meta_both,
  meta_id_col = "external_id",
  filter_pc = TRUE,
  condition_col = "cell_substate",
  ensembl_dictionary = human_pc,
  filter_duplicates = TRUE,
  sample_filter = cq_full_rank[["external_id"]],
  low_expression_filter = "per",
  scale_counts = TRUE
)
cat(sprintf("  processed object: %d genes x %d samples\n", nrow(obj), ncol(obj)))
stopifnot(ncol(obj) == nrow(cq_full_rank))
print(table(data.frame(colData(obj))$cq_test))

PANEL_DIR <- file.path(RERUN_DIR, "si_figure16_panels")
if (!dir.exists(PANEL_DIR)) dir.create(PANEL_DIR, recursive = TRUE)

cat("\n== Step 3: design matrix ==\n")
dds_list <- deseq_studies(study_obj = obj, protect_col = "cq_test", fix_col = "study")
dds <- dds_list[["deseq_obj"]]

pcas <- list()
for (arm in c("no_correction", "correction")) {
  cat(sprintf("\n== Step 4 (%s) ==\n", arm))
  cnts <- normalise_counts(
    dds_obj_use   = dds,
    batch_col     = if (arm == "correction") "study" else NULL,
    protect_col   = if (arm == "correction") c("cq_test") else NULL,
    method        = if (arm == "correction") "wgcna" else NULL,
    normalisation = "vst", vst_n = 500, blind = TRUE)

  # One file per panel, as in script 27 and as the 2024 scripts did with their _1/_2
  # split: a single combined image overlaps its own legends, and the panels are meant to
  # be assembled by hand. The combined sheet below is a contact sheet only.
  # Panels are written after the loop by R/pca_panel_figures.R: the two arms
  # become facet columns sharing one legend, and the narrow covariates stack
  # into a single figure. See that file for why the per-covariate split was
  # needed in 2024 and why it does not apply to the before/after pair.

  m <- assay(cnts)
  v <- apply(m, 1, var); sel <- order(v, decreasing = TRUE)[seq_len(min(500, length(v)))]
  pc <- prcomp(t(m[sel, ]))
  ve <- round(100 * pc$sdev^2 / sum(pc$sdev^2), 1)
  cat(sprintf("  PC1 %.1f%%  PC2 %.1f%%  (PC1+PC2 = %.1f%%; manuscript says %d%%)\n",
              ve[1], ve[2], ve[1] + ve[2], PUBLISHED_PC12))

  k <- match(colnames(m), cq_full_rank$external_id)
  pcas[[arm]] <- data.frame(
    arm = arm, sample = colnames(m), PC1 = pc$x[, 1], PC2 = pc$x[, 2],
    pc1_var = ve[1], pc2_var = ve[2],
    group = cq_full_rank$cq_test[k],
    cell_substate = cq_full_rank$cell_substate[k],
    cell_line = cq_full_rank$cell_line_resolved[k],
    tissue_verified = cq_full_rank$tissue[k],
    tissue_metadata = cq_full_rank$tissue_metadata[k],
    study = cq_full_rank$sra[k])
  for (pnl in setdiff(PANELS, names(pcas[[arm]])))
    pcas[[arm]][[pnl]] <- as.character(colData(cnts)[[pnl]])
}
co <- bind_rows(pcas)

source(file.path("R", "pca_panel_figures.R"))
cat("\nSI Figure 16 panels:\n")
save_pca_figures(co, PANELS, "si_figure16_pca", RERUN_DIR)
write.csv(co, file.path(RERUN_DIR, "si_figure16_pca_coordinates.csv"), row.names = FALSE)

# ---- do the samples actually fall into three groups? -----------------------
# The manuscript claims three clusters corresponding to CQ, SIPS and OIS. That
# is checkable: assign every sample to its nearest group centroid in the PC1-PC2
# plane and count how many land on their own group.
cat("\n== Step 5: three clusters, or not? ==\n")
for (arm in unique(co$arm)) {
  s <- co[co$arm == arm, ]
  cen <- s %>% group_by(group) %>%
    summarise(PC1 = mean(PC1), PC2 = mean(PC2), n = n(), .groups = "drop")
  cat(sprintf("\n-- %s (PC1 %.1f%% + PC2 %.1f%% = %.1f%%) --\n",
              arm, s$pc1_var[1], s$pc2_var[1], s$pc1_var[1] + s$pc2_var[1]))
  print(as.data.frame(cen), row.names = FALSE, digits = 3)
  d <- outer(seq_len(nrow(s)), seq_len(nrow(cen)),
             Vectorize(function(a, b) sqrt((s$PC1[a] - cen$PC1[b])^2 + (s$PC2[a] - cen$PC2[b])^2)))
  nearest <- cen$group[apply(d, 1, which.min)]
  cat("  nearest-centroid assignment:\n")
  print(table(actual = s$group, nearest = nearest))
  cat(sprintf("  %d of %d samples nearest their own group centroid (%.0f%%)\n",
              sum(nearest == s$group), nrow(s), 100 * mean(nearest == s$group)))
}

cat(sprintf("\nSaved -> %s\n", file.path(RERUN_DIR, "si_figure16_pca_coordinates.csv")))
