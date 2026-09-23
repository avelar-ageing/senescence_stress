# 27_si_figure1_pca.R
#
# SI Figure 1: PCA of the arrest samples before and after removing the study
# batch effect, with one panel per covariate.
#
# WHY THIS SCRIPT EXISTS. This figure is in the manuscript - it is SI Figure 1,
# and the "top two PCs accounted for 71% of sample variation, four separate
# clusters" sentence in 2.1.1 rests on it - but no script in this pipeline made
# it. The version in the submission was produced in 2024 by study_degs.R, which
# lives on the origin/main and origin/Scripts branches and was never carried into
# restructure-bulk-rna-seq, and by marian_variance.R, which is on no branch at
# all. So the first supplementary figure of the paper was not reproducible from
# the branch the paper is supposed to ship with.
#
# Nothing new is invented here: normalise_counts(), pca_rds() and save_p() are
# already in R/functions.R on this branch. This is the driver that calls them.
#
# ONE SUBSTANTIVE CHANGE FROM THE 2024 VERSION. The tissue panel now uses the
# VERIFIED tissue (meta_analysis/10, TISSUE_VERIFIED), not the hand-curated
# `tissue` column. That column called "Skin" the 6 samples now recorded as HCA2 (it
# (CVCL_E2UW) and the HDF series' own paper (Mitra 2018, Genome Biol) both say
# foreskin. The submitted figure therefore shows a three-level tissue panel whose
# third level is ten mislabelled samples. After correction there are two levels
# among these samples, and the only genuinely adult dermal strain, HDF161, is not
# in this object.
#
# It also writes the per-sample PC coordinates, so the claim that the mislabelled
# samples sit inside the foreskin cluster can be checked rather than eyeballed.
#
# Output: rerun_outputs/si_figure1_panels/si_figure1_pca_<arm>_<covariate>.png
#           -- one file per panel, for assembling the figure by hand. This is what the
#              2024 script did with its _1/_2 split, and for the same reason: a single
#              combined image overlaps its own legends.
#         rerun_outputs/si_figure1_pca_<arm>.png  -- combined contact sheet, for checking
#         rerun_outputs/si_figure1_pca_coordinates.csv

source("R/config.R")
source("R/functions.R")
suppressPackageStartupMessages({ library(DESeq2); library(dplyr); library(ggplot2) })

# The six covariates of the submitted SI Figure 1, in caption order: cell line,
# tissue, cell state, cell sub-state, sequencing platform, SRA accession.
PANELS <- c("cell_line", "tissue", "cell_state", "cell_substate",
            "sra.platform_model", "study")

PANEL_DIR <- file.path(RERUN_DIR, "si_figure1_panels")
if (!dir.exists(PANEL_DIR)) dir.create(PANEL_DIR, recursive = TRUE)

obj <- readRDS(file.path(RERUN_DIR, "cs_cq_all_study_processed_deseq2.rds"))
dds <- if (is.list(obj) && !is.null(obj$deseq_obj)) obj$deseq_obj else obj
cat(sprintf("deseq object: %d genes x %d samples\n", nrow(dds), ncol(dds)))

# ---- the corrected tissue ---------------------------------------------------
ann <- read.csv(file.path(RERUN_DIR, "immortalisation_annotation_corrected.csv"))
i <- match(colnames(dds), ann$external_id)
stopifnot(!any(is.na(i)))
old_tissue <- as.character(colData(dds)$tissue)
colData(dds)$tissue <- ann$tissue_verified[i]
# The cell_line panel must come from the verified annotation for the same reason.
# The deseq object was built before the per-sample HDF strain parse, so its
# colData still pools three strains under "HDF 10-2" and four under "HDF 10-5".
old_line <- as.character(colData(dds)$cell_line)
colData(dds)$cell_line <- ann$cell_line_resolved[i]
cat(sprintf("cell_line relabelled for %d of %d samples in this object\n",
            sum(old_line != colData(dds)$cell_line), ncol(dds)))
chg <- old_tissue != colData(dds)$tissue
cat(sprintf("tissue relabelled for %d of %d samples in this object\n", sum(chg), ncol(dds)))
if (any(chg)) print(table(old = old_tissue[chg], verified = colData(dds)$tissue[chg]))

pcas <- list()
for (arm in c("no_correction", "correction")) {
  cat(sprintf("\n== %s ==\n", arm))
  cnts <- normalise_counts(
    dds_obj_use   = dds,
    batch_col     = if (arm == "correction") "study" else NULL,
    protect_col   = if (arm == "correction") c("cell_substate") else NULL,
    method        = if (arm == "correction") "wgcna" else NULL,
    normalisation = "vst", vst_n = 500, blind = TRUE)
  # ONE FILE PER PANEL. This is how the submitted figure was made: the 2024 script wrote
  # cs_cq_no_batch_correction_1.png and _2.png separately so the panels could be assembled
  # by hand, because a single combined image overlaps its own legends. That is worse now,
  # not better -- resolving the strain labels took the cell_line panel from 7 coarse levels
  # to 16, and the SRA panel carries 34 -- so each covariate gets its own file at a size
  # that fits its legend, and the combined sheet is written only as a contact sheet for
  # reference. Coordinates are identical in every version; this is layout only.
  # Panels are written after the loop by R/pca_panel_figures.R, which pairs the
  # two arms as facet columns (one legend serves both) and stacks the narrow
  # covariates into a single figure. See that file for why the 2024 one-file-
  # per-covariate layout was needed and why it no longer applies to the arms.
  # per-sample coordinates, for checking rather than eyeballing
  m <- assay(cnts)
  v <- apply(m, 1, var); sel <- order(v, decreasing = TRUE)[seq_len(min(500, length(v)))]
  pc <- prcomp(t(m[sel, ]))
  ve <- round(100 * pc$sdev^2 / sum(pc$sdev^2), 1)
  cat(sprintf("  PC1 %.1f%%  PC2 %.1f%%  (PC1+PC2 = %.1f%%)\n", ve[1], ve[2], ve[1] + ve[2]))
  base <- data.frame(
    arm = arm, sample = colnames(m), PC1 = pc$x[, 1], PC2 = pc$x[, 2],
    pc1_var = ve[1], pc2_var = ve[2],
    cell_line = ann$cell_line_resolved[i], tissue_verified = colData(dds)$tissue,
    tissue_metadata = old_tissue, relabelled = chg,
    cell_substate = as.character(colData(dds)$cell_substate))
  # carry every covariate the figure plots, so the panels can be built after the loop
  for (pnl in setdiff(PANELS, names(base)))
    base[[pnl]] <- as.character(colData(cnts)[[pnl]])
  pcas[[arm]] <- base
}
co <- bind_rows(pcas)

source(file.path("R", "pca_panel_figures.R"))
cat("\nSI Figure 1 panels:\n")
save_pca_figures(co, PANELS, "si_figure1_pca", RERUN_DIR)
write.csv(co, file.path(RERUN_DIR, "si_figure1_pca_coordinates.csv"), row.names = FALSE)

# ---- does the relabelled group sit inside the foreskin cluster? -------------
cat("\n== are the relabelled samples inside the foreskin cluster? ==\n")
for (arm in unique(co$arm)) {
  s <- co[co$arm == arm, ]
  cen <- s %>% group_by(tissue_verified) %>%
    summarise(PC1 = mean(PC1), PC2 = mean(PC2), n = n(), .groups = "drop")
  cat(sprintf("\n-- %s --\n", arm)); print(as.data.frame(cen), row.names = FALSE, digits = 3)
  r <- s[s$relabelled, ]
  if (nrow(r)) {
    d <- function(row, t) {
      c1 <- cen[cen$tissue_verified == t, ]
      sqrt((row$PC1 - c1$PC1)^2 + (row$PC2 - c1$PC2)^2)
    }
    r$dist_foreskin <- vapply(seq_len(nrow(r)), function(k) d(r[k, ], "Foreskin"), numeric(1))
    r$dist_lung     <- vapply(seq_len(nrow(r)), function(k) d(r[k, ], "Lung"), numeric(1))
    r$nearer <- ifelse(r$dist_foreskin < r$dist_lung, "Foreskin", "Lung")
    print(r[, c("sample", "cell_line", "PC1", "PC2", "dist_foreskin", "dist_lung", "nearer")],
          row.names = FALSE, digits = 3)
    cat(sprintf("  %d of %d relabelled samples are nearer the foreskin centroid\n",
                sum(r$nearer == "Foreskin"), nrow(r)))
  }
}
cat(sprintf("\nSaved -> %s\n", file.path(RERUN_DIR, "si_figure1_pca_coordinates.csv")))
