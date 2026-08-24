# 04_tage_by_celltype.R
#
# tAge (transcriptomic age, Gladyshev-Lab/tAge EN elastic-net clock) for the
# ERP021140 temporal time course, per cell type (Fibroblast, Keratinocyte,
# Melanocyte) and per timepoint (none/4/10/20 days post-irradiation).
#
# NOTE ON SCALING (corrected -- see DISCREPANCY_REPORT for the walk-back):
# an earlier version of this script used scale_counts=TRUE here on the
# unverified assumption that the tAge EN clock was trained on recount3
# AUC-scaled input specifically. Checked the tAge package docs/README
# directly: there is no such requirement anywhere in them -- the package's
# own RLE_normalization() (edgeR-based) is its internal library-size/
# composition-bias correction step, applied regardless of what convention
# upstream counts used. Given that, and that this is a single-study dataset
# (same protocol throughout, so within-study normalization is what matters,
# same reasoning as 02_run_time_analysis.R's scale_counts=FALSE), this now
# matches the DEG-calling pipeline's own convention: scale_counts=FALSE.
#
# Each cell type is preprocessed (and control-subtracted against its own
# 'none' baseline) SEPARATELY, not pooled -- tAge_preprocessing's
# control_subtraction takes one global control group per call, and pooling
# would subtract a cross-cell-type average baseline against each sample,
# which is wrong when baseline expression differs hugely by cell type.

source("R/config.R")
source("R/functions.R")

suppressPackageStartupMessages({
  library(SummarizedExperiment)
  library(tAge)
  library(Biobase)
  library(dplyr)
  library(ggplot2)
})

SCRNA_DIR <- "/home/ro/APFS_copy/root/Backup/Documents/modules/scrna_seq/GSE226225/final_analysis"
MODEL_DIR <- file.path(SCRNA_DIR, "tAge_models")
Sys.setenv(RETICULATE_PYTHON = file.path(SCRNA_DIR, ".venv/bin/python"))
model_paths <- list(
  scaled_diff = file.path(MODEL_DIR, "EN_Chronoage_Multispecies_Multitissue_scaleddiff.pkl"),
  yugene_diff = file.path(MODEL_DIR, "EN_Chronoage_Multispecies_Multitissue_yugenediff.pkl")
)

recount_pheno <- readRDS(file.path(RERUN_DIR, "temporal_recount_pheno.rds"))
recount_pheno$time_after_treatment[is.na(recount_pheno$time_after_treatment)] <- "none"
# time_analysis() does this same space->underscore normalization internally
# (general_functions.R ~line 2590); replicated here since we're not calling
# time_analysis() itself for this script.
recount_pheno$time_after_treatment <- gsub(recount_pheno$time_after_treatment, pattern = " ", replacement = "_")

human_pc <- get_ensembl_release_pc()

cat("== Downloading ERP021140 (single study, all 3 cell types) ==\n")
erp_download <- download_studies(studies = "ERP021140", sra_organism = "human")
cat(sprintf("  Downloaded RSE: %d genes x %d samples\n", nrow(erp_download), ncol(erp_download)))

cell_types <- c("Fibroblast", "Keratinocyte", "Melanocyte")

run_one_celltype <- function(ct) {
  cat(sprintf("\n== %s ==\n", ct))
  samples_ct <- recount_pheno$sample_ID[recount_pheno$cell_type == ct]

  processed <- process_rse(
    rse = erp_download,
    meta_add = recount_pheno,
    meta_id_col = "sample_ID",
    filter_pc = TRUE,
    condition_col = "time_after_treatment",
    ensembl_dictionary = human_pc,
    filter_duplicates = TRUE,
    sample_filter = samples_ct,
    low_expression_filter = "all",
    scale_counts = FALSE  # see header note -- matches the DEG-calling pipeline's convention
  )
  cat(sprintf("  Processed: %d genes x %d samples\n", nrow(processed), ncol(processed)))

  eset <- ExpressionSet(
    assayData = as.matrix(assay(processed)),
    phenoData = AnnotatedDataFrame(as.data.frame(colData(processed)))
  )

  tAge_eset <- tAge_preprocessing(
    eset, species = "human", gene_mapping_type = "Gene.Symbol",
    control_group_column = "time_after_treatment", control_group_label = "none",
    verbose = TRUE, count_threshold = 10, percent_threshold = 20
  )
  result <- predict_tAge(tAge_eset, model_paths, species = "human", mode = "EN")
  result$cell_type <- ct
  result
}

tage_temporal <- do.call(rbind, lapply(cell_types, run_one_celltype))
tage_temporal$time_after_treatment <- factor(tage_temporal$time_after_treatment,
                                              levels = c("none", "4_days", "10_days", "20_days"))
write.csv(tage_temporal, file.path(RERUN_DIR, "tage_temporal_by_celltype.csv"), row.names = FALSE)

cat("\n== tAge spread by cell type x timepoint ==\n")
spread_stats <- tage_temporal %>%
  group_by(cell_type, time_after_treatment) %>%
  summarise(
    n = n(),
    scaled_diff_median = median(scaled_diff_EN_tAge), scaled_diff_sd = sd(scaled_diff_EN_tAge),
    yugene_diff_median = median(yugene_diff_EN_tAge), yugene_diff_sd = sd(yugene_diff_EN_tAge),
    .groups = "drop"
  )
print(spread_stats)
write.csv(spread_stats, file.path(RERUN_DIR, "tage_temporal_spread.csv"), row.names = FALSE)

# Significance testing lives in 07_tage_temporal_pairwise.R, which computes
# the canonical all-pairs family: every timepoint vs every other, per cell
# type, per model, BH-adjusted across 3 x 6 x 2 = 36 tests. That family
# subsumes the vs-baseline comparisons (none vs 4/10/20 days are 3 of the 6
# pairs) and additionally covers the between-timepoint contrasts the
# trajectory claims rest on, so it replaces the narrower vs-baseline-only
# family this script used to write out.

# Violin/boxplot: tAge trajectory over time, faceted by cell type
plot_df <- tage_temporal %>%
  dplyr::select(cell_type, time_after_treatment, scaled_diff_EN_tAge, yugene_diff_EN_tAge) %>%
  tidyr::pivot_longer(cols = c(scaled_diff_EN_tAge, yugene_diff_EN_tAge), names_to = "model", values_to = "tAge")
p <- ggplot(plot_df, aes(x = time_after_treatment, y = tAge, fill = time_after_treatment)) +
  geom_violin(alpha = 0.6, trim = FALSE) +
  geom_boxplot(width = 0.15, outlier.shape = NA, fill = "white") +
  facet_grid(model ~ cell_type, scales = "free_y") +
  theme_bw() + theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = NULL, y = "tAge", title = "Transcriptomic age over time post-irradiation, by cell type (ERP021140)")
ggsave(file.path(RERUN_DIR, "tage_temporal_by_celltype.png"), p, width = 10, height = 7, dpi = 300)

cat(sprintf("\nDone. Outputs -> %s\n", RERUN_DIR))
