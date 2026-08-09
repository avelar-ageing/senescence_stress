# 01_export_full_tage_matrices.R
#
# Step 1 of the partial-tAge pipeline (see PARTIAL_TAGE_METHOD.md). The EN
# tAge model is a linear ElasticNet (pred = intercept + sum(coef_i * z_i)),
# so a pathway's exact contribution to the prediction is sum(coef_i * z_i)
# over just that pathway's genes -- computed from the FULL,
# whole-transcriptome-normalized data (no re-normalization on a restricted
# gene subset, no forced imputation of missing pathway genes).
#
# This script runs tAge_preprocessing ONCE per group on the full gene set,
# and exports the scaled_diff/yugene_diff matrices that get fed to the model
# (samples x genes, mouse-ortholog-ID columns, NaN-padded to the 18,696
# reference list) plus group labels, for 04_partial_tage_decompose.py to
# consume.

source("R/config.R")
source("R/functions.R")
suppressPackageStartupMessages({
  library(tAge)
  library(Biobase)
})

OUT_DIR <- file.path(RERUN_DIR, "partial_tage")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

export_matrices <- function(eset, group_col, control_label, label) {
  cat(sprintf("-- %s --\n", label))
  tAge_eset <- suppressWarnings(tAge_preprocessing(
    eset, species = "human", gene_mapping_type = "Gene.Symbol",
    control_group_column = group_col, control_group_label = control_label,
    verbose = FALSE, count_threshold = 10, percent_threshold = 20
  ))

  for (variant in c("scaled_diff", "yugene_diff")) {
    m <- t(exprs(tAge_eset[[variant]]))  # samples x genes
    write.csv(data.frame(sample_id = rownames(m), m, check.names = FALSE),
              file.path(OUT_DIR, sprintf("%s_%s.csv", label, variant)), row.names = FALSE)
  }
  groups <- pData(tAge_eset$scaled_diff)[[group_col]]
  write.csv(data.frame(sample_id = rownames(pData(tAge_eset$scaled_diff)), group = groups),
            file.path(OUT_DIR, sprintf("%s_groups.csv", label)), row.names = FALSE)
  cat(sprintf("   exported %d samples x %d genes\n", nrow(m), ncol(m)))
}

# ── Meta-analysis: export once (all 6 conditions together) ─────────────────
cat("== Meta-analysis (all 6 conditions, one export) ==\n")
rse <- readRDS(file.path(RERUN_DIR, "cs_cq_all_study_processed.rds"))
eset <- ExpressionSet(assayData = as.matrix(assay(rse)), phenoData = AnnotatedDataFrame(as.data.frame(colData(rse))))
export_matrices(eset, "cell_substate", "Proliferating", "meta")

# ── Temporal: export once per cell type (pooled irradiated vs none) ────────
cat("\n== Temporal (per cell type) ==\n")
recount_pheno <- readRDS(file.path(RERUN_DIR, "temporal_recount_pheno.rds"))
recount_pheno$time_after_treatment[is.na(recount_pheno$time_after_treatment)] <- "none"
recount_pheno$time_after_treatment <- gsub(recount_pheno$time_after_treatment, pattern = " ", replacement = "_")
recount_pheno$irradiated <- ifelse(recount_pheno$time_after_treatment == "none", "none", "irradiated")

human_pc <- get_ensembl_release_pc()
erp_download <- download_studies(studies = "ERP021140", sra_organism = "human")

for (ct in c("Fibroblast", "Keratinocyte", "Melanocyte")) {
  samples_ct <- recount_pheno$sample_ID[recount_pheno$cell_type == ct]
  processed <- process_rse(
    rse = erp_download, meta_add = recount_pheno, meta_id_col = "sample_ID",
    filter_pc = TRUE, condition_col = "irradiated", ensembl_dictionary = human_pc,
    filter_duplicates = TRUE, sample_filter = samples_ct,
    low_expression_filter = "all", scale_counts = FALSE
  )
  eset_ct <- ExpressionSet(assayData = as.matrix(assay(processed)), phenoData = AnnotatedDataFrame(as.data.frame(colData(processed))))
  export_matrices(eset_ct, "irradiated", "none", ct)
}

cat(sprintf("\nDone. Exported to %s\n", OUT_DIR))
