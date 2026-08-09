# 14a2_export_temporal_bytimepoint_matrices.R
#
# Extends 14a to per-timepoint resolution (3 cell types x 3 timepoints vs
# each cell type's own 'none' baseline, 6v6 each) -- same exact-decomposition
# method, just at the finer temporal granularity used in script 10.

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
    m <- t(exprs(tAge_eset[[variant]]))
    write.csv(data.frame(sample_id = rownames(m), m, check.names = FALSE),
              file.path(OUT_DIR, sprintf("%s_%s.csv", label, variant)), row.names = FALSE)
  }
  groups <- pData(tAge_eset$scaled_diff)[[group_col]]
  write.csv(data.frame(sample_id = rownames(pData(tAge_eset$scaled_diff)), group = groups),
            file.path(OUT_DIR, sprintf("%s_groups.csv", label)), row.names = FALSE)
  cat(sprintf("   exported %d samples x %d genes\n", nrow(m), ncol(m)))
}

recount_pheno <- readRDS(file.path(RERUN_DIR, "temporal_recount_pheno.rds"))
recount_pheno$time_after_treatment[is.na(recount_pheno$time_after_treatment)] <- "none"
recount_pheno$time_after_treatment <- gsub(recount_pheno$time_after_treatment, pattern = " ", replacement = "_")

human_pc <- get_ensembl_release_pc()
erp_download <- download_studies(studies = "ERP021140", sra_organism = "human")

timepoints <- c("4_days", "10_days", "20_days")
for (ct in c("Fibroblast", "Keratinocyte", "Melanocyte")) {
  for (tp in timepoints) {
    samples_ct_tp <- recount_pheno$sample_ID[recount_pheno$cell_type == ct &
                                                recount_pheno$time_after_treatment %in% c("none", tp)]
    processed <- process_rse(
      rse = erp_download, meta_add = recount_pheno, meta_id_col = "sample_ID",
      filter_pc = TRUE, condition_col = "time_after_treatment", ensembl_dictionary = human_pc,
      filter_duplicates = TRUE, sample_filter = samples_ct_tp,
      low_expression_filter = "all", scale_counts = FALSE
    )
    eset_ct_tp <- ExpressionSet(assayData = as.matrix(assay(processed)), phenoData = AnnotatedDataFrame(as.data.frame(colData(processed))))
    export_matrices(eset_ct_tp, "time_after_treatment", "none", paste0(ct, "_", tp))
  }
}
cat(sprintf("\nDone. Exported to %s\n", OUT_DIR))
