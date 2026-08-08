# 10_pathway_tage_score_temporal_bytimepoint.R
#
# PILOT / EXPLORATORY. Expands 09_pathway_tage_score_all.R's pooled
# "irradiated vs none" temporal analysis into per-timepoint resolution:
# 3 cell types x 3 timepoints (4/10/20 days vs that cell type's own 'none'),
# 50 Hallmark pathways each = 450 pathway-restricted tAge runs.
#
# POWER CAVEAT (explicit, not hidden): each comparison here is only 6 vs 6
# samples. Wilcoxon's minimum achievable p-value at n=6 vs n=6 is
# 2/choose(12,6) = 2/924 ~= 0.00216 -- you cannot get a more significant
# result than that regardless of true effect size. BH correction is applied
# across all 450 x 2 models = 900 tests, which is a stricter bar than the
# pooled analysis's 300 tests. Treat "ns" here as INCONCLUSIVE at this sample
# size, not as evidence of no effect -- especially for pathways with fewer
# genes reaching the EN model's feature space.

source("R/config.R")
source("R/functions.R")
suppressPackageStartupMessages({
  library(msigdbr)
  library(tAge)
  library(Biobase)
  library(dplyr)
})

SCRNA_DIR <- "/home/ro/APFS_copy/root/Backup/Documents/modules/scrna_seq/GSE226225/final_analysis"
MODEL_DIR <- file.path(SCRNA_DIR, "tAge_models")
Sys.setenv(RETICULATE_PYTHON = file.path(SCRNA_DIR, ".venv/bin/python"))
model_paths <- list(
  scaled_diff = file.path(MODEL_DIR, "EN_Chronoage_Multispecies_Multitissue_scaleddiff.pkl"),
  yugene_diff = file.path(MODEL_DIR, "EN_Chronoage_Multispecies_Multitissue_yugenediff.pkl")
)

get_model_features <- function(pkl_path) {
  py <- file.path(SCRNA_DIR, ".venv/bin/python")
  tmp <- tempfile(fileext = ".txt")
  script_path <- tempfile(fileext = ".py")
  writeLines(sprintf(
    "import joblib, warnings\nwarnings.filterwarnings('ignore')\nm = joblib.load('%s')\nwith open('%s','w') as f:\n    f.write('\\n'.join(m.feature_names_in_))\n",
    pkl_path, tmp
  ), script_path)
  system2(py, script_path)
  readLines(tmp, warn = FALSE)
}
model_features <- get_model_features(model_paths$scaled_diff)

hallmark <- msigdbr(species = "Homo sapiens", collection = "H")
pathway_genes <- split(hallmark$gene_symbol, hallmark$gs_name)
cat(sprintf("EN model feature space: %d genes. Testing %d Hallmark pathways per (cell type x timepoint).\n",
            length(model_features), length(pathway_genes)))

cohens_d <- function(x, y) {
  nx <- length(x); ny <- length(y)
  pooled_sd <- sqrt(((nx - 1) * var(x) + (ny - 1) * var(y)) / (nx + ny - 2))
  (mean(x) - mean(y)) / pooled_sd
}

run_pathway_tage <- function(full_assay, full_pdata, group_col, test_label, control_label,
                              pathway_genes, model_features, model_paths, cell_type, timepoint, pw) {
  genes_pw <- unique(pathway_genes[[pw]])
  genes_in_data <- intersect(genes_pw, rownames(full_assay))

  res_row <- data.frame(
    cell_type = cell_type, timepoint = timepoint, pathway = pw,
    n_pathway_genes = length(genes_pw), n_in_expression_data = length(genes_in_data),
    n_overlap_EN_model_features = 0L, status = NA_character_,
    n_test = NA_integer_, n_control = NA_integer_,
    cohens_d_scaled = NA_real_, wilcox_p_scaled = NA_real_,
    cohens_d_yugene = NA_real_, wilcox_p_yugene = NA_real_,
    stringsAsFactors = FALSE
  )

  if (length(genes_in_data) < 10) {
    res_row$status <- "SKIPPED (<10 genes in data)"
    return(res_row)
  }

  n_model_overlap <- tryCatch({
    tmp_eset <- ExpressionSet(assayData = full_assay[genes_in_data, , drop = FALSE])
    mapped <- tAge:::map_genes(tmp_eset, "human", "Gene.Symbol", verbose = FALSE)
    length(intersect(unique(rownames(mapped)), model_features))
  }, error = function(e) 0L)
  res_row$n_overlap_EN_model_features <- n_model_overlap

  tryCatch({
    sub_assay <- full_assay[genes_in_data, , drop = FALSE]
    eset <- ExpressionSet(assayData = sub_assay, phenoData = AnnotatedDataFrame(full_pdata))

    tAge_eset <- suppressWarnings(tAge_preprocessing(
      eset, species = "human", gene_mapping_type = "Gene.Symbol",
      control_group_column = group_col, control_group_label = control_label,
      verbose = FALSE, count_threshold = 10, percent_threshold = 20
    ))
    pred <- suppressWarnings(predict_tAge(tAge_eset, model_paths, species = "human", mode = "EN"))

    test_rows <- pred[pred[[group_col]] == test_label, ]
    ctrl_rows <- pred[pred[[group_col]] == control_label, ]

    res_row$status <- "OK"
    res_row$n_test <- nrow(test_rows)
    res_row$n_control <- nrow(ctrl_rows)
    res_row$cohens_d_scaled <- cohens_d(test_rows$scaled_diff_EN_tAge, ctrl_rows$scaled_diff_EN_tAge)
    res_row$wilcox_p_scaled <- wilcox.test(test_rows$scaled_diff_EN_tAge, ctrl_rows$scaled_diff_EN_tAge)$p.value
    res_row$cohens_d_yugene <- cohens_d(test_rows$yugene_diff_EN_tAge, ctrl_rows$yugene_diff_EN_tAge)
    res_row$wilcox_p_yugene <- wilcox.test(test_rows$yugene_diff_EN_tAge, ctrl_rows$yugene_diff_EN_tAge)$p.value
  }, error = function(e) {
    res_row$status <<- paste("ERROR:", conditionMessage(e))
  })

  res_row
}

out_csv <- file.path(RERUN_DIR, "pathway_tage_score_temporal_bytimepoint.csv")
all_results <- list()
save_progress <- function() write.csv(do.call(rbind, all_results), out_csv, row.names = FALSE)

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
    full_assay <- as.matrix(assay(processed))
    full_pdata <- as.data.frame(colData(processed))
    cat(sprintf("-- %s, %s vs none (%d samples) --\n", ct, tp, ncol(processed)))

    for (pw in names(pathway_genes)) {
      row <- run_pathway_tage(full_assay, full_pdata, "time_after_treatment", tp, "none",
                               pathway_genes, model_features, model_paths, ct, tp, pw)
      all_results[[paste(ct, tp, pw)]] <- row
    }
    n_ok <- sum(sapply(all_results[grepl(paste0("^", ct, " ", tp, " "), names(all_results))],
                        function(r) r$status == "OK"))
    cat(sprintf("  done: %d/%d OK\n", n_ok, length(pathway_genes)))
    save_progress()
  }
}

final <- do.call(rbind, all_results)
final$p_adj_scaled <- p.adjust(final$wilcox_p_scaled, method = "BH")
final$p_adj_yugene <- p.adjust(final$wilcox_p_yugene, method = "BH")
final$efficiency_score_scaled <- abs(final$cohens_d_scaled) / log2(final$n_overlap_EN_model_features + 2)
final$efficiency_score_yugene <- abs(final$cohens_d_yugene) / log2(final$n_overlap_EN_model_features + 2)
write.csv(final, out_csv, row.names = FALSE)

cat(sprintf("\nDone. %d/%d rows OK. Saved -> %s\n", sum(final$status == "OK", na.rm = TRUE), nrow(final), out_csv))
