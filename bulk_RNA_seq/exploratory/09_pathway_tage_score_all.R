# 09_pathway_tage_score_all.R
#
# PILOT / EXPLORATORY (expands 08_pathway_restricted_tage_pilot.R to every
# meta-analysis condition AND the temporal data by cell type).
#
# For each (condition or cell-type) vs its control, loop over every MSigDB
# Hallmark pathway, restrict the expression matrix to just that pathway's
# genes, rerun tAge end to end, and score pathways by a size-adjusted effect:
#
#   cohens_d          = (mean_test - mean_control) / pooled_SD   [standardized
#                        effect, comparable across the two EN models and
#                        across conditions/cell types with different tAge scales]
#   genes_used        = n genes from that pathway actually in the EN model's
#                        10,487-feature mouse-ortholog space (the real driver
#                        of signal -- NOT the pathway's nominal size, most of
#                        which never reaches the model at all, per the SIPS
#                        pilot: even the biggest Hallmark set only ever
#                        touches ~1.5% of the model's feature space)
#   efficiency_score  = |cohens_d| / log2(genes_used + 2)
#                        -- "aging signal packed per unit of gene-set
#                        complexity"; log2 (not linear) so a pathway twice
#                        the size isn't automatically penalized 2x, only
#                        diminishing-returns-penalized.
#
# Meta-analysis: each of SSCQ/CICQ/RS/SIPS/OIS vs Proliferating (fibroblast,
# recount3 pooled).
# Temporal: each cell type's pooled irradiated samples (4+10+20 days) vs its
# own 'none' baseline -- the natural analog of "condition vs control" for a
# single-study time course, and gives 18 vs 6 samples per cell type (more
# power than any single timepoint alone).

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
cat(sprintf("EN model feature space: %d genes\n", length(model_features)))

hallmark <- msigdbr(species = "Homo sapiens", collection = "H")
pathway_genes <- split(hallmark$gene_symbol, hallmark$gs_name)
cat(sprintf("Testing %d Hallmark pathways per group\n", length(pathway_genes)))

cohens_d <- function(x, y) {
  nx <- length(x); ny <- length(y)
  pooled_sd <- sqrt(((nx - 1) * var(x) + (ny - 1) * var(y)) / (nx + ny - 2))
  (mean(x) - mean(y)) / pooled_sd
}

# Core runner: assay/pdata already subsetted to exactly the test+control rows.
run_pathway_tage <- function(full_assay, full_pdata, group_col, test_label, control_label,
                              pathway_genes, model_features, model_paths, analysis, group_name, pw) {
  genes_pw <- unique(pathway_genes[[pw]])
  genes_in_data <- intersect(genes_pw, rownames(full_assay))

  res_row <- data.frame(
    analysis = analysis, group = group_name, pathway = pw,
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

out_csv <- file.path(RERUN_DIR, "pathway_tage_score_all.csv")
all_results <- list()
save_progress <- function() write.csv(do.call(rbind, all_results), out_csv, row.names = FALSE)

# ── Meta-analysis: 5 conditions vs Proliferating ────────────────────────────
cat("\n== Meta-analysis conditions ==\n")
rse <- readRDS(file.path(RERUN_DIR, "cs_cq_all_study_processed.rds"))
meta_conditions <- c("Contact_inhibited CQ" = "CICQ", "Serum_starved CQ" = "SSCQ",
                     "Replicative CS" = "RS", "Stress-induced CS" = "SIPS",
                     "Oncogene-induced CS" = "OIS")

for (cond_value in names(meta_conditions)) {
  cond_label <- meta_conditions[[cond_value]]
  keep <- colData(rse)$cell_substate %in% c(cond_value, "Proliferating")
  rse_sub <- rse[, keep]
  full_assay <- as.matrix(assay(rse_sub))
  full_pdata <- as.data.frame(colData(rse_sub))
  cat(sprintf("-- %s vs Proliferating (%d samples) --\n", cond_label, ncol(rse_sub)))

  for (pw in names(pathway_genes)) {
    row <- run_pathway_tage(full_assay, full_pdata, "cell_substate", cond_value, "Proliferating",
                             pathway_genes, model_features, model_paths, "meta_analysis", cond_label, pw)
    all_results[[paste("meta", cond_label, pw)]] <- row
  }
  cat(sprintf("  done: %d/%d OK\n", sum(sapply(all_results[grepl(paste0("^meta ", cond_label, " "), names(all_results))],
                                                function(r) r$status == "OK")), length(pathway_genes)))
  save_progress()
}

# ── Temporal: each cell type, pooled irradiated (4+10+20 days) vs none ─────
cat("\n== Temporal cell types ==\n")
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
  full_assay <- as.matrix(assay(processed))
  full_pdata <- as.data.frame(colData(processed))
  cat(sprintf("-- %s: irradiated vs none (%d samples) --\n", ct, ncol(processed)))

  for (pw in names(pathway_genes)) {
    row <- run_pathway_tage(full_assay, full_pdata, "irradiated", "irradiated", "none",
                             pathway_genes, model_features, model_paths, "temporal", ct, pw)
    all_results[[paste("temporal", ct, pw)]] <- row
  }
  cat(sprintf("  done: %d/%d OK\n", sum(sapply(all_results[grepl(paste0("^temporal ", ct, " "), names(all_results))],
                                                function(r) r$status == "OK")), length(pathway_genes)))
  save_progress()
}

# ── Score + rank ─────────────────────────────────────────────────────────────
final <- do.call(rbind, all_results)
final$p_adj_scaled <- p.adjust(final$wilcox_p_scaled, method = "BH")
final$p_adj_yugene <- p.adjust(final$wilcox_p_yugene, method = "BH")
final$efficiency_score_scaled <- abs(final$cohens_d_scaled) / log2(final$n_overlap_EN_model_features + 2)
final$efficiency_score_yugene <- abs(final$cohens_d_yugene) / log2(final$n_overlap_EN_model_features + 2)
write.csv(final, out_csv, row.names = FALSE)

cat(sprintf("\nDone. %d/%d rows OK. Saved -> %s\n", sum(final$status == "OK", na.rm = TRUE), nrow(final), out_csv))
