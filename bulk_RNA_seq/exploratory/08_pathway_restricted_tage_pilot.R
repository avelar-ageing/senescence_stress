# 08_pathway_restricted_tage_pilot.R
#
# PILOT / EXPLORATORY (not part of the verified pipeline). For ONE condition
# (SIPS vs Proliferating, the strongest tAge signal in the full-gene-set
# analysis), loop over every MSigDB Hallmark pathway, restrict the input
# expression matrix to ONLY that pathway's genes, rerun tAge end to end on
# that restricted gene set, and report:
#   - how many of the pathway's genes actually survive into the EN model's
#     10,487-gene feature space (the rest get NA-imputed -- a pathway with
#     poor overlap is testing the imputer, not the pathway)
#   - whether the resulting pathway-restricted tAge still separates SIPS from
#     Proliferating (Wilcoxon), i.e. whether that pathway alone carries
#     enough age-like signal to be "powered" for this question at all.
#
# Sample n does NOT change across pathways (same 130 SIPS+Proliferating
# samples throughout) -- what varies is purely how much real (non-imputed)
# gene-level signal each pathway can feed the model.

source("R/config.R")
source("R/functions.R")
suppressPackageStartupMessages({
  library(msigdbr)
  library(tAge)
  library(Biobase)
})

SCRNA_DIR <- "/home/ro/APFS_copy/root/Backup/Documents/modules/scrna_seq/GSE226225/final_analysis"
MODEL_DIR <- file.path(SCRNA_DIR, "tAge_models")
Sys.setenv(RETICULATE_PYTHON = file.path(SCRNA_DIR, ".venv/bin/python"))
model_paths <- list(
  scaled_diff = file.path(MODEL_DIR, "EN_Chronoage_Multispecies_Multitissue_scaleddiff.pkl"),
  yugene_diff = file.path(MODEL_DIR, "EN_Chronoage_Multispecies_Multitissue_yugenediff.pkl")
)

# Model's actual feature space (mouse-ortholog gene symbols), for reporting
# per-pathway overlap -- extracted once via python/joblib, same approach as
# the earlier scRNA gene-dropout accounting.
get_model_features <- function(pkl_path) {
  py <- file.path(SCRNA_DIR, ".venv/bin/python")
  tmp <- tempfile(fileext = ".txt")
  script <- sprintf(
    "import joblib, warnings\nwarnings.filterwarnings('ignore')\nm = joblib.load('%s')\nwith open('%s','w') as f:\n    f.write('\\n'.join(m.feature_names_in_))\n",
    pkl_path, tmp
  )
  script_path <- tempfile(fileext = ".py")
  writeLines(script, script_path)
  system2(py, script_path)
  readLines(tmp)
}
model_features <- get_model_features(model_paths$scaled_diff)
cat(sprintf("EN model feature space: %d genes\n", length(model_features)))

# ── Subset to SIPS + Proliferating only ─────────────────────────────────────
rse <- readRDS(file.path(RERUN_DIR, "cs_cq_all_study_processed.rds"))
keep <- colData(rse)$cell_substate %in% c("Stress-induced CS", "Proliferating")
rse_sub <- rse[, keep]
cat(sprintf("SIPS + Proliferating subset: %d genes x %d samples\n", nrow(rse_sub), ncol(rse_sub)))

full_assay <- as.matrix(assay(rse_sub))
full_pdata <- as.data.frame(colData(rse_sub))

# ── Hallmark pathways ────────────────────────────────────────────────────────
hallmark <- msigdbr(species = "Homo sapiens", collection = "H")
pathway_genes <- split(hallmark$gene_symbol, hallmark$gs_name)
cat(sprintf("Testing %d Hallmark pathways\n", length(pathway_genes)))

results_list <- list()
out_csv <- file.path(RERUN_DIR, "pathway_restricted_tage_pilot_SIPS.csv")

for (pw in names(pathway_genes)) {
  genes_pw <- unique(pathway_genes[[pw]])
  genes_in_data <- intersect(genes_pw, rownames(full_assay))

  # model_features is the EN model's actual feature space, which is MOUSE
  # gene symbols (tAge maps human -> mouse orthologs internally before
  # prediction) -- map genes_in_data through the same ortholog table before
  # comparing, or every pathway looks like 0% overlap regardless of pathway.
  n_model_overlap <- 0L
  if (length(genes_in_data) >= 10) {
    mapped_mouse_genes <- tryCatch({
      tmp_eset <- ExpressionSet(assayData = full_assay[genes_in_data, , drop = FALSE])
      mapped <- tAge:::map_genes(tmp_eset, "human", "Gene.Symbol", verbose = FALSE)
      unique(rownames(mapped))
    }, error = function(e) character(0))
    n_model_overlap <- length(intersect(mapped_mouse_genes, model_features))
  }

  res_row <- data.frame(
    pathway = pw,
    n_pathway_genes = length(genes_pw),
    n_in_expression_data = length(genes_in_data),
    n_overlap_EN_model_features = n_model_overlap,
    pct_of_model_populated = round(100 * n_model_overlap / length(model_features), 2),
    status = NA_character_,
    n_sips = NA_integer_, n_prolif = NA_integer_,
    median_scaled_diff_SIPS = NA_real_, median_scaled_diff_Prolif = NA_real_,
    wilcox_p_scaled = NA_real_,
    median_yugene_diff_SIPS = NA_real_, median_yugene_diff_Prolif = NA_real_,
    wilcox_p_yugene = NA_real_,
    stringsAsFactors = FALSE
  )

  if (length(genes_in_data) < 10) {
    res_row$status <- "SKIPPED (<10 genes in data)"
    results_list[[pw]] <- res_row
    next
  }

  tryCatch({
    sub_assay <- full_assay[genes_in_data, , drop = FALSE]
    eset <- ExpressionSet(assayData = sub_assay, phenoData = AnnotatedDataFrame(full_pdata))

    tAge_eset <- suppressWarnings(tAge_preprocessing(
      eset, species = "human", gene_mapping_type = "Gene.Symbol",
      control_group_column = "cell_state", control_group_label = "Proliferating",
      verbose = FALSE, count_threshold = 10, percent_threshold = 20
    ))
    pred <- suppressWarnings(predict_tAge(tAge_eset, model_paths, species = "human", mode = "EN"))

    sips_rows <- pred[pred$cell_substate == "Stress-induced CS", ]
    prolif_rows <- pred[pred$cell_substate == "Proliferating", ]

    res_row$status <- "OK"
    res_row$n_sips <- nrow(sips_rows)
    res_row$n_prolif <- nrow(prolif_rows)
    res_row$median_scaled_diff_SIPS <- median(sips_rows$scaled_diff_EN_tAge, na.rm = TRUE)
    res_row$median_scaled_diff_Prolif <- median(prolif_rows$scaled_diff_EN_tAge, na.rm = TRUE)
    res_row$wilcox_p_scaled <- wilcox.test(sips_rows$scaled_diff_EN_tAge, prolif_rows$scaled_diff_EN_tAge)$p.value
    res_row$median_yugene_diff_SIPS <- median(sips_rows$yugene_diff_EN_tAge, na.rm = TRUE)
    res_row$median_yugene_diff_Prolif <- median(prolif_rows$yugene_diff_EN_tAge, na.rm = TRUE)
    res_row$wilcox_p_yugene <- wilcox.test(sips_rows$yugene_diff_EN_tAge, prolif_rows$yugene_diff_EN_tAge)$p.value
  }, error = function(e) {
    res_row$status <<- paste("ERROR:", conditionMessage(e))
  })

  results_list[[pw]] <- res_row
  cat(sprintf("  [%3d/%3d] %-45s genes=%3d/%3d model_overlap=%4d(%.1f%%) status=%s\n",
              length(results_list), length(pathway_genes), pw,
              res_row$n_in_expression_data, res_row$n_pathway_genes,
              res_row$n_overlap_EN_model_features, res_row$pct_of_model_populated, res_row$status))

  # write incrementally so partial progress survives if something crashes
  write.csv(do.call(rbind, results_list), out_csv, row.names = FALSE)
}

final <- do.call(rbind, results_list)
final$p_adj_scaled <- p.adjust(final$wilcox_p_scaled, method = "BH")
final$p_adj_yugene <- p.adjust(final$wilcox_p_yugene, method = "BH")
write.csv(final, out_csv, row.names = FALSE)

cat(sprintf("\nDone. %d/%d pathways completed OK.\n", sum(final$status == "OK", na.rm = TRUE), nrow(final)))
cat(sprintf("Saved -> %s\n", out_csv))
