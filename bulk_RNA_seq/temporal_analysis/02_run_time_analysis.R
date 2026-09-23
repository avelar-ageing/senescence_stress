# 02_run_time_analysis.R
#
# Runs the ERP021140 temporal DEG analysis for all three cell types
# (Fibroblast, Keratinocyte, Melanocyte), mirroring the three time_analysis()
# calls in the original Final/Scripts/for_github/temporal_deg_analysis.R
# (lines 34-111) exactly: same time_levels, same pval/log2FC cutoffs, same
# Keratinocyte batch-correction sample list, same simulation_n=10000.
#
# db/db_col are left NULL here: every original call used enrich_degs_test=FALSE,
# and db/db_col are only read inside time_analysis() when enrich_degs_test=TRUE
# (general_functions.R ~line 2810) -- so they are genuinely unused on this path.
# (The original script passed db=cellage, an object it never defines --
# a latent bug in the original that never mattered because it's dead code
# on the enrich_degs_test=FALSE path. Not reproduced here.)
#
# Requires 01_build_metadata.R to have already run.

source("R/config.R")
source("R/functions.R")

recount_pheno <- readRDS(file.path(RERUN_DIR, "temporal_recount_pheno.rds"))

cat("== Protein-coding + Entrez gene dictionary (Ensembl release 100, exact-pinned) ==\n")
human_pc <- get_ensembl_release_pc()
human_pc_entrez <- get_ensembl_release_pc_entrez()

search_1_time <- build_search(column = "sra", search_terms = c("ERP021140"))

# NB: time_analysis() builds its own subfolder path as
# paste0(results_dir, study_temp, sep) with NO separator inserted between
# results_dir and study_temp -- so results_dir must already end in "/", or
# this silently concatenates onto the parent dir name instead of nesting
# into it (that's what happened here on the previous run: it wrote to a
# sibling "rerun_outputsERP021140/" dir instead of "rerun_outputs/ERP021140/",
# then errored on dir.create() when a run tried to create that path a second
# time with a different parent).
temporal_results_dir <- paste0(RERUN_DIR, "/")
dir.create(temporal_results_dir, showWarnings = FALSE, recursive = TRUE)

run_one <- function(cell_type, fix_batch, batch_1 = NULL) {
  cat(sprintf("\n== time_analysis: %s ==\n", cell_type))
  time_analysis(
    cell_type = cell_type,
    recount_pheno = recount_pheno,
    ensembl_dictionary = human_pc,
    grab_all = FALSE,
    enrich_degs_test = FALSE,
    entrez_dictionary = human_pc_entrez,
    exclude_search = NULL,
    ensembl_gene_col = "ensembl_gene_id",
    treatment_time_col = "time_after_treatment",
    time_levels = c("none", "4_days", "10_days", "20_days"),
    pval_cutoff = 0.05,
    log2fc_cutoff = log2(1.5),
    independentFiltering = TRUE,
    db = NULL,
    db_col = NULL,
    search_term_study = search_1_time,
    sep = "/",
    results_dir = temporal_results_dir,
    fix_batch = fix_batch,
    batch_1 = batch_1,
    run_simulations = FALSE,  # skipped on request -- DEG counts don't depend on this;
                              # rerun with TRUE, simulation_n=10000 if the permutation
                              # test itself (not just the DEG counts) needs reproducing
    pca_plot_group = c("cell_state", "time_after_treatment"),
    vst_n = 500
  )
}

# Keratinocyte: known batch effect (two researchers processed it, main PDF
# p.21/p.37-38) -- same 12-sample batch_1 list as the original.
run_one("Keratinocyte", fix_batch = TRUE, batch_1 = c(
  "ERR1805235", "ERR1805236", "ERR1805238", "ERR1805230", "ERR1805231", "ERR1805224",
  "ERR1805223", "ERR1805222", "ERR1805239", "ERR1805240", "ERR1805241", "ERR1805229"
))
run_one("Melanocyte", fix_batch = FALSE)
run_one("Fibroblast", fix_batch = FALSE)

cat(sprintf("\nDone. Per-cell-type outputs under %s\n", temporal_results_dir))
