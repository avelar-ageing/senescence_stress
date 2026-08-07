# 01_build_metadata.R
#
# Reconstructs the ERP021140 temporal sample metadata table
# (`recount_filtered_quiescence_fixed` in the original
# Final/Scripts/for_github/temporal_deg_analysis.R, line 32) that fed
# time_analysis() for all three cell types. The original file lived at
# /Volumes/GoogleDrive/My Drive/PhD_copy/To Publish/CS Clusters/CS_studies/recount_final.csv
# -- a path from an even older project copy than the rest of the pipeline,
# and it is not present anywhere in this repo.
#
# What IS present: Final/ERP021140/{Fibroblast,Keratinocyte,Melanocyte}/sample_pheno.csv
# -- the final per-cell-type colData that time_analysis() itself saved out
# (general_functions.R's time_analysis(), ~line 2578: `save_csv(data =
# time_coldata, file_name = 'sample_pheno', ...)`). Row-binding those three
# files reconstructs the exact metadata time_analysis() needs, since they
# already carry every column (cell_type, cell_state, time_after_treatment,
# external_id, ...) merged in during the original run -- confirmed 24
# samples per cell type (6 samples x 4 timepoints), matching the paper
# (Methods 5.1: "For each cell type and condition, there were 6 samples").
#
# The one rename needed: time_analysis()'s internal process_rse() call hardcodes
# meta_id_col='sample_ID' (general_functions.R ~line 2566), but these CSVs use
# 'external_id' for the same recount3 sample accession -- renamed below.

source("R/config.R")

pheno_files <- c(
  Fibroblast   = file.path(DATA_DIR, "ERP021140", "Fibroblast",   "sample_pheno.csv"),
  Keratinocyte = file.path(DATA_DIR, "ERP021140", "Keratinocyte", "sample_pheno.csv"),
  Melanocyte   = file.path(DATA_DIR, "ERP021140", "Melanocyte",   "sample_pheno.csv")
)
stopifnot(all(file.exists(pheno_files)))

recount_pheno <- do.call(rbind, lapply(names(pheno_files), function(ct) {
  df <- read.csv(pheno_files[[ct]])
  df
}))
recount_pheno <- unique(recount_pheno)
names(recount_pheno)[names(recount_pheno) == "external_id"] <- "sample_ID"

cat(sprintf("Reconstructed recount_pheno: %d rows (expect ~72 = 3 cell types x 24 samples,\n", nrow(recount_pheno)))
cat("  fewer if some samples are shared/duplicated across the per-cell-type exports)\n")
print(table(recount_pheno$cell_type, recount_pheno$time_after_treatment))

saveRDS(recount_pheno, file.path(RERUN_DIR, "temporal_recount_pheno.rds"))
cat(sprintf("Saved to %s\n", file.path(RERUN_DIR, "temporal_recount_pheno.rds")))
