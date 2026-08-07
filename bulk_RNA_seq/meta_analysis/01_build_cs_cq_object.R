# 01_build_cs_cq_object.R
#
# Regenerates the pooled, recount3-scaled recount3 RangedSummarizedExperiment
# across ALL arrest/senescence conditions used in the paper's bulk RNA-seq
# meta-analysis: Proliferating, Contact-inhibited CQ, Serum-starved CQ,
# Replicative CS (RS), Oncogene-induced CS (OIS), Stress-induced CS (SIPS).
#
# This reconstructs the object the original pipeline saved as
# 'cs_cq_all_study_processed_degs.rds' (Final/Scripts/for_github/study_degs.R,
# ~line 466), which was never copied out of the author's local machine and so
# is not present in this repo. Every sample here comes from the SAME 34-study
# recount3 pool as the paper (Final/SI_tables/study_info_all.csv), which is
# 100% Fibroblast (Foreskin/Skin/Lung) by construction — no extra filtering
# needed to restrict to fibroblasts.
#
# Network + compute heavy: downloads gene-level counts for 34 SRA/ENA studies
# from recount3 and hits the Ensembl archive biomart for the protein-coding
# gene dictionary. Expect this to take a while on a fresh run.

# Run this from bulk_RNA_seq/ (e.g. `cd bulk_RNA_seq && Rscript meta_analysis/01_build_cs_cq_object.R`),
# or set BULK_RNASEQ_DIR to the systems_analysis_arrest/ project root first.
source("R/config.R")
source("R/functions.R")

cat("== Step 1: protein-coding gene dictionary (Ensembl release 100, exact-pinned) ==\n")
# Fetched from Ensembl's plain FTP GTF archive for release 100 (the release
# the original apr2020.archive.ensembl.org pipeline used) -- see
# get_ensembl_release_pc() in R/functions.R for why this, not live biomaRt.
human_pc <- get_ensembl_release_pc()
cat(sprintf("  %d protein-coding gene records (Ensembl release %d, exact-pinned)\n",
            nrow(human_pc), ENSEMBL_PINNED_RELEASE))

# Verification against the one object from the original pipeline we DO still
# have (cq_samples.rds, gene-symbol rownames already resolved through the
# original apr2020/Ensembl-100 dictionary): confirms the pin is exact.
old_cq_path <- file.path(SAVE_DIR_CSV, "cq_samples.rds")
if (file.exists(old_cq_path)) {
  old_genes <- rownames(readRDS(old_cq_path))
  overlap_n <- sum(old_genes %in% human_pc$external_gene_name)
  cat(sprintf("  Exact-pin verification: %d / %d genes from the ORIGINAL cq_samples.rds\n",
              overlap_n, length(old_genes)))
  cat(sprintf("  reproduced by this dictionary (%.1f%% agreement -- expect 100%%)\n",
              100 * overlap_n / length(old_genes)))
}

cat("\n== Step 2: sample/study metadata (34 studies, all Fibroblast) ==\n")
meta_both <- read.csv(file.path(SAVE_DIR_CSV, "study_info_all.csv"))
meta_both$cq_test <- meta_both$cell_substate
temp_label <- unique(meta_both$cq_test)[grepl(unique(meta_both$cq_test), pattern = "CQ")]
meta_both$cq_test <- ifelse(meta_both$cell_substate %in% temp_label, "CQ", meta_both$cell_substate)

cat(sprintf("  %d samples across %d studies\n", nrow(meta_both), length(unique(meta_both$study))))
print(table(meta_both$cell_substate))
stopifnot(all(meta_both$cell_type == "Fibroblast"))

cat("\n== Step 3: download gene-level counts from recount3 (34 studies) ==\n")
cs_cq_download <- download_studies(studies = unique(meta_both$study), sra_organism = "human")
cat(sprintf("  Downloaded RSE: %d genes x %d samples\n", nrow(cs_cq_download), ncol(cs_cq_download)))
saveRDS(cs_cq_download, file.path(RERUN_DIR, "cs_cq_download_raw.rds"))

cat("\n== Step 4: process_rse (protein-coding filter, ortholog/symbol rename,\n")
cat("   per-condition low-expression filter, recount3 AUC scale_counts=TRUE) ==\n")
cs_cq_all_study_processed <- process_rse(
  rse = cs_cq_download,
  meta_add = meta_both,
  meta_id_col = "external_id",
  filter_pc = TRUE,
  condition_col = "cell_substate",
  ensembl_dictionary = human_pc,
  filter_duplicates = TRUE,
  sample_filter = unique(meta_both[["external_id"]]),
  low_expression_filter = "per",
  scale_counts = TRUE
)
cat(sprintf("  Final object: %d genes x %d samples\n",
            nrow(cs_cq_all_study_processed), ncol(cs_cq_all_study_processed)))
print(table(data.frame(colData(cs_cq_all_study_processed))$cell_substate))

saveRDS(cs_cq_all_study_processed, file.path(RERUN_DIR, "cs_cq_all_study_processed.rds"))
cat(sprintf("\nSaved to %s\n", file.path(RERUN_DIR, "cs_cq_all_study_processed.rds")))

