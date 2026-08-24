# 02_call_degs.R
#
# DEG calling for the meta-analysis pipeline: each arrest/senescence
# condition vs its pooled Proliferating control, batch-corrected for `study`
# (DESeq2, ~cell_substate + study design). Mirrors the active (uncommented)
# block of the original Final/Scripts/for_github/study_degs.R (lines
# ~483-660) exactly -- same function calls, same arguments, same DEG
# thresholds (those live inside find_degs_new/deseq_studies in functions.R,
# unchanged). Produces arrest_degs_final.csv, the file every downstream
# meta-analysis script (DEG_analysis_recount3.R, cq_DEG_analysis_recount3.R,
# degs_v_stress.R, stress_responses.R, metabolism_script.R) consumes.
#
# Requires 01_build_cs_cq_object.R to have already run (reads its
# cs_cq_all_study_processed.rds from RERUN_DIR).

source("R/config.R")
source("R/functions.R")

cs_cq_all_study_processed <- readRDS(file.path(RERUN_DIR, "cs_cq_all_study_processed.rds"))

cell_counts_table <- data.frame(table(data.frame(colData(cs_cq_all_study_processed)) %>%
                                         dplyr::select(tissue, cell_substate)))
save_csv(cell_counts_table, file_name = "sample_subtype_counts_RERUN.csv", path = RERUN_DIR)

temp_meta <- data.frame(colData(cs_cq_all_study_processed))
save_csv(temp_meta, file_name = "sample_metadata_RERUN.csv", path = RERUN_DIR)

cat("== DESeq2 fit (~cell_substate + study) ==\n")
cs_cq_studies_processed_deseq <- deseq_studies(
  study_obj = cs_cq_all_study_processed,
  protect_col = "cell_substate",
  fix_col = "study"
)
saveRDS(cs_cq_studies_processed_deseq, file.path(RERUN_DIR, "cs_cq_all_study_processed_deseq2.rds"))

cat("\n== DEG calling: each condition vs Proliferating ==\n")
rs_degs    <- find_degs_new(rse_obj = cs_cq_all_study_processed, group_col = "cell_substate",
                             group_1 = "Replicative CS",       group_2 = "Proliferating", batch_col = "study")
ois_degs   <- find_degs_new(rse_obj = cs_cq_all_study_processed, group_col = "cell_substate",
                             group_1 = "Oncogene-induced CS",  group_2 = "Proliferating", batch_col = "study")
sips_degs  <- find_degs_new(rse_obj = cs_cq_all_study_processed, group_col = "cell_substate",
                             group_1 = "Stress-induced CS",    group_2 = "Proliferating", batch_col = "study")
cq_ci_degs <- find_degs_new(rse_obj = cs_cq_all_study_processed, group_col = "cell_substate",
                             group_1 = "Contact_inhibited CQ", group_2 = "Proliferating", batch_col = "study")
cq_ss_degs <- find_degs_new(rse_obj = cs_cq_all_study_processed, group_col = "cell_substate",
                             group_1 = "Serum_starved CQ",     group_2 = "Proliferating", batch_col = "study")

arrest_degs <- rbind(rs_degs$degs, ois_degs$degs, sips_degs$degs, cq_ci_degs$degs, cq_ss_degs$degs)
save_csv(arrest_degs, file_name = "arrest_degs_final_RERUN.csv", path = RERUN_DIR)

cat("\n== DEG counts by condition x direction (rerun) ==\n")
deg_counts_rerun <- arrest_degs %>%
  dplyr::filter(sig == "y") %>%
  dplyr::count(group_1, direction_1)
print(deg_counts_rerun)
save_csv(deg_counts_rerun, file_name = "deg_count_RERUN.csv", path = RERUN_DIR)

cat(sprintf("\nDone. arrest_degs_final_RERUN.csv -> %s\n", file.path(RERUN_DIR, "arrest_degs_final_RERUN.csv")))
