# 03_hallmark_enrichment.R
#
# Hallmark-pathway-vs-arrest-DEG overrepresentation, mirroring the main block
# of the original Final/Scripts/for_github/stress_responses.R (lines 1-191).
#
# PORTABILITY CHANGE (flagged, not silent): the original script read a
# manually-downloaded gmt file (`h.all.v2023.2.Hs.symbols.gmt`) from a
# personal Downloads folder. That file isn't in this repo and can't be
# fetched by a script. This version instead pulls the same MSigDB Hallmark
# (category "H") gene sets programmatically via msigdbr -- same source
# database, but likely a slightly different MSigDB build/version than
# "v2023.2" pinned in the original filename, which is a real (if usually
# tiny) source of drift -- see the discrepancy report for what msigdbr
# version this actually resolved to.
#
# SCOPE REDUCTION (flagged): the original also added a custom "Lysosomal
# Genes" set from a personal file (temp-autophagy.csv, filtered to
# class=='Lysosome') that isn't recoverable from anything in this repo
# (Final/SI_tables/autophagy_genes.csv exists but is a different, differently
# -schemed table with no 'class'/'Lysosome' column). Pure-Hallmark terms only
# here; "Lysosomal Genes" is dropped from pathway_of_interest.

source("R/config.R")
source("R/functions.R")
library(msigdbr)

cat("== MSigDB Hallmark gene sets (via msigdbr) ==\n")
hallmark <- msigdbr(species = "Homo sapiens", collection = "H")
cat(sprintf("  msigdbr package version: %s\n", as.character(packageVersion("msigdbr"))))
result <- data.frame(pathway = hallmark$gs_name, genes = hallmark$gene_symbol)
result$pathway <- gsub(result$pathway, pattern = "_V1", replacement = "")
result$pathway <- gsub(result$pathway, pattern = "_V2", replacement = "")
result$pathway <- gsub(result$pathway, pattern = "_", replacement = " ")
result <- unique(result)
cat(sprintf("  %d pathway x gene rows, %d unique pathways\n", nrow(result), length(unique(result$pathway))))

pathway_of_interest <- c(
  "HALLMARK TNFA SIGNALING VIA NFKB", "HALLMARK P53 PATHWAY", "HALLMARK MYC TARGETS",
  "HALLMARK MTORC1 SIGNALING", "HALLMARK MITOTIC SPINDLE", "HALLMARK INTERFERON GAMMA RESPONSE",
  "HALLMARK INTERFERON ALPHA RESPONSE", "HALLMARK INFLAMMATORY RESPONSE",
  "HALLMARK IL6 JAK STAT3 SIGNALING", "HALLMARK HYPOXIA", "HALLMARK G2M CHECKPOINT",
  "HALLMARK E2F TARGETS", "HALLMARK APOPTOSIS", "HALLMARK DNA REPAIR"
  # "Lysosomal Genes" dropped -- see header note.
)

save_csv(result, file_name = "stress_response_pathways_RERUN.csv", path = RERUN_DIR)

cat("\n== Overlap: Hallmark pathways vs arrest DEGs ==\n")
arrest_degs_merged <- read.csv(file.path(RERUN_DIR, "arrest_degs_final_RERUN.csv"))
arrest_degs_merged <- factor_column_and_modify(
  df = arrest_degs_merged, column = "group_1",
  old_list = c("Contact_inhibited_CQ", "Serum_starved_CQ", "Replicative_CS", "Stress_induced_CS", "Oncogene_induced_CS"),
  keyword = NULL
)
arrest_degs_merged_sig <- arrest_degs_merged[arrest_degs_merged$sig == "y", ]

ora_main_stress <- overlap_function(
  df_1 = result, df_2 = arrest_degs_merged_sig,
  gene_col_1 = "genes", gene_col_2 = "gene",
  group_col_2 = c("dir_accession"), group_col_1 = c("pathway"),
  carry_col_2 = c("group_1", "direction_1"),
  background = arrest_degs_merged$gene
)
save_csv(ora_main_stress, file_name = "hallmark_vs_arrest_degs_RERUN.csv", path = RERUN_DIR)

cat("\n== Significant (adj<0.05) Hallmark x condition x direction hits, pathway_of_interest only ==\n")
sig_hits <- ora_main_stress[ora_main_stress$pathway %in% pathway_of_interest & ora_main_stress$adj < 0.05, ]
print(sig_hits[order(sig_hits$pathway, sig_hits$group_1), c("pathway", "group_1", "direction_1", "actual", "odds", "adj")])
save_csv(sig_hits, file_name = "hallmark_vs_arrest_degs_SIGNIFICANT_RERUN.csv", path = RERUN_DIR)

cat(sprintf("\nDone. -> %s\n", file.path(RERUN_DIR, "hallmark_vs_arrest_degs_RERUN.csv")))
