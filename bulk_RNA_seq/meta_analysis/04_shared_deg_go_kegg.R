# 04_shared_deg_go_kegg.R
#
# GO + KEGG enrichment of the genes shared across all 5 arrest conditions
# (up-regulated set and down-regulated set separately), mirroring the
# corresponding block of the original
# Final/Scripts/for_github/DEG_analysis_recount3.R (~lines 424-505).
# This is the comparison point for the paper's SI Figure 5a/5b ("Cell cycle",
# "DNA replication", "Cellular senescence" as top downregulated-shared-DEG
# terms -- see the discrepancy report).

source("R/config.R")
source("R/functions.R")

arrest_degs_merged <- read.csv(file.path(RERUN_DIR, "arrest_degs_final_RERUN.csv"))

human_pc_entrez <- get_ensembl_release_pc_entrez()

arrest_degs_merged_up <- arrest_degs_merged[arrest_degs_merged$direction_1 == "up", ]
arrest_degs_merged_up_table <- table(arrest_degs_merged_up$gene)
arrest_degs_merged_up_background <- names(arrest_degs_merged_up_table)[arrest_degs_merged_up_table == 5]

arrest_degs_merged_down <- arrest_degs_merged[arrest_degs_merged$direction_1 == "down", ]
arrest_degs_merged_down_table <- table(arrest_degs_merged_down$gene)
arrest_degs_merged_down_background <- names(arrest_degs_merged_down_table)[arrest_degs_merged_down_table == 5]

arrest_degs_merged_up_sig <- arrest_degs_merged[arrest_degs_merged$direction_1 == "up" & arrest_degs_merged$sig == "y", ]
common_up <- find_common_genes(df = arrest_degs_merged_up_sig, group_col = "group_1_dir", gene_col = "gene")

arrest_degs_merged_down_sig <- arrest_degs_merged[arrest_degs_merged$direction_1 == "down" & arrest_degs_merged$sig == "y", ]
common_down <- find_common_genes(df = arrest_degs_merged_down_sig, group_col = "group_1_dir", gene_col = "gene")

cat(sprintf("Shared UP DEGs across all 5 conditions:   %d  (paper Table 2/SI Fig 3 figure: 101)\n", length(common_up)))
cat(sprintf("Shared DOWN DEGs across all 5 conditions: %d  (paper Table 2/SI Fig 3 figure: 316)\n", length(common_down)))
write.csv(data.frame(gene = common_up),   file.path(RERUN_DIR, "shared_up_genes_RERUN.csv"),   row.names = FALSE)
write.csv(data.frame(gene = common_down), file.path(RERUN_DIR, "shared_down_genes_RERUN.csv"), row.names = FALSE)

cat("\n== GO + KEGG enrichment: shared DOWN-regulated DEGs ==\n")
enrich_down <- enrich_genes(
  gene_list = common_down, background = arrest_degs_merged_down_background,
  gene_dictionary = human_pc_entrez, use_ensembl = FALSE
)
save_csv(enrich_down$enrichment, file_name = "common_deg_enrichment_RERUN.csv", path = RERUN_DIR)
top_kegg_down <- enrich_down$enrichment[enrich_down$enrichment$enrichment == "KEGG", ]
top_kegg_down <- top_kegg_down[order(top_kegg_down$p.adjust), ]
cat("Top 10 KEGG terms (shared down DEGs), by p.adjust (paper: Cell cycle top hit, incl. 'Cellular senescence' pathway):\n")
print(head(top_kegg_down[, c("Description", "GeneRatio", "pvalue", "p.adjust")], 10))

cat("\n== GO + KEGG enrichment: shared UP-regulated DEGs ==\n")
enrich_up <- enrich_genes(
  gene_list = common_up, background = arrest_degs_merged_up_background,
  gene_dictionary = human_pc_entrez, use_ensembl = FALSE
)
save_csv(enrich_up$enrichment, file_name = "common_up_deg_enrichment_RERUN.csv", path = RERUN_DIR)

cat(sprintf("\nDone. Enrichment CSVs -> %s\n", RERUN_DIR))
