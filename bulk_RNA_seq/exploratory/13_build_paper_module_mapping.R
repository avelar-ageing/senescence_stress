# 13_build_paper_module_mapping.R
#
# Builds a module -> clock-feature mapping from the paper's published WGCNA
# modules (Supplementary Table 5, extracted by 11_extract_paper_module_clocks.R),
# in the same format 04_partial_tage_decompose.py expects for Hallmark pathways
# (columns: pathway, mouse_gene_id). This lets the existing decomposition run
# unchanged with the paper's grouping unit substituted for MSigDB Hallmark.
#
# WHY THIS AND NOT THE PAPER'S MODULE CLOCKS. Evaluating the published module
# clocks directly fails validation (PARTIAL_TAGE_VS_PAPER.md section 6b): their
# coefficients are fitted with their own imputer/scaler, which is NOT published
# -- only coefficients and intercepts are -- so plugging them into our global
# model's feature space gives predictions that correlate r = 0.22-0.48 with the
# package's own tAge on the package's own example data, and that reverse the
# sign of the SIPS/OIS elevations on ours. What IS usable is module *membership*.
# Using the paper's modules as the grouping unit for our own (exact, verified)
# decomposition keeps the estimator we have validated while adopting the
# paper's data-derived grouping.
#
# WHY THE MULTISPECIES PANEL. It matches the Multispecies_Multitissue models
# both arms already use, and its coverage is complete: 14 modules, 1,248 genes,
# ZERO gene shared between modules (perfectly disjoint), and 1,248/1,248 present
# among the clock's 10,487 features. The rodent panel is 23 modules / 1,922
# genes but only 65% are in the clock's feature space, so it is not used here.
#
# Disjointness matters: it removes consequences (1) and (2) in
# PARTIAL_TAGE_VS_PAPER.md section 5 -- with non-overlapping modules the partial
# scores are far closer to a genuine decomposition and much less
# cross-correlated than the overlapping Hallmark sets.

source("R/config.R")

OUT_DIR <- file.path(RERUN_DIR, "partial_tage")
clocks_csv <- file.path(OUT_DIR, "paper_module_clocks.csv")
if (!file.exists(clocks_csv)) {
  stop(sprintf("%s not found -- run exploratory/11_extract_paper_module_clocks.R first.", clocks_csv))
}

d <- read.csv(clocks_csv)
s <- d[d$panel == "multispecies" &
         d$outcome == "Chronological age" &
         d$module != "All module genes", ]

# Readable label: module colour is the paper's identifier, the annotation is
# what it means. Both are kept so results can be traced back to the paper.
s$label <- sprintf("%s (%s)", s$annotation, s$module)

final <- unique(data.frame(pathway = s$label, mouse_gene_id = s$entrez_id))
final <- final[order(final$pathway, final$mouse_gene_id), ]

write.csv(final, file.path(OUT_DIR, "paper_module_mouse_ids.csv"), row.names = FALSE)

cat(sprintf("%d modules, %d gene assignments over %d unique genes (overlap = %d)\n",
            length(unique(final$pathway)), nrow(final),
            length(unique(final$mouse_gene_id)),
            nrow(final) - length(unique(final$mouse_gene_id))))
print(as.data.frame(table(final$pathway)), row.names = FALSE)
cat(sprintf("\nSaved -> %s\n", file.path(OUT_DIR, "paper_module_mouse_ids.csv")))
