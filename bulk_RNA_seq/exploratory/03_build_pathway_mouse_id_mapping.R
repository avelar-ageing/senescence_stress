# 03_build_pathway_mouse_id_mapping.R
#
# Builds the human-gene-symbol -> mouse-ortholog-ID mapping for every pathway
# in stress_response_pathways_RERUN.csv -- the exact pathway x gene table
# meta_analysis/03_hallmark_enrichment.R already builds and uses for the
# DEG-overlap enrichment (MSigDB Hallmark with V1/V2 merged, e.g.
# MYC_TARGETS_V1 + MYC_TARGETS_V2 -> one "HALLMARK MYC TARGETS" pathway, plus
# the custom "Lysosomal Genes" set, human_pc-filtered). Read from that single
# source of truth rather than re-deriving the merge here, so the tAge
# partial-decomposition pathways and the DEG-overlap enrichment pathways can
# never drift apart. This mapping to mouse orthologs is a fixed lookup
# (tAge's own ortholog table, not data-dependent), reused by every group's
# partial-tAge decomposition in 04_partial_tage_decompose.py.

source("R/config.R")
suppressPackageStartupMessages({
  library(tAge)
  library(Biobase)
})

OUT_DIR <- file.path(RERUN_DIR, "partial_tage")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

pathways_csv <- file.path(RERUN_DIR, "stress_response_pathways_RERUN.csv")
if (!file.exists(pathways_csv)) {
  stop(sprintf("%s not found -- run meta_analysis/03_hallmark_enrichment.R first.", pathways_csv))
}
all_pathway_genes <- read.csv(pathways_csv)
cat(sprintf("Loaded %d unique pathways, %d rows from %s\n",
            length(unique(all_pathway_genes$pathway)), nrow(all_pathway_genes), pathways_csv))

pathway_genes <- split(all_pathway_genes$genes, all_pathway_genes$pathway)

rows <- list()
for (pw in names(pathway_genes)) {
  genes <- unique(pathway_genes[[pw]])
  dummy_mat <- matrix(1, nrow = length(genes), ncol = 1, dimnames = list(genes, "s1"))
  eset <- ExpressionSet(assayData = dummy_mat)
  mapped <- tryCatch(tAge:::map_genes(eset, "human", "Gene.Symbol", verbose = FALSE), error = function(e) NULL)
  if (is.null(mapped) || nrow(mapped) == 0) next
  mouse_ids <- unique(rownames(mapped))
  rows[[pw]] <- data.frame(pathway = pw, mouse_gene_id = mouse_ids)
  cat(sprintf("%-45s %4d human genes -> %4d mouse IDs\n", pw, length(genes), length(mouse_ids)))
}
final <- do.call(rbind, rows)
write.csv(final, file.path(OUT_DIR, "hallmark_pathway_mouse_ids.csv"), row.names = FALSE)
cat(sprintf("\nSaved %d rows -> %s\n", nrow(final), file.path(OUT_DIR, "hallmark_pathway_mouse_ids.csv")))
