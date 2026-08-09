# 03_build_pathway_mouse_id_mapping.R
#
# Builds the human-gene-symbol -> mouse-ortholog-ID mapping for every MSigDB
# Hallmark pathway, once. This is a fixed lookup (tAge's own ortholog table,
# not data-dependent), reused by every group's partial-tAge decomposition in
# 04_partial_tage_decompose.py.

source("R/config.R")
suppressPackageStartupMessages({
  library(tAge)
  library(Biobase)
  library(msigdbr)
})

OUT_DIR <- file.path(RERUN_DIR, "partial_tage")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

hallmark <- msigdbr(species = "Homo sapiens", collection = "H")
pathway_genes <- split(hallmark$gene_symbol, hallmark$gs_name)

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
