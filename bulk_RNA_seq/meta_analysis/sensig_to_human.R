# sensig_to_human.R
#
# Maps the mouse "SenSig" senescence signature (Supplementary Table 1 of
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10651581/) onto human 1:1
# orthologs, for use in downstream overlap tests (macrophage_overlaps.R-style
# analyses).
#
# BUG FIX (flagged): the original script referenced `sensig_sig`, an object
# it never created (`sensig_sig$dir=...` after only ever building `sensig`) --
# it would have errored before reaching write.csv() if run top-to-bottom.
# Fixed to operate on `sensig` throughout.
#
# INPUT NOT IN REPO: the mouse SenSig table (11357_2023_785_MOESM2_ESM.csv,
# the paper's own supplementary file) isn't present anywhere in this repo --
# it needs to be re-downloaded from the PMC article above and pointed at via
# SENSIG_INPUT_CSV below.

source("R/config.R")

SENSIG_INPUT_CSV <- Sys.getenv("SENSIG_INPUT_CSV", file.path(RERUN_DIR, "11357_2023_785_MOESM2_ESM.csv"))
if (!file.exists(SENSIG_INPUT_CSV)) {
  stop(sprintf(
    "SenSig input not found at %s. Download Supplementary Table 1 from https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10651581/ and set SENSIG_INPUT_CSV.",
    SENSIG_INPUT_CSV
  ))
}

library(biomaRt)
mart_args <- list(biomart = "ENSEMBL_MART_ENSEMBL")
if (!is.null(ENSEMBL_ARCHIVE_HOST)) mart_args$host <- ENSEMBL_ARCHIVE_HOST

ensembl_mouse <- do.call(useMart, c(mart_args, list(dataset = "mmusculus_gene_ensembl")))
mouse_pc <- getBM(
  attributes = c("external_gene_name", "ensembl_gene_id",
                 "hsapiens_homolog_ensembl_gene", "hsapiens_homolog_associated_gene_name",
                 "hsapiens_homolog_orthology_type", "hsapiens_homolog_orthology_confidence"),
  mart = ensembl_mouse
)

# filter mouse protein-coding genes for confident 1:1 human orthologs
mouse_pc <- mouse_pc[mouse_pc$hsapiens_homolog_orthology_type == "ortholog_one2one", ]
mouse_pc <- mouse_pc[mouse_pc$hsapiens_homolog_orthology_confidence == 1, ]

sensig <- read.csv(SENSIG_INPUT_CSV)
sensig <- merge(sensig, mouse_pc, by.x = "gene_name", by.y = "external_gene_name")
sensig <- sensig[, colnames(sensig) %in% c("hsapiens_homolog_associated_gene_name",
                                           "hsapiens_homolog_ensembl_gene", "logFC", "FDR")]
sensig$dir <- ifelse(sensig$logFC < 0, "down", "up")

write.csv(sensig, file.path(RERUN_DIR, "sensig_homologues.csv"), row.names = FALSE)
cat(sprintf("Done. %d SenSig genes mapped to human orthologs -> %s\n",
            nrow(sensig), file.path(RERUN_DIR, "sensig_homologues.csv")))
