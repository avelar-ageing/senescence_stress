# 23_tpm_vs_counts_control.R
#
# Control for the one methodological compromise in the GSE175533 replication
# (meta_analysis/20-22): GEO deposits no count matrix for that series, so the
# analysis is run on SALMON TPM. This script tests, on OUR OWN count data where
# both versions can be produced, whether that substitution changes tAge.
#
# WHY IT SHOULD NOT, and what could still go wrong.
#   TPM_ij = c_ij / (len_i * L_j) * 1e6. tAge_preprocessing z-scores every gene
#   across samples (scale_eset) before subtracting the control group, and a z-score
#   is invariant to any gene-specific positive constant - so the len_i term cancels
#   EXACTLY, not approximately. The per-sample factor L_j is what RLE normalisation
#   removes in either case. The residual risk is not in the arithmetic but in GENE
#   SELECTION: filter_genes keeps genes with >= 10 counts in >= 20% of samples, and
#   that filter is applied to whatever matrix it is given, so dividing by gene
#   length reshuffles which genes clear the threshold. Short genes gain, long genes
#   lose. That is what this script measures.
#
# DESIGN. Take the meta-analysis count matrix, convert it to TPM using simulated
# transcript lengths (log-normal, median 2 kb - the real Ensembl lengths are not
# needed and arbitrary lengths are the stronger test, since invariance to ANY
# per-gene constant is the claim), rescale to a nominal 20M library exactly as
# script 20 does, and push both matrices through the identical preprocessing and
# the same three clocks. Repeated over several length draws so the answer cannot be
# an artefact of one seed.
#
# Output: rerun_outputs/gse175533/tpm_vs_counts_control.csv

source("R/config.R")
suppressPackageStartupMessages({
  library(tAge)
  library(Biobase)
  library(SummarizedExperiment)
})

OUT <- file.path(RERUN_DIR, "gse175533")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
SEEDS <- 1:5
LIB <- 20e6

rse <- readRDS(file.path(RERUN_DIR, "cs_cq_all_study_processed.rds"))
counts <- as.matrix(assay(rse))
pheno <- as.data.frame(colData(rse))

run_pipeline <- function(mat, label) {
  eset <- ExpressionSet(assayData = mat,
                        phenoData = AnnotatedDataFrame(pheno))
  pre <- suppressWarnings(tAge_preprocessing(
    eset, species = "human", gene_mapping_type = "Gene.Symbol",
    control_group_column = "cell_substate", control_group_label = "Proliferating",
    verbose = FALSE, count_threshold = 10, percent_threshold = 20))
  list(scaled = t(exprs(pre$scaled_diff)), yugene = t(exprs(pre$yugene_diff)),
       n_genes_kept = nrow(pre$RLE_normalized))
}

cat("== baseline: raw counts ==\n")
base <- run_pipeline(counts, "counts")
cat(sprintf("   genes retained: %d\n", base$n_genes_kept))

write_mat <- function(m, path) {
  write.csv(data.frame(sample_id = rownames(m), m, check.names = FALSE),
            path, row.names = FALSE)
}
write_mat(base$scaled, file.path(OUT, "control_counts_scaled_diff.csv"))
write_mat(base$yugene, file.path(OUT, "control_counts_yugene_diff.csv"))

res <- list()
for (s in SEEDS) {
  set.seed(s)
  len <- 10^rnorm(nrow(counts), mean = log10(2000), sd = 0.35)  # ~2 kb median
  rpk <- counts / len
  tpm <- sweep(rpk, 2, colSums(rpk), "/") * 1e6
  pseudo <- sweep(tpm, 2, colSums(tpm), "/") * LIB

  cat(sprintf("== seed %d: TPM-transformed ==\n", s))
  alt <- run_pipeline(pseudo, sprintf("tpm_seed%d", s))

  # gene-set overlap after the detection filter, the one real risk
  jac <- length(intersect(colnames(base$scaled)[colSums(!is.na(base$scaled)) > 0],
                          colnames(alt$scaled)[colSums(!is.na(alt$scaled)) > 0]))
  ub <- length(union(colnames(base$scaled)[colSums(!is.na(base$scaled)) > 0],
                     colnames(alt$scaled)[colSums(!is.na(alt$scaled)) > 0]))
  cat(sprintf("   genes retained: %d (counts %d); detected-gene Jaccard %.3f\n",
              alt$n_genes_kept, base$n_genes_kept, jac / ub))

  write_mat(alt$scaled, file.path(OUT, sprintf("control_tpm%d_scaled_diff.csv", s)))
  write_mat(alt$yugene, file.path(OUT, sprintf("control_tpm%d_yugene_diff.csv", s)))
  res[[length(res) + 1]] <- data.frame(
    seed = s, genes_counts = base$n_genes_kept, genes_tpm = alt$n_genes_kept,
    detected_gene_jaccard = jac / ub)
}
write.csv(do.call(rbind, res), file.path(OUT, "tpm_vs_counts_genes.csv"),
          row.names = FALSE)
cat(sprintf("\nMatrices written to %s; run 24_tpm_vs_counts_predict.py to compare tAge\n", OUT))
