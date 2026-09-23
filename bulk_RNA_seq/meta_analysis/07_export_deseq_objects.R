# 07_export_deseq_objects.R
#
# Exports the already-fitted DESeq2 objects (meta-analysis + all 3 temporal
# cell types) so a collaborator can run their OWN custom contrasts (e.g.
# OIS vs SIPS directly, or any other pairwise comparison) using
# DESeq2::results() on the pre-fit dispersions/coefficients, rather than
# refitting from scratch -- and a plain-CSV normalized-counts + metadata
# export for anyone who'd rather not load DESeq2 at all.

source("R/config.R")
suppressPackageStartupMessages(library(DESeq2))

EXPORT_DIR <- Sys.getenv("EXPORT_DIR", file.path(Sys.getenv("HOME"), "Downloads", "for cyril"))
META_OUT <- file.path(EXPORT_DIR, "meta_analysis")
TEMP_OUT <- file.path(EXPORT_DIR, "temporal")
dir.create(META_OUT, recursive = TRUE, showWarnings = FALSE)
dir.create(TEMP_OUT, recursive = TRUE, showWarnings = FALSE)

export_one <- function(deseq_rds_path, out_dir, label) {
  cat(sprintf("== %s ==\n", label))
  obj <- readRDS(deseq_rds_path)
  dds <- obj$deseq_obj

  # 1. The fitted DESeqDataSet itself -- has dispersions/coefficients already
  #    estimated, so results(dds, contrast=...) works directly without refitting.
  saveRDS(dds, file.path(out_dir, paste0(label, "_deseq_object.rds")))

  # 2. Plain CSV exports for anyone who doesn't want to touch DESeq2 directly.
  norm_counts <- counts(dds, normalized = TRUE)
  write.csv(data.frame(gene = rownames(norm_counts), norm_counts, check.names = FALSE),
            file.path(out_dir, paste0(label, "_normalized_counts.csv")), row.names = FALSE)

  meta <- as.data.frame(colData(dds))
  write.csv(data.frame(sample_id = rownames(meta), meta), file.path(out_dir, paste0(label, "_sample_metadata.csv")),
            row.names = FALSE)

  cat(sprintf("  dds: %d genes x %d samples, design %s\n",
              nrow(dds), ncol(dds), paste(deparse(design(dds)), collapse = "")))
  cat(sprintf("  available coefficients/levels for custom contrasts: %s\n",
              paste(resultsNames(dds), collapse = ", ")))
  cat(sprintf("  -> %s_deseq_object.rds / _normalized_counts.csv / _sample_metadata.csv\n\n", label))
}

export_one(file.path(RERUN_DIR, "cs_cq_all_study_processed_deseq2.rds"), META_OUT, "cs_cq")

for (ct in c("Fibroblast", "Keratinocyte", "Melanocyte")) {
  export_one(file.path(RERUN_DIR, "ERP021140", ct, "deseq.rds"), TEMP_OUT, ct)
}

cat(sprintf("Done. Exported to %s\n", EXPORT_DIR))
