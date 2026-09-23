# 06_compile_gene_tables_for_export.R
#
# Compiles, per senescence/quiescence condition (meta-analysis) AND per
# cell-type x timepoint (temporal), a full gene table (gene, log2FC, pvalue,
# padj, significant) plus the background gene universe that condition's/
# contrast's DESeq2 test was actually run against (i.e. every gene that
# passed that contrast's independent filtering, significant or not --
# reading the background off the DEG table itself, not re-deriving it, so
# it's guaranteed to match exactly what was tested).
#
# Output goes to EXPORT_DIR (default ~/Downloads/for cyril), in two
# subfolders: meta_analysis/ and temporal/.

source("R/config.R")

EXPORT_DIR <- Sys.getenv("EXPORT_DIR", file.path(Sys.getenv("HOME"), "Downloads", "for cyril"))
META_OUT <- file.path(EXPORT_DIR, "meta_analysis")
TEMP_OUT <- file.path(EXPORT_DIR, "temporal")
dir.create(META_OUT, recursive = TRUE, showWarnings = FALSE)
dir.create(TEMP_OUT, recursive = TRUE, showWarnings = FALSE)

clean_cols <- function(df, gene_col) {
  out <- data.frame(
    gene = df[[gene_col]],
    log2FoldChange = df$log2FoldChange,
    pvalue = df$pvalue,
    padj = df$padj,
    significant = ifelse(df$sig == "y", TRUE, FALSE)
  )
  out[order(out$padj), ]
}

# ── Meta-analysis: senescence (RS, OIS, SIPS) + quiescence (CICQ, SSCQ) ──────
cat("== Meta-analysis conditions ==\n")
arrest_degs <- read.csv(file.path(RERUN_DIR, "arrest_degs_final_RERUN.csv"))

condition_labels <- c(
  "Replicative_CS"       = "senescence_RS",
  "Oncogene_induced_CS"  = "senescence_OIS",
  "Stress_induced_CS"    = "senescence_SIPS",
  "Contact_inhibited_CQ" = "quiescence_CICQ",
  "Serum_starved_CQ"     = "quiescence_SSCQ"
)

for (grp in names(condition_labels)) {
  label <- condition_labels[[grp]]
  sub <- arrest_degs[arrest_degs$group_1 == grp, ]
  # each gene appears twice per group_1 in this table (once per direction_1
  # label, up/down) -- de-duplicate on gene, keeping the DESeq2 stats (which
  # are identical across the direction_1 duplicate rows; direction_1 is just
  # a label, not a different test)
  sub <- sub[!duplicated(sub$gene), ]

  gene_table <- clean_cols(sub, "gene")
  write.csv(gene_table, file.path(META_OUT, paste0(label, "_genes.csv")), row.names = FALSE)

  background <- data.frame(gene = sort(unique(sub$gene)))
  write.csv(background, file.path(META_OUT, paste0(label, "_background.csv")), row.names = FALSE)

  cat(sprintf("  %-18s -> %d genes tested, %d significant, background n=%d\n",
              label, nrow(gene_table), sum(gene_table$significant), nrow(background)))
}

# ── Temporal: 3 cell types x 3 timepoints (vs that cell type's own 'none') ──
cat("\n== Temporal conditions (per cell type x timepoint) ==\n")
cell_types <- c("Fibroblast", "Keratinocyte", "Melanocyte")
timepoints <- c("4_days", "10_days", "20_days")

for (ct in cell_types) {
  degs_path <- file.path(RERUN_DIR, "ERP021140", ct, "degs.csv")
  if (!file.exists(degs_path)) {
    message("  MISSING: ", degs_path, " -- skipping ", ct)
    next
  }
  degs <- read.csv(degs_path)
  for (tp in timepoints) {
    sub <- degs[degs$group_1 == tp, ]
    if (nrow(sub) == 0) next
    sub <- sub[!duplicated(sub$ensembl), ]

    label <- paste0(ct, "_", tp)
    gene_table <- clean_cols(sub, "ensembl")
    write.csv(gene_table, file.path(TEMP_OUT, paste0(label, "_genes.csv")), row.names = FALSE)

    background <- data.frame(gene = sort(unique(sub$ensembl)))
    write.csv(background, file.path(TEMP_OUT, paste0(label, "_background.csv")), row.names = FALSE)

    cat(sprintf("  %-24s -> %d genes tested, %d significant, background n=%d\n",
                label, nrow(gene_table), sum(gene_table$significant), nrow(background)))
  }
}

cat(sprintf("\nDone. Exported to %s\n", EXPORT_DIR))
