# 21_gse175533_tage.R
#
# Step 2 of the GSE175533 replication (see meta_analysis/20 for why this dataset
# and what is in it). Pushes the exported matrix through EXACTLY the same
# tAge_preprocessing call used for the meta-analysis and the temporal arm in
# exploratory/01, so the resulting scaled_diff / yugene_diff matrices are
# comparable in construction to everything else in 2.1.5.
#
# CONTROL GROUP. parental_dividing - parental WI-38 at PD 20-37, still dividing.
# This is the closest analogue of the "Proliferating" control used throughout, and
# it is what makes the resulting tAge values differences from a primary
# proliferating baseline rather than from an arbitrary group. The control choice
# fixes the zero point only; every contrast reported downstream is a difference
# between two groups, so it is invariant to that choice.
#
# Output: rerun_outputs/gse175533/gse175533_{scaled_diff,yugene_diff}.csv
#         rerun_outputs/gse175533/gse175533_groups.csv

source("R/config.R")
suppressPackageStartupMessages({
  library(tAge)
  library(Biobase)
})

DIR <- file.path(RERUN_DIR, "gse175533")

E <- read.csv(file.path(DIR, "gse175533_expression.csv"), row.names = 1,
              check.names = FALSE)
S <- read.csv(file.path(DIR, "gse175533_samples.csv"))
stopifnot(identical(colnames(E), S$sample_id))
rownames(S) <- S$sample_id

cat(sprintf("input: %d genes x %d samples\n", nrow(E), ncol(E)))
print(table(S$group))

eset <- ExpressionSet(assayData = as.matrix(E),
                      phenoData = AnnotatedDataFrame(S))

tAge_eset <- suppressWarnings(tAge_preprocessing(
  eset, species = "human", gene_mapping_type = "Gene.Symbol",
  control_group_column = "group", control_group_label = "parental_dividing",
  verbose = TRUE, count_threshold = 10, percent_threshold = 20
))

for (variant in c("scaled_diff", "yugene_diff")) {
  m <- t(exprs(tAge_eset[[variant]]))
  write.csv(data.frame(sample_id = rownames(m), m, check.names = FALSE),
            file.path(DIR, sprintf("gse175533_%s.csv", variant)), row.names = FALSE)
  cat(sprintf("%s: %d samples x %d genes\n", variant, nrow(m), ncol(m)))
}
write.csv(data.frame(sample_id = S$sample_id, group = S$group),
          file.path(DIR, "gse175533_groups.csv"), row.names = FALSE)
cat(sprintf("\nSaved -> %s\n", DIR))
