# 43_gene_annotation_coverage.R
#
# 2.1.5.3 (iii) says over 60% of the tAge movement sits in genes that no Hallmark set
# contains. An earlier draft of that sentence said instead that over 60% of the clock's
# "biological functionality" is unaccounted for. That is a different claim, and this
# script is what would be needed to make it: Hallmark is 50 coarse sets, and a gene's
# absence from all of them is no evidence that its function is unknown.
#
# THE TEST. For every clock feature, is there any GO biological-process term? Then the
# same two quantities as 42, on GO membership instead of Hallmark:
#   share of genes annotated
#   share of the total |contribution| carried by annotated genes, per condition
# If nearly every clock gene has a GO BP term, the "functionality" wording is wrong and
# the Hallmark wording in the manuscript is the only defensible one. If a real fraction
# has none, the stronger sentence can be written, with GO named as the source.
#
# Inputs: rerun_outputs/gene_contribution_table.csv (from 42_coverage_threshold_free.py)
# Output: rerun_outputs/gene_annotation_coverage.csv
source("R/config.R")
suppressMessages({library(dplyr); library(AnnotationDbi); library(org.Mm.eg.db); library(GO.db)})

g <- read.csv(file.path(RERUN_DIR, "gene_contribution_table.csv"))
g$gene <- as.character(g$gene)
for (cc in c("in_hallmark", "measured_in_arrest")) g[[cc]] <- as.logical(g[[cc]])
cat(sprintf("clock features: %d | in a Hallmark set: %d (%.1f%%)\n",
            nrow(g), sum(g$in_hallmark), 100 * mean(g$in_hallmark)))

# GO biological process, any evidence code, mouse Entrez ids
bp <- suppressMessages(AnnotationDbi::select(org.Mm.eg.db, keys = g$gene,
                                             keytype = "ENTREZID",
                                             columns = c("GOALL", "ONTOLOGYALL")))
bp <- bp[!is.na(bp$GOALL) & bp$ONTOLOGYALL == "BP", ]
n_bp <- tapply(bp$GOALL, bp$ENTREZID, function(x) length(unique(x)))
g$n_go_bp <- as.integer(ifelse(is.na(n_bp[g$gene]), 0, n_bp[g$gene]))
g$has_go_bp <- g$n_go_bp > 0
# a gene with only very general terms is not meaningfully characterised either; the
# median clock gene's count is reported so the reader can see where the cut would fall
cat(sprintf("with >= 1 GO BP term: %d (%.1f%%) | median terms per annotated gene: %d\n",
            sum(g$has_go_bp), 100 * mean(g$has_go_bp),
            as.integer(median(g$n_go_bp[g$has_go_bp]))))
cat(sprintf("in NO Hallmark set but WITH a GO BP term: %d of the %d unannotated by Hallmark (%.1f%%)\n",
            sum(!g$in_hallmark & g$has_go_bp), sum(!g$in_hallmark),
            100 * sum(!g$in_hallmark & g$has_go_bp) / sum(!g$in_hallmark)))

cols <- grep("^abs_contrib_", names(g), value = TRUE)
rows <- lapply(cols, function(cc) {
  v <- g[[cc]]
  keep <- !is.na(v)
  tot <- sum(v[keep])
  data.frame(group = sub("^abs_contrib_", "", cc),
             n_measured = sum(keep),
             pct_genes_hallmark = 100 * mean(g$in_hallmark[keep]),
             pct_genes_go_bp = 100 * mean(g$has_go_bp[keep]),
             pct_movement_hallmark = 100 * sum(v[keep & g$in_hallmark]) / tot,
             pct_movement_go_bp = 100 * sum(v[keep & g$has_go_bp]) / tot)
})
D <- bind_rows(rows)
print(D %>% mutate(across(where(is.numeric), ~round(.x, 1))), row.names = FALSE)
write.csv(D, file.path(RERUN_DIR, "gene_annotation_coverage.csv"), row.names = FALSE)
write.csv(g[, c("gene", "coef", "in_hallmark", "n_go_bp", "has_go_bp")],
          file.path(RERUN_DIR, "gene_go_annotation.csv"), row.names = FALSE)

# gene -> GO BP term mapping for 44_unaccounted_enrichment.py. GOALL is used, so
# ancestor terms are included: that is correct for enrichment, and term size is
# filtered downstream rather than here.
map <- unique(bp[, c("ENTREZID", "GOALL")])
names(map) <- c("gene", "go_id")
tn <- suppressMessages(AnnotationDbi::select(GO.db, keys = unique(map$go_id),
                                             keytype = "GOID", columns = "TERM"))
map$term <- tn$TERM[match(map$go_id, tn$GOID)]
# depth in the BP DAG, as the number of ancestors: a broad term like "metabolic
# process" has few, a specific one has many. Lets breadth be reported rather than
# assumed when the enriched terms are read.
anc <- as.list(GO.db::GOBPANCESTOR)
map$n_ancestors <- lengths(anc[map$go_id])
write.csv(map, file.path(RERUN_DIR, "gene_go_terms.csv"), row.names = FALSE)
cat(sprintf("GO BP mapping: %d gene-term pairs over %d terms\n",
            nrow(map), length(unique(map$go_id))))
cat("\nSaved -> gene_annotation_coverage.csv, gene_go_annotation.csv\n")
