# 11_extract_paper_module_clocks.R
#
# Extract the PAPER'S OWN module clocks (gene membership + trained coefficients + intercepts) from the
# published supplementary tables of Tyshkovskiy/Gladyshev et al. 2026, Nature.
#
# WHY THIS MATTERS. Our partial-tAge pipeline (exploratory/01-07) decomposes the fitted GLOBAL clock
# over MSigDB Hallmark pathways. The paper does something different: it TRAINS a separate elastic-net
# clock per WGCNA co-expression module (see DISCREPANCY_REPORT/PARTIAL_TAGE_VS_PAPER.md). We had
# assumed the module clocks were unobtainable -- they are absent from Zenodo record 18763485 (72 files,
# all full-transcriptome) and from the tAge R package, whose extdata carries only a module-to-function
# map with sizes, no membership.
#
# They are, however, published as Supplementary Table 5: sheet "(B) Module rodent clocks" and
# "(C) Module multi-species clocks" of MOESM7 give, per module clock, every gene's coefficient plus the
# intercept. That is the full linear model. So the paper's estimator IS reproducible without retraining:
#
#     module_tAge = intercept + sum_i coef_i * z_i    over that module's genes
#
# using the same tAge_preprocessing() output we already export. Genes with a non-zero coefficient in a
# module's column are that module's membership.
#
# WHICH SHEET TO USE. Human fibroblast data is scored with the Multispecies_Multitissue models, so
# sheet (C) is the relevant one: 14 modules plus an "All module genes" composite, for two outcomes
# (Chronological age, Mortality). Sheet (B) is the rodent panel: 23 modules, matching the paper's
# headline module count. Both are the "scaling" variant, i.e. they correspond to scaleddiff -- NO
# yugene module clocks appear to be published, which constrains any module analysis to one
# normalisation while our Hallmark decomposition runs both.
#
# INPUT. Download once (public, ~46 MB):
#   curl -o suppl.zip https://www.ebi.ac.uk/europepmc/webservices/rest/PMC13233323/supplementaryFiles
#   unzip suppl.zip -d suppl
# then point SUPPL_DIR at it.
#
# USAGE: SUPPL_DIR=/path/to/suppl Rscript exploratory/11_extract_paper_module_clocks.R
# OUT:   rerun_outputs/partial_tage/paper_module_clocks.csv   tidy: sheet, outcome, module, annotation,
#                                                             entrez_id, gene_symbol, coefficient
#        rerun_outputs/partial_tage/paper_module_sizes.csv    genes per module clock, vs the paper's
#                                                             stated module sizes
suppressPackageStartupMessages({ library(readxl); library(dplyr) })
SUPPL_DIR <- Sys.getenv('SUPPL_DIR', unset = NA)
if (is.na(SUPPL_DIR)) stop('set SUPPL_DIR to the unzipped supplementary directory')
XL  <- file.path(SUPPL_DIR, '41586_2026_10542_MOESM7_ESM.xlsx')
DICT<- file.path(SUPPL_DIR, '41586_2026_10542_MOESM8_ESM.xlsx')
OUT <- '/home/ro/APFS_copy/root/Backup/Documents/modules/systems_analysis_arrest/bulk_RNA_seq/rerun_outputs/partial_tage'
stopifnot(file.exists(XL))

# Sheet layout (verified): row1 = clock outcome, row2 = module ID, row3 = module annotation,
# row4 = column headers for the two ID columns, row5 = intercept, rows 6+ = per-gene coefficients.
parse_sheet <- function(sheet) {
  y <- suppressMessages(read_excel(XL, sheet = sheet, col_names = FALSE, .name_repair = 'minimal'))
  outcome <- as.character(unlist(y[1, ])); module <- as.character(unlist(y[2, ]))
  annot   <- as.character(unlist(y[3, ]))
  # carry the outcome label across its merged block
  for (i in seq_along(outcome)) if (is.na(outcome[i]) && i > 1) outcome[i] <- outcome[i-1]
  ids <- data.frame(entrez_id = as.character(unlist(y[-(1:5), 1])),
                    gene_symbol = as.character(unlist(y[-(1:5), 2])), stringsAsFactors = FALSE)
  out <- list()
  for (j in 3:ncol(y)) {
    if (is.na(module[j])) next
    co <- suppressWarnings(as.numeric(unlist(y[-(1:5), j])))
    ic <- suppressWarnings(as.numeric(y[[j]][5]))
    keep <- which(!is.na(co) & co != 0)
    if (!length(keep)) next
    out[[length(out)+1]] <- data.frame(sheet = sheet, outcome = outcome[j], module = module[j],
      annotation = annot[j], intercept = ic, ids[keep, ], coefficient = co[keep],
      stringsAsFactors = FALSE)
  }
  bind_rows(out)
}

res <- bind_rows(lapply(c('(B) Module rodent clocks','(C) Module multi-species clocks'), parse_sheet))
res$panel <- ifelse(grepl('rodent', res$sheet), 'rodent', 'multispecies')
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
write.csv(res, file.path(OUT, 'paper_module_clocks.csv'), row.names = FALSE)

sizes <- res %>% group_by(panel, outcome, module, annotation) %>%
  summarise(n_genes = n(), intercept = first(intercept), .groups = 'drop') %>% arrange(panel, outcome, desc(n_genes))
# cross-check against the published module dictionary (annotation + stated module size)
dict <- bind_rows(
  suppressMessages(read_excel(DICT, sheet = '(B) Rodent module dictionary')) %>% mutate(panel='rodent'),
  suppressMessages(read_excel(DICT, sheet = '(D) Multi-species module dict')) %>% mutate(panel='multispecies'))
names(dict)[1:3] <- c('module','annotation_dict','module_size_dict')
sizes <- left_join(sizes, dict[,c('panel','module','module_size_dict')], by = c('panel','module'))
write.csv(sizes, file.path(OUT, 'paper_module_sizes.csv'), row.names = FALSE)

cat('=== module clocks recovered ===\n')
print(as.data.frame(sizes %>% group_by(panel, outcome) %>%
  summarise(n_modules = n(), genes_min = min(n_genes), genes_max = max(n_genes),
            genes_total = sum(n_genes), .groups='drop')), row.names = FALSE)
cat('\n=== multispecies chronological-age modules (the panel relevant to human data) ===\n')
print(as.data.frame(sizes %>% filter(panel=='multispecies', grepl('hronological', outcome)) %>%
  select(module, annotation, n_genes, module_size_dict)), row.names = FALSE)
cat(sprintf('\nwrote %s/paper_module_clocks.csv (%d rows) and paper_module_sizes.csv\n', OUT, nrow(res)))
