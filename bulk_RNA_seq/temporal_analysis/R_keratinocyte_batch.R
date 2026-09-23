# R_keratinocyte_batch.R -- single source of truth for the keratinocyte batch.
#
# E-MTAB-5403's keratinocyte data was processed by two researchers (confirmed
# with the study's corresponding author; see main PDF p.21/p.37-38). The DEG
# arm already accounts for it: temporal_analysis/02_run_time_analysis.R calls
# run_one("Keratinocyte", fix_batch = TRUE, batch_1 = <these 12 ids>), and
# 03_temporal_overlaps.R repeats the same list. Fibroblasts and melanocytes are
# run with fix_batch = FALSE, i.e. Y = Time + e.
#
# The tAge arm did NOT account for it until now. The batch is balanced 3/3
# across all four timepoints, so it does not bias the timepoint estimates, but
# it adds a large between-batch offset to tAge (+5.8 YuGene, +17.1 scaled,
# +0.31 mortality units) that inflates the within-timepoint spread and costs
# power. Removing it is the two-step equivalent of the paper's
# Y = Time + Batch + e: subtract each batch's mean, then test timepoint as
# before. Point estimates are unchanged by construction; only the spread falls.
KERATINOCYTE_BATCH_1 <- c(
  "ERR1805235", "ERR1805236", "ERR1805238", "ERR1805230", "ERR1805231", "ERR1805224",
  "ERR1805223", "ERR1805222", "ERR1805239", "ERR1805240", "ERR1805241", "ERR1805229"
)

# Batch-centre the named tAge columns for keratinocytes only. Other cell types
# are returned untouched, matching fix_batch = FALSE in the DEG arm.
batch_centre_keratinocytes <- function(df, value_cols, id_col = "external_id",
                                       cell_col = "cell_type") {
  k <- df[[cell_col]] == "Keratinocyte"
  if (!any(k)) return(df)
  b <- ifelse(df[[id_col]][k] %in% KERATINOCYTE_BATCH_1, "b1", "b2")
  if (length(unique(b)) < 2L)
    stop("keratinocyte batch labels did not split the samples - check id_col")
  for (v in value_cols) {
    x <- df[[v]][k]
    df[[v]][k] <- x - ave(x, b) + mean(x)   # remove batch offset, keep the grand mean
  }
  attr(df, "keratinocyte_batch_centred") <- value_cols
  df
}
