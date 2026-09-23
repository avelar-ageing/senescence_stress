# pca_panel_figures.R -- consolidated before/after PCA figures for SI Figure 1
# (27_si_figure1_pca.R) and SI Figure 16 (28_si_figure16_pca.R).
#
# WHY THIS REPLACES ONE FILE PER COVARIATE PER ARM. The 2024 layout wrote each
# covariate separately because a single sheet of all six overlapped its own
# legends - cell_line carries 16 levels and study 34. That constraint is about
# LEGEND space, not panels, so it does not apply to the before/after split: the
# two arms plot the same samples with the same colour levels, so one legend
# serves both. Pairing the arms as facet columns therefore halves the file count
# for free. The low-cardinality covariates (<= 8 levels) additionally fit
# together as stacked rows of one figure; cell_line and study keep their own,
# which is what the legend constraint actually requires.
#
# Result per script: 3 files (one multi-covariate figure + one each for the two
# wide covariates) in place of 12 per-covariate panels and 2 contact sheets.
suppressMessages({library(ggplot2); library(patchwork)})

ARM_LABELS <- c(no_correction = "Before batch correction",
                correction    = "After batch correction")

# One covariate, both arms side by side, sharing a legend.
pca_arm_pair <- function(co, covariate, point_size = 2.2, legend_rows = NULL) {
  d <- co[!is.na(co[[covariate]]), ]
  d$arm <- factor(unname(ARM_LABELS[as.character(d$arm)]), levels = unname(ARM_LABELS))
  nlev <- length(unique(d[[covariate]]))
  if (is.null(legend_rows)) legend_rows <- if (nlev > 20) 9 else if (nlev > 8) 6 else 3
  vl <- tapply(d$pc1_var, d$arm, function(x) x[1])
  v2 <- tapply(d$pc2_var, d$arm, function(x) x[1])
  ggplot(d, aes(x = PC1, y = PC2, colour = .data[[covariate]])) +
    geom_point(size = point_size, alpha = 0.9) +
    facet_wrap(~arm, nrow = 1, scales = "free") +
    guides(colour = guide_legend(nrow = legend_rows, byrow = TRUE,
                                 override.aes = list(size = 3))) +
    theme_bw(base_size = 12) +
    theme(legend.position = "right", legend.title = element_text(face = "bold"),
          legend.text = element_text(size = 8), strip.background = element_blank(),
          strip.text = element_text(size = 11), panel.grid.minor = element_blank()) +
    labs(colour = covariate,
         x = sprintf("PC1 (%.1f%% / %.1f%%)", vl[1], vl[2]),
         y = sprintf("PC2 (%.1f%% / %.1f%%)", v2[1], v2[2]))
}

# Everything: narrow covariates stacked into one figure, wide ones separate.
save_pca_figures <- function(co, covariates, prefix, save_dir, wide_cutoff = 8) {
  present <- covariates[covariates %in% names(co)]
  nlev <- vapply(present, function(p) length(unique(co[[p]][!is.na(co[[p]])])), integer(1))
  present <- present[nlev >= 2]; nlev <- nlev[present]
  wide   <- present[nlev >  wide_cutoff]
  narrow <- present[nlev <= wide_cutoff]
  out <- character(0)
  if (length(narrow)) {
    ps <- lapply(narrow, function(p) pca_arm_pair(co, p))
    fig <- Reduce(`/`, ps) + plot_annotation(tag_levels = "a") &
      theme(plot.tag = element_text(face = "bold", size = 15))
    f <- file.path(save_dir, paste0(prefix, "_covariates.png"))
    ggsave(f, fig, width = 13, height = 3.6 * length(narrow), dpi = 300, limitsize = FALSE)
    cat(sprintf("  %s  (%s, before/after as columns)\n", basename(f),
                paste(narrow, collapse = ", ")))
    out <- c(out, f)
  }
  for (p in wide) {
    f <- file.path(save_dir, paste0(prefix, "_", p, ".png"))
    ggsave(f, pca_arm_pair(co, p), width = 14,
           height = if (nlev[[p]] > 20) 8.5 else 6.5, dpi = 300, limitsize = FALSE)
    cat(sprintf("  %s  (%d levels, own figure for legend space)\n", basename(f), nlev[[p]]))
    out <- c(out, f)
  }
  invisible(out)
}
