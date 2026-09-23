# 28_figure_withinset_cancellation.R
#
# Per-GENE-SET cancellation, for both arms, with the matched-null result marked.
#
# WHY. exploratory/24 sums over every measured clock gene and gives one bar per
# group. That shows a group's net is a small residual but says nothing about which
# sets it applies to, and the set-level sections are about individual sets. This
# plots the quantity those sections rest on: for each set, the summed positive and
# negative gene contributions and the net surviving between them.
#
# WHAT THE MARKING ADDS. A wide grey bar with a net near zero is a set whose
# reported contribution is a small difference of large opposing movements - but
# that on its own does not say whether the contribution is more than expected.
# Points are therefore outlined where the set beats its weight-matched random null
# (p_emp < 0.05) and labelled where it clears the threshold the text discusses
# (p < 0.004 cross-sectional, the eight strongest of 250). Reading the two together
# is the point: a set can survive heavy cancellation and still exceed its expected
# share, and a set with little cancellation can still be unremarkable.
#
# EXPRESSED AS A SURVIVING SHARE, NOT A FOLD-RATIO. "6x cancellation" is total
# movement over the net with both directions combined, so it is easily misread as
# per-side, and it is unstable: it reaches 719x for one set whose net is near zero,
# which is a property of the denominator. The share surviving, 1/ratio, is bounded
# and directly interpretable - a median of 17% within sets against about 5%
# transcriptome-wide.
#
# CLOCK: mortality only, as everywhere in the set-level analysis.
#
# Output: figure_withinset_cancellation.png            (arrest conditions)
#         figure_withinset_cancellation_temporal.png   (irradiation time course)
#         withinset_cancellation_summary.csv

source("R/config.R")
suppressPackageStartupMessages({ library(dplyr); library(ggplot2) })

COND <- c("CICQ", "SSCQ", "RS", "SIPS", "OIS")
TGROUP <- c(outer(c("Fibroblast", "Keratinocyte", "Melanocyte"),
                  c("4_days", "10_days", "20_days"), paste, sep = "_"))

d <- read.csv(file.path(RERUN_DIR, "withinset_cancellation.csv")) %>%
  mutate(set = gsub("^HALLMARK ", "", pathway), surviving = 1 / cancellation_ratio)
cat(sprintf("%d set x group rows (%s)\n", nrow(d),
            paste(names(table(d$arm)), table(d$arm), collapse = ", ")))

# ---- the matched-null result, which is the test the sections report ---------
ys <- read.csv(file.path(RERUN_DIR, "pathway_specificity_yardstick_mortality.csv")) %>%
  transmute(condition = label, pathway, p_emp)
yt <- read.csv(file.path(RERUN_DIR,
        "pathway_specificity_yardstick_mortality_temporal.csv")) %>%
  transmute(condition = label, pathway, p_emp)
d <- d %>% left_join(bind_rows(ys, yt), by = c("condition", "pathway"))
cat(sprintf("matched-null p joined for %d of %d rows\n", sum(!is.na(d$p_emp)), nrow(d)))

# FDR, NOT A RAW THRESHOLD (2026-08-27). This is the one analysis in the project
# where the simulation p-values are corrected, because it is a discovery family -
# 50 sets x 5 conditions, or x 9 groups, asking which sets are specific - so which
# ones we name has to carry a stated error rate. Benjamini-Hochberg within each
# arm's own family. Two levels are marked because the section is
# hypothesis-generating: FDR 10% is what it reports, FDR 5% is the conventional
# bar, and a reader should be able to see which results depend on the difference.
bh <- function(p) p.adjust(p, method = "BH")
# BH within each dataset's own family: 250 comparisons for the arrest conditions
# (50 sets x 5 conditions) and 450 for the time course (50 x 9 groups). The two
# are separate analyses of separate data reported in separate sections, so each
# carries its own error rate rather than being penalised for the other's tests.
d <- d %>% group_by(arm) %>% mutate(q = bh(p_emp)) %>% ungroup() %>%
  mutate(fdr05 = !is.na(q) & q < 0.05,
         fdr10 = !is.na(q) & q < 0.10 & !fdr05,
         any_fdr = fdr05 | fdr10)
cat("\nBH within each dataset's own family:\n")
print(d %>% group_by(arm) %>%
        summarise(n_tests = n(), q_lt_05 = sum(fdr05),
                  q_lt_10_only = sum(fdr10), .groups = "drop"))

plot_arm <- function(rows, levels_, file, width, height, xlab_note) {
  x <- rows %>% filter(condition %in% levels_) %>%
    # facet strips read "Fibroblast 4 days", not "Fibroblast_4_days"; the levels
    # themselves keep their underscores because they match the CSV keys.
    mutate(condition = factor(gsub("_", " ", condition),
                              levels = gsub("_", " ", levels_)))
  ord <- x %>% group_by(set) %>% summarise(m = median(surviving), .groups = "drop") %>%
    arrange(m)
  x$set <- factor(x$set, levels = ord$set)
  p <- ggplot(x, aes(y = set)) +
    geom_segment(aes(x = set_sum_negative, xend = set_sum_positive, yend = set),
                 colour = "grey75", linewidth = 1.4) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40") +
    geom_point(aes(x = set_net, fill = set_net > 0, colour = any_fdr,
                   size = any_fdr), shape = 21, stroke = 0.7) +
    geom_text(aes(x = set_net, label = ifelse(fdr05, "*", ifelse(fdr10, "\u2191", ""))),
              hjust = -0.4, vjust = 0.35, size = 4.4) +
    facet_wrap(~condition, nrow = if (length(levels_) > 5) 3 else 1) +
    scale_fill_manual(values = c(`TRUE` = "#B2182B", `FALSE` = "#2166AC"),
                      guide = "none") +
    scale_colour_manual(values = c(`TRUE` = "black", `FALSE` = "grey60"),
                        labels = c(`TRUE` = "FDR < 10%", `FALSE` = "not significant"),
                        name = NULL) +
    scale_size_manual(values = c(`TRUE` = 2.4, `FALSE` = 1.5), guide = "none") +
    theme_bw(base_size = 13) +
    theme(axis.text.y = element_text(size = 7.5), legend.position = "top",
          panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(),
          strip.background = element_blank()) +
    labs(x = xlab_note, y = NULL)
  ggsave(file.path(RERUN_DIR, file), p, width = width, height = height, dpi = 300)
  cat(sprintf("Saved -> %s\n", file))
}

XL <- function(nfam) paste0(
  "Summed gene contributions within each set (grey) and the net surviving (point), mortality-clock units\n",
  "sets ordered by share of movement surviving, least at the bottom\n",
  "against the weight-matched null, Benjamini-Hochberg within this family of ", nfam,
  " tests:   * FDR < 5%,   \u2191 FDR < 10%")

# One figure, both arms as facet columns/rows, rather than two files: the arms
# are the same measurement on two datasets and are read against each other.
# TWO FIGURES, not one. The merged version needed 20 x 26 inches to fit 14 facets
# of 50 gene sets and was unusable at page size, so each dataset gets its own
# file at a size that fits its own number of groups. The BH family is still
# POOLED across both (700 tests, see above), so a set named in either section
# carries the same error rate - only the rendering is split.
plot_arm(d %>% filter(arm == "cross_sectional"), COND,
         "figure_withinset_cancellation.png", 15, 11, XL(250))
plot_arm(d %>% filter(arm == "temporal"), TGROUP,
         "figure_withinset_cancellation_temporal.png", 15, 20, XL(450))

s <- d %>% group_by(arm) %>% summarise(
  n = n(), ratio_median = median(cancellation_ratio, na.rm = TRUE),
  ratio_max = max(cancellation_ratio, na.rm = TRUE),
  surviving_median = median(surviving, na.rm = TRUE),
  surviving_q25 = quantile(surviving, .25, na.rm = TRUE),
  surviving_q75 = quantile(surviving, .75, na.rm = TRUE),
  frac_genes_median = median(frac_genes_with_set_net),
  n_fdr10 = sum(any_fdr), .groups = "drop")
write.csv(s, file.path(RERUN_DIR, "withinset_cancellation_summary.csv"), row.names = FALSE)
print(s %>% mutate(across(c(starts_with("surviving"), starts_with("frac")),
                          ~round(100 * .x, 1)),
                   across(starts_with("ratio"), ~round(.x, 1))))

cat("\ndoes surviving share relate to beating the matched null?\n")
print(d %>% group_by(arm, any_fdr) %>%
        summarise(n = n(), surviving_pct = round(100 * median(surviving, na.rm = TRUE), 1),
                  .groups = "drop"))
