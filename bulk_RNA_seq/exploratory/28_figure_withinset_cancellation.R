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
# (p_emp < 0.05, the comparisons 2.2.4 counts as beating it) and labelled where it
# survives Benjamini-Hochberg within its family (* FDR < 5%, arrow FDR < 10%). Reading the two together
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
#         withinset_cancellation_by_null.csv     share of movement left after the up and
#             down contributions offset, and share of genes pushing the same way as the
#             set, for sets that beat their matched null vs the rest, at three thresholds

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
         any_fdr = fdr05 | fdr10,
         beats_null = !is.na(p_emp) & p_emp < 0.05)
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
    geom_point(aes(x = set_net, fill = set_net > 0, colour = beats_null,
                   size = beats_null), shape = 21, stroke = 0.7) +
    geom_text(aes(x = set_net, label = ifelse(fdr05, "*", ifelse(fdr10, "\u2191", ""))),
              hjust = -0.4, vjust = 0.35, size = 4.4) +
    facet_wrap(~condition, nrow = if (length(levels_) > 5) 3 else 1) +
    scale_fill_manual(values = c(`TRUE` = "#B2182B", `FALSE` = "#2166AC"),
                      guide = "none") +
    scale_colour_manual(values = c(`TRUE` = "black", `FALSE` = "grey60"),
                        labels = c(`TRUE` = "beats matched null, p < 0.05",
                                   `FALSE` = "does not beat it"),
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
  "black outline: beats the weight-matched null at p < 0.05;   Benjamini-Hochberg within this family of ", nfam,
  " tests:   * FDR < 5%,   \u2191 FDR < 10%")

# One figure, both arms as facet columns/rows, rather than two files: the arms
# are the same measurement on two datasets and are read against each other.
# TWO FIGURES, not one. The merged version needed 20 x 26 inches to fit 14 facets
# of 50 gene sets and was unusable at page size, so each dataset gets its own
# file at a size that fits its own number of groups. BH is applied within each
# dataset's own family (250 and 450 tests, see above).
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

# ---- sets that beat their matched null vs the rest --------------------------
# 2.2.4 compares the two groups on (i) the share of a set's gene movement left
# after its up and down contributions offset and (ii) the share of its genes
# pushing the same way as the set as a whole. Written at the threshold the text
# uses (p < 0.05, 74 of 450 time-course comparisons) and at both FDR levels, so
# the reader can see the comparison does not depend on the threshold.
by_null <- bind_rows(lapply(list(
    list(crit = "p_emp < 0.05", f = function(x) x$beats_null),
    list(crit = "BH q < 0.10",  f = function(x) x$any_fdr),
    list(crit = "BH q < 0.05",  f = function(x) x$fdr05)), function(k) {
  d %>% mutate(beats = k$f(d)) %>% group_by(arm, beats) %>%
    summarise(criterion = k$crit, n = n(),
              left_after_offset_median_pct = 100 * median(surviving, na.rm = TRUE),
              left_after_offset_q25_pct = 100 * quantile(surviving, .25, na.rm = TRUE),
              left_after_offset_q75_pct = 100 * quantile(surviving, .75, na.rm = TRUE),
              genes_same_way_as_set_median_pct = 100 * median(frac_genes_with_set_net),
              .groups = "drop")
  })) %>%
  mutate(group = ifelse(beats, "beats matched null", "rest")) %>%
  select(arm, criterion, group, n, everything(), -beats) %>%
  arrange(arm, criterion, desc(group == "beats matched null"))
write.csv(by_null, file.path(RERUN_DIR, "withinset_cancellation_by_null.csv"), row.names = FALSE)
cat("\nsets that beat their matched null vs the rest:\n")
print(by_null %>% mutate(across(where(is.numeric) & !n, ~round(.x, 1))), n = Inf, width = Inf)
