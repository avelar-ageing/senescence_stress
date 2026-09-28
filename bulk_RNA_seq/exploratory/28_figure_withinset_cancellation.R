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
#             set, for sets that beat their matched null vs the rest, at three thresholds;
#             also the share that would be left if every gene moved by the same amount (2f - 1)
#         withinset_cancellation_by_null_tests.csv   are the two groups different? rank-sum
#             test, and a permutation that shuffles the "beats its null" label only among
#             the conditions (or groups) of the same gene set, 20,000 times
#         figure_withinset_by_null.png            the two measures, sets that beat their null
#             vs the rest, both datasets

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
              # if every gene moved by the same amount, a set with a share f of its genes
              # pointing its way would keep 2f - 1 of its movement after the offset
              left_if_all_genes_equal_median_pct = 100 * median(2 * frac_genes_with_set_net - 1),
              .groups = "drop")
  })) %>%
  mutate(group = ifelse(beats, "beats matched null", "rest")) %>%
  select(arm, criterion, group, n, everything(), -beats) %>%
  arrange(arm, criterion, desc(group == "beats matched null"))
write.csv(by_null, file.path(RERUN_DIR, "withinset_cancellation_by_null.csv"), row.names = FALSE)
cat("\nsets that beat their matched null vs the rest:\n")
print(by_null %>% mutate(across(where(is.numeric) & !n, ~round(.x, 1))), n = Inf, width = Inf)

# ---- are the two groups different? (threshold used in the text: p < 0.05) -----
# The comparisons are not independent - the same gene set appears in every condition
# or group - so besides the rank-sum test the "beats its null" label is shuffled only
# among the conditions of the same set (20,000 times; p = (1 + k)/(B + 1)).
# The left-over share is partly selected for: a set beats its null on its net
# contribution, which is the left-over share times its movement.
set.seed(1); B <- 20000
# Wilcoxon rank-sum, normal approximation with tie and continuity corrections (as
# scipy.stats.mannwhitneyu). Written out because wilcox.test(exact = FALSE) in R 4.6
# returned 0 for the time-course left-over comparison, where p is about 2e-26.
rank_sum_p <- function(a, r) {
  n1 <- length(a); n2 <- length(r); N <- n1 + n2; rk <- rank(c(a, r))
  W <- sum(rk[seq_len(n1)]) - n1 * (n1 + 1) / 2; mu <- n1 * n2 / 2
  t <- table(rk); sig <- sqrt(n1 * n2 / 12 * ((N + 1) - sum(t^3 - t) / (N * (N - 1))))
  z <- (abs(W - mu) - 0.5) / sig
  2 * pnorm(-z)
}
perm_p <- function(x, col) {
  v <- x[[col]]; lab <- x$beats_null; by_set <- split(seq_len(nrow(x)), x$pathway)
  obs <- median(v[lab]) - median(v[!lab]); k <- 0L
  for (b in seq_len(B)) {
    L <- lab
    for (ii in by_set) L[ii] <- L[ii][sample.int(length(ii))]
    k <- k + (abs(median(v[L]) - median(v[!L])) >= abs(obs))
  }
  (1 + k) / (B + 1)
}
tests <- bind_rows(lapply(split(d, d$arm), function(x) {
  bind_rows(lapply(c(genes_same_way_as_set = "frac_genes_with_set_net",
                     left_after_offset = "surviving"), function(col) {
    a <- x[[col]][x$beats_null]; r <- x[[col]][!x$beats_null]
    data.frame(arm = x$arm[1], measure = col, n_beats = length(a), n_rest = length(r),
               median_beats_pct = 100 * median(a), median_rest_pct = 100 * median(r),
               p_rank_sum = rank_sum_p(a, r),
               p_within_set_permutation = perm_p(x, col), n_permutations = B)
  }), .id = "label")
}))
tests$measure <- tests$label; tests$label <- NULL
write.csv(tests, file.path(RERUN_DIR, "withinset_cancellation_by_null_tests.csv"), row.names = FALSE)
cat("\ngroup differences (p < 0.05 threshold):\n"); print(tests, digits = 3)

# ---- figure: the two measures, side by side -----------------------------------
fd <- bind_rows(
  d %>% transmute(arm, beats_null, measure = "genes pushing the set's way", value = 100 * frac_genes_with_set_net),
  d %>% transmute(arm, beats_null, measure = "movement left after up and down offset", value = 100 * surviving)) %>%
  mutate(dataset = ifelse(arm == "cross_sectional", "arrest conditions", "irradiation time course"),
         group = factor(ifelse(beats_null, "beats matched null", "rest"), levels = c("rest", "beats matched null")))
eq <- d %>% group_by(arm, beats_null) %>%
  summarise(value = 100 * median(2 * frac_genes_with_set_net - 1), .groups = "drop") %>%
  mutate(dataset = ifelse(arm == "cross_sectional", "arrest conditions", "irradiation time course"),
         group = factor(ifelse(beats_null, "beats matched null", "rest"), levels = c("rest", "beats matched null")),
         measure = "movement left after up and down offset")
pl <- tests %>% mutate(dataset = ifelse(arm == "cross_sectional", "arrest conditions", "irradiation time course"),
                       measure = ifelse(measure == "genes_same_way_as_set", "genes pushing the set's way",
                                        "movement left after up and down offset"),
                       lab = sprintf("rank-sum p = %.1e\nwithin-set shuffle p = %.1e%s", p_rank_sum, p_within_set_permutation,
                                     ifelse(p_within_set_permutation <= 1 / (B + 1) + 1e-12, " (floor)", "")))
pf <- ggplot(fd, aes(x = group, y = value)) +
  geom_boxplot(outlier.shape = NA, width = 0.55, fill = "grey92") +
  geom_jitter(width = 0.15, height = 0, size = 0.7, alpha = 0.45) +
  geom_point(data = eq, aes(x = group, y = value), shape = 95, size = 14, colour = "#B2182B") +
  geom_text(data = pl, aes(x = 1.5, y = Inf, label = lab), vjust = 1.3, size = 3.1, inherit.aes = FALSE) +
  facet_grid(measure ~ dataset, scales = "free_y") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.28))) +
  theme_bw(base_size = 12) + theme(strip.background = element_blank(), panel.grid.minor = element_blank()) +
  labs(x = NULL, y = "% (one point per gene set in one condition or group)",
       caption = "red dash: movement that would be left if every gene in the set moved by the same amount (median)")
ggsave(file.path(RERUN_DIR, "figure_withinset_by_null.png"), pf, width = 9, height = 8, dpi = 300)
cat("Saved -> withinset_cancellation_by_null_tests.csv, figure_withinset_by_null.png\n")
