# 31_figure_cancellation_sides.R
#
# The two results from 29_set_cancellation_structure.py and 30_set_concordance_null.py,
# neither of which had a figure.
#
# PANEL SET A - each direction tested separately (one file per dataset, as with
# 28_figure_withinset_cancellation.R, because 50 sets x 9 groups does not fit one
# sheet). Same layout as that figure: a bar from the summed negative to the summed
# positive contribution, with the net as a point. What is new is WHERE the markers
# sit. 28 marks the NET against its null; here each END of the bar is marked
# against the same weight-matched null, so a set whose genes move hard in both
# directions - invisible to a test of the net - is visible as a bar marked at both
# ends with an unmarked point in the middle. 42 of 250 comparisons in the arrest
# conditions and 88 of 450 in the time course are of that kind.
#
# PANEL B - the cross-condition concordance against its null. For each set, the
# median correlation between conditions of its per-gene contributions, with the
# random-set null (0.13) drawn as a reference line. Only E2F TARGETS and G2M
# CHECKPOINT clear the null in every condition pair.
#
# CLOCK: mortality only, as everywhere in the set-level analysis.
# FDR: q columns are BH within each dataset, computed in the python scripts.
source("R/config.R")
suppressMessages({library(dplyr); library(ggplot2); library(tidyr)})

sides <- read.csv(file.path(RERUN_DIR, "set_cancellation_sides.csv"))
conc  <- read.csv(file.path(RERUN_DIR, "set_concordance_null.csv"))

pretty_cond <- function(x) gsub("_", " ", x)
sides <- sides %>% mutate(set = gsub("^HALLMARK ", "", pathway),
                          cond = pretty_cond(condition),
                          mark_up = q_up < 0.05, mark_dn = q_down < 0.05,
                          mark_net = p_net < 0.05,
                          both = mark_up & mark_dn & !mark_net)

plot_sides <- function(d, file, width, height, nfam) {
  ord <- d %>% group_by(set) %>% summarise(m = median(up - down), .groups = "drop") %>% arrange(m)
  d$set <- factor(d$set, levels = ord$set)
  p <- ggplot(d, aes(y = set)) +
    geom_segment(aes(x = down, xend = up, yend = set), colour = "grey78", linewidth = 1.5) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40") +
    # each END marked against the null, which is what this figure adds
    geom_point(aes(x = up, colour = mark_up), size = 1.9) +
    geom_point(aes(x = down, colour = mark_dn), size = 1.9) +
    # the net, marked on its own (raw p, the pre-registered convention)
    geom_point(aes(x = net, shape = mark_net), fill = "white", colour = "black",
               size = 2.1, stroke = 0.6) +
    facet_wrap(~cond, nrow = if (length(unique(d$cond)) > 5) 3 else 1) +
    scale_colour_manual(values = c(`TRUE` = "#B2182B", `FALSE` = "grey55"),
                        labels = c(`TRUE` = "direction beats its null (FDR < 5%)",
                                   `FALSE` = "not significant"), name = NULL) +
    scale_shape_manual(values = c(`TRUE` = 23, `FALSE` = 21),
                       labels = c(`TRUE` = "net also significant", `FALSE` = "net not significant"),
                       name = NULL) +
    theme_bw(base_size = 12) +
    theme(axis.text.y = element_text(size = 6.5), legend.position = "top",
          legend.box = "vertical", panel.grid.minor = element_blank(),
          panel.grid.major.y = element_blank(), strip.background = element_blank()) +
    labs(y = NULL, x = paste0(
      "Summed gene contributions within each set, mortality-clock units: bar spans the ",
      "downward sum to the upward sum, hollow point is the net\n",
      "each direction tested against the weight-matched null, Benjamini-Hochberg within ",
      "this family of ", nfam, " tests\n",
      "sets ordered by total movement, least at the bottom"))
  ggsave(file.path(RERUN_DIR, file), p, width = width, height = height, dpi = 300, limitsize = FALSE)
  cat(sprintf("Saved -> %s  (%d of %d comparisons have a direction significant with the net not)\n",
              file, sum((d$mark_up | d$mark_dn) & !d$mark_net), nrow(d)))
}

a <- sides %>% filter(dataset == "arrest_conditions")
t <- sides %>% filter(dataset == "time_course")
plot_sides(a, "figure_cancellation_sides.png", 15, 11, 500)
plot_sides(t, "figure_cancellation_sides_temporal.png", 15, 20, 900)

# ---- concordance against its null -----------------------------------------
cc <- conc %>% filter(dataset == "arrest_conditions") %>%
  mutate(set = gsub("^HALLMARK ", "", pathway)) %>%
  group_by(set) %>%
  summarise(rho = median(rho), null = median(null_median),
            n_sig = sum(q_emp < 0.05), n_pairs = n(), .groups = "drop") %>%
  mutate(all_pairs = n_sig == n_pairs) %>%
  arrange(rho)
cc$set <- factor(cc$set, levels = cc$set)
NULL_LINE <- median(cc$null)
p2 <- ggplot(cc, aes(y = set, x = rho)) +
  geom_vline(xintercept = NULL_LINE, linetype = "dashed", colour = "#2166AC") +
  geom_segment(aes(x = NULL_LINE, xend = rho, yend = set), colour = "grey78", linewidth = 1.2) +
  geom_point(aes(colour = all_pairs, size = n_sig)) +
  scale_colour_manual(values = c(`TRUE` = "#B2182B", `FALSE` = "grey45"),
                      labels = c(`TRUE` = "beats its null in every condition pair",
                                 `FALSE` = "does not"), name = NULL) +
  scale_size_continuous(range = c(1.4, 4), name = "condition pairs\nwith FDR < 5%") +
  theme_bw(base_size = 12) +
  theme(axis.text.y = element_text(size = 7.5), legend.position = "right",
        panel.grid.minor = element_blank(), panel.grid.major.y = element_blank()) +
  labs(y = NULL, x = paste0(
    "Correlation between conditions of a set's per-gene contributions (median over the 10 pairs)\n",
    "dashed line: random gene groups of the same size (", sprintf("%.2f", NULL_LINE), ")"))
ggsave(file.path(RERUN_DIR, "figure_set_concordance.png"), p2,
       width = 11, height = 10, dpi = 300, limitsize = FALSE)
cat(sprintf("Saved -> figure_set_concordance.png  (%s beat the null in all pairs)\n",
            paste(cc$set[cc$all_pairs], collapse = ", ")))
