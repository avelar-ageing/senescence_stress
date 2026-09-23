# 35_figure_gene_level_reproducibility.R
#
# Two panels for the gene-level results.
#
# (a) HOW REPRODUCIBLE IS A SET'S PER-GENE CONTRIBUTION? Each of the 50 sets under
#     three settings, all a 3-versus-3 sample split with the controls split too
#     (34_ceiling_comparison.py). Referencing the six controls to all 91
#     proliferating samples, which is what the existing decomposition does, costs
#     0.36 of reproducibility against taking them from the same study. Matched
#     that way, the arrest conditions and the time course are indistinguishable.
#
# (b) WHERE DOES STRESSOR IDENTITY SIT? For each set, the cross-study ceiling
#     (32_set_reproducibility_ceiling.py) against its between-condition
#     correlation. Points on the diagonal are sets carried by the same genes in
#     every condition; points far below it are sets that recur while the genes
#     carrying them change. The cell-cycle sets sit on the diagonal and nearly
#     everything else falls below, which is the result.
#
# CLOCK: mortality only, as everywhere in the set-level analysis.
source("R/config.R")
suppressMessages({library(dplyr); library(ggplot2); library(patchwork); library(ggrepel)})

cmp  <- read.csv(file.path(RERUN_DIR, "set_ceiling_comparison.csv"))
ceil <- read.csv(file.path(RERUN_DIR, "set_reproducibility_ceiling.csv"))

LEV <- c("arrest_pooled_controls", "arrest_within_study", "time_course_within_study")
LAB <- c(arrest_pooled_controls = "Arrest conditions\ncontrols pooled across studies",
         arrest_within_study    = "Arrest conditions\ncontrols from the same study",
         time_course_within_study = "Irradiation time course\n(one study throughout)")
cmp <- cmp %>% mutate(setting = factor(setting, levels = LEV, labels = LAB[LEV]))

med <- cmp %>% group_by(setting) %>% summarise(m = median(ceiling_median), .groups = "drop")
pa <- ggplot(cmp, aes(x = setting, y = ceiling_median, fill = setting)) +
  geom_violin(alpha = 0.55, trim = FALSE, colour = "grey30") +
  geom_jitter(width = 0.09, size = 1.1, alpha = 0.55) +
  geom_text(data = med, aes(y = m, label = sprintf("%.2f", m)), x = c(1.36, 2.36, 3.36),
            inherit.aes = FALSE, size = 4.2, fontface = "bold") +
  scale_fill_manual(values = c("#2166AC", "#B2182B", "#009E73"), guide = "none") +
  scale_y_continuous(limits = c(0, 1)) +
  theme_bw(base_size = 13) +
  theme(axis.text.x = element_text(size = 10), panel.grid.minor = element_blank()) +
  labs(x = NULL, y = "reproducibility of a set's per-gene contributions\n(3 v 3 sample split, controls split too)",
       subtitle = "a  each point is one of the 50 gene sets")

lab_these <- c("HALLMARK E2F TARGETS", "HALLMARK G2M CHECKPOINT", "HALLMARK MYC TARGETS",
               "HALLMARK MITOTIC SPINDLE", "HALLMARK DNA REPAIR", "HALLMARK P53 PATHWAY",
               "HALLMARK EPITHELIAL MESENCHYMAL TRANSITION", "HALLMARK NOTCH SIGNALING",
               "HALLMARK TGF BETA SIGNALING", "HALLMARK HEDGEHOG SIGNALING",
               "HALLMARK UNFOLDED PROTEIN RESPONSE", "HALLMARK HYPOXIA")
ceil <- ceil %>% mutate(set = gsub("^HALLMARK ", "", pathway),
                        show = pathway %in% lab_these,
                        cellcycle = pathway %in% lab_these[1:5])
pb <- ggplot(ceil, aes(x = ceiling_median, y = between_median)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey45") +
  geom_point(aes(colour = cellcycle), size = 2.6, alpha = 0.9) +
  ggrepel::geom_text_repel(data = subset(ceil, show), aes(label = set), size = 3.2,
                           max.overlaps = 30, box.padding = 0.4, segment.colour = "grey60") +
  scale_colour_manual(values = c(`TRUE` = "#B2182B", `FALSE` = "grey40"),
                      labels = c(`TRUE` = "cell-cycle / DNA-repair sets", `FALSE` = "other sets"),
                      name = NULL) +
  coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
  theme_bw(base_size = 13) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank()) +
  labs(x = "reproducibility within a condition (cross-study ceiling)",
       y = "correlation between conditions",
       subtitle = "b  dashed line: the genes carrying a set are identical in every condition")

p <- pa + pb + plot_layout(widths = c(1, 1.15))
ggsave(file.path(RERUN_DIR, "figure_gene_level_reproducibility.png"), p,
       width = 15, height = 7.5, dpi = 300)
cat("Saved -> figure_gene_level_reproducibility.png\n")
cat(sprintf("  medians: %s\n", paste(sprintf("%s = %.2f", med$setting, med$m), collapse = "; ")))
cat(sprintf("  sets on/near the diagonal (gap < 0.15): %s\n",
            paste(ceil$set[ceil$ceiling_median - ceil$between_median < 0.15], collapse = ", ")))
