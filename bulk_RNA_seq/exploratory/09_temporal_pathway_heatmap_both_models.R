# 09_temporal_pathway_heatmap_both_models.R
#
# Temporal counterpart to exploratory/08_pathway_full_heatmap_both_models.R:
# comprehensive partial-tAge heatmap for the ERP021140 temporal time course,
# pooled irradiated vs none per cell type (Fibroblast/Keratinocyte/
# Melanocyte), all 50 pathways, both EN models shown side-by-side per cell
# type. Same geom_tile design as the meta-analysis figure (see that script's
# header for why pheatmap/clustering were dropped) -- alphabetical rows, no
# clustering, facets on top.

source("R/config.R")
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

df <- read.csv(file.path(RERUN_DIR, "partial_tage_ALL.csv"), check.names = FALSE)
pooled <- df[df$analysis == "temporal_pooled", ]
pooled$pathway_short <- gsub("^HALLMARK ", "", pooled$pathway)

sig_symbol <- function(p) ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", "")))

plot_df <- pooled %>%
  mutate(
    pathway_short = factor(pathway_short, levels = sort(unique(pathway_short), decreasing = TRUE)),
    label = factor(label, levels = c("Fibroblast", "Keratinocyte", "Melanocyte")),
    model_label = ifelse(model == "scaled", "Scaled", "YuGene"),
    model_label = factor(model_label, levels = c("Scaled", "YuGene")),
    stars = sig_symbol(p_adj)
  )

cap <- min(10, ceiling(quantile(abs(plot_df$cohens_d), 0.97, na.rm = TRUE)))

p <- ggplot(plot_df, aes(x = model_label, y = pathway_short, fill = cohens_d)) +
  geom_tile(colour = "grey80") +
  geom_text(aes(label = stars), size = 4) +
  facet_grid(~label) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0,
                       limits = c(-cap, cap), oob = scales::squish, name = "Cohen's d") +
  theme_minimal(base_size = 16) +
  theme(
    axis.text.y = element_text(size = 11),
    axis.text.x = element_text(size = 14, angle = 30, hjust = 1),
    axis.ticks.y = element_line(colour = "grey40"),
    axis.title = element_blank(),
    strip.text = element_text(size = 16, face = "plain"),
    legend.text = element_text(size = 13),
    legend.title = element_text(size = 15),
    panel.grid = element_blank(),
    panel.spacing = unit(0.4, "lines")
  )

ggsave(file.path(RERUN_DIR, "figure_temporal_pathway_heatmap.png"), p, width = 10, height = 16, dpi = 300)
cat(sprintf("Saved -> %s\n", file.path(RERUN_DIR, "figure_temporal_pathway_heatmap.png")))

n_sig <- pooled %>% group_by(pathway_short) %>% summarise(n_sig = sum(p_adj < 0.05))
cat("\nDistribution of n significant tests (of 6: 3 cell types x 2 models), all 50 pathways:\n")
print(table(n_sig$n_sig))
