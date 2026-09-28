# 48_figure_set_driver_genes.R
#
# WHY. exploratory/47 finds, for every gene set that beats its matched null, how much of
# the set's excess each gene carries. This figure shows which genes carry which sets,
# and which genes carry several sets at once.
#
# LAYOUT, one figure per dataset:
#   columns  gene sets that beat their matched null (p < 0.05) in at least one group
#            (condition, or cell type x timepoint); the number of such groups in brackets
#   rows     each set's largest contributor in any group, plus every gene among the three
#            largest contributors of two or more sets; genes carrying most sets at the top
#   dot      the gene is among the set's three largest contributors in at least one group
#            size   = number of groups in which it is
#            colour = its median share of the set's excess over the matched null, in %
#            black outline = in at least one group the set no longer beats its null
#                            without this gene
#   grey dot the gene belongs to the set but is not among its three largest there
#
# A gene's excess in a group is the same in every set that contains it (exploratory/47),
# so a gene shared by several sets carries each of them by the same amount.
#
# USAGE: Rscript exploratory/48_figure_set_driver_genes.R      (from bulk_RNA_seq/)
# OUT:   rerun_outputs/figure_set_driver_genes.png            (arrest conditions)
#        rerun_outputs/figure_set_driver_genes_temporal.png   (irradiation time course)

source("R/config.R")
suppressPackageStartupMessages({ library(dplyr); library(ggplot2) })

G <- read.csv(file.path(RERUN_DIR, "set_driver_genes.csv"))
S <- read.csv(file.path(RERUN_DIR, "set_driver_summary.csv"))
G$set <- sub("^HALLMARK ", "", G$pathway)
G$set_depends_on_gene <- as.logical(G$set_depends_on_gene)   # pandas writes True/False
stopifnot(!anyNA(G$set_depends_on_gene))
S$set <- sub("^HALLMARK ", "", S$pathway)

plot_arm <- function(arm_, file, n_groups_label) {
  g <- G %>% filter(arm == arm_)
  s <- S %>% filter(arm == arm_)
  top3 <- g %>% filter(rank_in_set <= 3)
  keep <- union(g$gene[g$rank_in_set == 1],
                top3 %>% group_by(gene) %>% summarise(n = n_distinct(set), .groups = "drop") %>%
                  filter(n >= 2) %>% pull(gene))
  cells <- top3 %>% filter(gene %in% keep) %>% group_by(gene, set) %>%
    summarise(n_groups = n(), share_pct = 100 * median(share_of_set_excess),
              depends = any(set_depends_on_gene), .groups = "drop")
  member <- g %>% filter(gene %in% keep) %>% distinct(gene, set) %>%
    anti_join(cells, by = c("gene", "set"))
  set_n <- s %>% count(set, name = "n_beat") %>% arrange(desc(n_beat), set)
  gene_ord <- cells %>% group_by(gene) %>%
    summarise(n_sets = n_distinct(set), n = sum(n_groups), .groups = "drop") %>%
    arrange(n_sets, n, desc(gene))
  lab_set <- setNames(sprintf("%s (%d)", set_n$set, set_n$n_beat), set_n$set)
  lab_gene <- setNames(sprintf("%s (%d)", gene_ord$gene, gene_ord$n_sets), gene_ord$gene)
  lv <- function(x, y) factor(x, levels = y)
  cells$set <- lv(cells$set, set_n$set); member$set <- lv(member$set, set_n$set)
  cells$gene <- lv(cells$gene, gene_ord$gene); member$gene <- lv(member$gene, gene_ord$gene)

  p <- ggplot() +
    geom_point(data = member, aes(x = set, y = gene), colour = "grey75", size = 0.9) +
    geom_point(data = cells, aes(x = set, y = gene, size = n_groups, fill = share_pct,
                                 colour = depends), shape = 21, stroke = 0.8) +
    scale_fill_gradient(low = "#FDDBC7", high = "#B2182B", limits = c(0, NA),
                        oob = scales::squish, name = "share of the set's\nexcess (%)") +
    scale_colour_manual(values = c(`TRUE` = "black", `FALSE` = "grey60"),
                        labels = c(`TRUE` = "set no longer beats its null without this gene",
                                   `FALSE` = "set still beats its null without it"),
                        name = NULL) +
    scale_size_continuous(range = c(2, 7), breaks = seq_len(max(cells$n_groups)),
                          name = paste0("number of ", n_groups_label)) +
    # limits fix the order: with two layers, ggplot otherwise orders by first appearance
    scale_x_discrete(limits = set_n$set, labels = lab_set) +
    scale_y_discrete(limits = gene_ord$gene, labels = lab_gene) +
    theme_bw(base_size = 12) +
    theme(axis.text.x = element_text(angle = 50, hjust = 1, size = 9),
          axis.text.y = element_text(face = "italic", size = 9),
          panel.grid.minor = element_blank(), legend.position = "right") +
    guides(colour = guide_legend(order = 1, override.aes = list(size = 4, fill = "white")),
           size = guide_legend(order = 2), fill = guide_colourbar(order = 3)) +
    labs(x = paste0("gene sets that beat their matched null (number of ", n_groups_label, ")"),
         y = "gene (number of sets it helps carry)")
  h <- 2 + 0.24 * length(keep); w <- 4.5 + 0.42 * nrow(set_n)
  ggsave(file.path(RERUN_DIR, file), p, width = max(w, 10), height = max(h, 6), dpi = 300)
  cat(sprintf("Saved -> %s (%d genes x %d sets)\n", file, length(keep), nrow(set_n)))
}

plot_arm("arrest", "figure_set_driver_genes.png", "conditions")
plot_arm("time_course", "figure_set_driver_genes_temporal.png", "groups")
