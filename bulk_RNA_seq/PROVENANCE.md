# Provenance of the tAge results (sections 2.1.5 - 2.2.4)

Every result stated in the manuscript, the file that holds it, and the script
that writes that file. Regenerate with `exploratory/46_build_provenance.py`,
which fails if a citation, an output or a write call goes missing.

Paths are relative to `bulk_RNA_seq/`. Scripts take the output directory and,
where a clock is needed, the model directory as arguments; see `R/config.R`
and each script's header.

| Section | Result | Output | Script | Lang | Needs first |
|---|---|---|---|---|---|
| 2.1.5 | WI-38 parental vs hTERT tAge trajectory and contrasts | `rerun_outputs/gse175533/gse175533_tage.csv` | `meta_analysis/22_gse175533_htert_contrasts.py` | Py | meta_analysis/20_gse175533_export.py; 21_gse175533_tage.R |
|  |  | `rerun_outputs/gse175533/gse175533_contrasts.csv` | `meta_analysis/22_gse175533_htert_contrasts.py` | Py |  |
| 2.1.5 | WI-38 figure | `rerun_outputs/figure_immortalisation_htert.png` | `meta_analysis/25_gse175533_figures.R` | R | meta_analysis/22_gse175533_htert_contrasts.py |
| 2.1.5.1 | proliferating-control heterogeneity between studies and strains | `rerun_outputs/proliferating_baseline_heterogeneity.csv` | `meta_analysis/11_proliferating_baseline_heterogeneity.R` | R | meta_analysis/05_tage_all_conditions.R |
| 2.1.5.1 | baseline heterogeneity figure | `rerun_outputs/figure_proliferating_baseline_heterogeneity.png` | `meta_analysis/17_universal_tage_figures.R` | R | meta_analysis/11_proliferating_baseline_heterogeneity.R |
| 2.1.5.1 | within-study effect and permutation p for each condition, all three clocks | `rerun_outputs/section_2_1_5_summary_table.csv` | `meta_analysis/26_section_summary_table.py` | Py | meta_analysis/13_condition_within_study.R; 18_mortality_clock.py |
| 2.1.5.1 / 2.1.5.2 | control SD, study-median range and immortalisation gap per clock | `rerun_outputs/clock_comparison.csv` | `meta_analysis/29_clock_comparison.R` | R | meta_analysis/05_tage_all_conditions.R; 18_mortality_clock.py |
| 2.1.5.2 | curated immortalisation status and cell line per sample, with GEO evidence | `rerun_outputs/immortalisation_annotation_corrected.csv` | `meta_analysis/10_immortalisation_annotation_audit.R` | R | meta_analysis/05_tage_all_conditions.R |
| 2.1.5.2 | immortalisation and tissue-background contrasts (tests A-F) | `rerun_outputs/immortalisation_contrasts.csv` | `meta_analysis/16_immortalisation_contrasts.R` | R | rerun_outputs/immortalisation_annotation_corrected.csv |
| 2.1.5.3 | clock representation per set; five-gene dominance per clock | `rerun_outputs/pathway_representation.csv` | `exploratory/14_pathway_representation.py` | Py | exploratory/01_export_full_tage_matrices.R |
|  |  | `rerun_outputs/pathway_representation_mortality.csv` | `exploratory/14_pathway_representation.py` | Py |  |
| 2.1.5.3 / 2.2.4 | transcriptome-wide cancellation, significant-set counts, separation vs magnitude | `rerun_outputs/decomposition_diagnostics.csv` | `exploratory/23_decomposition_diagnostics.py` | Py | exploratory/01_export_full_tage_matrices.R; 02_export_temporal_bytimepoint_matrices.R |
| 2.1.5.3.2 (i) | subsampling every condition to 11v11, 10,000 draws | `rerun_outputs/decomposition_subsampling.csv` | `exploratory/23_decomposition_diagnostics.py` | Py | as above |
| 2.1.5.3.2 (ii) / 2.2.4 | gene-set coverage of the clock, weight-matched null, n_half, n_eff | `rerun_outputs/coverage_threshold_free.csv` | `exploratory/42_coverage_threshold_free.py` | Py | exploratory/01_...R; 02_...R |
|  |  | `rerun_outputs/clock_unaccounted_shares.csv` | `exploratory/42_coverage_threshold_free.py` | Py |  |
|  |  | `rerun_outputs/gene_contribution_table.csv` | `exploratory/42_coverage_threshold_free.py` | Py |  |
| 2.1.5.3.2 (ii) / 2.2.4 | coverage figure (panels a and b) | `rerun_outputs/figure_signal_concentration.png` | `exploratory/40_figure_signal_concentration.R` | R | exploratory/42_coverage_threshold_free.py |
| 2.1.5.4 | recurrent top genes, CS/CQ signatures, matched-random nulls | `rerun_outputs/signal_concentration.csv` | `exploratory/39_signal_concentration.py` | Py | exploratory/01_...R; 02_...R |
|  |  | `rerun_outputs/signal_concentration_top_genes.csv` | `exploratory/39_signal_concentration.py` | Py |  |
|  |  | `rerun_outputs/signal_concentration_classes.csv` | `exploratory/39_signal_concentration.py` | Py |  |
|  |  | `rerun_outputs/signal_concentration_directions.csv` | `exploratory/39_signal_concentration.py` | Py |  |
| 2.1.5.1 | arrest signal with cell-cycle features removed, plus size control | `rerun_outputs/arrest_minus_cellcycle.csv` | `exploratory/37_arrest_minus_cellcycle.py` | Py | exploratory/01_export_full_tage_matrices.R |
| 2.1.5.3.2 (iii) | split-half reproducibility ceiling, within-study vs pooled controls | `rerun_outputs/set_ceiling_comparison.csv` | `exploratory/34_ceiling_comparison.py` | Py | exploratory/01_export_full_tage_matrices.R; 03_build_pathway_mouse_id_mapping.R; meta_analysis/10_immortalisation_annotation_audit.R |
|  |  | `rerun_outputs/set_ceiling_paired_tests.csv` | `exploratory/34_ceiling_comparison.py` | Py |  |
| 2.1.5.3.1 | matched-null yardstick, 250 arrest comparisons | `rerun_outputs/pathway_specificity_yardstick_mortality.csv` | `exploratory/20_pathway_specificity_yardstick.py` | Py | run with --mortality; exploratory/01_export_full_tage_matrices.R; 14_pathway_representation.py --mortality |
| 2.1.5.3.1 | condition effect within a shared cell line (IMR90) | `rerun_outputs/condition_within_cellline.csv` | `meta_analysis/12_condition_within_cellline.R` | R | meta_analysis/05_tage_all_conditions.R |
| 2.2.3 | every tAge value and FDR in the time course, all three clocks | `rerun_outputs/tage_temporal_tests.csv` | `temporal_analysis/07_tage_temporal_tests.R` | R | temporal_analysis/04_tage_by_celltype.R; 08_mortality_temporal.py; R_keratinocyte_batch.R |
| 2.2.4 | per-set cancellation, arrest and time course | `rerun_outputs/withinset_cancellation.csv` | `exploratory/23_decomposition_diagnostics.py` | Py | exploratory/01_...R; 02_...R |
| 2.2.4 | cancellation summary per dataset (retention, same-direction fraction) | `rerun_outputs/withinset_cancellation_summary.csv` | `exploratory/28_figure_withinset_cancellation.R` | R | exploratory/23_decomposition_diagnostics.py |
| 2.2.4 | within-set cancellation figure, time course | `rerun_outputs/figure_withinset_cancellation_temporal.png` | `exploratory/28_figure_withinset_cancellation.R` | R | exploratory/23_decomposition_diagnostics.py |
| 2.2.4 | up/down side test figure, time course | `rerun_outputs/figure_cancellation_sides_temporal.png` | `exploratory/31_figure_cancellation_sides.R` | R | exploratory/29_set_cancellation_structure.py |
| 2.2.4 | matched-null yardstick, 450 time-course comparisons | `rerun_outputs/pathway_specificity_yardstick_mortality_temporal.csv` | `exploratory/22_mortality_pathway_decomposition.py` | Py | exploratory/02_export_temporal_bytimepoint_matrices.R |
| 2.2.4 | set recurrence across groups and its exact probability | `rerun_outputs/yardstick_recurrence.csv` | `exploratory/25_yardstick_calibration.py` | Py | exploratory/22_mortality_pathway_decomposition.py |
| 2.2.4 | within- vs between-cell-type profile similarity, PCA separation | `rerun_outputs/tage_modularity_mortality.csv` | `exploratory/21_tage_modularity.R` | R | exploratory/22_mortality_pathway_decomposition.py |
| 2.2.4 | gene-shuffle nulls: within/between gap, PC variance, sign changes, melanocyte dissent | `rerun_outputs/composition_vs_magnitude.csv` | `exploratory/27_composition_vs_magnitude.py` | Py | exploratory/02_...R |

27 results, 35 output files, 24 scripts (12 R, 12 Python).
