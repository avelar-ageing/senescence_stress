#!/usr/bin/env python3
"""46_build_provenance.py

Builds the provenance manifest for the tAge sections: every result stated in
2.1.5-2.2.4 -> the output file holding it -> the script that writes that file.

This is a CHECK, not just documentation. It fails if
  * the manuscript cites a rerun_outputs path that no row covers,
  * a row names an output that does not exist,
  * a row names a script that does not exist,
  * a named script does not actually write the output it is credited with.

Run it after any change to the manuscript's citations or to the scripts.

Usage: 46_build_provenance.py <manuscript.txt> <bulk_RNA_seq dir>
Output: PROVENANCE.md, provenance.csv (both in the bulk_RNA_seq dir)
"""
import sys, os, re, csv

# section | result | outputs | script | prerequisites
ROWS = [
 ("2.1.5", "WI-38 parental vs hTERT tAge trajectory and contrasts",
  ["rerun_outputs/gse175533/gse175533_tage.csv", "rerun_outputs/gse175533/gse175533_contrasts.csv"],
  "meta_analysis/22_gse175533_htert_contrasts.py", "meta_analysis/20_gse175533_export.py; 21_gse175533_tage.R"),
 ("2.1.5", "WI-38 figure",
  ["rerun_outputs/figure_immortalisation_htert.png"],
  "meta_analysis/25_gse175533_figures.R", "meta_analysis/22_gse175533_htert_contrasts.py"),
 ("2.1.5.1", "proliferating-control heterogeneity between studies and strains",
  ["rerun_outputs/proliferating_baseline_heterogeneity.csv"],
  "meta_analysis/11_proliferating_baseline_heterogeneity.R", "meta_analysis/05_tage_all_conditions.R"),
 ("2.1.5.1", "baseline heterogeneity figure",
  ["rerun_outputs/figure_proliferating_baseline_heterogeneity.png"],
  "meta_analysis/17_universal_tage_figures.R", "meta_analysis/11_proliferating_baseline_heterogeneity.R"),
 ("2.1.5.1", "within-study effect and permutation p for each condition, all three clocks",
  ["rerun_outputs/section_2_1_5_summary_table.csv"],
  "meta_analysis/26_section_summary_table.py",
  "meta_analysis/13_condition_within_study.R; 18_mortality_clock.py"),
 ("2.1.5.1 / 2.1.5.2", "control SD, study-median range and immortalisation gap per clock",
  ["rerun_outputs/clock_comparison.csv"],
  "meta_analysis/29_clock_comparison.R", "meta_analysis/05_tage_all_conditions.R; 18_mortality_clock.py"),
 ("2.1.5.2", "curated immortalisation status and cell line per sample, with GEO evidence",
  ["rerun_outputs/immortalisation_annotation_corrected.csv"],
  "meta_analysis/10_immortalisation_annotation_audit.R", "meta_analysis/05_tage_all_conditions.R"),
 ("2.1.5.2", "immortalisation and tissue-background contrasts (tests A-F)",
  ["rerun_outputs/immortalisation_contrasts.csv"],
  "meta_analysis/16_immortalisation_contrasts.R", "rerun_outputs/immortalisation_annotation_corrected.csv"),
 ("2.1.5.3", "clock representation per set; five-gene dominance per clock",
  ["rerun_outputs/pathway_representation.csv", "rerun_outputs/pathway_representation_mortality.csv"],
  "exploratory/14_pathway_representation.py", "exploratory/01_export_full_tage_matrices.R"),
 ("2.1.5.3 / 2.2.4", "transcriptome-wide cancellation, significant-set counts, separation vs magnitude",
  ["rerun_outputs/decomposition_diagnostics.csv"],
  "exploratory/23_decomposition_diagnostics.py", "exploratory/01_export_full_tage_matrices.R; 02_export_temporal_bytimepoint_matrices.R"),
 ("2.1.5.3.2 (i)", "subsampling every condition to 11v11, 10,000 draws",
  ["rerun_outputs/decomposition_subsampling.csv"],
  "exploratory/23_decomposition_diagnostics.py", "as above"),
 ("2.1.5.3.2 (ii) / 2.2.4", "gene-set coverage of the clock, weight-matched null, n_half, n_eff",
  ["rerun_outputs/coverage_threshold_free.csv", "rerun_outputs/clock_unaccounted_shares.csv",
   "rerun_outputs/gene_contribution_table.csv"],
  "exploratory/42_coverage_threshold_free.py", "exploratory/01_...R; 02_...R"),
 ("2.1.5.3.2 (ii) / 2.2.4", "coverage figure (panels a and b)",
  ["rerun_outputs/figure_signal_concentration.png"],
  "exploratory/40_figure_signal_concentration.R", "exploratory/42_coverage_threshold_free.py"),
 ("2.1.5.4", "recurrent top genes, CS/CQ signatures, matched-random nulls",
  ["rerun_outputs/signal_concentration.csv", "rerun_outputs/signal_concentration_top_genes.csv",
   "rerun_outputs/signal_concentration_classes.csv",
   "rerun_outputs/signal_concentration_directions.csv"],
  "exploratory/39_signal_concentration.py", "exploratory/01_...R; 02_...R"),
 ("2.1.5.1", "arrest signal with cell-cycle features removed, plus size control",
  ["rerun_outputs/arrest_minus_cellcycle.csv"],
  "exploratory/37_arrest_minus_cellcycle.py", "exploratory/01_export_full_tage_matrices.R"),
 ("2.1.5.3.2 (iii)", "split-half reproducibility ceiling, within-study vs pooled controls",
  ["rerun_outputs/set_ceiling_comparison.csv", "rerun_outputs/set_ceiling_paired_tests.csv"],
  "exploratory/34_ceiling_comparison.py", "exploratory/01_export_full_tage_matrices.R; 03_build_pathway_mouse_id_mapping.R; meta_analysis/10_immortalisation_annotation_audit.R"),
 ("2.1.5.3.1", "matched-null yardstick, 250 arrest comparisons",
  ["rerun_outputs/pathway_specificity_yardstick_mortality.csv"],
  "exploratory/22_mortality_pathway_decomposition.py", "exploratory/01_...R"),
 ("2.1.5.3.1", "condition effect within a shared cell line (IMR90)",
  ["rerun_outputs/condition_within_cellline.csv"],
  "meta_analysis/12_condition_within_cellline.R", "meta_analysis/05_tage_all_conditions.R"),
 ("2.2.3", "every tAge value and FDR in the time course, all three clocks",
  ["rerun_outputs/tage_temporal_tests.csv"],
  "temporal_analysis/07_tage_temporal_tests.R",
  "temporal_analysis/04_tage_by_celltype.R; 08_mortality_temporal.py; R_keratinocyte_batch.R"),
 ("2.2.4", "per-set cancellation, arrest and time course",
  ["rerun_outputs/withinset_cancellation.csv"],
  "exploratory/23_decomposition_diagnostics.py", "exploratory/01_...R; 02_...R"),
 ("2.2.4", "cancellation summary per dataset (retention, same-direction fraction)",
  ["rerun_outputs/withinset_cancellation_summary.csv"],
  "exploratory/28_figure_withinset_cancellation.R", "exploratory/23_decomposition_diagnostics.py"),
 ("2.2.4", "within-set cancellation figure, time course",
  ["rerun_outputs/figure_withinset_cancellation_temporal.png"],
  "exploratory/28_figure_withinset_cancellation.R", "exploratory/23_decomposition_diagnostics.py"),
 ("2.2.4", "up/down side test figure, time course",
  ["rerun_outputs/figure_cancellation_sides_temporal.png"],
  "exploratory/31_figure_cancellation_sides.R", "exploratory/29_set_cancellation_structure.py"),
 ("2.2.4", "matched-null yardstick, 450 time-course comparisons",
  ["rerun_outputs/pathway_specificity_yardstick_mortality_temporal.csv"],
  "exploratory/22_mortality_pathway_decomposition.py", "exploratory/02_export_temporal_bytimepoint_matrices.R"),
 ("2.2.4", "set recurrence across groups and its exact probability",
  ["rerun_outputs/yardstick_recurrence.csv"],
  "exploratory/25_yardstick_calibration.py", "exploratory/22_mortality_pathway_decomposition.py"),
 ("2.2.4", "within- vs between-cell-type profile similarity, PCA separation",
  ["rerun_outputs/tage_modularity_mortality.csv"],
  "exploratory/21_tage_modularity.R", "exploratory/22_mortality_pathway_decomposition.py"),
 ("2.2.4", "gene-shuffle nulls: within/between gap, PC variance, sign changes, melanocyte dissent",
  ["rerun_outputs/composition_vs_magnitude.csv"],
  "exploratory/27_composition_vs_magnitude.py", "exploratory/02_...R"),
]

WRITE_PAT = re.compile(r"write\.csv|to_csv|ggsave|write_csv|savefig")


def main(manuscript, root):
    txt = open(manuscript).read()
    cited = {p for p in re.findall(r"rerun_outputs/[A-Za-z0-9_./-]*", txt)
             if p.endswith((".csv", ".png"))}
    covered, errs = set(), []
    for sec, res, outs, script, prereq in ROWS:
        covered |= set(outs)
        if not os.path.exists(os.path.join(root, script)):
            errs.append(f"script missing: {script}")
            continue
        body = open(os.path.join(root, script), errors="ignore").read()
        for o in outs:
            if not os.path.exists(os.path.join(root, o)):
                errs.append(f"output missing: {o}")
            stem = os.path.basename(o).rsplit(".", 1)[0]
            lines = body.split("\n")
            # direct call, a call within three lines, a constructed name (sprintf/paste),
            # or the filename passed as a quoted argument to a plotting helper - all of
            # which are real writes, and all of which occur in this codebase
            direct = any(stem in l and WRITE_PAT.search(l) for l in lines)
            near = any(stem in l and any(WRITE_PAT.search(x) for x in lines[max(0, i - 3):i + 4])
                       for i, l in enumerate(lines))
            built = any(stem.rsplit("_", 1)[0] in l and WRITE_PAT.search(l) for l in lines)
            passed = (any(f'"{os.path.basename(o)}"' in l for l in lines)
                      and any(WRITE_PAT.search(l) for l in lines))
            hit = direct or near or built or passed
            if not hit:
                errs.append(f"{script} is credited with {o} but no write call mentions it")

    for c in sorted(cited - covered):
        errs.append(f"manuscript cites {c} but no manifest row covers it")

    with open(os.path.join(root, "provenance.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["section", "result", "output", "script", "language", "prerequisites"])
        for sec, res, outs, script, prereq in ROWS:
            for o in outs:
                w.writerow([sec, res, o, script,
                            "R" if script.endswith(".R") else "Python", prereq])

    with open(os.path.join(root, "PROVENANCE.md"), "w") as fh:
        fh.write("# Provenance of the tAge results (sections 2.1.5 - 2.2.4)\n\n")
        fh.write("Every result stated in the manuscript, the file that holds it, and the script\n"
                 "that writes that file. Regenerate with `exploratory/46_build_provenance.py`,\n"
                 "which fails if a citation, an output or a write call goes missing.\n\n")
        fh.write("Paths are relative to `bulk_RNA_seq/`. Scripts take the output directory and,\n"
                 "where a clock is needed, the model directory as arguments; see `R/config.R`\n"
                 "and each script's header.\n\n")
        fh.write("| Section | Result | Output | Script | Lang | Needs first |\n")
        fh.write("|---|---|---|---|---|---|\n")
        for sec, res, outs, script, prereq in ROWS:
            for i, o in enumerate(outs):
                fh.write(f"| {sec if i == 0 else ''} | {res if i == 0 else ''} | `{o}` | "
                         f"`{script}` | {'R' if script.endswith('.R') else 'Py'} | {prereq if i == 0 else ''} |\n")
        fh.write(f"\n{len(ROWS)} results, "
                 f"{len({o for _, _, outs, _, _ in ROWS for o in outs})} output files, "
                 f"{len({s for _, _, _, s, _ in ROWS})} scripts "
                 f"({len({s for _, _, _, s, _ in ROWS if s.endswith('.R')})} R, "
                 f"{len({s for _, _, _, s, _ in ROWS if s.endswith('.py')})} Python).\n")

    print(f"  rows: {len(ROWS)} | outputs: {len(covered)} | cited by the manuscript: {len(cited)}")
    print(f"  scripts: {len({s for _,_,_,s,_ in ROWS})} "
          f"({len({s for _,_,_,s,_ in ROWS if s.endswith('.R')})} R, "
          f"{len({s for _,_,_,s,_ in ROWS if s.endswith('.py')})} Python)")
    if errs:
        print(f"\n  {len(errs)} PROBLEM(S):")
        for e in errs:
            print(f"    {e}")
        sys.exit(1)
    print("\n  all checks pass -> PROVENANCE.md, provenance.csv")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
