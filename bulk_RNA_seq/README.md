# Bulk RNA-seq analysis

This is the bulk RNA-seq arm of the paper (companion to the GSE226225 scRNA-seq analysis).
It has two independent pipelines, pooling recount3 data in different ways:

- **`meta_analysis/`** — cross-study recount3 meta-analysis. Pools 34 SRA/ENA studies to
  build one object spanning Proliferating, Contact-inhibited CQ, Serum-starved CQ,
  Replicative CS (RS), Stress-induced CS (SIPS), and Oncogene-induced CS (OIS) fibroblast
  samples, calls DEGs for each condition vs. Proliferating, and runs Hallmark/GO/KEGG
  enrichment. **This is the pipeline that needs recount3's AUC-based count scaling**
  (`transform_counts()`, applied in `process_rse()`) because it spans studies of different
  sequencing depths/read lengths.
- **`temporal_analysis/`** — the single-study ERP021140 (Hernandez-Segura et al.) time
  course: fibroblast/keratinocyte/melanocyte at none/4/10/20 days post-irradiation.
  A single study doesn't need cross-study count scaling, so this pipeline does not apply it
  (DESeq2's own size factors handle within-study depth differences) — see
  `../DISCREPANCY_REPORT/` for the verification of that scaling decision.

## Setup and reproduction

`REPRODUCE.md`: software versions (`env/renv.lock`, `env/requirements.txt`), the inputs
that are not in the repository, and how to run and check the pipeline.
`bash run_pipeline.sh` runs the 41 scripts in `pipeline.tsv` that produce every output the
tAge sections cite; `PROVENANCE.md` maps each result to its output and script.
`CLAUDE.md`: orientation for agents working in this directory.

Run scripts from inside `bulk_RNA_seq/`: they `source("R/config.R")` by relative path.
`BULK_RNASEQ_DIR` overrides the directory that holds `bulk_RNA_seq/` and `Final/`
(default: the parent of the working directory).

## Run order

The complete order for the tAge results is `pipeline.tsv`. The notes below cover the
DEG, enrichment and SI-figure scripts.

**meta_analysis/**
1. `01_build_cs_cq_object.R` — downloads the 34 studies from recount3, builds the pooled,
   protein-coding-filtered, AUC-scaled RSE (`rerun_outputs/cs_cq_all_study_processed.rds`).
2. `02_call_degs.R` — DESeq2 fit (`~study + cell_substate`) + DEG calling per condition vs
   Proliferating (`rerun_outputs/arrest_degs_final_RERUN.csv`).
3. `03_hallmark_enrichment.R` — MSigDB Hallmark pathway overlap with arrest DEGs.
4. `04_shared_deg_go_kegg.R` — GO/KEGG enrichment of DEGs shared across all 5 conditions.
5. `05_tage_all_conditions.R` — transcriptomic age (Gladyshev-Lab/tAge) per condition.

Scripts `06`-`28` were added in the annotation/confound pass. Three of them are
load-bearing for the supplementary figures and should be run in this order after `01`:

- `10_immortalisation_annotation_audit.R` — audits `cell_line`, `immortalised` and
  `tissue`; writes `rerun_outputs/immortalisation_annotation_corrected.csv`, which
  carries `tissue_verified` (a per-strain map with a source for every entry, and a hard
  stop on any unmapped strain).
- `15_patch_si_table_annotation.py` — applies those corrections **in place** to the SI
  tables, preserving the submitted labels as `cell_line_assubmitted` and
  `tissue_assubmitted`. Idempotent. Fixes `tissue` for 10 samples in two strains
  (HCA2-hTERT, HDF 12-3: "Skin" -> "Foreskin"), so the tissue panels of the PCA figures
  are right at the source rather than only inside the figure scripts.
- `27_si_figure1_pca.R` and `28_si_figure16_pca.R` — SI Figure 1 (all arrest samples) and
  SI Figure 16 (SIPS/OIS/CQ), each before and after removing the study batch effect.
  Neither existed on this branch: the submitted versions came from `study_degs.R` and
  `marian_variance.R`, which are on origin/main and on no branch respectively. Script 28
  also needed `filter_dataframe()`, which was defined inline in `study_degs.R` and is now
  in `R/functions.R`. Both reproduce the published variance figures (71% / 70.7% and
  83% / 83.4%).

  Note on numbering: SI Figure 16 is the number in the current supplementary file
  (`Final/final/SI_figures_final.pdf`, 30 Sep 2024) and in `modules_october_2024.docx`.
  The same figure is SI Figure 13 in `Final/nature_ageing/SI_figures.pdf` and is cited as
  "SI Figure 12" in the older `modules_final.docx`.

The remaining meta_analysis/ scripts (`DEG_analysis_recount3.R`, `cq_DEG_analysis_recount3.R`,
`degs_v_stress.R`, `stress_responses.R`, `enrich_cellage.R`, `kegg_subset.R`,
`macrophage_overlaps.R`, `post_preprint.R`, `metabolism_script.R`, `sensig_to_human.R`) are the
original downstream figure/overlap scripts, **portabilized (paths fixed) but not
independently re-executed and re-verified line-by-line** in this pass — see the "ported,
not reverified" note in `../DISCREPANCY_REPORT/DISCREPANCY_REPORT.md`. Several still read
personal-machine input files that aren't recoverable from anything in this repo (curated
gene lists like `DNA_repair.csv`, `SCAPs.csv`, WikiPathways text exports, etc.) — each is
commented inline with what's missing.

**temporal_analysis/**
1. `01_build_metadata.R` — reconstructs the ERP021140 sample metadata (3 cell types) from
   the per-cell-type `sample_pheno.csv` files already in `Final/ERP021140/`.
2. `02_run_time_analysis.R` — DEG calling + 10,000-simulation overlap test per cell type.
3. `03_temporal_overlaps.R` — cross-cell-type/cross-timepoint overlap analysis and
   CellAge/metabolism-pathway comparisons (ported, not independently re-verified).

## What's fully rerun-and-verified vs. ported-only

See `DISCREPANCY_REPORT/DISCREPANCY_REPORT.md` for the full verification: sample counts
match the paper exactly (6/6 conditions); the Ensembl gene dictionary is pinned to the
release-100 GTF (100% agreement with the original `cq_samples.rds`). MSigDB Hallmark sets
come from the msigdbr release in `env/renv.lock` (26.1.0, MSigDB 2026.1.Hs), not the
original version.

## Where the original (unrestructured) code lives

The complete original codebase — every script, including stale/superseded/exploratory
ones, and all the manuscript drafts/figures — is archived unchanged at
`../archive_original/` (gitignored; not part of this repo).
