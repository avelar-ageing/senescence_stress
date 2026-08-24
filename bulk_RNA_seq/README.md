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

## Setup

```r
# R packages: recount3, biomaRt, DESeq2, edgeR, WGCNA, tibble, pheatmap, gplots, tidyr,
# dplyr, GeneOverlap, EnhancedVolcano, RColorBrewer, clusterProfiler, stringi, rrvgo,
# org.Hs.eg.db, msigdbr, ComplexUpset, ggpubr, ggplot2
```

Set `BULK_RNASEQ_DIR` to this repo's `systems_analysis_arrest/`-equivalent root if you're
not running scripts from inside `bulk_RNA_seq/` itself (see `R/config.R`).

## Run order

**meta_analysis/**
1. `01_build_cs_cq_object.R` — downloads the 34 studies from recount3, builds the pooled,
   protein-coding-filtered, AUC-scaled RSE (`rerun_outputs/cs_cq_all_study_processed.rds`).
2. `02_call_degs.R` — DESeq2 fit (`~study + cell_substate`) + DEG calling per condition vs
   Proliferating (`rerun_outputs/arrest_degs_final_RERUN.csv`).
3. `03_hallmark_enrichment.R` — MSigDB Hallmark pathway overlap with arrest DEGs.
4. `04_shared_deg_go_kegg.R` — GO/KEGG enrichment of DEGs shared across all 5 conditions.
5. `05_tage_all_conditions.R` — transcriptomic age (Gladyshev-Lab/tAge) per condition.

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

See `../DISCREPANCY_REPORT/DISCREPANCY_REPORT.md` for the full verification: sample counts
match the paper exactly (6/6 conditions); the Ensembl gene dictionary and MSigDB Hallmark
gene sets could not be pinned to the exact original versions (external archives retired
since the original run) and that gap is measured and disclosed, not hidden.

## Where the original (unrestructured) code lives

The complete original codebase — every script, including stale/superseded/exploratory
ones, and all the manuscript drafts/figures — is archived unchanged at
`../archive_original/` (gitignored; not part of this repo).
