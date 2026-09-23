# Reproducing the bulk tAge results

Every result in the tAge sections (2.1.5–2.2.4) comes from a file in `rerun_outputs/`
written by a script listed in `pipeline.tsv`; `PROVENANCE.md` maps each result to its file
and script. This page lists what to install, which inputs to supply, and how to run and
check the pipeline.

## 1. Software

| | Version used | Pinned in |
|---|---|---|
| R | 4.6.0 | `env/renv.lock` (259 packages; Bioconductor 3.23) |
| tAge R package | 1.0.1, GitHub `Gladyshev-Lab/tAge` @ `3f1b4914c683` | `env/renv.lock` |
| msigdbr | 26.1.0 (MSigDB 2026.1.Hs) | `env/renv.lock` |
| Python | 3.12.3 | `env/requirements.txt` (scikit-learn 1.9.0, numpy 2.5.3, pandas 3.0.5, scipy 1.18.1, joblib 1.6.0) |

macOS (Apple silicon), once, from `bulk_RNA_seq/`:

```
xcode-select --install                   # compilers for packages built from source
# install R 4.6 arm64 from https://cloud.r-project.org/bin/macosx/
# install gfortran from https://mac.r-project.org/tools/
brew install python@3.12
bash env/setup_macos.sh                  # .venv, R packages from renv.lock, input check
```

On Linux, the same script works with R 4.6 and `python3.12` on the path.

`setup_macos.sh` installs the R packages into `bulk_RNA_seq/.Rlib`, a library used only by
this project; `run_pipeline.sh` puts it first and runs `env/check_versions.sh`, which stops
if any installed R or Python package differs from the lockfiles.

R reaches Python through reticulate. `R/config.R` points it at `.venv/bin/python`
(override with `TAGE_PYTHON`), so the R and Python scripts use the same scikit-learn.

## 2. Inputs

`bash env/check_inputs.sh` checks each against the md5 in `env/external_inputs.tsv`.

Not in this repository (supply them yourself):

| File | Where it comes from |
|---|---|
| `external_inputs/models/EN_Chronoage_Multispecies_Multitissue_scaleddiff.pkl` | tAge clock (Tyshkovskiy et al.), from the authors; not redistributed here |
| `external_inputs/models/EN_Chronoage_Multispecies_Multitissue_yugenediff.pkl` | as above |
| `external_inputs/models/EN_Mortality_Multispecies_Multitissue_scaleddiff.pkl` | as above |
| `external_inputs/geo/GSE175533_TPM.xlsx` | GEO GSE175533 supplementary file `GSE175533_hTERT.RS.RIS.CD.TPM_table.xlsx` (46,113,681 bytes), renamed; strict-OOXML workbook read by `meta_analysis/xlsx_strict.py` |

In this repository, at the repository root (tracked although `Final/` is otherwise ignored):

| File | Contents |
|---|---|
| `Final/SI_tables/study_info_all.pre_annotation_fix.csv` | sample sheet of the 34 studies, as submitted (SI table). Path set by `STUDY_INFO_CSV` in `R/config.R` |
| `Final/SI_tables/lyso_genes.csv` | lysosomal gene set (SI table) |
| `Final/SI_tables/enrichment_background.csv` | enrichment background (SI table) |
| `Final/ERP021140/{Fibroblast,Keratinocyte,Melanocyte}/sample_pheno.csv` | time-course sample sheets |

`R/config.R` reads `Final/` from the directory above `bulk_RNA_seq/`: the repository root in
a clone, `systems_analysis_arrest/` in the original working directory. `meta_analysis/01`
stops if it is given the corrected sample sheet (`Final/SI_tables/study_info_all.csv`,
rewritten in place by `meta_analysis/15` on 2026-08-31): built from it, `meta_analysis/10`
records the corrected labels as the submitted ones, and the audit columns of
`immortalisation_annotation_corrected.csv` change. `Final/SI_tables/cq_samples.rds` is
optional; when present, `meta_analysis/01` checks the gene dictionary against it.

Downloaded by the scripts (network needed on the first run only):

| Script | Source | Cached as |
|---|---|---|
| `meta_analysis/01` | Ensembl release-100 GTF (47 MB) | `rerun_outputs/Homo_sapiens.GRCh38.100.gtf.gz`, `human_pc_ensembl100_exact.rds` |
| `meta_analysis/01` | recount3, 34 studies | `rerun_outputs/cs_cq_download_raw.rds` |
| `exploratory/01, 02`, `temporal_analysis/04` | recount3, ERP021140 | `rerun_outputs/erp021140_download_raw.rds` |
| `meta_analysis/14` | NCBI GEO sample and series records | `.geo_cache/` |

GEO answers some networks with a reCAPTCHA page. `meta_analysis/14` then stops without
caching anything. Its output, `rerun_outputs/sample_level_line_verification.csv`, can be
copied from an earlier run; `run_pipeline.sh` skips a script whose outputs are present.

## 3. Run

```
bash run_pipeline.sh --dry-run    # what will run
bash run_pipeline.sh              # everything whose outputs are missing
```

Logs go to `rerun_outputs/logs/`. The driver stops at the first script that fails or that
does not write the outputs listed for it in `pipeline.tsv`. After changing a script, re-run
it and everything downstream with `--from <script>`.

Stages: A builds the pooled objects (recount3 download, gene dictionary, GSE175533 export);
B computes tAge and the per-gene matrices; C runs the analyses, nulls and figures.

## 4. Check

```
.venv/bin/python env/compare_outputs.py <reference rerun_outputs> rerun_outputs
```

compares every CSV listed in `pipeline.tsv` (numeric columns within
|a − b| ≤ 1e-9 + 1e-6·|b|, other columns exactly) and reports figures and RDS files by md5
only.

## 5. Clean-checkout test, 2026-09-23

Fresh clone of `restructure-bulk-rna-seq`, fresh `.venv` from `env/requirements.txt`, only
the inputs in section 2 (the six `Final/` tables then outside git), R packages at the `env/renv.lock` versions (Linux x86_64, R 4.6.0,
8 cores). Seeded with two files, to avoid the network step and the GEO block:
`cs_cq_download_raw.rds` (recount3, 34 studies) and
`sample_level_line_verification.csv` (`meta_analysis/14`).

| | |
|---|---|
| Result | all 63 outputs listed in `pipeline.tsv` identical to the working tree's `rerun_outputs/` (`env/compare_outputs.py`, CSVs within 1e-6 relative; figures byte-identical) |
| Gene dictionary | 13,681 / 13,681 genes of the original `cq_samples.rds` |
| Run time | stage A 3 min, stage B 79 min (of which `meta_analysis/02`, six DESeq2 fits, 65–70 min), stage C 19 min |

Found and fixed during the test:

- `meta_analysis/01` read the corrected sample sheet; it now reads the submitted one
  (section 2).
- `pipeline.tsv` credited three outputs to the wrong invocation (`exploratory/14`,
  `exploratory/20 --mortality`, `meta_analysis/21` vs `22`).
- A GTF download cut off by R's 60 s timeout left a partial file that the next run would
  have read as complete.
- Another project replaced `msigdbr` 26.1.0 with 25.1.1 in the shared R library during
  the run (2 more Hallmark rows; 19 set-level CSVs changed). Hence `.Rlib` and
  `env/check_versions.sh`.
- Six working-tree outputs predated later script or annotation changes
  (`mortality_tage.csv`, `mortality_within_study.csv`, `section_2_1_5_summary_table.csv`,
  `yardstick_recurrence.csv` row order, `gse175533_contrasts.csv`,
  `pathway_specificity_yardstick_mortality_temporal.csv` `p_floor`). They were
  regenerated in place; the previous versions are kept as `*.PRE_REFRESH_20260923`.

Not covered: the recount3 and GEO downloads themselves, and macOS. numpy wheels for macOS
arm64 link Apple Accelerate instead of OpenBLAS, which can change trailing digits;
`compare_outputs.py` reports any difference above 1e-6.
