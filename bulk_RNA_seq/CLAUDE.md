# bulk_RNA_seq — orientation for agents

Bulk RNA-seq side of the senescence/quiescence paper: transcriptomic age (tAge) of
arrested fibroblasts. The single-cell analysis (GSE226225) is a separate project directory
and is not in this repository.

## Where things are

| What | Where |
|---|---|
| Working tree (edit here) | the local project directory `systems_analysis_arrest/bulk_RNA_seq/` — **not a git repository**. Machine-specific paths: `../CLAUDE.md` (local only, not in git) |
| Git clone | a separate clone of this repository, branch `restructure-bulk-rna-seq` |
| Copy working tree → clone | `../sync_to_repo.sh` (dry run), `../sync_to_repo.sh --apply` (local tool). Additive; runs the checks below |
| Manuscript (Results 2.1.5–2.2.4) | outside the repository; never commit it |
| Outputs | `rerun_outputs/` — gitignored, about 1 GB |
| Inputs | not in git: `external_inputs/` (clocks, GSE175533 table). In git: six tables under `Final/` at the repository root (`git add -f`; `Final/` is otherwise ignored). md5s in `env/external_inputs.tsv` |

Commit only in the clone. The GitHub repository is **public**: push only after the owner
says yes to that specific push. Never commit `rerun_outputs/`, `*.pkl` (third-party
clocks), manuscript drafts, edit logs or review notes. `sync_to_repo.sh` holds back
`DISCREPANCY_REPORT/METHODS_TAGE.md` and `DISCREPANCY_REPORT/RESULTS_2_1_5_REVISION.md`
until the owner decides. Do not delete or move files in the working tree; back up before
editing (`cp f f.PRE_<TAG>`).

## Running

```
bash env/setup_macos.sh          # once: .venv + R packages from env/renv.lock + input check
bash run_pipeline.sh --dry-run   # what would run
bash run_pipeline.sh             # runs every script in pipeline.tsv whose outputs are missing
bash run_pipeline.sh --from exploratory/39_signal_concentration.py   # after editing 39
```

- `pipeline.tsv`: 42 invocations of the 41 scripts behind every cited result, in execution
  order, with the outputs each writes. Scripts not listed are legacy or uncited (see README).
- `PROVENANCE.md` / `provenance.csv`: manuscript result → output file → script. Rebuild with
  `.venv/bin/python exploratory/46_build_provenance.py <manuscript.txt> .`; it fails if a
  cited file, output or write call is missing. Every result in the manuscript must cite a
  CSV, and every CSV and figure must come from a script in `pipeline.tsv`.
- Run every script from `bulk_RNA_seq/`: R scripts `source("R/config.R")` relative to it.
  Python scripts take `<rerun_dir> <model_dir>` as arguments.
- `env/compare_outputs.py <ref_dir> <new_dir>`: CSV-by-CSV comparison with tolerance.

## Data

- Arrest conditions: 34 recount3 studies, 230 samples, 91 proliferating controls, five
  conditions — contact-inhibited quiescence (CICQ), serum-starved quiescence (SSCQ),
  replicative senescence (RS), stress-induced (SIPS), oncogene-induced (OIS). Built by
  `meta_analysis/01`.
- Time course: ERP021140 / E-MTAB-5403, fibroblast, keratinocyte, melanocyte, 0/4/10/20
  days after irradiation.
- GSE175533: parental vs hTERT WI-38 over population doublings (`meta_analysis/20–25`).
- Clocks: three elastic nets from the tAge package authors — chronological scaled
  difference (1,839 non-zero genes), chronological YuGene (1,938), mortality (10,487, all
  non-zero). YuGene is a separately fitted model, not a renormalised copy. Report all
  three.
- Gene sets: 49 Hallmark sets + "Lysosomal Genes" = 50. MSigDB has 50 Hallmark sets;
  `exploratory/03` merges MYC_TARGETS_V1 and V2 into "HALLMARK MYC TARGETS".

## Statistical conventions

- Empirical p = (1 + k)/(B + 1), floor 1/(B + 1). No BH on simulation/permutation nulls:
  report the raw p and the floor.
- A contrast whose smallest attainable p is above 0.05 is **not testable**: report group
  sizes and means, no p. Rank-sum floor 2/C(n1+n2, min): 3v3 → 0.100, 7v2 → 0.056,
  3v6 → 0.024, 6v6 → 0.0022. Paired signed-rank floor 2/2^n.
- B = 20,000 draws for the set-level nulls: the largest family is 450 tests, and
  0.05/450 = 1.1e-4; at 10,000 draws the floor (1e-4) would sit 1.1× below it.
- Condition effects are estimated within study (labels permuted within study,
  `meta_analysis/13`, weight n_t·n_c/(n_t+n_c)). Pooling across studies mixes in the
  spread between studies' proliferating controls (study medians span 87.8 scaled-difference
  units, larger than most condition effects).
- Immortalisation and cell line: read `rerun_outputs/immortalisation_annotation_corrected.csv`.
  The `immortalised` and `cell_line` columns of `sample_metadata_RERUN.csv` are wrong
  (20 immortalised there vs 46 from GEO). IMR90 ER:RAS is not immortalised.
- Manuscript edits: new text wrapped in `** **`; no nested markers.

## Pitfalls already found

- `Final/SI_tables/study_info_all.csv` was corrected in place by `meta_analysis/15`
  (2026-08-31). `meta_analysis/01` builds from the submitted version,
  `study_info_all.pre_annotation_fix.csv` (`STUDY_INFO_CSV` in `R/config.R`), and stops if
  handed the corrected one. From the corrected table, counts, tAge and DEGs are identical
  but the "curated" audit columns of `immortalisation_annotation_corrected.csv` change
  (clean-checkout test, 2026-09-23).
- The R user library is shared with other projects. On 2026-09-23 one of them replaced
  msigdbr 26.1.0 with 25.1.1 mid-run and 19 set-level CSVs changed. The pipeline now uses
  `.Rlib` first, and `run_pipeline.sh` stops if `env/check_versions.sh` finds any package
  off the lockfile.
- A script that runs in several modes has one `pipeline.tsv` row per mode:
  `exploratory/14` (chronological, `--mortality`), `exploratory/20 --mortality` (the
  250-comparison arrest yardstick), `exploratory/21 mortality`.
- `download_studies()` (R/functions.R) skips a study that fails to download with only a
  message. Pipeline scripts call `download_studies_cached()`, which stops unless all
  studies arrived and caches the raw download in `rerun_outputs/*_download_raw.rds`.
- GEO serves a reCAPTCHA page (HTTP 200) to some clients. `.geo_cache/` in the Linux
  working tree holds 200 such pages and no records. `meta_analysis/14` now rejects them;
  on a blocked network it stops, so its output `sample_level_line_verification.csv`
  (built earlier from real GEO text: 46 immortalised, 194/200 agree) has to be carried
  over. One SOFT record per study (33) is in `DISCREPANCY_REPORT/evidence/geo_soft/` on
  the branch.
- `exploratory/34`, setting A ("pooled_controls") pools treated and control samples across
  studies, not only the controls.
- `exploratory/05, 08, 09, 10, 15, 18` read `partial_scores` files that only the
  single-cell code wrote. They are not on the manuscript path; `exploratory/21 mortality`
  reads only `mortality_partial_tage_ALL.csv`.
- The three `.pkl` clocks were pickled with scikit-learn 1.3.2. The environment runs 1.9.0;
  the tAge package patches `SimpleImputer` on load (`tAge/python/tage_predict.py`).
  `meta_analysis/05` and `temporal_analysis/04` previously ran through the single-cell
  project's Python (scikit-learn 1.8.0); they now use `.venv`.
- Legacy scripts in `meta_analysis/` (`metabolism_script.R`, `degs_v_stress.R`,
  `stress_responses.R`, …) read files from a previous machine (`/Users/ravelarvargas/…`)
  and are not part of the pipeline.
- `exploratory/41` is superseded by `42`; `exploratory/43` (GO annotation) and `45`
  (time-course monotonicity) are not cited.
