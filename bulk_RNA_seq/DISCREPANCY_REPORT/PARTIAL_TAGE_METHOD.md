# Pathway-level tAge: corrected to an exact decomposition method

## What changed and why

The original pathway-level analyses (`exploratory/08`, `09`, `10`) worked by **restricting the
expression matrix to one pathway's genes and rerunning `tAge_preprocessing()` +
`predict_tAge()` from scratch** on that tiny subset. This forced the EN model to run on
input where **only 0.1–1.5% of its 10,487 features were real** — the rest silently
imputed with a training-set constant (`SimpleImputer`). That constant is identical across
every sample in a comparison, so it doesn't bias a *group difference* (Cohen's d), but the
approach still means each pathway's preprocessing (RLE + YuGene normalization) was
re-derived from a handful of genes instead of the whole transcriptome — a real,
un-quantified source of instability, on top of running the model far outside how it was
validated.

**Checked whether the tAge authors solved this properly first** (per instruction, before
inventing a fix): the [tAge R package](https://github.com/Gladyshev-Lab/tAge) is
inference-only — no training or gene-subset-clock code. But the associated paper
(Tyshkovskiy, Gladyshev et al. 2026, *Nature*, ["Universal transcriptomic hallmarks of
mammalian ageing and mortality"](https://pmc.ncbi.nlm.nih.gov/articles/PMC13233323/))
describes exactly this problem's solution: **"partial tAge differences predicted using
only genes from the respective module."** Implemented that instead of a noise simulation
on the flawed approach.

## The method (exact, not approximate)

The EN model is a linear `ElasticNet` inside a `Pipeline` (impute → scale → select → predict), so:

```
prediction = intercept + Σ_i (coef_i × z_i)
```

where `z_i` is the imputed+standardized value of feature `i`. A pathway's exact
contribution to that sum is `Σ_i (coef_i × z_i)` over just the genes in that pathway — computed
**once**, from the correctly, whole-transcriptome-normalized data (the same
`tAge_preprocessing()` output already used for the validated, full-gene-set tAge scores),
with no re-imputation and no re-normalization on a restricted gene set.

Pipeline (`exploratory/14a*` → `14b` → `14c`/`14d`):
1. `14a_export_full_tage_matrices.R` / `14a2_..._bytimepoint.R` — run `tAge_preprocessing()`
   ONCE per group (5 meta-analysis conditions, 3 pooled temporal cell types, 9 per-timepoint
   temporal comparisons), export the `scaled_diff`/`yugene_diff` matrices that actually feed
   the model (samples × 18,696 mouse-ortholog-ID columns).
2. `14b_partial_tage_decompose.py` — loads each EN model, applies its own
   imputer+scaler, multiplies by `coef_`, and sums per Hallmark pathway (using a
   human-symbol→mouse-ID mapping built once via `tAge:::map_genes`, `exploratory/14a`'s
   `hallmark_pathway_mouse_ids.csv`).
3. **Verified exact, not just plausible**: reconstructing the full prediction (sum over
   *all* features + intercept) reproduces (a) a direct `model.predict()` call and (b) the
   original R `predict_tAge()` output already reported earlier in this conversation, to
   **floating-point precision (max diff = 0.0)**, across all 17 group comparisons × 2
   models.
4. `14c`/`14d` — Cohen's d / Wilcoxon per pathway per group from these exact partial
   scores; consolidated into `rerun_outputs/partial_tage_ALL.csv` (1,700 rows: 17 group
   comparisons × 50 Hallmark pathways × 2 EN models).

## Old (pathway-restricted rerun) vs. New (exact decomposition) — what held up, what didn't

| Group | Pathway | OLD d (padj) | NEW d (padj) | Verdict |
|---|---|---|---|---|
| Melanocyte | INTERFERON_ALPHA_RESPONSE | −7.29 (6.3e-5) | **−10.50** (5.3e-5) | Holds, even stronger |
| Melanocyte | DNA_REPAIR | 3.98 (6.3e-5) | 3.41 (5.3e-5) | Holds |
| OIS | NOTCH_SIGNALING | −1.75 (8.4e-12) | −1.74 (9.7e-12) | Holds, near-identical |
| RS | MITOTIC_SPINDLE | 2.77 (1.4e-6) | 1.98 (1.7e-5) | Holds, smaller but still strong |
| **SIPS** | **PI3K_AKT_MTOR_SIGNALING** | **1.13 (3.4e-10)** | **0.30 (2.1e-2)** | **Was substantially overstated** — real effect is much weaker, though still nominally significant |

**Practical implication**: most headline findings from the earlier exploratory work
survive being redone properly (some even strengthen). The one clear exception —
SIPS's PI3K-AKT-mTOR-signaling claim — should not be cited at its original effect size;
it's real but far smaller than first reported.

## Cross-analysis correlation, redone (Fibroblast/Keratinocyte/Melanocyte vs. meta-analysis conditions)

| | CICQ | SSCQ | RS | SIPS | OIS |
|---|---|---|---|---|---|
| Fibroblast | 0.37 | 0.54 | 0.41 | 0.71 | 0.60 |
| Keratinocyte | 0.41 | 0.39 | 0.09 | 0.59 | 0.47 |
| Melanocyte | −0.13 | 0.06 | −0.01 | 0.30 | 0.14 |

Same qualitative story as before the correction: temporal Fibroblast resembles the
meta-analysis fibroblast conditions (strongest to SIPS), Melanocyte resembles none of
them. Magnitudes shifted modestly but the conclusion is unchanged and, if anything, more
trustworthy now.

## Canonical files (supersede the `08`/`09`/`10`/`11`/`12`/`13` outputs)

- `rerun_outputs/partial_tage_ALL.csv` — the master table (all 17 groups × 50 pathways × 2 models).
- `rerun_outputs/pathway_effect_heatmap.png` — cross-analysis heatmap (`exploratory/15`).
- `rerun_outputs/pathway_divergence_meta_conditions.png` — within-meta-analysis divergence (`exploratory/16`).

The earlier `08`–`13` scripts/CSVs/figures are kept for the historical record (and because
most of their conclusions held up), but are no longer the source of truth — use
`partial_tage_ALL.csv` and the two regenerated figures above going forward.
