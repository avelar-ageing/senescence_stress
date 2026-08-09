# Pathway-level tAge: the "partial tAge" decomposition method

## What this is

Per-pathway tAge effects for the senescence/quiescence meta-analysis (CICQ/SSCQ/RS/SIPS/OIS
vs. pooled Proliferating) and the ERP021140 temporal time course (Fibroblast/Keratinocyte/
Melanocyte, pooled and per-timepoint), using the method described in the paper behind the
tAge package (Tyshkovskiy, Gladyshev et al. 2026, *Nature*,
["Universal transcriptomic hallmarks of mammalian ageing and mortality"](https://pmc.ncbi.nlm.nih.gov/articles/PMC13233323/)):
**"partial tAge differences predicted using only genes from the respective module."**

## The method

The EN tAge model is a linear `ElasticNet` inside a `Pipeline` (impute → scale → select →
predict):

```
prediction = intercept + Σ_i (coef_i × z_i)
```

where `z_i` is the imputed+standardized value of feature `i`. A pathway's exact
contribution to that sum is `Σ_i (coef_i × z_i)` restricted to just the genes in that
pathway — computed from the correctly, whole-transcriptome-normalized data (the same
`tAge_preprocessing()` output the full-gene-set tAge scores use), so every real,
measured gene contributes its actual value and nothing needs to be imputed beyond what
the original full analysis already required.

**Verified exact**: reconstructing the full prediction (sum over *all* features +
intercept) reproduces a direct `model.predict()` call to floating-point precision
(max diff = 0.0), across all 17 group comparisons × 2 EN models (`scaled_diff`,
`yugene_diff`) × 50 MSigDB Hallmark pathways.

## Pipeline (`exploratory/01`–`07`)

1. **`01_export_full_tage_matrices.R`** — runs `tAge_preprocessing()` once per group
   (5 meta-analysis conditions together, 3 pooled temporal cell types), exports the
   `scaled_diff`/`yugene_diff` matrices actually fed to the model (samples × 18,696
   mouse-ortholog-ID columns) plus group labels.
2. **`02_export_temporal_bytimepoint_matrices.R`** — same, at per-timepoint resolution
   (3 cell types × 3 timepoints, 6 vs. 6 samples each, vs. that cell type's own baseline).
3. **`03_build_pathway_mouse_id_mapping.R`** — maps every MSigDB Hallmark pathway's human
   gene symbols to the mouse-ortholog IDs the model actually uses (via `tAge:::map_genes`,
   the same fixed lookup table `tAge_preprocessing()` uses internally). Computed once,
   reused by every group.
4. **`04_partial_tage_decompose.py`** — for each EN model, loads its fitted pipeline,
   applies its own imputer + scaler to the exported matrix, multiplies by the model's
   coefficients, and sums per pathway. Outputs one partial score per sample per pathway,
   plus the exact-reconstruction sanity check.
5. **`05_consolidate_partial_scores.R`** — Cohen's d / Wilcoxon (BH-adjusted within each
   analysis × model) per pathway per group, across all 17 group comparisons. Output:
   `rerun_outputs/partial_tage_ALL.csv` (1,700 rows: 17 groups × 50 pathways × 2 models).
6. **`06_pathway_effect_heatmap.R`** — cross-analysis figure (meta-analysis conditions +
   temporal cell types side by side). Pathway selection: significant (padj<0.05) in ≥6 of
   the 8 groups — a topic-blind recurrence filter, not a hand-picked list; rows
   hierarchically clustered.
7. **`07_pathway_divergence_meta_conditions.R`** — which pathways most separate the 5
   meta-analysis conditions *from each other* (not vs. Proliferating in general): top 20
   ranked by range of Cohen's d across the 5 conditions.

## Headline results

**Cross-analysis correlation** (Spearman, Cohen's d, all 50 pathways, yugene_diff model —
does each temporal cell type's pathway-effect profile resemble any meta-analysis
condition's?):

| | CICQ | SSCQ | RS | SIPS | OIS |
|---|---|---|---|---|---|
| Fibroblast | 0.37 | 0.54 | 0.41 | **0.71** | 0.60 |
| Keratinocyte | 0.41 | 0.39 | 0.09 | 0.59 | 0.47 |
| Melanocyte | −0.13 | 0.06 | −0.01 | 0.30 | 0.14 |

Temporal Fibroblast's pathway-effect signature most resembles the meta-analysis's own
fibroblast conditions — strongest to SIPS (both acute-stress-induced), weakest to CICQ
(a mechanistically different trigger, contact inhibition). Melanocyte resembles none of
the five meta-analysis conditions.

**Melanocyte's distinct signature**: `HALLMARK_INTERFERON_ALPHA_RESPONSE` is its strongest
single pathway effect (Cohen's d ≈ −10.5, yugene_diff model, pooled irradiated vs. none),
present already at 4 days and sustained through 20 days (per-timepoint breakdown in
`partial_tage_ALL.csv`, `analysis == "temporal_bytimepoint"`) — absent from Fibroblast's
and Keratinocyte's top pathways in either direction.

**mTOR signaling across meta-analysis conditions**: `HALLMARK_PI3K_AKT_MTOR_SIGNALING` is
strongest in SIPS (Cohen's d ≈ 0.30, yugene_diff — nominally significant but a modest
effect, not a dominant one), present but weaker in OIS, and not detectable in RS, CICQ,
or SSCQ.

## Canonical files

- `rerun_outputs/partial_tage_ALL.csv` — the master table (all 17 groups × 50 pathways × 2 models).
- `rerun_outputs/pathway_effect_heatmap.png` — cross-analysis heatmap.
- `rerun_outputs/pathway_divergence_meta_conditions.png` — within-meta-analysis divergence.
- `rerun_outputs/partial_tage/` — intermediate exported matrices, pathway-mapping table,
  and per-sample partial scores (inputs to the two figures above).
