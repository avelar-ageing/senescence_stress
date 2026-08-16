# Pathway-level tAge: the "partial tAge" decomposition method

## What this is

Per-pathway tAge effects for the senescence/quiescence meta-analysis (CICQ/SSCQ/RS/SIPS/OIS
vs. pooled Proliferating) and the ERP021140 temporal time course (Fibroblast/Keratinocyte/
Melanocyte, pooled and per-timepoint), **adapted from** the module-level analysis in the paper
behind the tAge package (Tyshkovskiy, Gladyshev et al. 2026, *Nature*,
["Universal transcriptomic hallmarks of mammalian ageing and mortality"](https://pmc.ncbi.nlm.nih.gov/articles/PMC13233323/)).

We apply the paper's **gene-contribution definition** (coefficient x differential expression;
their Fig. 3g "logFC x clock coefficient") but **aggregate over MSigDB Hallmark pathways using
the fitted global clock**, rather than training a separate clock per WGCNA module as the paper
does. The paper's phrase "partial tAge differences predicted using only genes from the
respective module" refers to a module clock's own prediction, not a post-hoc decomposition of
the global model. This is a defensible adaptation, not a reproduction -- see
`PARTIAL_TAGE_VS_PAPER.md` for the full comparison and for three consequences of using
overlapping curated sets instead of near-disjoint data-derived modules.

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

**Pathway list (step 3)**: read directly from `rerun_outputs/stress_response_pathways_RERUN.csv`
— the exact pathway × gene table `meta_analysis/03_hallmark_enrichment.R` already builds for the
DEG-overlap enrichment (MSigDB Hallmark with `_V1`/`_V2` merged, e.g. `MYC_TARGETS_V1` +
`MYC_TARGETS_V2` → one `HALLMARK MYC TARGETS` pathway, plus the custom **"Lysosomal Genes"** set
recovered from `Final/SI_tables/lyso_genes.csv`, 191 genes, `human_pc`-filtered). 49 Hallmark
terms + Lysosomal Genes = 50 pathways. This is read from that single source of truth rather than
re-derived here, so the tAge partial-decomposition pathways and the DEG-overlap enrichment
pathways can never drift apart.
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
| Fibroblast | 0.31 | 0.51 | 0.45 | **0.76** | 0.65 |
| Keratinocyte | 0.36 | 0.38 | 0.06 | 0.56 | 0.45 |
| Melanocyte | −0.13 | 0.03 | 0.00 | 0.31 | 0.15 |

Temporal Fibroblast's pathway-effect signature most resembles the meta-analysis's own
fibroblast conditions — strongest to SIPS (both acute-stress-induced), weakest to CICQ
(a mechanistically different trigger, contact inhibition). Melanocyte resembles none of
the five meta-analysis conditions.

**Melanocyte's distinct signature**: `HALLMARK INTERFERON ALPHA RESPONSE` is its strongest
single pathway effect (Cohen's d ≈ −10.5, yugene_diff model, pooled irradiated vs. none),
present already at 4 days and sustained through 20 days (per-timepoint breakdown in
`partial_tage_ALL.csv`, `analysis == "temporal_bytimepoint"`) — absent from Fibroblast's
and Keratinocyte's top pathways in either direction.

**mTOR signaling across meta-analysis conditions**: `HALLMARK PI3K AKT MTOR SIGNALING` is
strongest in SIPS (Cohen's d ≈ 0.30, yugene_diff — nominally significant but a modest
effect, not a dominant one), present but weaker in OIS, and not detectable in RS, CICQ,
or SSCQ.

**Lysosomal Genes** (recovered custom pathway, not in default MSigDB Hallmark): one of the
strongest and most consistent effects in the whole panel — significant (padj<0.05) in 6/8
of the primary meta-analysis + pooled-temporal groups, yugene_diff model. Effect size is
particularly large in temporal Fibroblast, and grows with time (Cohen's d = 2.99 at 4 days,
3.04 at 10 days, 5.75 at 20 days, irradiated vs. none) — consistent with progressive
lysosomal/autophagic dysfunction accumulating over the senescence time course. Also
significant across all five meta-analysis conditions (CICQ/SSCQ/RS/SIPS/OIS, d = 1.11–1.88)
and in pooled temporal Fibroblast (d = 2.92), but not significant in Keratinocyte or
Melanocyte (pooled or per-timepoint, except a small negative effect at Melanocyte 10 days).

## Canonical files

- `rerun_outputs/partial_tage_ALL.csv` — the master table (all 17 groups × 50 pathways × 2 models).
- `rerun_outputs/pathway_effect_heatmap.png` — cross-analysis heatmap.
- `rerun_outputs/pathway_divergence_meta_conditions.png` — within-meta-analysis divergence.
- `rerun_outputs/partial_tage/` — intermediate exported matrices, pathway-mapping table,
  and per-sample partial scores (inputs to the two figures above).
