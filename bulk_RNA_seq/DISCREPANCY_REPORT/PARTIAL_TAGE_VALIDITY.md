# Is the Hallmark partial-tAge decomposition valid? Three tests

Written 2026-08-17 to answer the question directly, rather than by argument from method
description. Prompted by the scRNA arm noting our approach does not match the paper's.

**Short answer: the quantity is exactly what we say it is, but it does not support
pathway-level biological claims. It supports gene-level ones.**

---

## What is not in doubt

1. **The arithmetic is an identity, not an approximation.** `pred = intercept + Σ coef_i·z_i`, so
   `Σ_{i∈pathway} coef_i·z_i` is literally the part of that sample's predicted age contributed by
   those genes. Reconstruction over all features reproduces `model.predict()` to max diff 0.0 across
   every group, model and resolution.
2. **The gene-level quantity is the paper's own.** Their Fig. 3g is "logFC × clock coefficient" and
   Fig. 4c "slope × clock coefficient". We compute coefficient × control-subtracted expression. Same
   measure. We aggregate it over a gene set; TACO's own "Pathway Contribution" tab indicates the
   authors aggregate it over pathways too (though whether by summation is UNVERIFIED).
3. **It differs from the paper's estimator.** They *train* a clock per WGCNA module; we *restrict* a
   fitted global clock. Different estimators, and they can disagree. This is a real difference but
   not by itself an error — see below for what actually limits us.

## Test 1 — size-matched null: 9% survive

1,000 random equal-size draws from the clock's own 10,487 features, per pathway per group per
model. Only **72/800 (9%)** reach p_emp < 0.05; **15/400** group–pathway pairs clear it on both
models. Caveat that cuts the other way: the draw pool is genes already selected for
age-association, so this asks "does this set beat an arbitrary slice of the age clock", which is
deliberately harsh, and non-survival is not evidence of no effect.

## Test 2 — convergent validity against DEG enrichment: weak

The same project already has Hallmark overrepresentation of each condition's DEGs
(`hallmark_vs_arrest_degs_RERUN.csv`) — an independent instrument, same data, same gene sets. If a
pathway's partial-tAge effect tracks real differential expression of that pathway, the two should
agree. Across 250 condition × pathway pairs:

| | scaled | yugene |
|---|---|---|
| \|partial d\| vs DEG −log10(padj), Spearman | +0.090 (p = 0.16) | +0.127 (p = 0.044) |
| \|partial d\| vs max odds ratio | +0.012 | +0.037 |
| median \|d\|, DEG-enriched vs not | 0.68 vs 0.59 (p = 0.068) | 0.73 vs 0.53 (p = 0.0095) |

Directionally correct — DEG-enriched pathways do have larger partial effects — but the association
is weak. Partly expected, since partial tAge depends on coefficient magnitude and not only on
expression change. Still, ρ ≈ 0.1 is not the convergence you would want before making
pathway-level claims.

## Test 3 — where the signal actually comes from: a handful of genes

This is the finding that decides the question.

| | scaled | yugene |
|---|---|---|
| clock genes per pathway (median) | 125 | 125 |
| of those, **non-zero coefficient** | **20 (21%)** | **22 (24%)** |
| share of pathway \|contribution\| from its **single largest gene** (median) | **26.9%** | **23.4%** |
| same, maximum across pathways | 73.4% | 62.3% |
| share from **top 5 genes** (median) | **72.9%** | **64.4%** |
| pathways where one gene is >25% of signal | 28/50 | 24/50 |

So ~79% of a pathway's genes contribute nothing (elastic net zeroed them), and a median of five
genes carry roughly two thirds to three quarters of what remains. **"Pathway X contributes N to
tAge" is in practice "these five genes contribute N".** This is the elastic-net redundancy problem
made concrete: with correlated predictors, EN assigns weight semi-arbitrarily among them, so which
genes carry a pathway's score is partly an artefact of fitting, and a pathway whose signal is
redundant with a non-pathway gene can be silenced entirely.

Note this is one thing a module-*trained* clock would partly avoid, since refitting within a module
spreads weight across that module's genes instead of letting global competition concentrate it.
That is the substantive methodological point behind the scRNA arm's observation — not the
gene-contribution definition, which is fine.

## The same test on our most robust finding

Melanocyte `HALLMARK INTERFERON ALPHA RESPONSE`, the one effect that survives the size-matched null
(d = −10.50, null q95 5.26, p_emp = 0.005). Between-group contribution difference, irradiated vs none:

| | scaled | yugene |
|---|---|---|
| clock genes / non-zero | 66 / 13 | 66 / 15 |
| top gene share | **38.2%** | **50.9%** |
| top 5 share | **82.1%** | **91.4%** |

Leading genes, **consistent across both models**: `Isg15` (largest in both), `Herc6`, `Usp18`,
plus `Helz2`, `Cd74`, `Lgals3bp` (scaled) and `Csf1`, `Irf9` (yugene).

Every leading gene is a canonical interferon-stimulated gene. So the biology is coherent and the
effect is genuinely interferon-specific — but it is carried by a small ISG set led by `Isg15`, not
by "the interferon pathway" as a whole.

## Conclusion and recommendation

- **Valid as:** a descriptive decomposition of the clock's own prediction, and a
  hypothesis-generating screen that points at specific genes.
- **Not valid as:** evidence that a pathway drives transcriptomic ageing, nor a variance partition,
  nor a substitute for the paper's module clocks.
- **Recommended framing:** demote pathway-level results to supplementary/hypothesis-generating.
  Where a finding is kept in the main text, name the genes carrying it and report the top-gene
  share, so the reader can see the claim's actual support. The melanocyte result is strong enough to
  keep on those terms — reframed as an ISG-led effect led by `Isg15`, which is a more specific and
  more checkable claim than the pathway-level version.
- **Always:** "the clock's prediction restricted to these genes moves this way", never "this pathway
  activates/suppresses". Direction of gene expression has not been checked for any pathway.

## Reproduce

Tests 2 and 3 were run ad hoc; the numbers above are from `rerun_outputs/partial_tage_ALL.csv`,
`hallmark_vs_arrest_degs_RERUN.csv`, `partial_tage/hallmark_pathway_mouse_ids.csv` and the two
`.pkl` models. Test 1 is `exploratory/12_partial_tage_size_matched_null.py` →
`partial_tage_size_matched_null.csv`. Worth promoting tests 2 and 3 to a numbered script if these
figures are cited.
