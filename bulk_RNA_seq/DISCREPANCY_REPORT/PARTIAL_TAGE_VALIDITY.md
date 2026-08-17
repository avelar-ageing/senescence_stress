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

## Test 1 — size-matched null: 10% survive

1,000 random equal-size draws from the clock's own 10,487 features, per set per group per model.
**Coverage completed 2026-08-17** — the first run covered the meta-analysis and pooled temporal
comparisons only (800 tests); the 9 per-timepoint comparisons, which the temporal write-up
actually leans on, were missing and have now been added. Full coverage is 1,700 tests
(meta 500, temporal pooled 300, temporal per-timepoint 900).

**171/1700 (10.1%)** reach p_emp < 0.05. Both-model survivors: 7 in the meta-analysis, 8 pooled
temporal, 17 per-timepoint — **32 of ~850 group x set combinations**.

Caveat that cuts the other way: the draw pool is genes already selected for age-association, so
this asks "does this set beat an arbitrary slice of the age clock", which is deliberately harsh,
and non-survival is not evidence of no effect.

### The null and the representation filter are in tension

**Every gene set that survives the null in 2 or more comparisons is POOR-tier:**

| set | comparisons survived | eff_n | tier |
|---|---|---|---|
| MITOTIC SPINDLE | 5 | 3.4 | POOR |
| **INTERFERON ALPHA RESPONSE** | **4** | 5.9 | POOR |
| APICAL JUNCTION | 3 | 6.3 | POOR |
| E2F TARGETS | 3 | 3.9 | POOR |
| G2M CHECKPOINT | 3 | 3.4 | POOR |
| HEME METABOLISM | 3 | 5.7 | POOR |
| ADIPOGENESIS | 2 | 7.1 | POOR |
| MYC TARGETS | 2 | 3.0 | POOR |
| UV RESPONSE UP | 2 | 5.0 | POOR |

No WELL or MODERATE set survives in more than one comparison. **Only 3 group x set combinations in
the entire analysis pass both filters**, and they share no pattern: meta OIS x XENOBIOTIC
METABOLISM, pooled Melanocyte x P53 PATHWAY, Melanocyte 4-day x COMPLEMENT.

*Mechanism UNVERIFIED.* Concentration does not simply inflate effect size (eff_n vs max |d|,
rho = -0.09, p = 0.53), so the tension is not explained by that. Whatever its cause, the practical
consequence is that **the two filters cannot both be satisfied**, and a set passing the null should
be read as reproducible-and-concentrated rather than as validated.

### The interferon case specifically

`HALLMARK INTERFERON ALPHA RESPONSE` in melanocytes survives the null in **all 8 melanocyte tests**
(pooled + 3 timepoints, x 2 models), p_emp 0.001-0.033, d = -7.2 to -22.6 against null 95th
percentiles of 5.3-9.9. By reproducibility it is the strongest result in the analysis after
MITOTIC SPINDLE. By representation it is POOR (eff_n 5.9; `Isg15` alone is 38-51% of the effect).

The two tests therefore disagree on this set more sharply than on any other, and the resolution is
the level of claim: **it does not support "the interferon pathway", but it does support a specific,
reproducible statement about `Isg15`, `Herc6` and `Usp18`.** Dropping it entirely would discard the
most consistently reproducible signal in the decomposition; keeping it as a pathway claim would
overstate a five-gene effect.

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

## Test 3b — which gene sets ARE well represented? (added 2026-08-17)

Test 3's medians hide a real gradient, so representation was computed per gene set:
`exploratory/14_pathway_representation.py` -> `rerun_outputs/pathway_representation.csv`.
Metric: **effective number of contributing genes**, eff_n = 1/sum(share^2) on mean |coef*z|
shares (inverse-Simpson). eff_n = 10 means the set behaves like ~10 equally weighted genes.
Tier = min(eff_n) across both models.

**8 WELL (eff_n >= 14), 7 MODERATE (>= 9), 35 POOR.** So most sets cannot support pathway-level
claims, but a real minority can.

| tier | gene set | eff_n | non-zero genes (s/y) | both-model-sig groups (of 17) | survives size null |
|---|---|---|---|---|---|
| WELL | EPITHELIAL MESENCHYMAL TRANSITION | 20.9 | 55/57 | 2 | 0 |
| WELL | ESTROGEN RESPONSE EARLY | 17.7 | 50/48 | 3 | 0 |
| WELL | **TNFA SIGNALING VIA NFKB** | 16.3 | 50/47 | **10** | 0 |
| WELL | ESTROGEN RESPONSE LATE | 16.2 | 44/37 | 2 | 0 |
| WELL | IL2 STAT5 SIGNALING | 16.2 | 35/47 | 1 | 0 |
| WELL | **MYOGENESIS** | 16.0 | 52/44 | **11** | 0 |
| WELL | **P53 PATHWAY** | 15.0 | 36/36 | **11** | **1** |
| WELL | HYPOXIA | 14.4 | 32/46 | 3 | 0 |
| MOD | **GLYCOLYSIS** | 14.0 | 30/34 | 5 | 0 |
| MOD | **XENOBIOTIC METABOLISM** | 12.2 | 38/43 | 6 | **1** |
| MOD | **Lysosomal Genes** | 12.0 | 22/34 | 5 | 0 |
| MOD | **FATTY ACID METABOLISM** | 10.5 | 24/28 | 6 | 0 |
| MOD | **APOPTOSIS** | 9.9 | 25/22 | 8 | 0 |
| MOD | **COMPLEMENT** | 9.9 | 35/35 | **13** | 0 |
| MOD | MTORC1 SIGNALING | 9.7 | 28/30 | 4 | 0 |

**The lysosomal set is defensible** — eff_n 12.0, 22/34 non-zero genes, significant on both models
in 5 of 17 group comparisons. It is better represented than most Hallmark sets.

**Representation and effect size are uncorrelated** (eff_n vs n both-model-sig groups,
Spearman rho = -0.11, p = 0.44; eff_n vs max |d|, rho = -0.09, p = 0.53). An earlier impression
that poorly-represented sets produce *more* findings was **not** supported when tested. The two
filters are simply independent.

**Consequence: the two quality filters select almost disjoint sets.** Of the size-matched-null
survivors, nearly all are POOR-tier: MITOTIC SPINDLE (eff_n 3.4, 3 survivors), G2M CHECKPOINT
(3.4, 2), HEME METABOLISM (5.7, 2), E2F TARGETS (3.9, 1), MYC TARGETS (3.0, 1), ADIPOGENESIS
(7.1, 1), APICAL JUNCTION (6.3, 1), INTERFERON ALPHA RESPONSE (5.9, 1), PROTEIN SECRETION
(2.4, 1). **Exactly two gene sets pass both filters: P53 PATHWAY (eff_n 15.0) and XENOBIOTIC
METABOLISM (eff_n 12.2).**

Sets to avoid entirely for pathway-level statements, being both poorly represented and previously
quoted as findings: PI3K AKT MTOR SIGNALING (eff_n **1.8**, and the basis of the proposed
"DNA damage drives a shared mTOR component" claim), PROTEIN SECRETION (2.4, top-5 share 100%),
REACTIVE OXYGEN SPECIES PATHWAY (2.4, 100%), WNT BETA CATENIN SIGNALING (2.7, only 4 non-zero
genes, 100%), HEDGEHOG SIGNALING (2.9, 96%), MYC TARGETS (3.0), MITOTIC SPINDLE (3.4),
G2M CHECKPOINT (3.4), SPERMATOGENESIS (3.4, 96%), E2F TARGETS (3.9).

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
