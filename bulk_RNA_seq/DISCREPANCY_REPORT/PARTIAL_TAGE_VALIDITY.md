# Is the Hallmark partial-tAge decomposition valid? Four tests

Written 2026-08-17 to answer the question directly, rather than by argument from method
description. Prompted by the scRNA arm noting our approach does not match the paper's.

**Short answer: the quantity is exactly what we say it is. Whether it supports a pathway-level
claim depends on the gene set — 15 of 50 sets are adequately represented in the clock, the other 35
are carried by <=5 genes and support only gene-level statements.**

*Note: Test 1 below (size-matched null) was subsequently EXCLUDED — see its header for why. The
conclusions here rest on Tests 2, 3 and 3b.*

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

## Test 1 — size-matched null: EXCLUDED (decision 2026-08-17)

> **This test is excluded from the analysis and should not be cited.** It was run to completion
> (1,700 tests, `exploratory/12_partial_tage_size_matched_null.py`,
> `rerun_outputs/partial_tage_size_matched_null.csv`) and is kept for provenance only. Reasons:
>
> 1. **The draw pool is not a neutral background.** Draws come from the clock's own 10,487
>    features — genes selected precisely *because* they change with age across mammals. The test
>    therefore asks "is this set more informative than an arbitrary slice of an age clock", which
>    is a specificity question, not the question any Results claim makes.
> 2. **The null is not centred on no effect.** Because every clock feature carries age signal and
>    the groups genuinely differ, random draws already separate them: median null mean |d| was
>    0.55-0.80 (meta) and 1.0-2.1 (temporal). A "null" distribution centred on a real effect does
>    not test what a null is supposed to test.
> 3. **Empirically it selects for the wrong property.** Every gene set surviving in >=2 comparisons
>    is POOR-represented (MITOTIC SPINDLE eff_n 3.4 with 5 survivals, INTERFERON ALPHA 5.9 with 4,
>    APICAL JUNCTION / E2F / G2M / HEME 3 each). **No WELL or MODERATE set survives more than
>    once**, and only 3 of ~850 group x set combinations pass both this test and the representation
>    filter. A test that anti-correlates with the property we actually care about is not a quality
>    filter. (Mechanism unverified; concentration does not simply inflate effect size, eff_n vs
>    max |d| rho = -0.09, p = 0.53.)
>
> **What is lost by excluding it:** no guard remains against the size dependence that motivated it
> (signed d vs set size, rho = 0.325). Test 3/3b's representation filter partly covers this — it
> directly measures whether a set's signal is carried by few genes — but it does not ask whether an
> equally-sized arbitrary set would do as well.
>
> **The replacement that does answer the Results question is a label-permutation null**: hold the
> gene set fixed (removing size and concentration confounds entirely) and permute group labels,
> testing whether the set separates *these groups* more than chance. It requires no choice of
> background gene pool. **Proposed, not yet run — this is the top outstanding methodological item.**
>
> **Knock-on effect on conclusions:** the retraction of the "universal / shared core" framing
> originally rested on this null. That retraction still stands, but now on representation grounds
> instead — of the seven sets in that framing, only P53 PATHWAY (eff_n 15.0) and COMPLEMENT (9.9)
> clear the representation filter; E2F TARGETS (3.9), PI3K AKT MTOR SIGNALING (1.8), SPERMATOGENESIS
> (3.4), HEDGEHOG SIGNALING (2.9) and UV RESPONSE UP (5.0) do not.

## Test 1 (excluded). Original text: size-matched null, 10% survive

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

## Test 4 — label-permutation null: fair, but adds nothing (run 2026-08-17)

`exploratory/15_partial_tage_label_permutation.py` -> `rerun_outputs/partial_tage_label_permutation.csv`.
2,000 permutations of the group labels per set per group per model, gene set held fixed, so there is
no background pool to choose and none of Test 1's unfairness. 1,700 tests.

**It is redundant with what we already had.** The Wilcoxon rank-sum test *is* an exact permutation
test, differing only in statistic (ranks vs means), so the two share a null:

| comparison vs the existing Wilcoxon p | result |
|---|---|
| Spearman rho(p_permutation, uncorrected wilcox_p) | **0.956** |
| same call at 0.05 (uncorrected) | **94.8%** |
| same call vs BH-adjusted p_adj | 93.4% |

**And it is blind to the property that actually matters.** Pass rates are identical across
representation tiers — WELL 64.0%, MODERATE 65.5%, POOR 63.8% — because the set is held fixed, so
neither size nor concentration can influence the result. It therefore does **not** substitute for
what Test 1 was attempting.

**What it is worth:** a confirmation that the Wilcoxon-based findings are not artefacts of the rank
statistic or of any asymptotic approximation. That is a real if modest contribution, and it means
the BH-adjusted Wilcoxon results can be reported as-is without a distribution-free caveat.

Under this test the well-represented sets perform strongly — both-model survivors out of 17
comparisons: COMPLEMENT 16, P53 PATHWAY 12, TNFA SIGNALING VIA NFKB 12, APOPTOSIS 11,
MYOGENESIS 11, FATTY ACID METABOLISM 7, GLYCOLYSIS 6, MTORC1 SIGNALING 6, XENOBIOTIC METABOLISM 6,
HYPOXIA 5, Lysosomal Genes 5.

### Consequence for the interferon decision

Melanocyte `HALLMARK INTERFERON ALPHA RESPONSE` passes decisively under the fair test:
p_permutation = 0.0005 pooled (both models) and 0.0011-0.0043 at every individual timepoint, against
a floor of 0.0022 for 6v6. So **the only argument against it is representation** (eff_n 5.9; `Isg15`
alone 38-51% of the effect), which is an argument about the *level* of the claim, not about whether
the effect is real. Recommended resolution unchanged: report it as a gene-level finding
(`Isg15`, `Herc6`, `Usp18`), not as "the interferon pathway".

### Net position on filters

Neither null discriminates usefully: Test 1 is unfair and selects *for* gene concentration; Test 4
is fair but redundant and blind to concentration. **The representation filter (Tests 3/3b) is the
only quality filter that discriminates**, so pathway-level claims should rest on: both-model
agreement in sign, BH-adjusted Wilcoxon significance, and eff_n tier. Nothing else earns its place.

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
