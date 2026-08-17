# Bulk tAge work: what was done and what was found

Session handoff, written 2026-08-17. Branch `restructure-bulk-rna-seq`, commits
`2e8bdfe`..`0152457` (21 commits). Scope: the bulk RNA-seq tAge arm — universal
(whole-transcriptome) tAge and the pathway-level decomposition, meta-analysis and
ERP021140 time course. Written to be compared against a parallel write-up from another
instance, so it separates **what changed in the code**, **what was measured**, and **what
is still open** — and flags where my own earlier statements were wrong.

---

## 1. A real bug found and fixed first (`2e8bdfe`)

The pathway-level pipeline was using a **different gene-set list** from the DEG-enrichment
arm of the same project. `exploratory/03_build_pathway_mouse_id_mapping.R` pulled raw
msigdbr Hallmark, which meant:

- `MYC_TARGETS_V1` and `MYC_TARGETS_V2` stayed as two separate pathways, where the
  enrichment arm merges them into one `HALLMARK MYC TARGETS`
- a curated **191-gene lysosomal set was absent entirely** (it had been written off as
  unrecoverable; it was in fact sitting in `Final/SI_tables/lyso_genes.csv`)

Fix: `03_...R` now reads the pathway x gene table from
`rerun_outputs/stress_response_pathways_RERUN.csv` — the file
`meta_analysis/03_hallmark_enrichment.R` already produces — so the two arms cannot drift.
Also fixed a knock-on bug this exposed: with space-separated pathway names,
`read.csv()` silently mangles `"HALLMARK MYC TARGETS"` to `"HALLMARK.MYC.TARGETS"`, which
broke the label regex in scripts 06/07 (`check.names = FALSE` + regex fixed).

Whole `exploratory/01–07` pipeline re-run. **Lysosomal Genes turned out to be one of the
strongest and most consistent effects in the Hallmark panel** — significant in 6/8 primary
groups, and in temporal fibroblasts it grows with time (Cohen's d = 2.99 / 3.04 / 5.75 at
4 / 10 / 20 days). Entirely invisible before the fix. Cross-analysis correlations shifted
slightly (e.g. Fibroblast–SIPS 0.71 → 0.76).

**Note for comparison:** any number produced from `partial_tage_ALL.csv` *before*
`2e8bdfe` is stale.

## 2. Statistical families were mixed; now separated

Two things were wrong in how p-values were being quoted:

- **Meta-analysis.** Two BH families exist: each condition vs Proliferating (5 x 2 models
  = 10 tests, `tage_wilcoxon_vs_proliferating.csv`) and all-pairs between conditions
  (15 pairs x 2 = 30 tests, `tage_pairwise_all_conditions.csv`). The same underlying
  Wilcoxon p appears in both with **different** FDRs (SSCQ vs Proliferating on
  scaled_diff: raw p = 0.206 → FDR 0.206 in the 10-test family, 0.309 in the 30-test
  family). I had quoted 0.31 in one paragraph and 0.21 in another. Convention now:
  vs-Proliferating claims come from the 10-test family, condition-vs-condition claims from
  the 30-test family, never mixed.
- **Time course.** Only vs-baseline tests existed (18-test family). Between-timepoint
  contrasts — which every trajectory claim depends on — were **absent**. Added
  `temporal_analysis/07_tage_temporal_pairwise.R`: all 6 timepoint pairs x 3 cell types
  x 2 models = 36 tests, one BH family. This **subsumes** the vs-baseline comparisons
  (none-vs-4/10/20d are 3 of the 6 pairs), so the 18-test family was deleted along with
  the figure script built on it (`05_tage_spread_figure.R`, superseded by `06`); its
  Kruskal–Wallis omnibus was carried into `07` so nothing was lost. FDRs shift slightly
  under the wider correction (Fibroblast 4d scaled: 0.010 → 0.014). Commits `5fc3245`,
  `26bbef9`.

## 3. Figures built

All in `rerun_outputs/`, mirrored to `DISCREPANCY_REPORT/evidence/figures/`.

| figure | script | content |
|---|---|---|
| `figure_universal_tage_differences.png` | `meta_analysis/08_...R` | tAge by condition, both models, vs-Proliferating + SSCQ-pairwise brackets |
| `figure_pathway_heatmap_all_both_models.png` | `exploratory/08_...R` | all 50 pathways x 5 conditions x 2 models |
| `figure_temporal_universal_tage.png` | `temporal_analysis/06_...R` | tAge trajectory, 3 cell types x 2 models, median trend line |
| `figure_temporal_pathway_heatmap.png` | `exploratory/09_...R` | 50 pathways x 3 cell types x 2 models (pooled) |
| `figure_temporal_pathway_heatmap_bytimepoint.png` | `exploratory/10_...R` | 50 pathways x 3 cell types x 3 timepoints x 2 models |

Two decisions worth recording. **pheatmap was dropped for `geom_tile`** — with clustering
off it contributed nothing but its annotation strips, which `facet_grid` reproduces with
better control. **Row clustering was removed in favour of alphabetical order**: these
figures are comprehensive lookup references (nothing filtered out), so being able to find
a named pathway beats a data-driven row order, and with 50 mostly-correlated columns the
cluster order was not stable enough to carry weight.

## 4. A hard statistical limit at 6v6

In the per-timepoint analysis every comparison is 6 vs 6. The smallest attainable
two-sided Wilcoxon p is `2/C(12,6) = 2.16e-3` **even under complete rank separation**, so
after BH across the 900-test family nothing can fall below FDR ~4.6e-3 and only **19
distinct raw p-values** exist across all 900 tests. `***` is impossible there by
construction. Effect size is the only thing carrying resolution at that level, and stars
in that figure are **not** comparable to the pooled figures where n is doubled. Documented
in the script header.

## 5. The size-matched null, and what it broke (`8e2376d`)

`PARTIAL_TAGE_VS_PAPER.md` §5 flagged that partial-tAge effect size partly tracks how many
clock features a set contains (confirmed: Spearman rho = **0.325**). Implemented
`exploratory/12_partial_tage_size_matched_null.py` — 1,000 random equal-size draws from the
clock's own 10,487 features, per set per group per model.

**Result: only 72/800 tests (9.0%) reach p_emp < 0.05, and only 15 group–pathway pairs
clear it on both models.** Survivors:

| group | surviving both models |
|---|---|
| meta OIS | APICAL JUNCTION, MYC TARGETS, XENOBIOTIC METABOLISM |
| meta RS | MITOTIC SPINDLE, PROTEIN SECRETION |
| meta SIPS | HEME METABOLISM, MITOTIC SPINDLE |
| temporal Fibroblast | G2M CHECKPOINT, HEME METABOLISM, MITOTIC SPINDLE |
| temporal Keratinocyte | ADIPOGENESIS, E2F TARGETS, G2M CHECKPOINT |
| temporal Melanocyte | **INTERFERON ALPHA RESPONSE**, P53 PATHWAY |

**Consequence — this retracts something I had already drafted.** I had written a
"shared / universal core" framing for §2.1.5.2 and §2.1.5.4 (COMPLEMENT, UV RESPONSE UP,
P53 PATHWAY, E2F TARGETS, PI3K AKT MTOR SIGNALING, SPERMATOGENESIS, HEDGEHOG SIGNALING as
cell-type-independent). **It does not survive size matching** and must be withdrawn or
downgraded to "consistent in direction across cell types, but not exceeding a size-matched
null".

**The melanocyte result strengthens instead.** Melanocyte `INTERFERON ALPHA RESPONSE`
is d = **-10.50** against a null 95th percentile of 5.26 (p_emp = 0.005) — the clearest
pathway-specific signal anywhere in the dataset. Melanocyte P53 PATHWAY also survives
(d = 5.79 vs null q95 5.41, p = 0.039).

**Important framing correction (the user pushed back on this and was right).** The draw
pool is the clock's own features — genes *already selected for age-association across
mammals*. So this is not a random-gene background; it asks "does this set beat an arbitrary
slice of the age clock itself", which is deliberately harsh. Structurally there is no
alternative pool, since a non-clock gene has coefficient exactly 0 and contributes nothing.
**Non-survival is not evidence of no effect.** Both facts must be reported together.

## 6. Two numbers in `PARTIAL_TAGE_VS_PAPER.md` §5 did not reproduce

Checked all three "measured" consequences against `rerun_outputs/`:

| claim | stated | measured | verdict |
|---|---|---|---|
| partials sum vs full tAge | 80.6 vs 33.7, ~2.4x | ratio unstable: median 0.96 scaled, 1.76 yugene; only 4,251/10,487 clock features are in any Hallmark set | **corrected** |
| inter-pathway median \|rho\| | 0.545, 34% > 0.7 | **0.20** per-sample (scaled 0.212, yugene 0.201), 1.5–2.3% > 0.7. Highest alternative reading (effect-size profiles across the 17 comparisons) = 0.35–0.37, 14–16% | **corrected** |
| \|d\| vs set size | rho = 0.33 | rho = **0.325** | confirmed |

Qualitative conclusions all stand; the dependence is milder than stated, so BH across 50 is
mildly optimistic rather than badly so.

## 7. Overstated alignment with the paper, now reworded

`PARTIAL_TAGE_METHOD.md` claimed we use "the method described in the paper". That is wrong:
the paper **trains a separate elastic-net clock per WGCNA module**; we **post-hoc restrict a
fitted global clock** over curated overlapping Hallmark sets. What we genuinely share is the
*gene-level contribution definition* (coefficient x differential expression, their Fig. 3g).
Reworded to "adapted from", with the difference stated explicitly.

## 8. The paper's modules: both routes tried, neither usable

Full detail in `MODULE_GROUPING_FINDINGS.md`. Summary:

**Route A — the paper's module CLOCKS. Blocked.** Supplementary Table 5 publishes per-gene
coefficients and intercepts, which looks like the complete linear model, but the fitted
imputer/scaler is **not** published, and EN coefficients are meaningless without the
standardisation they were fitted against. Evaluating them fails a positive control on the
**tAge package's own example mouse data** (r = 0.48 rodent, 0.22 multispecies against the
package's own `predict_tAge`), and on our data reverses the sign of the SIPS and OIS
elevations. Not our setup, not an input convention (raw vs scaled agree to 4 d.p.). Needs
the authors' scaler statistics.

**Route B — module MEMBERSHIPS as a grouping for our own verified decomposition. Ran
cleanly (26 decompositions, all exact 0.0) and performed worse than Hallmark:**

| | paper modules | Hallmark |
|---|---|---|
| passes size-matched null | 17/224 (7.6%) | 72/800 (9.0%) |
| **survives on both models** | **1/112 (0.9%)** | **15/400 (3.8%)** |
| median \|d\| | 0.86 | 0.87 |
| median \|d\| / null q95 | 0.40 | 0.39 |
| cross-model sign agreement | **0.55** | 0.775 |
| inter-set median \|rho\| | 0.20 | 0.201 |

Raw effect sizes are indistinguishable; the entire deficit is **reproducibility** — 4x
fewer both-model survivors, sign agreement down from 78% to 55%, with some modules flipping
sign while both models are significant (Fibroblast "Fatty acid metabolism/Peroxisome":
d = -6.18 scaled vs +10.20 yugene). Probable cause: modules span only 12% of the clock's
features, so model-specific coefficient differences average down far less.

**Two findings from this that matter beyond the module question:**

1. **Disjointness bought nothing.** The modules are *perfectly* disjoint (zero shared
   genes) yet show the same inter-set median |rho| (0.20) as overlapping Hallmark. So
   set-to-set dependence is **not caused by gene overlap** — it is shared sample-level
   structure. §5 consequence (2) is misattributed.
2. **The module annotations do not describe their gene content.** `HALLMARK INTERFERON
   ALPHA RESPONSE` and the paper's `darkmagenta / Interferon signaling` module share
   **zero** genes. Canonical ISGs sit elsewhere: **Isg15 and Irf7 in "Chromatin
   modification"**, Rsad2 in "mRNA splicing", Usp18/Ifih1 in "Mitochondrial translation".
   No module is enriched for its own annotation's markers above background. Verified not
   an extraction bug — IDs resolve (1,248/1,248), entrez↔symbol 98.5% matches the package
   table (rest are synonyms), sheet layout read directly from MOESM7 and matches the
   parser's assumption, labels/sizes agree with the official dictionary (MOESM8 sheet D).
   These are co-expression modules whose annotation is a top-enrichment label, not a
   membership criterion.

**Therefore the melanocyte interferon finding can be neither confirmed nor refuted by the
modules**, and must not be cited as corroborated by them.

## 9. Still open

1. **Label-permutation null — NOT RUN.** Permute group labels, recompute d for the real
   set. This asks the question the results sections actually make ("does this set separate
   the groups more than chance") rather than the much harsher size-matched question, and
   would be a fairer complement for both arms. This is the top outstanding item.
2. **§2.1.5.2 and §2.1.5.4 need rewriting** with the null integrated and the
   universal-core claims withdrawn/downgraded. Drafts exist but are superseded by §5 above.
3. **Rodent module panel not examined** for the same annotation/content mismatch (only
   multispecies was).
4. **Route A remains open only if** someone obtains the module clocks' scaler statistics
   from the authors.
5. **Interpretation guard, applies throughout:** partial-tAge scores are a set's
   *model-weighted contribution to the age prediction*, not pathway activity. A positive
   score can come from up-regulated positively-weighted genes, down-regulated
   negatively-weighted genes, or both. All text says "raises/lowers the tAge contribution",
   never "activates/suppresses". Verifying an activation claim needs the raw expression
   direction, which has not been done for any pathway.

## 10. Canonical files

**Universal tAge:** `tage_all_conditions.csv`, `tage_wilcoxon_vs_proliferating.csv`
(10-test family), `tage_pairwise_all_conditions.csv` (30-test family),
`tage_spread_by_condition.csv`, `tage_temporal_by_celltype.csv`,
`tage_temporal_pairwise_all_timepoints.csv` (36-test family, canonical for the time
course), `tage_temporal_kruskal.csv`, `tage_temporal_spread.csv`.

**Decomposed:** `stress_response_pathways_RERUN.csv` (single source of truth for the gene
sets), `partial_tage/hallmark_pathway_mouse_ids.csv`, `partial_tage_ALL.csv` (1,700 rows),
`partial_tage_size_matched_null.csv`.

**Module arm, provenance only — do not report as pathway results:**
`partial_tage/paper_module_clocks.csv`, `partial_tage/paper_module_mouse_ids.csv`,
`module_tage_ALL.csv`, `module_size_matched_null.csv`.

**Deleted:** `tage_temporal_wilcoxon_vs_none.csv` and `temporal_analysis/05_tage_spread_figure.R`
(18-test family, superseded by the 36-test family and script `06`).
