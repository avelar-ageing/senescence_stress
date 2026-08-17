# The paper's WGCNA modules as a grouping unit: what we found

**Bottom line: do not use them.** They were run end to end, correctly, and they do not
outperform MSigDB Hallmark on any axis we can measure. They are markedly *less* reproducible
across the two EN normalisations, and they cannot be used to test our headline interferon result
because the gene list published for the paper's interferon module contains none of the interferon
genes — a defect in the published supplement, not in the modules themselves
(`RECONCILIATION_MODULE_CLOCKS.md` §3). **Read that document alongside this one:** it retracts the
route-A diagnosis in §2 below and corrects the interpretation in §5.4.

Written 2026-08-16. Self-contained handoff. Every number below was computed on
`rerun_outputs/` and is reproducible from the scripts named. Anything not verified is marked.

---

## 1. Why we tried this

`PARTIAL_TAGE_VS_PAPER.md` §5 lists three consequences of decomposing a global clock over
*overlapping curated* Hallmark sets rather than the paper's *near-disjoint data-derived* WGCNA
modules: the partials don't partition the prediction, the pathways aren't independent, and
effect size partly tracks set size. Adopting the paper's modules was supposed to fix the first
two. It was ranked the highest-value outstanding action.

## 2. Two routes, and why route A is dead

**Route A — run the paper's module CLOCKS. Blocked, but NOT for the reason first given.**

> **This paragraph's original diagnosis was wrong and is retracted (2026-08-17).** It claimed the
> fitted imputer/scaler is unpublished and the coefficients therefore unusable. In fact MOESM7
> sheet (A) col 27's published coefficients equal the fitted pkl's to 9.9e-17 and reconstruct
> `model.predict()` to 2.6e-07, so published coefficients DO live in the pkl's standardised space.
> The "positive control failure" cited below was also a bad test (24 near-identical-age Klotho
> mice). Route A is blocked instead by a gene-to-module assignment defect in MOESM7 sheets (B)/(C)
> that contradicts the paper's own enrichment tables. Full evidence:
> `RECONCILIATION_MODULE_CLOCKS.md` §1 and §3.

Original text: Supplementary Table 5 (MOESM7) publishes every module clock's per-gene coefficients
and intercept, which looks like the full linear model. It isn't: the fitted `SimpleImputer` and
`StandardScaler` are **not** published [RETRACTED — they are recoverable]. Evaluating them in our
global model's feature space fails a positive control on the tAge package's own example mouse data
(r = 0.48 rodent panel, 0.22 multispecies) [RETRACTED — invalid control], and on our data it
reverses the sign of the SIPS and OIS tAge elevations [this concerns the all-module composite;
individual module clocks behave coherently].

**Route B — use the module MEMBERSHIPS as a grouping for our own verified decomposition.**
This keeps the estimator we have validated (exact linear decomposition of the fitted global
clock, reconstruction max diff 0.0) and swaps only the gene grouping. This route works
mechanically and is what everything below reports.

## 3. What was run

| step | script | output |
|---|---|---|
| extract module clocks from supplement | `exploratory/11_extract_paper_module_clocks.R` | `partial_tage/paper_module_clocks.csv` |
| build module → clock-feature map | `exploratory/13_build_paper_module_mapping.R` | `partial_tage/paper_module_mouse_ids.csv` |
| decompose (reused **unchanged**) | `exploratory/04_partial_tage_decompose.py` | `partial_tage/*_modulescores_*.csv` |
| consolidate stats | 05-equivalent | `rerun_outputs/module_tage_ALL.csv` |
| size-matched null | `exploratory/12_partial_tage_size_matched_null.py` | `rerun_outputs/module_size_matched_null.csv` |

Panel used: **multispecies, Chronological age** (matches the `Multispecies_Multitissue` models
both arms already use). 26 decompositions, all exact (max diff 0.0). 17 group comparisons x 14
modules x 2 models = 476 rows.

The rodent panel was **not** used: 23 modules, 1,922 genes, but only 65% are in the clock's
feature space. The multispecies panel is complete (1,248/1,248).

## 4. On paper this should have worked

- 14 modules, 1,248 genes
- **Zero** genes shared between modules — perfectly disjoint, exactly the property Hallmark lacks
- **1,248/1,248** present among the clock's 10,487 features — no coverage loss
- Module sizes reproduce the published dictionary (MOESM8 sheet D) exactly for 9 of 14, with
  small shortfalls (75 vs 78, 73 vs 79, 65 vs 68, 23 vs 24) where elastic net zeroed coefficients

## 5. It did not work. Four measured results

### 5.1 Disjointness bought almost nothing
Inter-set correlation of partial scores across samples (meta-analysis, yugene):

| | median \|Spearman rho\| | frac pairs > 0.7 |
|---|---|---|
| paper modules (disjoint) | **0.20** | **0%** |
| Hallmark (overlapping) | 0.201 | 1.5% |

Only the extreme tail improved. The dependence between sets is driven by shared sample-level
structure, not by shared genes — so §5 consequence (2) is not actually caused by overlap, and
removing overlap does not remove it.

### 5.2 Modules are LESS reproducible across normalisations
Agreement between the `scaled` and `yugene` EN models across all 17 group comparisons:

| | sign agreement | both significant AND same sign |
|---|---|---|
| paper modules | **0.55** | **0.235** |
| Hallmark | 0.775 | 0.382 |

Several modules flip sign with *both* models significant — e.g. Fibroblast "Fatty acid
metabolism/Peroxisome" d = **-6.18** (scaled) vs **+10.20** (yugene). Likely mechanism: modules
span only 12% of the clock's features, so model-specific coefficient differences are averaged
down far less than in the broader Hallmark sets.

### 5.3 They do not beat a size-matched null any better than Hallmark
1,000 equal-size draws from the clock's own 10,487 features, per set per group per model:

| | p_emp < 0.05 | survivors on BOTH models | median \|d\| | median \|d\| / null q95 |
|---|---|---|---|---|
| paper modules | 17/224 (7.6%) | **1 / 112 (0.9%)** | 0.86 | 0.40 |
| Hallmark | 72/800 (9.0%) | **15 / 400 (3.8%)** | 0.87 | 0.39 |

Raw effect sizes are indistinguishable (0.86 vs 0.87; ratio to null 0.40 vs 0.39). The
difference is entirely in reproducibility: modules yield **4x fewer** both-model survivors.
The single module survivor is Keratinocyte "Cell cycle / DNA replication (green)".

*Caveat on this null, which applies equally to both arms:* the draw pool is the clock's own
features, i.e. genes already selected for age-association across mammals. It is therefore not a
neutral background — it asks "does this set beat an arbitrary slice of the age clock itself",
a deliberately hard question. Structurally there is no alternative pool: a non-clock gene has
coefficient exactly 0 and contributes nothing. Non-survival is not evidence of no effect.

### 5.4 The module gene lists do not match their annotations

> **INTERPRETATION CORRECTED 2026-08-17.** The mismatch below is real, but the explanation given
> at the end of this section (that WGCNA annotations are loose "top-enrichment labels" and
> co-expression need not group canonical pathway members) is **wrong**. The paper's own enrichment
> table gives darkmagenta an odds ratio of 594 for Hallmark Interferon Alpha Response at
> padj = 8e-58 — it genuinely IS an interferon module. The gene lists published in MOESM7 sheets
> (B)/(C) are what is wrong, not the labels. See `RECONCILIATION_MODULE_CLOCKS.md` §3 for the
> evidence and two failed repair attempts.
This is the most important finding for anyone hoping to interpret module results.

`HALLMARK INTERFERON ALPHA RESPONSE` (85 clock genes) and the paper's
`darkmagenta / Interferon signaling` module (45 clock genes) share **zero genes**.

Canonical interferon-stimulated genes that *are* in the sheet's 1,284-gene universe sit in
entirely different modules:

| gene | assigned module | its annotation |
|---|---|---|
| Isg15 | orange | Chromatin modification |
| Irf7 | orange | Chromatin modification |
| Rsad2 | darkred | mRNA splicing |
| Usp18 | pink | Mitochondrial translation / OxPhos |
| Ifih1 | pink | Mitochondrial translation / OxPhos |

No module is enriched for its own annotation's canonical markers above background (e.g.
"mRNA splicing", 79 genes, contains 0 of 9 canonical splicing factors; "Adaptive immunity /
T cell signaling", 73 genes, contains 0 of 10 canonical T-cell markers).

**This was checked hard and is NOT an extraction bug:**
- Both ID sets are mouse Entrez and resolve in the package gene table (1,248/1,248 modules;
  4,234/4,251 Hallmark)
- Extracted entrez↔symbol pairing is **98.5%** identical to the package's `Gene_table_mouse.csv`;
  the 1.5% are symbol-vintage synonyms (`Atp5mc3`/`Atp5g3`, `Bbln`/`1110008P14Rik`)
- The MOESM7 sheet layout was read directly and matches `11_extract_paper_module_clocks.R`'s
  assumption: columns 1–2 are a shared gene index, row 2 = module colour, row 3 = annotation,
  row 5 = intercept, rows 6+ = coefficients
- Module colours, annotations and sizes agree with the official dictionary
  (MOESM8 sheet D: darkmagenta = Interferon signaling, orange = Chromatin modification)

**Interpretation.** These are WGCNA *co-expression* modules derived from rodent ageing data. The
functional annotation is a top-enrichment label for the module as a whole, not a membership
criterion — and rodent co-expression structure need not group canonical pathway members
together, still less transfer to irradiated human skin cells. The labels are therefore not
usable as pathway identities in our setting.

## 6. Consequences

1. **Do not report module-level results as pathway-level findings.** The outputs
   (`module_tage_ALL.csv`, `*_modulescores_*.csv`) are kept for provenance only.
2. **The melanocyte interferon result can be neither confirmed nor refuted by the modules**,
   because the module labelled "Interferon signaling" does not contain the interferon genes
   that drive it. It is not independent corroboration and must not be cited as such.
3. **Hallmark stays the interpretable arm** — as `PARTIAL_TAGE_VS_PAPER.md` §4 already argued
   on other grounds. This work strengthens rather than weakens that position: the paper's own
   grouping was tried fairly and performed worse.
4. **§5 consequence (2) needs rewording.** Set-to-set dependence is not caused by gene overlap;
   perfectly disjoint sets show the same median |rho| (0.20). It is sample-level structure.
5. **Route A stays open only if someone obtains the module clocks' scaler statistics** from the
   authors. Without them the published coefficients cannot be applied.

## 7. Open / not done

- **UNVERIFIED:** whether the module annotations are similarly non-descriptive in the *rodent*
  panel (23 modules). Only the multispecies panel was examined in detail.
- A **label-permutation null** (permute group labels, recompute d for the real set) has been
  proposed but not run. It asks the question the results sections actually make — "does this set
  separate the groups more than chance" — rather than the harder size-matched question, and
  would be a fairer complement for both arms.
- No figures were produced for the module arm, deliberately.
