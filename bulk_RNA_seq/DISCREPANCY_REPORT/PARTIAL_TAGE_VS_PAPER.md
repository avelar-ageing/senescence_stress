# Partial tAge: how our implementation differs from the tAge paper

**Status:** the implementation is arithmetically correct and its gene-level definition matches the
paper's. What differs is the *grouping unit* (overlapping curated Hallmark sets vs the paper's
non-overlapping WGCNA modules) and the *estimator* (post-hoc restriction of the fitted global clock vs
a separately trained clock per module). `PARTIAL_TAGE_METHOD.md` claims we use "the method described in
the paper" — that overstates it and should be reworded (§4). **Nothing already computed needs redoing.**

**Headline update (§6):** the paper's module clocks turned out to be fully published as Supplementary
Table 5 — every gene coefficient and intercept, for 14 multi-species and 23 rodent modules. So the
paper's own estimator can now be run directly on our data, with no retraining. Extraction is
implemented. This is the single most useful thing in this document.

Written 2026-08-16 after checking the paper, the package, the Zenodo model record and the TACO web
tool. Sources at the bottom. Anything I could not verify is marked **UNVERIFIED** rather than guessed.

---

## 1. What we implement

`exploratory/01–07` + `DISCREPANCY_REPORT/PARTIAL_TAGE_METHOD.md`.

The EN clock is linear, `pred = intercept + Σ_i coef_i · z_i`, where `z_i` is the imputed+scaled
feature value. We restrict that sum to the genes of an MSigDB Hallmark pathway:

```
partial_tAge(pathway) = Σ_{i ∈ pathway} coef_i · z_i        (× species adjustment)
```

computed from the full whole-transcriptome `tAge_preprocessing()` output — no re-normalisation of a
gene subset, no extra imputation. Reconstructing over *all* features reproduces `model.predict()` to
floating-point precision (max diff 0.0). Grouping unit: **50 sets = 49 Hallmark + a 191-gene
lysosomal list**, read from `stress_response_pathways_RERUN.csv` so the tAge and DEG-enrichment arms
cannot drift.

## 2. What the paper does

Tyshkovskiy, Gladyshev et al. 2026, *Nature*.

- **Modules are data-derived, not curated.** WGCNA on their rodent meta-dataset ("relative expression
  centred within each dataset, tissue and sex") gave 28 robust modules, 23 retained, 30–630 genes
  each, with "largely distinct pathway enrichments with **minimal overlap**" (Supplementary Table
  6a,b). 2,141 genes across all modules.
- **Module clocks were TRAINED per module**, not restricted post hoc:
  > "To develop a compact panel of biomarkers based on these transcriptomic network components, we
  > **trained** multi-tissue clocks of chronological age and expected mortality for each module"

  So the Extended Data Fig. 6h caption — "partial tAge differences predicted using only genes from the
  respective module" — refers to *the module clock's own prediction*, not a decomposition of the
  global model.
- **Gene-level contribution is defined the way we compute it.** Fig. 3g: "Values are logFC × clock
  coefficient"; Fig. 4c uses "slope × clock coefficient". Our `coef_i · z_i` on control-subtracted
  (`scaled_diff` / `yugene_diff`) input is the same quantity — coefficient × a differential expression
  value. **This is the part where we are aligned with the paper.**
- **Partial tAge is never applied to single-cell data in the paper.** Single-cell gets composite tAge
  only. Reported accuracy: ~1M reads per metacell → r = 0.90 (95% of maximum); ~250K reads → ~89% of
  maximum.

## 3. Side by side

| | paper | our pipeline |
|---|---|---|
| grouping unit | 23 WGCNA co-expression modules, data-derived | 49 Hallmark + lysosomal, curated |
| module overlap | "minimal" | substantial — Hallmark sets share genes |
| how the score is produced | separate EN clock **trained** per module | post-hoc restriction of the fitted global clock |
| gene-level contribution | logFC × coefficient | coef × z on control-subtracted input — **same idea** |
| single-cell / pseudobulk | composite tAge only | (our scRNA extension goes further — see §7) |

## 4. Is our method wrong?

**No.** Three separate claims, kept apart:

1. **The arithmetic is exact and verified.** The decomposition is an identity, not an approximation.
2. **The gene-level quantity matches the paper's own definition** (coefficient × differential
   expression). We are not inventing a contribution measure.
3. **The grouping unit and the training step differ.** We aggregate contributions of a *fitted global
   model* over *overlapping curated* pathways; the paper trains a *separate clock* per
   *non-overlapping data-derived* module. These are different estimators and can disagree.

So the method is a defensible adaptation, not an error. But the sentence in `PARTIAL_TAGE_METHOD.md`
that we use "the method described in the paper behind the tAge package" should become something like:

> Adapted from the module-level analysis in Tyshkovskiy/Gladyshev et al. 2026: we apply the paper's
> gene-contribution definition (coefficient × differential expression) but aggregate over MSigDB
> Hallmark pathways using the fitted global clock, rather than training a separate clock per WGCNA
> module as the paper does.

## 5. Three consequences of the overlap difference

These are properties of *our* variant that the paper's non-overlapping modules largely avoid.

*Re-verified 2026-08-16 against `rerun_outputs/`. Two of the three original figures did not
reproduce and are corrected below; the qualitative conclusions stand.*

1. **The partials do not partition the prediction.** Only **4,251 of the clock's 10,487 features**
   sit in any Hallmark set, and the sets themselves overlap, so the 50 partials neither sum to the
   prediction nor keep a fixed relation to it. The ratio (Σ partials / full tAge) is unstable across
   samples and models — median 0.96 on `scaled`, 1.76 on `yugene`, individual samples past 3.5x.
   Never present these as "pathway X contributes N of the total".
   *(The original 80.6-vs-33.7 / 2.4x example did not reproduce on either model.)*
2. **The pathways are not independent tests.** ~~median |Spearman rho| 0.545, 34% above 0.7~~ —
   **CORRECTED**, this did not reproduce under any interpretation tried. Per-sample across the
   meta-analysis cohort, median |Spearman rho| between pathway partial scores is **0.20** (scaled
   0.212, yugene 0.201), with **1.5-2.3%** of pairs above 0.7. Correlating pathway *effect-size
   profiles* across the 17 group comparisons instead — the highest reading found — gives median
   |rho| **0.35-0.37** with 14-16% above 0.7. Dependence is real but moderate; BH across 50 is
   mildly optimistic rather than badly so. Still report an axis, not a winning pathway.
3. **Effect size partly tracks pathway size.** Signed Cohen's d vs `n_genes_in_model` is ρ = **0.325**
   (confirmed; gene counts span 19-164). The size-matched null is **now implemented** —
   `exploratory/11_partial_tage_size_matched_null.py`, 1,000 random equal-size draws from the clock's
   own 10,487 features per pathway per group per model, output
   `rerun_outputs/partial_tage_size_matched_null.csv`. See §5b.

A fourth caveat is intrinsic to elastic net and applies to the paper too: EN distributes weight
arbitrarily among correlated genes, so which genes carry a pathway's score is partly an artefact of
fitting. The honest claim is "the clock's prediction, restricted to this pathway's genes, moves this
way" — not "this pathway drives ageing".

## 5b. What the size-matched null does to the results (added 2026-08-16)

**It is severe, and it should be.** Every clock feature carries age signal and the groups genuinely
differ, so a *random* equal-sized slice of the clock already separates them: median null mean |d| is
0.55-0.80 in the meta-analysis and 1.0-2.1 in the time course, with null 95th percentiles reaching
|d| ~5.9. The null therefore asks "is this pathway more informative than an arbitrary equal-sized
slice of the same clock?" — the right bar for a *pathway-specific* claim, but not for "does this
pathway move tAge".

Across 800 pathway x group x model tests only **72 (9%)** reach p_empirical < 0.05, and only **15
group-pathway combinations** clear it on *both* models with consistent sign:

| group | surviving on both models |
|---|---|
| meta OIS | APICAL JUNCTION, MYC TARGETS, XENOBIOTIC METABOLISM |
| meta RS | MITOTIC SPINDLE, PROTEIN SECRETION |
| meta SIPS | HEME METABOLISM, MITOTIC SPINDLE |
| temporal Fibroblast | G2M CHECKPOINT, HEME METABOLISM, MITOTIC SPINDLE |
| temporal Keratinocyte | ADIPOGENESIS, E2F TARGETS, G2M CHECKPOINT |
| temporal Melanocyte | INTERFERON ALPHA RESPONSE, P53 PATHWAY |

**Consequences for the write-up.** The "shared / universal core" framing (COMPLEMENT, UV RESPONSE UP,
P53 PATHWAY, E2F TARGETS, PI3K AKT MTOR SIGNALING, SPERMATOGENESIS, HEDGEHOG SIGNALING as
cell-type-independent) does **not** survive size matching in most groups. It must be withdrawn, or
explicitly downgraded to "consistent in direction across cell types, but not exceeding a
size-matched null".

The **melanocyte result strengthens**: INTERFERON ALPHA RESPONSE in melanocytes is d = -10.50 against
a null 95th percentile of 5.26 (p_emp = 0.005) — the clearest pathway-specific signal anywhere in the
dataset. Melanocyte P53 PATHWAY also survives (d = 5.79 vs null q95 5.41, p = 0.039), as do MITOTIC
SPINDLE in RS/SIPS/Fibroblast and MYC TARGETS in OIS (d = -2.52, p = 0.013).

Non-survival is **not** evidence of no effect: it means the effect is no larger than a random
equal-sized slice of an age clock applied to genuinely age-discordant groups. Report both facts
together. Per the project convention, empirical null p-values are reported raw, with the p-floor
(1/1001) and effect size, and are not BH-corrected — a pre-registered robustness check is not a
discovery family.

## 6. The paper's module clocks ARE reproducible — resolved 2026-08-16

Earlier drafts of this document said the module clocks were unobtainable. **That was wrong**, and the
correction matters more than anything else here.

- **Not on Zenodo, not in the package.** Zenodo record 18763485 holds 72 model files, all
  full-transcriptome (`[BR|EN]_[Chronoage|Mortality|NormalizedAge]_[Species]_[Tissue]_[scaleddiff|yugenediff].pkl`).
  The installed `tAge` package ships only `extdata/metadata/Module_to_function_map.csv` — 28 modules
  with names, functional labels, gene counts and a +/- direction, but **no gene membership**.
- **They are published as Supplementary Table 5.** `MOESM7_ESM.xlsx`, sheets "(B) Module rodent clocks"
  and "(C) Module multi-species clocks", give **every gene's coefficient plus the intercept for each
  module clock**. That is the complete linear model, so the paper's estimator can be evaluated directly:

  ```
  module_tAge = intercept + sum_i coef_i * z_i      over that module's genes
  ```

  on the same `tAge_preprocessing()` output we already export. **No retraining required.**
- **Extraction is implemented**: `exploratory/11_extract_paper_module_clocks.R` →
  `rerun_outputs/partial_tage/paper_module_clocks.csv` (tidy: panel, outcome, module, annotation,
  intercept, entrez_id, gene_symbol, coefficient) and `paper_module_sizes.csv`. Download the
  supplement once with
  `curl -o suppl.zip https://www.ebi.ac.uk/europepmc/webservices/rest/PMC13233323/supplementaryFiles`.

**What was recovered**

| panel | outcome | modules | genes/module | total |
|---|---|---|---|---|
| multispecies | Chronological age | 14 + composite | 20–238 | 1,284 composite |
| multispecies | Mortality | 14 + composite | 20–238 | |
| rodent | Chronological Age | 23 + composite | 12–2,141 | 2,141 composite |
| rodent | Mortality | 23 + composite | 12–2,141 | |

Recovered sizes cross-check against the published module dictionary (MOESM8, sheets B/D): 238/238,
186/186, 169/169, 117/117, 92/92, 79/79, 45/45, 21/21, 20/20 exact, with small shortfalls (75 vs 78,
73 vs 79, 65 vs 68, 23 vs 24) where elastic net zeroed a few coefficients. Expected, worth a footnote.

**Use the multispecies panel for human data** — it matches the `Multispecies_Multitissue` models both
arms already use. Its 14 modules: innate immunity/inflammation (238), mitochondrial translation/OxPhos
(186), chromatin modification (169), cell cycle/DNA replication (117), muscle contraction/cytoskeleton/
glycolysis (92), mRNA splicing (79), ECM organization/EMT (75), adaptive immunity/T cell signaling (73),
OxPhos/heme metabolism (65), interferon signaling (45), VEGF signaling (45), protein processing in
ER/UPR (23), protein folding/translation (21), fatty acid metabolism/peroxisome (20).

## 6b. VALIDATION FAILURE — do not run the module clocks yet (added 2026-08-16)

§6 says the module clocks "can be evaluated directly ... **No retraining required**". **That is not
yet safe.** The naive evaluation fails a basic sanity check, so no module-level results should be
produced or reported until this is resolved.

**What was tested.** The multispecies `All module genes` composite chronological-age clock
(1,284 genes, intercept 0.032118972, from `paper_module_clocks.csv`) evaluated as
`intercept + Σ coef_i · z_i` over its genes, where `z` is our exported `meta_scaled_diff` matrix
after the global clock's own imputer + scaler, × the 122.5 human species adjustment — i.e. exactly
the recipe §6 proposes, on the matrices we already export.

**Three results, in increasing severity:**

1. **Gene IDs are fine.** Both panels use mouse symbols (Title case) and mouse Entrez IDs; all
   1,284 composite genes are present among the clock's 10,487 features. Not a mapping bug.
2. **It barely tracks the global clock.** Pearson r = **0.263** against the
   `Multispecies_Multitissue` prediction on the same 230 samples. Two chronological-age clocks on
   identical input should agree far better. (Standardisation convention is not the culprit:
   using raw imputed values instead of scaler `z` gives r = 0.263 to four decimals, because the
   `scaled_diff` input is already approximately standardised, so the StandardScaler is near-identity.)
3. **It reverses the study's most robust finding.** Against Proliferating controls:

   | condition | module composite d (p) | global clock d (p) |
   |---|---|---|
   | CICQ | +0.66 (0.033) | +1.21 (3.2×10⁻⁵) |
   | SSCQ | +0.27 (0.46) | +0.37 (0.21) |
   | RS | +0.57 (0.075) | +1.02 (2.9×10⁻³) |
   | **SIPS** | **−0.20 (0.33)** | **+1.13 (3.6×10⁻⁷)** |
   | **OIS** | **−0.10 (0.081)** | **+0.96 (1.5×10⁻⁶)** |

   SIPS and OIS are the two conditions with the largest, most model-consistent tAge elevation in
   the entire study (2.1.5.1). The module composite makes them *negative and non-significant*.

**Also note:** the composite's intercept is bit-identical to our global clock's
(0.03211897213960977) — but so is the `yugenediff` model's, so that is a shared training-target
constant, **not** evidence the models match. Coefficients are uncorrelated between them
(r = −0.009 over the 1,284 shared genes; only 234 are among our clock's 1,839 non-zero features).

**Most likely cause: an input-convention mismatch we have not identified.** The paper trained module
clocks on "relative expression centred within each dataset, tissue and sex", whereas our matrices are
control-subtracted against a designated control group; these are not obviously the same
transformation. Other candidates: a different species adjustment for the module panel, or a
normalisation step applied before their scaler that we do not replicate.

**Required before any module analysis is trusted:** a positive control. Evaluate the module clocks on
data with known ages where the published prediction is recoverable — e.g. the package's own
`Exprs_example.csv`/`Metadata_example.csv`, or a public ageing dataset the paper reports — and confirm
the composite reproduces sensible chronological-age estimates. If it cannot be made to, the module
route is closed and the Hallmark decomposition stands alone (which §4 already establishes as
defensible). **Until then, §8 action 3 is blocked, not merely pending.**

**One constraint (if it is unblocked):** the published module clocks are the "scaling" variant only — i.e. `scaleddiff`.
There appear to be **no yugene module clocks**. Both arms currently report `scaled` and `yugene`
throughout, so any module analysis will be single-normalisation and that asymmetry must be stated.

- **TACO** (https://app.gladyshevlab.org/TACO/) exposes "Gene contribution", "Pathway Contribution" and
  "Module clocks", with enrichment across "KEGG, Hallmark, and Module" categories — so scoring Hallmark
  pathways against clock contributions is a supported use of this clock family, which independently
  supports our grouping choice. **UNVERIFIED:** whether TACO's "Pathway Contribution" sums gene
  contributions (our method) or runs enrichment on top-contributing genes. It is a Shiny app, so this
  cannot be settled by fetching; it needs an interactive session. Low priority now that the module
  clocks are in hand.

## 7. The scRNA-seq extension (context for whoever picks this up)

The same decomposition has been applied to the GSE226225 pseudobulks
(`scrna_seq/GSE226225/final_analysis/partial_tage_0{1,2,3,4}*`), at three resolutions: cluster,
per-replicate substate, and replicate-pooled condition patch; both `scaled` and `yugene` models.
`partial_tage_02_decompose.py` is this repo's `04_partial_tage_decompose.py` copied **unchanged**, and
the pathway→mouse-ID mapping is read from this repo's `hallmark_pathway_mouse_ids.csv`, so the two
arms cannot drift.

Status there: all six decompositions exact (0.0); reconstructed full tAge correlates r = 0.9992 with
the previously published cluster-level tAge. Clock-gene coverage per pathway halves relative to bulk
(median ~57–60 vs 125), and two pathways become unusable — `PANCREAS BETA CELLS` (8 genes) and
`WNT BETA CATENIN SIGNALING` (12 genes, zero-variance partials).

**Note this goes beyond the paper**, which applies partial tAge only to bulk. It also runs at
`coverage_threshold = 2e5` per pseudobulk, against the package README's single-cell example of `1e7`
and the paper's ~250K-reads-per-metacell figure (~89% of maximum accuracy). That threshold was chosen
from a plateau sweep on our own data and is defensible, but it needs stating in Methods with the
sweep as justification.

## 8. Recommended actions, ranked

1. ~~Reword the "method described in the paper" sentence in `PARTIAL_TAGE_METHOD.md`~~ — **DONE
   2026-08-16.** The file now says "adapted from", states what we borrow (the gene-contribution
   definition) and what we change (global clock + curated overlapping sets vs per-module trained
   clocks + near-disjoint data-derived modules), and points here.
2. ~~Add the size-matched null~~ — **DONE 2026-08-16**,
   `exploratory/12_partial_tage_size_matched_null.py` → `partial_tage_size_matched_null.csv`.
   Consequences in §5b: only 9% of pathway effects clear it, the "universal core" framing does not
   survive, the melanocyte interferon result does and strengthens.
3. **BLOCKED — see §6b.** Running the paper's module clocks was attempted and **fails validation**:
   the multispecies composite correlates r = 0.26 with the global clock on the same samples and
   reverses the sign of the SIPS and OIS tAge elevation. Do not report module-level results until a
   positive control on known-age data succeeds. The coefficients in `paper_module_clocks.csv` are
   correctly extracted; the problem is the input convention for applying them.
4. Note the normalisation asymmetry: module clocks exist for `scaleddiff` only, while both arms report
   `scaled` and `yugene` everywhere else.
5. Cross-check one dataset against TACO's "Pathway Contribution" output to learn what it computes
   (low priority now).
6. Record consequences (1)–(3) wherever the Hallmark partial-tAge results are reported, so no reader
   takes the partials as a variance partition. These largely do not apply to the module clocks, whose
   modules barely overlap and which are genuine standalone predictions rather than decomposed terms.

## Sources

- Tyshkovskiy, Gladyshev et al. (2026) *Universal transcriptomic hallmarks of mammalian ageing and
  mortality*, Nature — https://pmc.ncbi.nlm.nih.gov/articles/PMC13233323/
- https://github.com/Gladyshev-Lab/tAge
- https://zenodo.org/records/18763485 (72 model files, all full-transcriptome)
- https://app.gladyshevlab.org/TACO/
