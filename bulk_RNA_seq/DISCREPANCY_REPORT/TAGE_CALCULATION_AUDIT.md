# How tAge is computed, and whether the reported quantities are sound

Written 2026-08-18. Reproducible via `exploratory/17_tage_calculation_audit.py`.
Read alongside `PARTIAL_TAGE_VALIDITY.md` (which audits the pathway decomposition)
and `PARTIAL_TAGE_VS_PAPER.md` (which audits the method against the source paper).

**Verdict: within a preprocessing run the reported differences are sound, and are
provably insensitive to both imputation and control-reference noise. Across runs,
magnitudes are not comparable. Absolute tAge is not interpretable.**

---

## 1. The chain, and where it depends on the batch

From `tAge_preprocessing()`, read directly from the installed package:

| step | what it does | batch-dependent? |
|---|---|---|
| `filter_genes` | keep genes with >=10 counts in >=20% of samples | **YES** — threshold applied within the run |
| `map_genes` | human -> mouse orthologs, fixed lookup table | no |
| `RLE_normalization` | edgeR size factors | **YES** — reference built from the run's samples |
| `log_transform` | per value | no |
| `scale_eset` | `scale(exprs)` on a genes x samples matrix, so R standardises **columns** = each SAMPLE across genes | **no** — sample-intrinsic |
| `YuGene` | `apply(counts, 2, ...)`, also per sample | **no** — sample-intrinsic |
| `.align_to_gene_list` | align to the 18,696-gene reference; absent genes become NA | no |
| `control_subtraction` | subtract the per-gene **median of the control samples** | **YES** — by design; defines the zero |
| model `SimpleImputer(median)` | fill NA with training medians | no (fixed) |
| model `StandardScaler` | training-set statistics | no (fixed) |
| `ElasticNet` -> x122.5 | linear prediction, then species adjustment | no (fixed) |

Note the two normalisation steps that matter most, `scale_eset` and `YuGene`, are
**per sample** and therefore do not depend on batch composition. An earlier draft of
this audit asserted the opposite (that `scale_eset` standardised genes across
samples); that was wrong — R's `scale()` operates on columns, and the matrix is
genes x samples.

## 2. Units — reconciling two statements that sounded contradictory

Both of these are true and they describe different things:

**(a) The scale is years-dimensioned.** The clocks were trained on chronological age
"divided by species maximum recorded lifespan", so the model predicts a *fraction of
maximum lifespan*. `tage_predict.py` then multiplies by
`PREDICTIONS_SPECIES_ADJ = {"human": 122.5, "mouse": 48, "rat": 50.4, "monkey": 39}`.
For human, 122.5 is the maximum recorded human lifespan in years. So **122.5 units =
one maximum human lifespan**, and a value of +54.3 is 44% of that. Describing the
scale in these terms is legitimate.

**(b) The quantity is not an age, and its magnitude is not calibrated.** Two separate
reasons, which an earlier draft wrongly merged into a blanket "not years":

- **Reference point.** Because the input is control-subtracted, the model scores each
  sample's *deviation from the control group*, not its age. Control samples therefore
  score at the intercept: the Proliferating group's own median tAge is 2.61 (scaled) /
  0.96 (YuGene) against an intercept x 122.5 of 3.93. Zero is the control group, not
  birth. "+54.3" is not "54 years old".
- **Calibration.** Whether +54.3 corresponds to 54 years of real ageing is unverified
  here: the clock was trained on tissues and applied to cultured cells; ~30% of its
  coefficient weight sits on genes this experiment never measured, so the response is
  attenuated (section 3); no donor ages exist in either dataset to check against; and
  the two normalisations differ roughly two-fold on the same samples (SIPS raw 56.9 vs
  30.9).

**Practical convention:** report values as tAge units, define the scale once
(122.5 units = one maximum human lifespan), and state that they are differences from
the named control group. Do not quote them as years of ageing. Note the scRNA arm
currently reports years ("clusters differ by 7-12 years"); the two arms need one
agreed convention.

## 3. Coverage: how much of the clock is actually operating

The model has 10,487 features, 1,839 with non-zero coefficients. Because
`filter_genes` runs within each batch, coverage differs by run:

| run | n | features observed | % | non-zero observed | **% of non-zero weight measured** |
|---|---|---|---|---|---|
| meta-analysis | 230 | 9,105 | 86.8% | 1,410 | **80.5%** |
| Fibroblast pooled | 24 | 8,529 | 81.3% | 1,217 | 70.2% |
| Keratinocyte pooled | 24 | 8,546 | 81.5% | 1,234 | 70.0% |
| Melanocyte pooled | 24 | 8,711 | 83.1% | 1,266 | 72.8% |
| Fibroblast 4 days | 12 | 8,497 | 81.0% | 1,204 | **69.8%** |
| Melanocyte 20 days | 12 | 8,696 | 82.9% | 1,262 | 73.0% |

So the clock runs on 70-80% of its weight, and **the fraction differs by ~10
percentage points between the meta-analysis and the temporal runs**. This attenuates
the response and is an independent reason cross-run magnitudes are not comparable. It
does **not** bias anything within a run — see next section.

## 4. Two concerns tested and dismissed

### 4.1 Imputation is difference-neutral

Unmeasured features are filled with the training median, so they take the **same value
in every sample of a run** and cancel from any between-group difference. Verified:

| run | imputed features | max \|contribution difference\| among them | observed-only total | whole-transcriptome |
|---|---|---|---|---|
| Fibroblast 20 days | 1,898 | 0.000e+00 | +21.0091 | +21.0091 |
| meta (SIPS) | 1,382 | 1.5e-16 | +31.8260 | +31.8260 |

An earlier draft of this audit described the prediction as "part real signal, part
fixed offset", implying imputation corrupts the estimates. It does not. It costs
**sensitivity, not accuracy**.

### 4.2 Control-reference noise is difference-neutral

`control_subtraction` uses a per-gene median estimated from the control samples — only
**6** in the temporal runs, where a median's sampling error is roughly 0.5 sigma. The
worry was that this noisy reference is *shared* by every sample in the run and could
therefore move a whole gene set's contribution coherently, in a way the Wilcoxon test
could not detect.

**It cannot.** Subtracting a per-gene constant `c` shifts the test and control means
equally, so `mean_test - mean_ctrl` is invariant; and since every per-sample value
shifts by the same amount, ranks are unchanged and the Wilcoxon test is invariant too.
Bootstrapping the control reference (B = 200) confirms it:

| run | n control | whole-transcriptome observed | bootstrap SD | max set-contribution SD |
|---|---|---|---|---|
| Fibroblast 20 days | 6 | +21.0091 | 8.5e-15 | 2.0e-14 |
| Melanocyte 20 days | 6 | +11.3432 | 6.4e-15 | 1.8e-14 |

Zero to floating-point precision. The small control group affects the *precision of the
Wilcoxon test* through the usual sample-size route, but introduces no reference-estimation
error into the reported effects.

## 5. What remains: cross-run comparability

This is the one substantive problem, and it has three independent sources, all
quantified above:

1. **Gene retention** differs by run (12,236 genes in meta vs 11,298-11,614 temporal).
2. **Measured coefficient weight** differs (80.5% meta vs ~70% temporal).
3. **Input scale** differs (median \|z\| 0.165 meta vs 0.073-0.118 temporal; SD 0.427 vs
   0.228-0.328), which passes linearly into the prediction.

Plus a fourth, conceptual: each run is referenced to a **different control group**
(pooled Proliferating fibroblasts for the meta-analysis; that cell type's own untreated
samples for each temporal run). Because control subtraction forces each control group
to ~0, the real between-cell-type baseline difference is removed and cannot be recovered
from these outputs.

**Permitted**
- All within-run comparisons. The six meta-analysis conditions share one run, so that
  entire analysis, including every condition-vs-condition contrast, is sound.
- Rank and correlation analyses across runs, which are scale-invariant.
- Within the temporal arm, magnitude comparisons that substantially exceed the scale
  drift (roughly 1.5-fold; drift is 8-21%).

**Forbidden**
- Magnitude comparison between the meta-analysis and the time course.
- Absolute tAge values, or any statement in years of ageing.
- Comparing baselines between cell types.

**If cross-arm comparison is needed**, the fix is a single joint preprocessing run with
one common control group. That was deliberately avoided in
`temporal_analysis/04_tage_by_celltype.R` because subtracting a cross-cell-type average
baseline is inappropriate for measuring within-cell-type irradiation effects — so a
joint run would have to be an *additional* analysis, not a replacement.

## 6. Consequences for claims already drafted

- The modularity argument rests on three correlational results (within-cell-type
  profile similarity 0.62 vs between-cell-type 0.12; three PCA axes for 80% of
  variance), all scale-invariant and therefore unaffected.
- The fourth modularity line — that groups with near-identical *aggregate* tAge have
  unrelated compositions (rho = -0.07), illustrated by keratinocyte 4 days (+19.0) and
  melanocyte 20 days (+19.5) — **is compromised**, because it compares aggregates
  across runs. Drop it.
- Every effect reported within the meta-analysis stands.
- Statements of the form "tAge rose by X in cell type A and Y in cell type B" are safe
  only where X/Y substantially exceeds 1.5.
