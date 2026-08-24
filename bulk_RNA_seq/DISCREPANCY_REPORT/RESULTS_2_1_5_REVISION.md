# 2.1.5.x: section-by-section revision after the baseline-confound work

Written 2026-08-18. Supersedes the drafts held in conversation only; the prose had
never been committed to a file, so the pre-revision text of each affected claim is
quoted here before its replacement.

Driving results: `meta_analysis/11` (controls are not exchangeable),
`meta_analysis/13` (within-study condition effects),
`exploratory/18` (within-study pathway contributions), `meta_analysis/10`
(GEO-verified cell line and hTERT status).

---

## Summary of what changes

| section | verdict | why |
|---|---|---|
| 2.1.5.1 universal meta tAge | **rewrite** | pooled contrasts confounded by control baseline; CICQ and SSCQ do not survive on scaled_diff |
| 2.1.5.2 pathway meta tAge | **rewrite** | same confound, inherited; every contribution recomputed within study |
| 2.1.5.3 universal temporal tAge | **no change on this account** | single study (ERP021140, all 72 samples); every contrast is already within study and within cell type |
| 2.1.5.4 pathway temporal tAge | **two sentences only** | the temporal numbers stand; the cross-arm comparability claims are wrong and must go |
| Methods | **two additions, one deletion** | add the within-study design; delete the cross-arm comparability sentence |

The single largest change is that **`scaled_diff` fails and `yugene_diff` holds.**
That is now supported three independent ways: scaled_diff responds to
immortalisation status among controls (script 09), its control baselines are twice
as dispersed as yugene's (script 11: study medians span 87.8 vs 72.0 units, IQR
47.4 vs 15.9), and two of its five condition effects vanish under within-study
control while none of yugene's do (script 13).

---

## 2.1.5.1 Universal Transcriptomic Age Differences

### What was wrong

The section reported pooled condition-versus-Proliferating contrasts. Pooling is
not valid here: among the 91 untreated Proliferating controls alone, study medians
span 87.8 tAge units on scaled_diff and cell-line medians 70.7 (IMR90 -15.3 vs
HCA2/BJ primary +45.4). The condition effects being interpreted were 8-61 units,
i.e. the same size as the variation between untreated control groups. Condition is
also confounded with cell line - CICQ contains no IMR90 at all, OIS is 28 of 48
IMR90-derived, while Proliferating is 33 of 91 IMR90 - so a pooled CICQ contrast
compares a high-baseline strain against an IMR90-dominated control group.

Two specific claims do not survive:

- > "The shift was clearest in SIPS (+54.3 / +29.9 units ...) ... and, intriguingly, CICQ (+47.6 / +22.1 ...)"

  CICQ's scaled_diff effect falls from +47.6 to **+5.9 (p = 0.14, not
  significant)** within study. The "intriguingly, CICQ" framing was an artefact.
  Earlier in this work I attributed it to immortalisation; that was wrong in
  attribution though right in direction - it is baseline strain composition.

- > "SIPS > CICQ > RS > OIS > SSCQ still holds, and that's the claim the data supports."

  This ordering is not supportable. Within study on yugene the order is RS > SIPS >
  OIS > CICQ > SSCQ, and RS's estimate rests on 6 strata totalling 11 arrested and
  9 control samples. The conditions cannot be confidently ranked; only SSCQ being
  smallest is stable.

### Replacement text

> Applying the tAge elastic-net clock, every arrest condition showed a significantly
> higher tAge score than proliferating controls (Table x). Because tAge is computed
> from a control-subtracted expression matrix, scores are differences from the
> proliferating group on the clock's lifespan-relative scale, where 122.5 units
> corresponds to one maximum human lifespan; they are not estimates of donor age,
> and the zero point is the control group rather than birth.
>
> Untreated proliferating fibroblasts were found not to be interchangeable between
> studies: among the 91 proliferating samples alone, median tAge differed by up to
> 87.8 units between studies and by 70.7 units between cell strains (IMR90 -15.3
> versus primary HCA2/BJ +45.4 on the scaled-difference model), and cell strain,
> tissue of origin, study and hTERT-immortalisation status all predicted control
> tAge (Kruskal-Wallis, FDR < 0.05 on both models). Since these baseline
> differences are as large as the effects of interest, and since arrest condition is
> confounded with cell strain across studies, each arrested sample was compared only
> with proliferating controls from its own study. Each study here uses a single cell
> line, so this holds cell strain, immortalisation status, tissue and batch constant
> by construction. Per-study differences were combined with precision weighting and
> tested by permuting the condition label within study (20,000 permutations); 32 of
> 34 studies contain internal controls, covering 223 of 230 samples.
>
> On the YuGene model all five conditions remained significant, with effects close
> to their pooled values and consistent in direction across every contributing
> study: replicative senescence +32.5 units (6/6 studies, FDR = 1.9x10-4),
> stress-induced senescence +27.8 (4/4, FDR = 8.3x10-5), oncogene-induced
> senescence +26.9 (15/15, FDR = 8.3x10-5), contact-inhibited quiescence +19.1
> (6/6, FDR = 8.3x10-5) and serum-starved quiescence +7.8 (6/6, FDR = 2.8x10-3).
> On the scaled-difference model only the three senescence subtypes survived - RS
> +42.3 (FDR = 5.0x10-4), OIS +33.3 (FDR = 1.3x10-4) and SIPS +23.3 (FDR =
> 1.3x10-4) - while contact-inhibited quiescence fell from +47.6 pooled to +5.9
> (p = 0.14) and serum-starved quiescence from +21.7 to -3.4 (p = 0.11). The
> scaled-difference model was also the more sensitive of the two to control
> composition, differing by +31.1 units between immortalised and primary
> proliferating samples where YuGene differed by +6.6, and showing twice YuGene's
> between-study dispersion among controls. Subsequent magnitude statements
> therefore use the YuGene model.
>
> Serum-starved quiescence gave the smallest elevation on both models. The
> remaining four conditions could not be reliably ordered: the within-study
> ranking differs from the pooled one, and the replicative-senescence estimate
> rests on six studies contributing 11 arrested and 9 control samples in total.
> Pairwise contrasts between conditions were correspondingly unable to separate the
> three senescence subtypes from contact inhibition (all FDR > 0.13).

### Also required

- Delete any "primary-only" sensitivity table (the `09` exclusion analysis). It is
  superseded: the within-study design removes the immortalisation confound without
  dropping samples, and 09's own framing mis-attributed the CICQ effect.
- Replicative senescence cannot be compared with an hTERT control, because hTERT
  maintains telomeres and so structurally precludes replicative senescence. Keep
  this as a stated limitation, not a gap to be filled.

### Immortalisation: one Results sentence, one Discussion sentence

Decision (2026-08-21): the immortalisation work is NOT written up as a result. The
effect is not estimable - hTERT status is perfectly nested within study (0 of 34
studies contain both), and the immortalised OIS studies additionally differ in
induction agent (BRAFV600E rather than HRASG12V in the largest) and exposure
duration (10-28 days versus 6). Scripts 09, 12 and 16 are retained as provenance.
It gets one sentence in each of Results and Discussion, no table and no figure.

NOTE ON WORDING (corrected 2026-08-21): an earlier draft of this sentence said OIS
was "the only arrest condition containing both". That is WRONG. Immortalised
samples occur in four conditions - Proliferating 21/91, CICQ 3/19, SSCQ 5/22 and
OIS 17/48 - and are absent only from RS (0/11) and SIPS (0/39). OIS is the only
SENESCENCE SUBTYPE with both, which is the defensible version of the claim.

Two further constraints:
  - "Older" is specific to the scaled-difference model. On YuGene nothing is
    significant and the OIS difference is NEGATIVE (-4.5).
  - Even on scaled_diff the direction is inconsistent across conditions:
    Proliferating +31.1 (FDR 6.5e-4) and OIS +28.2 (FDR 0.029) are significant,
    CICQ is +16.3 (p 0.36) and SSCQ is -13.8 (p 0.70). So the sentence must not
    imply a uniform elevation.

Results, after the model-preference sentence in 2.1.5.1:

> Immortalisation status was unevenly distributed across conditions: hTERT-
> immortalised samples were present among proliferating controls (21 of 91), in both
> quiescence conditions (3 of 19 and 5 of 22) and in oncogene-induced senescence (17
> of 48), but absent from stress-induced senescence. On the
> scaled-difference model immortalised samples scored higher than primary ones among
> proliferating controls (+31.1 units, FDR = 6.5x10-4) and in oncogene-induced
> senescence (+28.2, FDR = 0.029), though not in either quiescence condition (+16.3
> and -13.8, both n.s.), and no difference was significant on YuGene. Because
> immortalisation status is perfectly nested within study, none of these differences
> can be separated from the laboratories and protocols that use immortalised lines,
> and we do not interpret them as effects of immortalisation.

Discussion:

> The apparent elevation of transcriptomic age in hTERT-immortalised cells warrants
> dedicated study. It cannot be resolved in a meta-analysis, because immortalisation
> status is confounded with study; it requires parental and immortalised cells of the
> same strain arrested and profiled together. Doing so would also establish whether
> the transcriptomic age of senescent cells depends on immortalisation status, which
> matters for interpretation well beyond this dataset, since immortalised fibroblasts
> are widely used in oncogene-induced senescence work - 5 of the 15 such studies
> analysed here - and replicative senescence cannot be modelled in them at all.

---

## 2.1.5.2 Pathway-level Transcriptomic Age Differences

### What was wrong

Every contribution in this section came from the same pooled contrasts, and a set's
contribution is a weighted sum over its genes of the same control-subtracted
matrix, so nothing insulated it. Recomputed within study (`exploratory/18`,
same design as above, BH across 50 sets x 5 conditions within each model), of the
17 interpretable sets: 15 to 16 of 50 sets change sign in CICQ, and 8 to 11 in RS.
Of the sets the section reported, these specific claims fail:

| claim in draft | within study |
|---|---|
| KRAS SIGNALING UP "largest single effect", CICQ +8.5/+4.8 | holds, +7.2/+4.8, 6/6 studies both models |
| TNFa "reverses sign between quiescence and senescence" | **holds and is now cleaner** - negative in CICQ and SSCQ, positive in SIPS and OIS |
| GLYCOLYSIS "positive in CICQ (+1.8/+3.3)" | **fails** - not significant in CICQ within study |
| APOPTOSIS "significant effects fall on exactly the three senescence subtypes" | **fails** - significant in RS, SIPS, OIS but no longer the clean split, and OIS agreement is only 9/15 studies |
| P53 PATHWAY "most strongly in SIPS and SSCQ (+4.0/+4.4)" | **SSCQ fails**; RS gains (+5.0/+6.0, 6/6 both models) |
| EMT and XENOBIOTIC METABOLISM "restricted to OIS" | holds for both; but XENOBIOTIC in CICQ now disagrees in sign between models and is excluded |
| ESTROGEN RESPONSE LATE in SSCQ (-0.9/-0.7) | **fails** - models disagree in sign |
| FATTY ACID METABOLISM in CICQ | **fails** in CICQ; holds in SIPS and OIS |
| HYPOXIA "single-condition, SIPS" | now significant in RS, SIPS and SSCQ |

Note also that the within-study permutation test is **more powerful** than the
pooled Wilcoxon, because it uses means and pools information across studies. So the
raw count of significant sets goes up (25-41 of 50 pooled, 26-43 within study) and
that increase is **not** evidence of robustness. What carries evidential weight is
sign stability and per-study consistency, which the new output reports.

### Reporting rule used

A set is reported only if it is interpretable (five largest genes carry <= 65% of
its contribution under both models, 17 of 50 sets), significant within study at
FDR < 0.05 on **both** models, and **agrees in sign between the two models**. The
last requirement is new and removes three otherwise-significant results where the
two normalisations disagree in direction (CICQ XENOBIOTIC METABOLISM, SSCQ EMT,
SIPS ESTROGEN RESPONSE LATE); when the direction depends on the normalisation it
is not identifiable. Per-study direction agreement is reported as a descriptive
count rather than used as a threshold, deliberately - an agreement cutoff cannot be
justified at a particular value and behaves differently at n = 4 and n = 15 strata.

Counts passing: CICQ 8, SIPS 9, OIS 13, RS 5, SSCQ 5 of 17.

### Replacement text

> To identify which programmes carry these shifts, we decomposed each condition's
> tAge prediction into per-gene-set contributions using the EN model's own
> coefficients (see xx in Methods). Because the clock is sparse, a gene set can only
> support a set-level claim if its contribution is spread across many genes rather
> than dominated by a few; we therefore interpreted a set at set level only where
> its five largest-contributing genes accounted for no more than 65% of its total
> contribution under both models, retaining 17 of 50 sets (remainder in SI Table x).
> Contributions were estimated within study exactly as for whole-transcriptome tAge
> above, are reported as the difference in mean contribution between groups in the
> same species-adjusted tAge units, and are called significant only where both
> models agree in sign at FDR < 0.05.
>
> KRAS SIGNALING UP (155 clock genes, 30% and 28% contributing, top five carrying
> 58%) was the largest single effect, raising tAge contribution by +7.2/+4.8 units
> in contact-inhibited quiescence, +5.1/+2.4 in serum-starved quiescence and
> +6.4/+4.6 in oncogene-induced senescence, in every contributing study in all
> three, while not reaching significance in replicative or stress-induced
> senescence - a pattern grouping the quiescence conditions with OIS rather than
> separating quiescence from senescence. MYOGENESIS (top five 46%) was significant
> in the same OIS direction (+6.9/+3.1) plus SIPS (+3.9/+1.4) and CICQ (+2.6/+3.8).
>
> TNFa SIGNALING VIA NFkB, the second least gene-dominated set tested (164 genes,
> top five 44%), reverses sign between quiescence and senescence: it lowers tAge
> contribution in contact-inhibited (-2.2/-2.9) and serum-starved (-1.4/-1.5)
> quiescence and raises it in stress-induced (+2.1/+0.7) and oncogene-induced
> (+3.5/+1.6) senescence. This is the clearest quiescence-versus-senescence
> dissociation in the panel, it is not attributable to a single gene, and it is
> consistent in direction in every contributing study for three of the four
> conditions.
>
> Three sets contributed in four of the five conditions. COMPLEMENT (127 genes,
> 28%/28% contributing, top five 60%) rose in CICQ (+3.6/+1.9), SSCQ (+2.2/+1.5),
> SIPS (+3.0/+3.5) and OIS (+3.1/+3.8); the curated lysosomal set (138 genes,
> 16%/25% contributing, top five 55%), the only non-Hallmark set tested, rose in
> CICQ (+2.9/+2.7), SSCQ (+2.0/+2.6) and OIS (+2.1/+2.5), and reached significance
> in SIPS as well but weakly and in only 2 of its 4 studies (+0.6/+1.6). Together
> with TNFa these are the only sets contributing across four conditions, and no set
> contributes in all five.
>
> Replicative senescence was distinguished by a damage-response and
> structural-remodelling signature that was consistent in every contributing study:
> P53 PATHWAY (+5.0/+6.0), EPITHELIAL MESENCHYMAL TRANSITION (+4.1/+2.0),
> APOPTOSIS (+3.8/+2.1) and HYPOXIA (+2.3/+1.2), with GLYCOLYSIS lowered
> (-1.9/-2.1). P53 PATHWAY was also the strongest stress-induced-senescence
> contribution (+3.6/+5.4, 4/4 studies both models), alongside APOPTOSIS
> (+2.1/+1.6) and MYOGENESIS (+3.9/+1.4).
>
> Two large effects were restricted to oncogene-induced senescence, both on
> well-spread sets: EPITHELIAL MESENCHYMAL TRANSITION, the least gene-dominated set
> tested (153 genes, 36%/37% contributing, top five 41%; +7.2/+5.7 in all 15
> studies), and XENOBIOTIC METABOLISM (top five 52%; +3.1/+7.4, also all 15).
> GLYCOLYSIS was lowered in OIS (-2.0/-1.8, 14 of 15 studies) as in replicative
> senescence. UV RESPONSE DN met the interpretability criterion but showed no
> effect significant on both models in any condition.
>
> These results describe how the clock's prediction is composed rather than which
> programmes drive ageing. Because elastic-net regularisation distributes weight
> semi-arbitrarily among correlated genes, the genes carrying a set's contribution
> are partly a property of model fitting, and the partial contributions neither sum
> to the total prediction nor partition it, so they cannot be read as fractions of a
> condition's overall tAge shift.

---

## 2.1.5.3 Universal Transcriptomic Age Across the Time Course

**No change required on account of the baseline confound.** The time course is a
single study (ERP021140, all 72 samples), so every comparison is already within
study, and each cell type is compared with its own untreated controls. The
between-study baseline problem cannot arise.

One disclosure to add, because it affects both arms:

> The time-course dataset also contributes 30 samples to the cross-sectional
> meta-analysis above (6 proliferating, 6 serum-starved quiescent and 18
> stress-induced senescent), i.e. 18 of the 39 stress-induced senescence samples.
> The two analyses are therefore not independent.

---

## 2.1.5.4 Pathway-level Transcriptomic Age Across the Time Course

The temporal numbers all stand. Two claims must be deleted, both asserting
cross-arm comparability, which `TAGE_CALCULATION_AUDIT.md` §5 lists as forbidden:
the runs differ in gene retention (12,236 vs 11,298-11,614), in measured
coefficient weight (80.5% vs ~70%) and in input scale (SD 0.427 vs 0.228-0.328),
and each is referenced to a different control group.

- Delete from the opening: > "and on the same tAge-unit scale, so temporal and cross-sectional contributions are directly comparable."

  Replace with: > "and on the same tAge-unit scale. Because each preprocessing run
  retains a different gene set and is referenced to its own control group,
  magnitudes are compared within an arm only, not between the time course and the
  meta-analysis."

- Delete from the Methods addition: > "and comparable between the meta-analysis and the time course."

  The rest of that Methods paragraph, on why Cohen's d was not used, is unaffected
  and should stay - though note the example it gives (COMPLEMENT d = 18.3 vs 1.6
  for +4.8 vs +3.8 units) is itself a cross-arm magnitude comparison; it is
  legitimate there only because the point being made is that d is *not* comparable.

Everything else in the section - the melanocyte divergence, the Isg15 gene-level
result, the fibroblast broadening - rests on within-cell-type contrasts inside one
study and on rank/direction rather than cross-arm magnitude, and is unaffected.

---

## Methods

Add:

> **Within-study estimation.** Because untreated proliferating controls differed
> systematically between studies and arrest condition was confounded with cell
> strain, condition effects and gene-set contributions were estimated within study.
> For a study s contributing n_t arrested and n_c control samples, d_s was the
> difference in group means; these were combined as sum(w_s d_s)/sum(w_s) with
> w_s = n_t n_c/(n_t + n_c). Null distributions were generated by permuting the
> condition label within each study (20,000 permutations for whole-transcriptome
> tAge, 10,000 for gene sets, the latter shared across sets within a condition and
> model so that between-set correlation is preserved), which holds every
> study-level baseline fixed. Each study contributes a single cell line, so cell
> strain, hTERT-immortalisation status, tissue of origin and batch are held constant
> by construction. Studies without internal proliferating controls (2 of 34, 7 of
> 230 samples) were excluded from these tests.

> **Cell line and immortalisation annotation.** hTERT-immortalisation status and
> cell strain were taken from the GEO SOFT records of all contributing studies
> (sample- and series-level) rather than from the aggregated recount3 metadata,
> which carries neither the growth- nor treatment-protocol fields and mislabels
> both attributes. 46 of 230 samples across 8 studies are hTERT-immortalised.
> Inducible oncogene constructs were not counted as immortalisation. Per-study
> annotations and the establishing quotations are in SI Table x.

Delete the cross-arm comparability clause noted under 2.1.5.4.

---

## Still open

1. **The pathway panel is close to saturated.** 26-43 of 50 sets reach FDR < 0.05
   in a given condition. Because a set's contribution largely tracks the global
   shift, "significant in at least one condition" carries little information. The
   interpretability gate and the both-models sign rule reduce this, but a
   set-specific null (the outstanding label-permutation item in
   `SESSION_HANDOFF_BULK_TAGE.md` §9) is still the right fix and is still not run.
2. **Figures not yet regenerated.** `figure_universal_tage_differences.png` and
   `figure_pathway_heatmap_all_both_models.png` both display pooled contrasts and
   are now inconsistent with the text.
3. **The 26 primary studies rest on absence of an immortalisation statement**
   rather than a positive declaration of primary status.
4. **Induction protocol still varies within condition across studies**; within-study
   estimation removes baseline confounding, not protocol heterogeneity.

---

## 2.1.5.2 SECOND REVISION: the yardstick result (2026-08-21)

`exploratory/20_pathway_specificity_yardstick.py` supplies what the section was
missing - a scale against which a set's contribution can be called large or small.
The result substantially demotes the section, so it is recorded separately from the
first revision above rather than folded into it.

### The yardstick

A set's contribution is a sum of coef_i x z_i, so what it should contribute if it
were unremarkable is set by how much of the clock's coefficient weight it holds:

    expected_s = (sum |coef| in set / sum |coef| overall) x whole-transcriptome shift

Significance comes from random sets matched to the real set on gene count AND on
the distribution of |coef| (stratified by coefficient decile, zero-coefficient
genes forming their own stratum), drawn from the clock's own features, B = 20,000,
all computed on the within-study contrasts. Validation: the per-gene decomposition
sums exactly to the whole-transcriptome within-study effects of script 13
(CICQ +5.88, SSCQ -3.43, RS +42.33, SIPS +23.29, OIS +33.30 on scaled_diff), and
the matched null's mean tracks the weight-share expectation (r = 0.93).

Why this null is legitimate where the earlier size-matched one was not: that null
was asked "does this set have an effect", for which a pool of age-selected genes is
unfair. This one is asked "does this set carry more of the shift than an arbitrary
equally-weighted slice of the same clock". For a question about disproportion,
other clock genes are the correct comparator, and no alternative exists - a
non-clock gene has coefficient exactly 0 and cannot contribute.

### The result: almost nothing is disproportionate

Of 85 interpretable set x condition tests per model, **two** survive BH, both on
the scaled-difference model and both the same set: KRAS SIGNALING UP in
contact-inhibited quiescence (observed +7.17 against +0.14 expected, z = 4.24, FDR
= 0.030) and in serum-starved quiescence (+5.06, z = 3.85, FDR = 0.045). **Nothing
survives on both models.**

Across all interpretable tests the median |z| is 0.22 (scaled) and 0.26 (yugene),
and only 11% and 8% exceed |z| = 2 - close to what chance produces. The 17
interpretable sets hold 40% of the clock's coefficient weight and carry between
43% and 126% of each condition's shift.

Strongest sub-threshold cases, for completeness: OIS XENOBIOTIC METABOLISM
(yugene z = 3.97, FDR 0.127), OIS KRAS SIGNALING UP (scaled z = 3.12, FDR 0.095),
OIS MYOGENESIS (scaled z = 3.03, FDR 0.095), SIPS P53 PATHWAY (yugene z = 3.09,
FDR 0.293).

### What this means for the text

The contributions reported in the first revision are real - non-zero,
direction-consistent across studies - but they are **proportional to each set's
weight in the clock**. The pathway decomposition therefore describes how the
prediction is COMPOSED; it does not identify programmes that are specifically
implicated. Three consequences:

1. Remove "the largest single effect" and any comparative ranking of sets. A set's
   contribution being large is mostly a statement about its size and weight.
2. KRAS SIGNALING UP is the one set that exceeds its share, and only on
   scaled_diff - the model this work otherwise argues against. It should be
   reported as the single exception with that caveat, not as a headline.
3. The honest framing of the section is that arrest raises transcriptomic age
   broadly across the clock rather than through identifiable programmes. That is a
   result, and it is the opposite of what the original draft implied.

### Reporting caveat

observed/expected divides by (weight share x total shift), so where a condition's
shift is near zero the ratio explodes and flips sign meaninglessly: CICQ on
scaled_diff (+5.9 units) yields ratios up to 52 and SSCQ on scaled_diff (-3.4) down
to -63. Report z, not ratio. The output carries a `ratio_reliable` flag
(|total shift| >= 10 units), false for CICQ-scaled, SSCQ-scaled and SSCQ-yugene.

---

## Modularity: now scripted, and re-checked (2026-08-21)

`exploratory/21_tage_modularity.R`. Previously computed ad hoc in conversation, so
neither reproducible nor checkable. All statistics are rank-based, which is what
makes them safe across preprocessing runs where magnitudes are not.

| original claim | status | corrected value |
|---|---|---|
| within-cell-type profile similarity 0.62 vs between 0.12 | **holds, numbers wrong** | 0.840 vs 0.272 (scaled), 0.761 vs 0.253 (yugene) |
| three PCA axes for 80% of variance | **partly** | 2 axes (scaled), 3 (yugene) |
| COMPLEMENT the only set significant in all three cell types | **holds** | confirmed, under the interpretable + both-models rule |
| groups with equal aggregate tAge have unrelated composition | **dropped** | cross-run comparison, not recoverable |

### 1. Composition is cell-type specific

Profiles of the same cell type at different timepoints correlate at rho = 0.840
(scaled) and 0.761 (yugene); profiles from different cell types at 0.272 and 0.253.
The gap (0.568, 0.508) is significant against a null that shuffles which cell type
each group belongs to (10,000 permutations, p = 0.0037 and 0.0030). This is the
central modularity result and it now has a test, which the ad hoc version did not.
Per cell type: fibroblast 0.928/0.857, melanocyte 0.867/0.791, keratinocyte
0.724/0.636 - keratinocytes are the least internally consistent over time.

### 2. The dominant axes ARE cell type

More informative than the variance-explained count: of the spread on PC1, 96.6%
(scaled) and 97.8% (yugene) lies between cell types rather than within them, and on
PC2 87.8% and 69.1%. Melanocytes sit apart on PC1 (mean score +5.6 against
keratinocyte -3.8 and fibroblast -1.8), which is the same divergence 2.1.5.4
describes set by set. PC3 is not cell-type organised (12.7%, 29.7%).

Note the PCA is run on per-group standardised profiles, so the shared magnitude of
the tAge shift is removed before decomposition; PC1 is therefore compositional, not
the global shift.

### 3. Breadth

Under the rule the text uses - interpretable sets, significant on both models -
COMPLEMENT is the only set reaching significance in all three cell types, as
claimed. Eight sets qualify before the interpretability gate (ADIPOGENESIS, APICAL
JUNCTION, APICAL SURFACE, COMPLEMENT, E2F TARGETS, HEDGEHOG SIGNALING, PI3K AKT
MTOR SIGNALING, UV RESPONSE UP) and only COMPLEMENT survives it. Counting loosely -
all 50 sets, either model - gives 14 and 11, so the claim must always be stated
with its qualifiers or it looks wrong.

### 4. New: the meta-analysis conditions are modular too

Condition profiles form two blocks. The quiescence conditions correlate with each
other (CICQ-SSCQ rho = 0.84 scaled, 0.73 yugene) and the senescence subtypes with
each other (RS-SIPS-OIS 0.65-0.85), while across the blocks correlation is low
(CICQ-RS 0.26, SSCQ-RS 0.22-0.34). So composition separates quiescence from
senescence even where aggregate tAge does not - CICQ's aggregate elevation is
comparable to the senescence subtypes on yugene, but what carries it is not.

### Reconciling this with the yardstick result

These look contradictory and are not. The yardstick asks whether an individual set
carries more than a weight-matched slice of the clock WITHIN one comparison, and
almost none does. Modularity asks whether the pattern across sets is the same in
different cell types, and it is not. Were contributions purely proportional to
weight share, every group's profile would be identical and between-cell-type
correlation would approach 1; it is 0.25. So no single set stands out, yet the
ensemble differs by cell type. Both statements should appear together, because
either alone is misleading.

---

## Mortality-only set-level sections: gate, exclusions, and the two tests (2026-08-22)

Three decisions, recorded because each changes reported numbers.

### 1. The interpretability gate is clock-specific

The 65% top-5 gate was originally computed on the chronological clocks and wrongly
applied to mortality-clock results. Gene domination is a property of the model's
sparsity:

| clock | non-zero genes per set | median top-5 share | sets passing |
|---|---|---|---|
| chronological scaled | 20 | 0.729 | 18/50 |
| chronological YuGene | 22 | 0.644 | 26/50 |
| mortality | 125 | 0.246 | 49/50 |

Only HALLMARK PANCREAS BETA CELLS fails on the mortality clock. The old gate
discarded 32 well-represented sets and set the BH family to 85 tests where it
should have been 245 (meta) and 441 (temporal). `exploratory/14 --mortality` now
emits `pathway_representation_mortality.csv`; scripts 20 and 22 read the gate
matching their clock.

### 2. Proliferation-associated sets are excluded

E2F targets, MYC targets, G2M checkpoint and mitotic spindle are dropped from
interpretation, leaving 45 sets. Every contrast in this work is arrested against
dividing cells, so these sets must move by construction and their movement is
confirmation that arrest occurred, not a finding. This removes E2F TARGETS in SIPS
from the cross-sectional hits (3 -> 2) and E2F TARGETS in fibroblasts and
keratinocytes plus MYC TARGETS in fibroblasts from the temporal hits (12 -> 9).

### 3. Two different tests, and where they disagree

The Wilcoxon asks whether a set's contribution differs from zero; the matched-null
yardstick asks whether it exceeds what the set's weight in the clock predicts. The
first is nearly always significant, because the whole transcriptome shifts. Only
the second is evidence of specificity, and the two can disagree.

HYPOXIA is the case that matters. It passes the Wilcoxon decisively in all three
cell types (FDR 1e-4 to 1e-3) AND runs about elevenfold above its expected share
in all nine cell type x timepoint tests, reaching nominal significance in seven,
but no single test survives BH over 441. P53 PATHWAY is roughly twice as large per
test (median obs/exp 19 versus 11, median z 3.4 versus 2.1) and does survive. So
the difference is multiplicity, not effect size, and hypoxia is reported as a
consistent but sub-threshold excess following the same pattern as p53 rather than
as an independent result. In the cross-sectional arm hypoxia is likewise the next
strongest set in exactly the two conditions where p53 survives (RS z = 3.20,
SIPS z = 2.86; nominal p = 0.003 and 0.006).

An earlier note in this file described hypoxia as "consistently non-zero, not
disproportionate". That was wrong: it is disproportionate, by about elevenfold,
and fails only on correction.

### Current hit lists

Cross-sectional (245 tests, 45 sets x 5 conditions, cell-cycle sets excluded):
P53 PATHWAY in SIPS (z 4.61, FDR 0.012) and RS (4.19, 0.018). Nothing else.

Temporal (441 tests): P53 PATHWAY in melanocytes at 4, 10 and 20 days (4.01, 4.31,
5.35) and fibroblasts at 4 days (4.03); ANGIOGENESIS in keratinocytes at 10 and 20
days (4.11, 4.60); ADIPOGENESIS in keratinocytes at 10 days (3.45); and negative
excesses in fibroblasts for IL6 JAK STAT3 SIGNALING (-3.93) and TNFA SIGNALING VIA
NFKB (-3.31) at 20 days. Keratinocyte 4-day and fibroblast 10-day P53 fall just
outside (both FDR 0.061).

### REVERSED same day: proliferation sets are NOT excluded

The exclusion recorded immediately above was a misreading of an instruction to
remove chronological-clock material; it was not an instruction to drop cell-cycle
sets. They are restored. The reasoning for keeping them is sound on the clock being
used: E2F targets, MYC targets, G2M checkpoint and mitotic spindle are well
represented in the dense mortality clock (all four inside the gate), and the
yardstick asks whether a set carries MORE than its weight share, which arrest alone
does not guarantee. A set moving because the cells stopped dividing would move by
about its weight share and fail the test; exceeding it is a real observation.

Restored hits, and the gate stays at 49 sets, so no BH value changes - the earlier
removal was a reporting filter applied after the fact, not a re-run:

  cross-sectional, 3 of 245: P53 PATHWAY in SIPS (z 4.61) and RS (4.19), and
  E2F TARGETS in SIPS (4.12).

  temporal, 12 of 441: adds E2F TARGETS in fibroblasts (4.10) and keratinocytes
  (4.16) at 4 days, and MYC TARGETS contributing significantly LESS than its weight
  predicts in fibroblasts at 4 days (-3.74).

E2F TARGETS appearing in stress-induced senescence and in both fibroblasts and
keratinocytes at 4 days is coherent with the p53 result, since the two are
mechanistically linked through p21. Note the direction guard still applies: a
positive contribution can arise from downregulated genes with negative
coefficients, which is what cell-cycle genes do during arrest, so these remain
statements about the contribution to the prediction rather than about pathway
activity.

### Gate dropped entirely for the mortality sections (2026-08-24)

The 65% top-5 gate is removed from both set-level sections and all 50 sets are
reported. Reasoning: on the mortality clock the gate excluded one set,
HALLMARK PANCREAS BETA CELLS, and that set reaches nothing on the yardstick - 0 of
5 cross-sectional tests and 1 of 9 temporal at nominal p, none after correction.
Including it changes no result: hits stay at 3 and 12 with families of 250 and 450
instead of 245 and 441.

The quantity the gate measures is a chronological-clock problem that does not exist
on the mortality clock. A typical Hallmark set has 125 genes in the clock's feature
space on all three models, but only about 20 carry non-zero weight on the
chronological clocks against all 125 on mortality, so five genes are 75% of a
20-gene effective set there and 25% of a 125-gene one here. Reporting a criterion
that filters nothing invites the reader to think it did work.

The gate remains relevant to any chronological set-level analysis, where it
excludes 33 of 50 sets, so the criterion stays documented in Methods and in
exploratory/14.

Text now reads: "between 19 and 38 of the 50" sets significant per condition;
"Three of 250" and "Twelve of 450" tests surviving; "Of the 50 sets, only four
contribute significantly in all three cell types"; "Thirty-one of the 50 differ in
sign between cell types". Hypoxia's family reference updated to 450.

NUMBER TO WATCH. Four medians exist for the top-5 concentration and they are easily
confused: scaled difference 0.729, YuGene 0.644, the per-set MAXIMUM across those
two 0.750, and mortality 0.246. The 0.750 figure is a composite that describes
neither chronological clock and is meaningful only as the gate criterion, since the
gate required a set to pass on both. Prose comparing "the chronological clocks"
must quote 73% and 64%, not 75%. An earlier draft of the replacement sentence used
75% and was corrected.

### "Most sets reach significance" was wrong, and the fix is stronger

The claim overstated: only two of five conditions reach a majority of the 50 sets
(SIPS 38, OIS 36), with CICQ at 25, SSCQ 21 and RS 19. Checking what the count
actually tracks gives a better version of the same argument:

  condition  sig sets  n samples  whole-transcriptome effect
  RS             19        11              +1.183
  SIPS           38        39              +0.763
  OIS            36        48              +0.634
  CICQ           25        19              +0.413
  SSCQ           21        22              +0.167

The count correlates with SAMPLE SIZE at rho = 0.80 and with the whole-transcriptome
effect at rho = 0.00. Replicative senescence has the largest effect of any condition
and the fewest significant sets, because it has the fewest samples. So the Wilcoxon
count measures power, not biology, which is a sharper reason to disregard it than
"most sets are significant" ever was - and it is the direct motivation for the
weight-matched yardstick.

### The same caveat holds in the temporal arm, for a different reason

Asked whether the power explanation carries over. It cannot: every temporal
comparison is 6 versus 6, so n is constant. With n fixed the count instead tracks
the whole-transcriptome shift, but only loosely (rho = 0.49 over the nine cell type
x timepoint groups, against 0.80 for sample size in the cross-sectional arm):

  Fibroblast    4d 28/50 (shift 0.481)   10d 30/50 (0.618)   20d 41/50 (0.710)
  Keratinocyte  4d 17/50 (0.486)         10d 22/50 (0.326)   20d 20/50 (0.192)
  Melanocyte    4d 28/50 (0.709)         10d 29/50 (0.256)   20d 21/50 (0.413)

Fibroblasts at 20 days give 41 of 50 on a shift of 0.71 while keratinocytes at 4
days give 17 on a comparable 0.49, so the relationship is weak. The conclusion is
unchanged - the Wilcoxon count is not a readout of which sets matter - but 2.2.4
now states its own reason rather than inheriting the cross-sectional one, which
would be wrong there.

### Two assertions about the significance count, tested (2026-08-24)

Asked for proof of two things the text had been asserting.

CLAIM 1: "the whole transcriptome shifts, so most sets follow". FALSE as written.
Gene-level contributions to the mortality tAge difference are large in BOTH
directions and largely cancel; the net shift is the residual:

  cond  sum(+)  sum(-)     net   cancellation
  CICQ   +7.26   -6.85   +0.413      34x
  SSCQ   +6.01   -5.84   +0.167      71x
  RS     +7.46   -6.28   +1.183      12x
  SIPS   +5.02   -4.26   +0.763      12x
  OIS    +6.98   -6.35   +0.634      21x

Only about 50% of measured clock genes move in the net direction - indistinguishable
from chance - and 3 to 20 genes account for half the net shift. There IS a
transcriptome-wide change (14-29% of 19,439 tested genes are differentially
expressed, 2,667-5,567 per condition) but it is not a coherent shift that carries
every set with it. The sentence has been removed.

CLAIM 2: "the count tracks power rather than effect size". PARTLY true; the clean
version was also an assertion. Correlation cannot settle it because n and effect
size are themselves correlated across the five conditions (rho = 0.60). Subsampling
every condition to 11 versus 11, 200 draws, does settle it:

  cond  full n  count(full)  count(11v11)
  CICQ      19           24     6 [1-13]
  SSCQ      22           20     2 [0-9]
  RS        11           17    11 [8-14]
  SIPS      39           40    15 [9-23]
  OIS       48           36    16 [10-22]

Matching n compresses the spread from 23 to 14 but does not remove it, and it
REORDERS the conditions: RS rises from last to mid-rank, so its apparently low count
is indeed a power artefact and per sample it is among the strongest; SSCQ falls to
last, consistent with its small shift. So both factors operate and neither alone
explains the counts, which is what the text now says.

---

## Decomposition diagnostics, scripted and tested (2026-08-24)

`exploratory/23_decomposition_diagnostics.py` and `24_figure_cancellation.R`. These
replace ad hoc checks and two claims that had been asserted. B = 10,000 draws, chosen
from the Monte Carlo error rather than by habit: MCSE of a median is 1.253 SD/sqrt(B),
so about +-0.06 sets here against +-0.4 at the B = 200 an earlier version used. That
200 had no statistical justification - the Mann-Whitney call simply had not been
vectorised, and scipy computes all 50 sets in one call across an axis.

### Cancellation: measured, and it applies to BOTH arms

Summed positive and negative per-gene contributions against the net:

  cross-sectional   CICQ 34x   SSCQ 71x   RS 12x   SIPS 12x   OIS 21x
  temporal          17x to 87x across the nine cell type x timepoint groups

In every one of the fourteen groups, contributions in either direction reach 4 to 9
units while the net is 0.17 to 1.18, and the fraction of measured clock genes moving
in the net direction is 49.8% to 51.0% - indistinguishable from chance. So the
question "does the temporal arm behave the same way" is answered yes, and if anything
more strongly: keratinocytes at 20 days show the largest cancellation of any group
(87x on a net of +0.19).

### A test I designed that could not work

The temporal arm was initially subsampled to 3 v 3 to mirror the cross-sectional
test. That is void: the smallest attainable two-sided Mann-Whitney p at 3 v 3 is
2/C(6,3) = 0.100, so nothing can reach significance after correction whatever the
data, and it returned nine zeros that were arithmetic rather than biology. The script
now refuses to subsample where the floor exceeds 0.05 and says so. For reference the
floor is 2.8e-6 at 11 v 11 (cross-sectional, usable) and 2.2e-3 at 6 v 6.

### What the counts track, tested separately per arm

  cross-sectional, n varies 11 to 48: count vs n rho = 0.80, count vs effect size
  rho = 0.60 - but n and effect size themselves correlate at 0.60, so correlation
  cannot separate them. Subsampling to 11 v 11 can, and does: counts fall from
  24/20/17/40/36 to 6/2/11/15/16, compressing the spread from 23 to 14 and reordering
  the conditions. Replicative senescence rises from last to mid-rank, so its low count
  is a power artefact and per sample it is among the strongest; serum-starved
  quiescence falls to last, consistent with its small shift.

  temporal, n constant at 6 v 6: no power explanation is available by construction, so
  the counts are compared against effect size directly. Count vs median separation
  rho = 0.92; count vs net shift rho = 0.44. Range 12 to 42 of 50.

Note the temporal counts here (12-42) use BH within each group over 50 tests, whereas
the 17-41 quoted earlier used BH within the whole temporal analysis over 450. Quote
the per-group figures, which are what script 23 reproduces.

### The claim that IS supported

Differential expression, against 19,439 tested genes: CICQ 5,567 (28.6%), SSCQ 4,263
(21.9%), OIS 4,185 (21.5%), RS 3,775 (19.4%), SIPS 2,667 (13.7%). So a
transcriptome-wide change is real; what is not supported is that it is a coherent
shift carrying every set with it.

### BH family audit across all scripts (2026-08-24)

Checked every p.adjust/BH call in the arm after noticing the two sections quoted
counts from different families.

CONVENTION, used consistently by exploratory/05, 18 and 22: correct WITHIN AN
ANALYSIS. That gives 250 tests for the five arrest conditions (50 sets x 5), 450 for
the nine temporal groups, 150 for the three pooled temporal comparisons. The
whole-transcriptome scripts follow the analogous rule, correcting within model:
meta_analysis/05 uses 10 tests vs-Proliferating and 30 for all pairs; 09, 11, 12, 13
and 16 all correct within model; temporal_analysis/07 uses one 36-test family and 08
one 18-test family.

THE INCONSISTENCY. exploratory/23 corrected within each GROUP over 50 tests, and
2.2.4 had been quoting its 12-42 range while 2.1.5.2 quoted 19-38 from the
analysis-wide family. Same kind of statement, two different families.

RESOLVED, with ONE family and no dual reporting. Script 23 uses the analysis-wide
family throughout, including inside the subsampling loop: for each draw the resampled
condition's 50 p-values are combined with the other four conditions' unchanged
p-values to form the full 250, BH is applied to that, and the count is taken among the
resampled condition's sets. The other conditions supply the rest of the family exactly
as they do in the real analysis, so observed and subsampled counts sit on the same
footing and match the figures quoted in the Results. Verified equal to the established
convention in all 14 groups.

Under the single family the subsampling conclusion is unchanged but the numbers move:
observed 19/21/25/36/38 (RS, SSCQ, CICQ, OIS, SIPS) fall to 15/8/10/18/17, compressing
the spread from 19 to 10. RS still rises from last to mid-rank and SSCQ still falls to
last. MCSE 0.027-0.053 sets at B = 10,000, so the draw count is ample.

2.2.4 updated to 17-41 and rho 0.94/0.49. The meta count range is no longer quoted in
2.1.5.2, the rewritten paragraph having dropped it - that avoids the earlier "most
sets" error and loses nothing, since the subsampling result is what carries the point.

### DEFINITION TO GET RIGHT: the cancellation ratio is BOTH sides, not one

`cancellation_ratio` in both diagnostics files is (sum of positive - sum of negative)
/ |net|, i.e. the TOTAL movement in either direction divided by what survives it.
Per side the figure is roughly half that:

                            combined ratio    each side alone
  within a gene set          6.0 (median)      3.1 up, 3.0 down
  whole transcriptome        12 to 87          5 to 44

Prose saying "each side is N times the net" must use the per-side column, not the
ratio. A draft said "what is left is about six times smaller than either side" and
"the two sides being twelve to eighty-seven times larger than the difference", both
wrong by a factor of two for exactly this reason. The text now avoids the ratio
altogether and gives a worked example instead: in CICQ, +7.3 up, -6.9 down, +0.4 net.

### No BH on the simulation nulls (2026-08-24)

House convention applied: Monte Carlo nulls are reported as raw empirical p with the
p-floor and an effect size, not BH-adjusted. A simulation null asks whether one set
beats its own matched background; it is not a draw from a discovery family, and
adjusting it conflates the two. Scripts 20 and 22 now print the excess over chance
instead, and their p_emp_adj column is renamed p_emp_adj_DEPRECATED and kept only for
provenance.

What replaces correction, in the text:

  excess over chance   cross-sectional, of 250 comparisons: 40 at p<0.05 against 12
                       expected, 17 at p<0.01 against 2, 3 at p<0.001 against none.
                       Temporal, of 450: 74, 31 and 10 against 22, 4 and none.
  effect size          z against the matched null, quoted with every result.
  recurrence           in the temporal arm HYPOXIA exceeds expectation in 7 of the 9
                       groups and P53 PATHWAY in 6, which is stronger evidence than
                       any single comparison.

The interpreted results are unchanged, because the three cross-sectional comparisons
at p < 0.001 are exactly the three that had survived BH: SIPS P53 (z 4.61), RS P53
(4.19), SIPS E2F (4.12). The framing changes from "three survive correction" to "the
excess over chance is concentrated in three strong results", which is both more honest
and more informative.

NUMBER CORRECTED: a draft gave HYPOXIA as z = 3.20 and 3.12 in RS and SIPS. The SIPS
value is 2.86; 3.12 is SIPS MYOGENESIS. Temporal maxima are HYPOXIA 2.58 and P53 5.35,
quoted as 2.6 and 5.4.

### The interpretation rule, replacing correction (2026-08-24)

Dropping BH left the question of which comparisons to interpret. Reporting "40 reach
p < 0.05 but we discuss 3" is selection after the fact and would rightly be read as
cherry-picking. The rule is therefore fixed by the size of the family, not chosen
after looking:

  interpret a comparison only at p < 1/N, N being the number of comparisons, so that
  fewer than one is expected to pass by chance.

That is p < 0.004 for the 250 cross-sectional comparisons, giving 8, and p < 0.0022
for the 450 temporal ones, giving 14. Both scripts now compute an `interpreted`
column so the rule lives in code rather than in prose, and print the conventional
thresholds alongside for context.

This is more permissive than BH was and reinstates results BH had discarded:

  cross-sectional, 8: P53 PATHWAY in SIPS (z 4.61) and RS (4.19); E2F TARGETS (4.12)
  and MYOGENESIS (3.12) in SIPS; HYPOXIA in RS (3.20); and three depletions - EMT in
  both quiescence conditions (-3.30, -3.08) and MYC TARGETS in SIPS (-3.30).

  temporal, 14: P53 PATHWAY in six of the nine groups, now including keratinocytes at
  4 days and fibroblasts at 10 days, which BH had put "just outside" at FDR 0.061.

The headline is unchanged and slightly strengthened: everything elevated in the
cross-sectional arm is confined to the two damage-triggered subtypes, and p53 now
appears in every cell type of the time course rather than two of three.
