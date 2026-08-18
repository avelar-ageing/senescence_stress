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
