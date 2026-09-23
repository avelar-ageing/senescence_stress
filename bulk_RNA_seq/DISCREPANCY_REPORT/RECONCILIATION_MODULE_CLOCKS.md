# Reconciliation: bulk vs scRNA arm on the paper's module clocks

Written 2026-08-17 after reading `scrna_seq/GSE226225/final_analysis/PARTIAL_TAGE_SCRNA_REPORT.md`
and re-testing both arms' claims. **Both arms had a wrong claim. Neither was wrong in the way it
thought.** Net effect: one of my conclusions is retracted, one of theirs is invalidated, and a
**defect in the published supplement** is identified that neither arm had found.

---

## 1. Their validation is correct — I reproduced it, and it retracts my §6b diagnosis

`PARTIAL_TAGE_SCRNA_REPORT.md` §2.3 claims the published coefficients live in the same standardised
space as the fitted `.pkl`, validating to `max|diff| = 0.000000, r = 1.000000`. **Confirmed
independently:**

MOESM7 sheet **(A)**, column 27 (`Chronological Age, Multi-species / Multi-tissue, Scaling`):

| check | result |
|---|---|
| non-zero coefficients | **1,839** — exactly the pkl's count |
| max \|published coef − pkl coef\| | **9.9e-17** |
| intercept difference | 2.1e-09 |
| max \|reconstruction − `model.predict()`\| | **2.6e-07** |
| r | 0.9999999999999997 |

So published coefficients **are** the pkl's coefficients, and applying them with the pkl's own
imputer + scaler is the correct procedure.

**Therefore my `PARTIAL_TAGE_VS_PAPER.md` §6b diagnosis — "the fitted imputer/scaler is not
published, so the coefficients are unusable" — is WRONG and is retracted.** The scaler is
recoverable, because `StandardScaler` statistics are per-feature and so identical for any gene
regardless of which other genes are in the model.

My supporting evidence was also weak and I withdraw it:

- **The "positive control failure" was a bad test.** I evaluated module composites on the tAge
  package's example data and got r = 0.22–0.48 against `predict_tAge`. That data is 24 Klotho mice
  of near-identical age — with almost no true age variance, the correlation between two different
  clocks measures noise. It was not a valid control.
- **Low r against the global clock is expected, not diagnostic.** A 20–238-gene component clock
  should not track a 10,487-feature composite. Demanding that it does was the wrong bar.
- **Individual module clocks in fact behave coherently on bulk data.** Cohen's d vs Proliferating
  across all five arrest conditions: `ECM organization/EMT` **+1.00, +1.81, +0.78, +1.41, +1.10**
  (consistently positive); `Protein processing in ER/UPR` **−1.02, −2.07, −2.32, −2.65, −1.93**
  (consistently negative); `Interferon signaling` +0.42 to +1.36. The machinery works.

The one real oddity: the `All module genes` composite (1,284 genes) reverses SIPS (−0.20 vs
global +1.13) and OIS (−0.10 vs +0.96). That is worth noting but it is the composite, not the
individual module clocks.

## 2. But their validation does not cover what they used it for

**Sheet (A) is the composite full-transcriptome clocks. The module clocks are sheets (B) and (C).**
Validating (A) shows the *file's coefficient convention* is right; it says nothing about the
*gene↔module assignment* in (B)/(C). Their report generalises the (A) result to the module clocks
and proceeds to §5–§6 findings on that basis. That step is unsupported — and, as below, wrong.

## 3. The published supplement is internally inconsistent (new finding, neither arm had it)

`MOESM7` sheet (C) and `MOESM8` sheet (C) contradict each other about which genes are in which
module.

**The paper's own enrichment table** (MOESM8, "(C) Multi-species module enrich") for `darkmagenta`:

| term | odds ratio | padj |
|---|---|---|
| Interferon Gamma Response | 488.9 | 2.2e-63 |
| Interferon Alpha Response | 594.4 | 8.1e-58 |
| Defense Response to Virus | 103.2 | 7.2e-27 |
| Interferon Signaling (Reactome) | 91.4 | 3.6e-26 |

OR ≈ 594 at padj 1e-58 means essentially the whole 45-gene module is interferon-stimulated genes.

**MOESM7 sheet (C)'s `darkmagenta` column** has exactly 45 non-zero coefficients — the right count —
but they are `Abca2, Actn2, Ankrd10, Arhgdib, Bard1, Ccar1, Cdc25c, Cenpn, Cenpu, Cep250, Ckap2, …`
— cell-cycle and cytoskeletal genes, and **zero** overlap with Hallmark Interferon Alpha Response.

The interferon genes *are* in sheet (C)'s 1,284-gene universe (34 of the 85 mouse-mapped Hallmark
IFN-α genes) but are scattered across **12 different modules**, with darkmagenta receiving none:

| gene | assigned module in MOESM7 (C) |
|---|---|
| Irf7 | orange (Chromatin modification) |
| Ifit2 | orange |
| Irf1 | pink (Mito translation/OxPhos) |
| Cxcl10 | darkred (mRNA splicing) |
| B2m, Eif2ak2, Lgals3bp | turquoise (Innate immunity) |
| Cd74 | blue (Muscle contraction) |
| Ifi35 | ivory (Fatty acid metabolism) |

Scattered roughly in proportion to module size — the signature of a **scrambled row assignment**.

### This is not our extraction

Everything checkable about `11_extract_paper_module_clocks.R` is right:

- Sheet layout read directly from the file and matches the parser (row 2 = module, row 3 =
  annotation, row 5 = intercept, rows 6+ = coefficients, cols 1–2 = gene index)
- Module sizes match the dictionary: 9 of 14 **exact**, rest −1 to −6 (EN zeroing)
- Labels confirmed by **three** independent sources: MOESM7 row 3, MOESM8 sheet D, and the tAge
  package's own `extdata/metadata/Module_to_function_map_54.csv` — all say darkmagenta =
  Interferon signaling, 45 genes
- Union of the 14 module gene sets (1,248) is a clean subset of the composite (1,284)
- Extracted entrez↔symbol pairing is 98.5% identical to the package's `Gene_table_mouse.csv`
  (residual 1.5% are symbol-vintage synonyms, e.g. `Atp5mc3`/`Atp5g3`)

### Two repair attempts, both failed

- **Constant row offset** (−3…+3): darkmagenta IFN-α hits stay 0–3 of 45. No offset works.
- **Sheet (A) gene order** — sheet (C)'s ID columns are near-alphabetical (`Abca2, Acadvl, Acads,
  Acta1, Actn2, Adam15…`) whereas sheet (A)'s are not, suggesting the ID columns were sorted
  independently of the coefficient block. Re-assigning sheet (C)'s coefficient rows to sheet (A)'s
  order restricted to those 1,284 genes gives **3 of 45** — better than 0, still nowhere near the
  ~40 that OR = 594 implies. Rejected.

**Conclusion: the gene↔module assignment in MOESM7 sheets (B)/(C) cannot be reconciled with the
paper's own enrichment tables, and the correct mapping is not recoverable from the published
files.** Sheet (A) is unaffected and verified exact.

### A self-correction

My `MODULE_GROUPING_FINDINGS.md` §5.4 reported the same mismatch but **interpreted it wrongly** —
I concluded the annotations were merely loose "top-enrichment labels" and that WGCNA co-expression
need not group canonical pathway members. The paper's own enrichment table refutes that: darkmagenta
*is* an interferon module. The gene lists in MOESM7 are wrong, not the labels. That section needs
rewriting.

## 4. What this means for each arm

**scRNA arm — module-clock results are not safely interpretable.** Their §5 (between-cluster module
effects: "mRNA splicing −1.82", "mito translation/OxPhos −2.05") and §6 within-cluster table
("mRNA splicing ρ 0.60", "ECM/EMT ρ 0.35") are computed from these gene lists, so the *numbers* are
arithmetically fine but the *labels* are unreliable. In particular:

- §6's row "inflammatory … module clocks **flat** (ρ 0.04)" is attributed to turquoise losing ~60%
  of its genes. That explanation is not needed — the label may simply not describe the gene set.
- The §6 "genuine tension" between Hallmark inflammatory rising and module inflammatory flat is
  **not a biological tension**; it is most likely a mislabelled gene set.
- §5's calibration observation (tens to >100 nominal years vs 2–12 composite gaps) stands
  independently and is a second, separate reason for caution.

**bulk arm — my route-B module-membership analysis is affected the same way.** The label-free
statistics survive (they don't depend on which module is which): 17/224 pass the size-matched null;
**1/112** both-model survivors vs Hallmark's 15/400; cross-model sign agreement **0.55** vs 0.775;
inter-set median |ρ| 0.20 despite perfect disjointness. Those still say the module grouping performs
worse. But no module-labelled result from either arm should be reported.

## 5. Corrections to specific claims in the scRNA report

| their §  | claim | status |
|---|---|---|
| §2.3 | published coefficients validate exactly | **correct, reproduced** — but covers sheet (A) only, not the module clocks |
| §7.1 | "in the bulk data the 50 partials overshoot the total ~2.4×" | **stale.** Did not reproduce. Ratio is unstable: median 0.96 (`scaled`), 1.76 (`yugene`); only 4,251/10,487 clock features sit in any Hallmark set |
| §7.2 | "bulk median \|ρ\| between pathway partials is 0.545, a third of pairs above 0.7" | **stale.** Measured 0.20 per-sample (scaled 0.212, yugene 0.201), 1.5–2.3% above 0.7. Highest alternative reading (effect-size profiles across 17 comparisons) is 0.35–0.37 |
| §7.2 | "much less true of the modules, which barely overlap — a reason to prefer them" | **refuted.** Perfectly disjoint modules give the *same* median \|ρ\| (0.20). Dependence is sample-level structure, not gene overlap |
| §7.4 | "No size-matched null yet … not implemented in either arm" | implemented in bulk, then **excluded on 2026-08-17** — the draw pool is the clock's own age-selected features, so it is not a neutral background, and it selects gene-concentrated sets over well-represented ones. Do not adopt it in the scRNA arm. A **label-permutation null** is the neutral replacement and is outstanding in both arms. See `PARTIAL_TAGE_VALIDITY.md` Test 1 |
| §7.7 / §8 | no yugene module clocks | **confirmed.** Sheet (A) has both normalisations; (B) and (C) are headed "scaling" only |
| §8 | "Module clocks are fully published … earlier belief that they were unavailable was wrong" | **half right.** The coefficients are published and the convention is verified; the gene assignment is broken. Not usable as published |
| §9.3 | both arms should add module clocks | **contraindicated** pending §3 |

Their **incidental finding** that x.2.3.2's numbers come from `scaled_diff` (not yugene), and that
the models diverge substantially at cluster level (ETO cl3: 10.4 vs 44.3 yr), is useful and
consistent with the bulk arm's decision to report both models everywhere.

## 6. One cross-arm claim that needs care, not assertion

Their §9.1 asks whether the arms agree on lysosomal decline. **They are not yet comparable:**

- bulk: Cohen's d **+1.11 to +1.88** (meta-analysis, arrested vs proliferating) and **+2.99 → +5.75**
  (temporal fibroblast, 4→20 d). Positive = lysosomal genes *raise* the tAge contribution in
  arrested/irradiated cells.
- scRNA: Spearman ρ **−0.51** (within-cluster, lysosomal partial score vs full tAge across substates).
  Negative = as a substate's tAge rises, the lysosomal contribution *falls*.

These are different estimands (between-group difference vs within-group gradient) and are not in
conflict, but "both arms find lysosomal decline" would be wrong. A joint statement must name the
comparison in each arm. **Note also** that lysosomal only exists as a gene set in this project
because it was restored in bulk commit `2e8bdfe`; anything computed before that lacks it.

## 7. Recommended actions

1. **Do not report module-clock or module-membership results in either arm** until §3 is resolved.
   Both arms' Hallmark decompositions are unaffected.
2. **Contact the authors** with the §3 evidence — this is a supplement defect worth reporting, and
   the correct membership would make the module analysis genuinely valuable in both arms.
   Alternative source: TACO's "Module clocks" selector may expose membership directly.
3. **scRNA arm: refresh the three stale bulk figures** quoted in §7 (2.4× → unstable ratio;
   0.545/34% → 0.20/1.5–2.3%; size-matched null now exists) and drop the §7.2 parenthetical
   preferring modules.
4. **bulk arm: rewrite** `MODULE_GROUPING_FINDINGS.md` §5.4 (wrong interpretation, see §3) and
   `PARTIAL_TAGE_VS_PAPER.md` §6b (wrong diagnosis, see §1).
5. **Still outstanding in both arms:** a label-permutation null. Proposed, not run.
6. **Sheet (A) is a genuine asset** either way: it means published composite coefficients can be
   used to verify any clock evaluation independently of the `.pkl` files.
