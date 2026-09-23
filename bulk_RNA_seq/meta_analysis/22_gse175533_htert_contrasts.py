#!/usr/bin/env python
"""22_gse175533_htert_contrasts.py

Step 3 of the GSE175533 replication. Predicts all three clocks on the matrices
from script 21 and runs the immortalisation contrasts our own data cannot support.

THE QUESTION. Script 16 finds immortalised proliferating controls +31.1 scaled_diff
units above primary ones (FDR 6.5e-4), null on yugene (+6.6, p=0.14) - but 0 of 34
studies contain both, so that gap is inseparable from a study effect, and study
baselines among proliferating controls span 87.8 units (script 11). Tyshkovskiy et
al. report that hTERT "abolished the increase in tAge" in this dataset, which reads
as the opposite result. It is not the same quantity: theirs is a SLOPE over time
within one line, ours is a LEVEL between lines. This script measures both in one
laboratory, one strain, so the two can be compared directly.

FOUR CONTRASTS, because no single one is unconfounded:

  A dividing-vs-dividing   hTERT (PD 46-109) vs parental (PD 20-37)
      The analogue of our own claim - immortalised against primary cells, both
      proliferating. Confounded by population doubling: the arms cannot be
      PD-matched while both divide, because the hTERT arm starts at the parental
      replicative limit.

  B limit-matched          hTERT (PD 46, 51) vs parental (PD 45-53)
      Matched on doublings instead, which forces the parental arm to be at or
      approaching senescence. Confounded in the opposite direction from A: it
      compares immortalised dividing cells with primary cells that are exhausted.
      A and B bracket the truth rather than either one settling it.

  C PD-46 exact            hTERT PD 46 (n=3) vs parental PD 46 (n=3)
      The only exact PD match in the dataset. Reported for its effect size only:
      the smallest attainable two-sided rank p at 3 vs 3 is 2/C(6,3) = 0.100, so
      this comparison CANNOT return a significant result and no p-value from it is
      interpretable. Stated here so it is not read as a null.

  D1 within-arm slopes     tAge ~ time, each arm separately, ALL of its samples
      The reproduction of the authors' own test. Their Fig. 4b regresses
      batch-adjusted tAge residuals on time in culture within each arm and asks
      whether the rise present in the parental arm is absent in the hTERT arm. All
      samples of each arm are used, as they did: restricting to timepoints shared
      by both arms would discard parental TP8-TP10 (PD 50, 52, 53), which is
      precisely where the parental rise occurs, and would report a spurious null.

  D2 arm effect            tAge ~ timepoint + arm + interaction, shared TPs only
      Our question, with time held constant. Only the shared timepoints can enter
      here, because a model spanning TPs present in one arm alone would attribute
      that arm's unmatched time range to the arm term. The cost is the one D1
      avoids: the parental senescence transition is excluded, so this model must
      not be read as evidence about slopes. D1 and D2 answer different questions
      on deliberately different sample sets.

PSEUDOREPLICATION, AND WHY NO P-VALUE FROM THIS SCRIPT IS INFERENTIAL FOR THE
BETWEEN-ARM LEVEL CONTRAST (added 2026-08-26). Each timepoint contributes three
replicate libraries from the same culture at the same doubling, so the 18 and 30
samples are not 18 and 30 independent observations. Collapsing to timepoint means
(6 hTERT, 10 parental) fixes that layer but NOT the one below it: GSE175533's growth
protocol is a single WI-38 stock from Coriell at PDL 15, propagated continuously and
split every four to five days, so each arm is ONE lineage sampled repeatedly. The
timepoints within an arm are therefore serially dependent, and for a comparison of
LEVEL between arms the independent unit is the lineage - one against one. No valid
p-value exists for that contrast, and the p = 0.0043 reported below is the rank-test
floor under an independence assumption the design does not satisfy. It is retained
as a descriptor of complete separation, not as evidence of significance; the
manuscript quotes the direction and the separation, never the p.
  The authors' own genotype (intercept) term has the same limitation. Their SLOPE
test is better founded, because a regression on time within one lineage uses the
serial structure rather than assuming it away - though autocorrelation still makes
its nominal p optimistic.
  The independent unit for the timepoint-level contrasts is the timepoint (6 hTERT,
10 parental). Sample-level rank tests here are therefore anticonservative in exactly
the way the between-condition tAge tests were, and every group contrast is repeated
on timepoint means (suffix _bytimepoint) with the sample-level version kept only for
comparison. The timepoint-level test is the one to quote. Its cost is power: at 6
versus 5 timepoints the smallest attainable two-sided rank p is 2/C(11,6) = 0.0043,
which is reported alongside so a null is never mistaken for an absence of effect.

TIMEPOINT, NOT DOUBLING, IS THE TIME AXIS in D, because the shared TP schedule is
calendar time in parallel culture whereas the arms' PD ranges hardly overlap; a
model in PD would extrapolate each arm into the other's unobserved range. The
authors likewise report tAge against time and note it tracks time more closely than
doublings or passage.

UNITS. Chronological clocks are multiplied by 122.5, the maximum human lifespan the
model was trained against, exactly as in meta_analysis/05, so the values are on the
same scale as the +31.1 they are being compared with. The mortality clock is left in
native units (see script 18): its target is not a lifespan fraction.

Usage: 22_gse175533_htert_contrasts.py <rerun_dir> <model_dir>
Output: <rerun_dir>/gse175533/gse175533_tage.csv
        <rerun_dir>/gse175533/gse175533_contrasts.csv
"""
import itertools
import math
import os
import sys
import warnings

import joblib
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, t as tdist

warnings.filterwarnings("ignore")

MAX_LIFESPAN = 122.5
SEED, NPERM = 1, 20000
MODELS = {
    "mortality":          ("EN_Mortality_Multispecies_Multitissue_scaleddiff.pkl",  "scaled_diff", 1.0),
    "chrono_scaled_diff": ("EN_Chronoage_Multispecies_Multitissue_scaleddiff.pkl",  "scaled_diff", MAX_LIFESPAN),
    "chrono_yugene_diff": ("EN_Chronoage_Multispecies_Multitissue_yugenediff.pkl",  "yugene_diff", MAX_LIFESPAN),
}


def _patch(imp):
    if not hasattr(imp, "_fill_dtype"):
        imp._fill_dtype = imp.statistics_.dtype if hasattr(imp, "statistics_") else np.float64


def predict(m, path, scale):
    feats = list(map(str, m.feature_names_in_))
    e = pd.read_csv(path)
    sid = e["sample_id"].values
    e = e.drop(columns=["sample_id"])
    e.columns = e.columns.map(str)
    for g in [g for g in feats if g not in e.columns]:
        e[g] = np.nan
    X = e.loc[:, feats]
    Z = m.named_steps["scaler"].transform(m.named_steps["imputation"].transform(X))
    est = m.named_steps["estimator"]
    return sid, (est.intercept_ + Z @ est.coef_) * scale


def min_two_sided_rank_p(n1, n2):
    """Smallest attainable two-sided Mann-Whitney p: complete separation.

    For the between-arm level contrasts this floor is the ONLY thing the p-value
    conveys - see the pseudoreplication note in the module docstring. Do not quote
    it as significance.
    """
    return 2.0 / math.comb(n1 + n2, n1)


def compare(a, b, label, clock, note=""):
    a, b = np.asarray(a, float), np.asarray(b, float)
    p = mannwhitneyu(a, b, alternative="two-sided").pvalue
    floor = min_two_sided_rank_p(len(a), len(b))
    return dict(contrast=label, clock=clock, n_a=len(a), n_b=len(b),
                median_a=np.median(a), median_b=np.median(b),
                diff=np.median(a) - np.median(b),
                mean_diff=a.mean() - b.mean(),
                p=p, p_floor=floor,
                interpretable=bool(floor <= 0.05), note=note)


def ols(X, y):
    """Returns (beta, se, p) with t-distributed p-values; X includes an intercept."""
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    resid = y - X @ beta
    dof = len(y) - np.linalg.matrix_rank(X)
    s2 = resid @ resid / dof
    XtXi = np.linalg.pinv(X.T @ X)
    se = np.sqrt(np.diag(XtXi) * s2)
    tt = beta / se
    return beta, se, 2 * tdist.sf(np.abs(tt), dof), dof


def main(rerun_dir, model_dir):
    D = os.path.join(rerun_dir, "gse175533")
    rng = np.random.default_rng(SEED)
    S = pd.read_csv(os.path.join(D, "gse175533_samples.csv")).set_index("sample_id")

    tage = pd.DataFrame(index=S.index)
    for name, (pkl, variant, scale) in MODELS.items():
        m = joblib.load(os.path.join(model_dir, pkl))
        _patch(m.named_steps["imputation"])
        sid, pred = predict(m, os.path.join(D, f"gse175533_{variant}.csv"), scale)
        tage[name] = pd.Series(pred, index=sid).reindex(tage.index).values
    d = S.join(tage)
    d.to_csv(os.path.join(D, "gse175533_tage.csv"))

    print("=" * 78)
    print("GSE175533: tAge by arm and population doubling")
    print("=" * 78)
    print(d.groupby(["arm", "population_doublings"])[list(MODELS)]
          .median().round(3).to_string())

    TPi = d.timepoint.str.replace("TP", "").astype(int)
    d = d.assign(tp_index=TPi)

    rows = []
    for clock in MODELS:
        v = d[clock]
        # --- A: the analogue of our own claim --------------------------------
        rows.append(compare(v[(d.arm == "hTERT")],
                            v[(d.arm == "parental") & (d.state == "dividing")],
                            "A_dividing_vs_dividing", clock,
                            "hTERT PD 46-109 vs parental PD 20-37; PD-confounded"))
        # --- B: matched on doublings instead --------------------------------
        rows.append(compare(v[(d.arm == "hTERT") & (d.population_doublings <= 51)],
                            v[(d.arm == "parental") & (d.state == "late")],
                            "B_limit_matched", clock,
                            "hTERT PD 46,51 vs parental PD 45-53; parental arm senescing"))
        # --- C: the only exact PD match -------------------------------------
        rows.append(compare(v[(d.arm == "hTERT") & (d.population_doublings == 46)],
                            v[(d.arm == "parental") & (d.population_doublings == 46)],
                            "C_PD46_exact", clock,
                            "3v3: minimum attainable p is 0.100, effect size only"))
    # --- the same contrasts on timepoint means, the independent unit --------
    tp_mean = d.groupby(["arm", "state", "timepoint", "population_doublings"],
                        as_index=False)[list(MODELS)].mean()
    for clock in MODELS:
        v = tp_mean[clock]
        rows.append(compare(v[tp_mean.arm == "hTERT"],
                            v[(tp_mean.arm == "parental") & (tp_mean.state == "dividing")],
                            "A_dividing_vs_dividing_bytimepoint", clock,
                            "timepoint means; the independent unit"))
        rows.append(compare(v[(tp_mean.arm == "hTERT") & (tp_mean.population_doublings <= 51)],
                            v[(tp_mean.arm == "parental") & (tp_mean.state == "late")],
                            "B_limit_matched_bytimepoint", clock,
                            "timepoint means; the independent unit"))
    res = pd.DataFrame(rows)

    print("\n" + "=" * 78)
    print("Contrasts A-C  (positive diff = immortalised arm HIGHER, our direction)")
    print("=" * 78)
    for c in ["A_dividing_vs_dividing", "A_dividing_vs_dividing_bytimepoint",
              "B_limit_matched", "B_limit_matched_bytimepoint", "C_PD46_exact"]:
        s = res[res.contrast == c]
        print(f"\n-- {c} --   {s.note.iloc[0]}")
        print(s[["clock", "n_a", "n_b", "median_a", "median_b", "diff", "p",
                 "p_floor", "interpretable"]].round(4).to_string(index=False))

    # --- D1: the authors' test - within-arm slope, all of each arm's samples ---
    print("\n" + "=" * 78)
    print("Contrast D1  within-arm slope over time, ALL samples of each arm")
    print("            (the authors' test: does the parental rise disappear in hTERT?)")
    print("=" * 78)
    d1rows = []
    for clock in MODELS:
        for xname in ("tp_index", "population_doublings"):
            for a in ("parental", "hTERT"):
                g = d[d.arm == a]
                y = g[clock].values
                x = g[xname].values.astype(float)
                b, se, p, dof = ols(np.column_stack([np.ones_like(x), x]), y)
                d1rows.append(dict(clock=clock, x=xname, arm=a, n=len(y),
                                   x_range=f"{x.min():.0f}-{x.max():.0f}",
                                   slope=b[1], se=se[1], p=p[1], dof=dof))
    d1 = pd.DataFrame(d1rows)
    for xname in ("tp_index", "population_doublings"):
        print(f"\n-- time axis: {xname} --")
        print(d1[d1.x == xname][["clock", "arm", "n", "x_range", "slope", "se", "p"]]
              .round(4).to_string(index=False))
    print("\n  authors' claim reproduced if parental slope > 0 and significant"
          " while hTERT slope is not")
    for clock in MODELS:
        for xname in ("tp_index", "population_doublings"):
            g = d1[(d1.clock == clock) & (d1.x == xname)].set_index("arm")
            par_up = g.loc["parental", "slope"] > 0 and g.loc["parental", "p"] < 0.05
            ht_up = g.loc["hTERT", "slope"] > 0 and g.loc["hTERT", "p"] < 0.05
            verdict = ("REPRODUCED" if par_up and not ht_up else
                       "both rise" if par_up and ht_up else
                       "no parental rise to abolish")
            print(f"  {clock:<19} x={xname:<21} {verdict}")

    # --- D2: our question, time held constant ------------------------------
    print("\n" + "=" * 78)
    print("Contrast D2  tAge ~ timepoint + arm + timepoint:arm, shared timepoints only")
    print("=" * 78)
    shared = sorted(set(d.tp_index[d.arm == "hTERT"]) & set(d.tp_index[d.arm == "parental"]))
    sub = d[d.tp_index.isin(shared)].copy()
    print(f"shared timepoints: {shared}  "
          f"(hTERT n={sum(sub.arm=='hTERT')}, parental n={sum(sub.arm=='parental')}); "
          f"parental TP8-10 (PD 50/52/53) are EXCLUDED - see D1")
    drows = []
    for clock in MODELS:
        y = sub[clock].values
        tp = sub.tp_index.values.astype(float)
        arm = (sub.arm == "hTERT").values.astype(float)
        beta, se, p, dof = ols(np.column_stack([np.ones_like(y), tp, arm, tp * arm]), y)
        # Permutation null on the arm label. Permuting SAMPLE labels would break
        # the replicate structure and count each library as independent, so the
        # label is permuted at the level of the (arm, timepoint) culture and
        # applied to all three of its replicates together.
        obs = beta[2]
        cell = (sub.arm + "_" + sub.timepoint).values
        cells = np.unique(cell)
        cell_arm = np.array([arm[cell == c][0] for c in cells])
        cnt = 0
        for _ in range(NPERM):
            perm_map = dict(zip(cells, rng.permutation(cell_arm)))
            pa = np.array([perm_map[c] for c in cell])
            if pa.std() == 0:
                continue
            b2, *_ = ols(np.column_stack([np.ones_like(y), tp, pa, tp * pa]), y)
            cnt += abs(b2[2]) >= abs(obs)
        drows.append(dict(clock=clock, arm_effect=beta[2], se_arm=se[2],
                          p_arm_effect=p[2], p_arm_perm=(1 + cnt) / (NPERM + 1),
                          slope_difference=beta[3], p_slope_difference=p[3], dof=dof))
    dres = pd.DataFrame(drows)
    print(dres.round(4).to_string(index=False))
    print("\n  arm_effect: hTERT minus parental with timepoint held constant."
          "\n              NEGATIVE means the immortalised arm scores LOWER.")

    # --- E: put our own confounded estimate next to it ---------------------
    print("\n" + "=" * 78)
    print("E  our own immortalisation estimate vs this one, same units")
    print("=" * 78)
    try:
        ours = pd.read_csv(os.path.join(rerun_dir, "immortalisation_contrasts.csv"))
        o = ours[(ours.test == "A_naive_within_condition_CONFOUNDED") &
                 (ours.stratum == "Proliferating")]
        for _, r in o.iterrows():
            key = "chrono_scaled_diff" if r.model == "scaled_diff" else "chrono_yugene_diff"
            here = res[(res.contrast == "A_dividing_vs_dividing_bytimepoint") &
                       (res.clock == key)].iloc[0]
            print(f"  {r.model:<12} 34-study meta-analysis: {r['diff']:+8.2f} "
                  f"(p={r.p:.2g}, {int(r.n_studies_imm)} vs {int(int(r.n_studies_prim))} studies, "
                  f"study-confounded)")
            print(f"  {'':<12} GSE175533 one lab   : {here['diff']:+8.2f} "
                  f"(p={here.p:.2g}, floor {here.p_floor:.4f}, one strain, "
                  f"study held constant, timepoint-level)")
            print(f"  {'':<12} -> {'SAME' if np.sign(r['diff']) == np.sign(here['diff']) else 'OPPOSITE'} sign\n")
    except FileNotFoundError:
        print("  immortalisation_contrasts.csv not found; run meta_analysis/16 first")

    out = os.path.join(D, "gse175533_contrasts.csv")
    pd.concat([res, d1.assign(contrast="D1_within_arm_slope"),
           dres.assign(contrast="D2_arm_effect")]).to_csv(out, index=False)
    print(f"\nSaved -> {os.path.join(D, 'gse175533_tage.csv')}\n         {out}")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
