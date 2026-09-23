#!/usr/bin/env python
"""25_yardstick_calibration.py

Two checks on the matched-null yardstick, answering the question that multiplicity
correction would otherwise be invoked for: of the comparisons that beat their null,
which are real?

A  CALIBRATION. Are the empirical p-values valid in the first place? The gene-to-set
   assignment is shuffled, which destroys any real set structure while leaving the
   contribution vector, the coefficient strata and the matching untouched. Under that
   null the p-values should be uniform, so about 5% of sets should fall below 0.05. If
   many more do, the null is mis-specified; if far fewer, it is conservative.

   The number of draws does not enter this. Draws control how precisely each p is
   estimated, not the rate at which a valid 5% test fires under the null - 20,000
   draws give a precise p, not a conservative one.

B  RECURRENCE. This is what identifies real results without any adjustment. Each
   group (condition, or cell type x timepoint) is an independent test of the same set,
   so under the null the number of groups in which a set beats its own null is
   Binomial(n_groups, 0.05). A set clearing 6 of 9 groups is expected 1.2e-6 of the
   time, i.e. 0.0001 times among 50 sets - which is a far stronger statement than any
   single comparison at p < 0.001, and needs no correction because the multiplicity is
   built into the statistic.

   This is why the reported results lean on sets that recur rather than on the most
   extreme single comparison.

Output: rerun_outputs/yardstick_calibration.csv
        rerun_outputs/yardstick_recurrence.csv

Usage: 25_yardstick_calibration.py <rerun_dir> <model_dir> [n_reps] [n_draws]
"""
import sys, warnings
import joblib, numpy as np, pandas as pd
from scipy.stats import binom
warnings.filterwarnings("ignore")

MORT = "EN_Mortality_Multispecies_Multitissue_scaleddiff.pkl"
SEED = 7
CONDS = {"Contact_inhibited CQ": "CICQ", "Serum_starved CQ": "SSCQ",
         "Replicative CS": "RS", "Stress-induced CS": "SIPS",
         "Oncogene-induced CS": "OIS"}


def main(rerun_dir, model_dir, n_reps=10, n_draws=4000):
    PT = f"{rerun_dir}/partial_tage"
    rng = np.random.default_rng(SEED)
    m = joblib.load(f"{model_dir}/{MORT}")
    imp = m.named_steps["imputation"]
    if not hasattr(imp, "_fill_dtype"):
        imp._fill_dtype = imp.statistics_.dtype
    feats = list(map(str, m.feature_names_in_))
    idx = {g: i for i, g in enumerate(feats)}
    coef = m.named_steps["estimator"].coef_; aco = np.abs(coef)
    e = pd.read_csv(f"{PT}/meta_scaled_diff.csv"); sid = e["sample_id"].values
    e = e.drop(columns=["sample_id"]); e.columns = e.columns.map(str)
    for g in [g for g in feats if g not in e.columns]:
        e[g] = np.nan
    Z = m.named_steps["scaler"].transform(imp.transform(e.loc[:, feats]))
    C = Z * coef[None, :]
    grp = pd.read_csv(f"{PT}/meta_groups.csv").set_index("sample_id")["group"]
    md = pd.read_csv(f"{rerun_dir}/sample_metadata_RERUN.csv").set_index("external_id")
    g = np.array([grp.get(s) for s in sid])
    st = np.array([md.loc[s, "study"] for s in sid])
    pw = pd.read_csv(f"{PT}/hallmark_pathway_mouse_ids.csv")
    sets = {n: np.array([idx[x] for x in s.mouse_gene_id.astype(str) if x in idx])
            for n, s in pw.groupby("pathway")}
    sets = {k: v for k, v in sets.items() if len(v)}
    q = np.quantile(aco, np.linspace(0, 1, 11)[1:-1]); strat = np.searchsorted(q, aco)
    by_strat = {s: np.where(strat == s)[0] for s in np.unique(strat)}

    def within_study_d(cond):
        acc = np.zeros(C.shape[1]); w = 0.0
        for s in np.unique(st):
            ti = np.where((st == s) & (g == cond))[0]
            ci = np.where((st == s) & (g == "Proliferating"))[0]
            if len(ti) and len(ci):
                ww = len(ti) * len(ci) / (len(ti) + len(ci))
                acc += ww * (C[ti].mean(0) - C[ci].mean(0)); w += ww
        return acc / w

    def p_emp(d, cols, B):
        counts = {s: int((strat[cols] == s).sum()) for s in np.unique(strat[cols])}
        draws = np.zeros(B)
        for s_, k in counts.items():
            draws += rng.choice(d[by_strat[s_]], size=(B, k), replace=True).sum(1)
        o = d[cols].sum()
        return (1 + np.sum(np.abs(draws - draws.mean()) >= abs(o - draws.mean()))) / (B + 1)

    # ---- A calibration ----
    rows = []
    print(f"A  CALIBRATION: gene-to-set assignment shuffled, {n_reps} replicates,"
          f" {n_draws} draws per comparison")
    for cond_value, lab in CONDS.items():
        d = within_study_d(cond_value)
        real = np.array([p_emp(d, c, n_draws) for c in sets.values()])
        shuf = []
        for _ in range(n_reps):
            dp = d[rng.permutation(len(d))]
            ps = np.array([p_emp(dp, c, n_draws) for c in sets.values()])
            shuf.append((ps < 0.05).mean())
        rows.append(dict(condition=lab, n_sets=len(sets),
                         frac_p05_real=float((real < 0.05).mean()),
                         frac_p05_shuffled_mean=float(np.mean(shuf)),
                         frac_p05_shuffled_min=float(np.min(shuf)),
                         frac_p05_shuffled_max=float(np.max(shuf)),
                         n_reps=n_reps, n_draws=n_draws))
        print(f"   {lab:<5} real {100*(real<0.05).mean():>5.1f}%   shuffled"
              f" {100*np.mean(shuf):>5.1f}% (range {100*np.min(shuf):.1f}-{100*np.max(shuf):.1f}%)"
              f"   nominal 5%")
    cal = pd.DataFrame(rows)
    cal.to_csv(f"{rerun_dir}/yardstick_calibration.csv", index=False)

    # ---- B recurrence ----
    print("\nB  RECURRENCE: groups in which a set beats its own null, vs Binomial(n, 0.05)")
    out = []
    for f, arm, ngrp in [("pathway_specificity_yardstick_mortality.csv", "cross_sectional", 5),
                         ("pathway_specificity_yardstick_mortality_temporal.csv", "temporal", 9)]:
        y = pd.read_csv(f"{rerun_dir}/{f}")
        hits = y[y.p_emp < 0.05].pathway.value_counts()
        for pwy, k in hits.items():
            p1 = float(binom.sf(k - 1, ngrp, 0.05))
            out.append(dict(arm=arm, pathway=pwy, n_groups=ngrp, n_hit=int(k),
                            p_recurrence=p1, expected_among_50_sets=50 * p1))
    rec = pd.DataFrame(out).sort_values(["arm", "p_recurrence"])
    rec.to_csv(f"{rerun_dir}/yardstick_recurrence.csv", index=False)
    for arm in rec.arm.unique():
        s = rec[(rec.arm == arm) & (rec.expected_among_50_sets < 0.05)]
        print(f"   {arm}: sets whose recurrence is expected < 0.05 times among 50")
        for _, r in s.iterrows():
            print(f"     {r.pathway.replace('HALLMARK ', ''):<32} {r.n_hit}/{r.n_groups}"
                  f"   P = {r.p_recurrence:.1e}   expected {r.expected_among_50_sets:.4f}")
    print(f"\nSaved -> yardstick_calibration.csv, yardstick_recurrence.csv")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2],
         int(sys.argv[3]) if len(sys.argv) > 3 else 10,
         int(sys.argv[4]) if len(sys.argv) > 4 else 4000)
