#!/usr/bin/env python
"""08_mortality_temporal.py

Full mortality-clock analysis of the time course, matching the statistical
treatment the chronological arm received in 07_tage_temporal_pairwise.R.

WHY. The mortality clock was initially run only as vs-baseline corroboration, which
is not enough to support trajectory claims: statements like "peaks at 4 days then
partially reverses" are between-timepoint claims and need between-timepoint tests.

PREPROCESSING UNIT. Values come from the per-CELL-TYPE runs, in which all four
timepoints of one cell type are preprocessed together, so timepoints are directly
comparable within a cell type. The per-timepoint runs used elsewhere cannot support
between-timepoint contrasts, because each is separately preprocessed and control-
subtracted, so 4-day and 10-day values from those runs are on different scales.

FAMILY. All 6 timepoint pairs x 3 cell types = 18 tests, one BH family, which
subsumes the vs-baseline comparisons (none-vs-4/10/20 are 3 of the 6 pairs). A
Kruskal-Wallis omnibus per cell type is reported alongside. Only one mortality
normalisation exists, so there is no second model to corroborate against.

Output: rerun_outputs/mortality_temporal_pairwise.csv
Usage: 08_mortality_temporal.py <rerun_dir> <model_dir>
"""
import sys, warnings, itertools
import joblib, numpy as np, pandas as pd
from scipy.stats import mannwhitneyu, kruskal
warnings.filterwarnings("ignore")

MORT = "EN_Mortality_Multispecies_Multitissue_scaleddiff.pkl"
CELLS = ["Fibroblast", "Keratinocyte", "Melanocyte"]
ORDER = ["none", "4_days", "10_days", "20_days"]


def bh(p):
    p = np.asarray(p, float); n = len(p); o = np.argsort(p); a = np.empty(n)
    a[o] = np.minimum.accumulate((p[o] * n / (np.arange(n) + 1))[::-1])[::-1]
    return np.clip(a, 0, 1)


def main(rerun_dir, model_dir):
    PT = f"{rerun_dir}/partial_tage"
    m = joblib.load(f"{model_dir}/{MORT}")
    imp = m.named_steps["imputation"]
    if not hasattr(imp, "_fill_dtype"):
        imp._fill_dtype = imp.statistics_.dtype
    feats = list(map(str, m.feature_names_in_))
    coef = m.named_steps["estimator"].coef_
    md = pd.read_csv(f"{rerun_dir}/tage_temporal_by_celltype.csv").set_index("external_id")

    per_sample, rows, omni = [], [], []
    for ct in CELLS:
        e = pd.read_csv(f"{PT}/{ct}_scaled_diff.csv")
        sid = e["sample_id"].values
        e = e.drop(columns=["sample_id"]); e.columns = e.columns.map(str)
        for g in [g for g in feats if g not in e.columns]:
            e[g] = np.nan
        Z = m.named_steps["scaler"].transform(imp.transform(e.loc[:, feats]))
        pred = m.named_steps["estimator"].intercept_ + Z @ coef
        tp = md.loc[sid, "time_after_treatment"].values
        per_sample.append(pd.DataFrame(dict(external_id=sid, cell_type=ct,
                                            timepoint=tp, mortality_tAge=pred)))
        groups = [pred[tp == t] for t in ORDER]
        omni.append(dict(test="kruskal", cell_type=ct,
                         statistic=kruskal(*groups).statistic,
                         p=kruskal(*groups).pvalue))
        for a, b in itertools.combinations(ORDER, 2):
            x, y = pred[tp == b], pred[tp == a]
            rows.append(dict(test="pairwise_timepoints", cell_type=ct,
                             group_1=a, group_2=b, n_1=len(y), n_2=len(x),
                             median_1=np.median(y), median_2=np.median(x),
                             diff=np.median(x) - np.median(y),
                             p=mannwhitneyu(x, y).pvalue))
    ps = pd.concat(per_sample, ignore_index=True)
    ps.to_csv(f"{rerun_dir}/mortality_temporal_by_celltype.csv", index=False)

    r = pd.DataFrame(rows); r["p_adj"] = bh(r.p.values)
    o = pd.DataFrame(omni); o["p_adj"] = bh(o.p.values)

    print("=== median mortality tAge per cell type x timepoint (per-cell-type runs) ===")
    print(ps.pivot_table(index="cell_type", columns="timepoint",
                         values="mortality_tAge", aggfunc="median")[ORDER].round(3).to_string())
    print("\n=== Kruskal-Wallis across all four timepoints ===")
    print(o.round(5).to_string(index=False))
    print("\n=== all 6 timepoint pairs x 3 cell types (18-test BH family) ===")
    for ct in CELLS:
        print(f"\n-- {ct} --")
        print(r[r.cell_type == ct][["group_1", "group_2", "diff", "p", "p_adj"]]
              .round(4).to_string(index=False))
    print("\n=== trajectory summary (vs baseline, and the late contrasts) ===")
    for ct in CELLS:
        s = r[r.cell_type == ct].set_index(["group_1", "group_2"])
        vb = [s.loc[("none", t), "diff"] for t in ORDER[1:]]
        sig = ["*" if s.loc[("none", t), "p_adj"] < 0.05 else "" for t in ORDER[1:]]
        late = s.loc[("4_days", "20_days")]
        mid = s.loc[("4_days", "10_days")]
        print(f"  {ct:<13} vs baseline {vb[0]:+.3f}{sig[0]:<1} {vb[1]:+.3f}{sig[1]:<1} "
              f"{vb[2]:+.3f}{sig[2]:<1} | 4->10 {mid['diff']:+.3f} "
              f"(BH {mid['p_adj']:.3f}) | 4->20 {late['diff']:+.3f} (BH {late['p_adj']:.3f})")
    pd.concat([r, o]).to_csv(f"{rerun_dir}/mortality_temporal_pairwise.csv", index=False)
    print(f"\nSaved -> {rerun_dir}/mortality_temporal_pairwise.csv and "
          f"mortality_temporal_by_celltype.csv")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
