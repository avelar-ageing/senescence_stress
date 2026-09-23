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

Output: rerun_outputs/mortality_temporal_by_celltype.csv (per-sample predictions only;
testing moved to 07_tage_temporal_tests.R so all three clocks share one BH family)
Usage: 08_mortality_temporal.py <rerun_dir> <model_dir>
"""
import sys, warnings, itertools
import joblib, numpy as np, pandas as pd
from scipy.stats import mannwhitneyu, kruskal

# Keratinocyte batch -- mirrors temporal_analysis/R_keratinocyte_batch.R and the
# fix_batch = TRUE list in 02_run_time_analysis.R (that cell type was processed
# by two researchers). Balanced 3/3 across timepoints, so the timepoint
# estimates are unchanged; centring only removes the batch offset from the
# within-timepoint spread. Fibroblasts and melanocytes untouched.
KERATINOCYTE_BATCH_1 = {
    "ERR1805235", "ERR1805236", "ERR1805238", "ERR1805230", "ERR1805231", "ERR1805224",
    "ERR1805223", "ERR1805222", "ERR1805239", "ERR1805240", "ERR1805241", "ERR1805229",
}


def batch_centre(pred, sample_ids, cell_type):
    """Subtract each batch's mean, keeping the grand mean. Keratinocytes only."""
    if cell_type != "Keratinocyte":
        return pred
    b = np.array([s in KERATINOCYTE_BATCH_1 for s in sample_ids])
    if b.all() or not b.any():
        raise ValueError("keratinocyte batch labels did not split the samples")
    out = pred.astype(float).copy()
    grand = out.mean()
    out[b] = out[b] - out[b].mean() + grand
    out[~b] = out[~b] - out[~b].mean() + grand
    return out
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
        pred_raw = m.named_steps["estimator"].intercept_ + Z @ coef
        # mortality_tAge is the batch-centred value used by every test and
        # figure; mortality_tAge_raw keeps the uncorrected prediction so the
        # correction stays reversible and auditable.
        pred = batch_centre(pred_raw, sid, ct)
        tp = md.loc[sid, "time_after_treatment"].values
        per_sample.append(pd.DataFrame(dict(external_id=sid, cell_type=ct,
                                            timepoint=tp, mortality_tAge=pred,
                                            mortality_tAge_raw=pred_raw)))
    ps = pd.concat(per_sample, ignore_index=True)
    ps.to_csv(f"{rerun_dir}/mortality_temporal_by_celltype.csv", index=False)

    print("\nTesting is done by temporal_analysis/07_tage_temporal_tests.R,\n"
          "which pools all three clocks into one BH family.")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
