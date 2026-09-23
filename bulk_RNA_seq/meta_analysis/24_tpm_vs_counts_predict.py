#!/usr/bin/env python
"""24_tpm_vs_counts_predict.py

Second half of the TPM-versus-counts control (meta_analysis/23). Predicts all three
clocks on the count-derived and TPM-derived matrices of OUR OWN data and reports
how much the substitution moves the numbers we report.

WHAT WOULD INVALIDATE THE GSE175533 REPLICATION. Not a shift in absolute tAge - all
reported quantities are differences between groups, so a common offset is harmless.
What matters is whether the CONTRASTS change: if the immortalised-versus-primary
difference or the condition effects moved appreciably between the two inputs, then
the sign reversal found in GSE175533 could be an artefact of feeding TPM to a
pipeline built for counts. So the comparison reported here is on the contrasts, not
on the per-sample values, with the per-sample agreement given only as context.

Usage: 24_tpm_vs_counts_predict.py <rerun_dir> <model_dir>
Output: <rerun_dir>/gse175533/tpm_vs_counts_control.csv
"""
import os
import sys
import warnings

import joblib
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, pearsonr, spearmanr

warnings.filterwarnings("ignore")

MAX_LIFESPAN = 122.5
MODELS = {
    "mortality":          ("EN_Mortality_Multispecies_Multitissue_scaleddiff.pkl", "scaled_diff", 1.0),
    "chrono_scaled_diff": ("EN_Chronoage_Multispecies_Multitissue_scaleddiff.pkl", "scaled_diff", MAX_LIFESPAN),
    "chrono_yugene_diff": ("EN_Chronoage_Multispecies_Multitissue_yugenediff.pkl", "yugene_diff", MAX_LIFESPAN),
}
CONDS = ["CICQ", "SSCQ", "RS", "SIPS", "OIS"]
SUBSTATE = {"Contact_inhibited CQ": "CICQ", "Serum_starved CQ": "SSCQ",
            "Replicative CS": "RS", "Stress-induced CS": "SIPS",
            "Oncogene-induced CS": "OIS"}


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
    Z = m.named_steps["scaler"].transform(
        m.named_steps["imputation"].transform(e.loc[:, feats]))
    est = m.named_steps["estimator"]
    return pd.Series((est.intercept_ + Z @ est.coef_) * scale, index=sid)


def main(rerun_dir, model_dir):
    D = os.path.join(rerun_dir, "gse175533")
    meta = pd.read_csv(os.path.join(rerun_dir, "sample_metadata_RERUN.csv"))
    ann = pd.read_csv(os.path.join(rerun_dir,
                                   "immortalisation_annotation_corrected.csv")).set_index("external_id")
    seeds = sorted(int(f.split("tpm")[1].split("_")[0])
                   for f in os.listdir(D) if f.startswith("control_tpm")
                   and f.endswith("scaled_diff.csv"))
    print(f"TPM length draws available: {seeds}")

    loaded = {}
    for name, (pkl, variant, scale) in MODELS.items():
        m = joblib.load(os.path.join(model_dir, pkl))
        _patch(m.named_steps["imputation"])
        loaded[name] = (m, variant, scale)

    def all_clocks(tag):
        out = {}
        for name, (m, variant, scale) in loaded.items():
            out[name] = predict(m, os.path.join(D, f"control_{tag}_{variant}.csv"), scale)
        return pd.DataFrame(out)

    base = all_clocks("counts")
    info = meta.set_index("external_id").reindex(base.index)
    cond = info.cell_substate.map(SUBSTATE).fillna("Proliferating")
    imm = ann.reindex(base.index).immortalised

    def contrasts(P):
        r = {}
        for c in CONDS:
            for clock in MODELS:
                # within-study, the estimator used throughout 2.1.5
                num = den = 0.0
                for st, g in P.groupby(info.study.values):
                    cc = cond.loc[g.index]
                    x = g[clock][cc == c].values
                    y = g[clock][cc == "Proliferating"].values
                    if len(x) and len(y):
                        w = len(x) * len(y) / (len(x) + len(y))
                        num += w * (x.mean() - y.mean())
                        den += w
                r[(f"within_study_{c}", clock)] = num / den if den else np.nan
        for clock in MODELS:
            p = cond == "Proliferating"
            a = P[clock][p & (imm == "yes")]
            b = P[clock][p & (imm == "no")]
            r[("immortalised_vs_primary_proliferating", clock)] = a.median() - b.median()
        return r

    cb = contrasts(base)
    rows = []
    for s in seeds:
        alt = all_clocks(f"tpm{s}")
        common = base.index.intersection(alt.index)
        for clock in MODELS:
            rows.append(dict(seed=s, clock=clock, level="per_sample",
                             quantity="tAge",
                             counts=np.nan, tpm=np.nan,
                             pearson=pearsonr(base.loc[common, clock],
                                              alt.loc[common, clock])[0],
                             spearman=spearmanr(base.loc[common, clock],
                                                alt.loc[common, clock])[0],
                             abs_diff_median=float(np.median(np.abs(
                                 base.loc[common, clock] - alt.loc[common, clock])))))
        ca = contrasts(alt)
        for (q, clock), v in cb.items():
            rows.append(dict(seed=s, clock=clock, level="contrast", quantity=q,
                             counts=v, tpm=ca[(q, clock)],
                             pearson=np.nan, spearman=np.nan,
                             abs_diff_median=abs(v - ca[(q, clock)])))
    R = pd.DataFrame(rows)

    print("\n" + "=" * 78)
    print("per-sample agreement between counts- and TPM-derived tAge")
    print("=" * 78)
    print(R[R.level == "per_sample"].groupby("clock")[
        ["pearson", "spearman", "abs_diff_median"]].agg(["min", "max"]).round(4).to_string())

    print("\n" + "=" * 78)
    print("the quantities we actually report: counts vs TPM")
    print("=" * 78)
    c = R[R.level == "contrast"]
    piv = c.pivot_table(index=["quantity", "clock"], values=["counts", "tpm"],
                        aggfunc={"counts": "first", "tpm": ["min", "max"]})
    piv.columns = ["counts", "tpm_max", "tpm_min"]
    piv["worst_shift"] = (piv[["tpm_min", "tpm_max"]].sub(piv.counts, axis=0)
                          .abs().max(axis=1))
    piv["sign_stable"] = (np.sign(piv.tpm_min) == np.sign(piv.counts)) & \
                         (np.sign(piv.tpm_max) == np.sign(piv.counts))
    print(piv[["counts", "tpm_min", "tpm_max", "worst_shift", "sign_stable"]]
          .round(3).to_string())

    worst = piv.worst_shift.max()
    print(f"\nlargest movement of any reported contrast across {len(seeds)} length draws: "
          f"{worst:.3f} units")
    print(f"sign preserved for every contrast: {bool(piv.sign_stable.all())}")
    print("\nInterpretation: the gene-length term cancels in the per-gene z-score, so"
          "\nthe only route by which TPM could matter is the detection filter changing"
          "\nwhich genes are kept (see 23: Jaccard ~0.97). If the shifts above are small"
          "\nrelative to the effects reported, feeding TPM to this pipeline is safe.")

    out = os.path.join(D, "tpm_vs_counts_control.csv")
    R.to_csv(out, index=False)
    print(f"\nSaved -> {out}")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
