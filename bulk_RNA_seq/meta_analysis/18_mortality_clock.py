#!/usr/bin/env python
"""18_mortality_clock.py

Computes MORTALITY tAge alongside the chronological tAge used elsewhere, and runs
the same within-study analysis on it.

WHY. Tyshkovskiy et al. publish two clocks. Everything in 2.1.5 so far uses the
chronological one. Their own hTERT experiment - transducing WI-38 fibroblasts and
culturing them against untransduced controls - is reported with the MORTALITY
clock, so that finding cannot be compared with ours until this clock is run. The
mortality clock is also a structurally different instrument, which matters for the
pathway question:

              features   non-zero   alpha    l1_ratio
  chronoage     10,487      1,839   0.001    0.20   (sparse elastic net)
  mortality     10,487     10,487   0.056    0.00   (pure RIDGE, fully dense)

Both share the same feature list; their coefficients correlate 0.74 over the 1,839
features non-zero in both. Because ridge shrinks correlated genes together rather
than selecting arbitrarily among them, and because every gene carries weight, the
two objections that limit set-level interpretation of the chronological clock -
25-57 effective genes per set, and semi-arbitrary allocation among co-expressed
genes - are much weaker here.

UNITS. No species adjustment is applied. The chronological clock is multiplied by
122.5 because it was trained on age divided by maximum lifespan, so the product is
years. The mortality target is not a lifespan fraction, so scaling it by a lifespan
would be meaningless. Values are therefore in the model's native units and are NOT
comparable in magnitude with the chronological values. Every quantity reported is a
difference from controls, so the choice of scale does not affect any test.

LIMITATION, stated up front: only a scaleddiff mortality clock exists. There is no
YuGene mortality model, so the two-model agreement criterion used throughout the
chronological analysis is unavailable, and the one model available is the
normalisation shown to be the more sensitive of the two to control composition.

Output: rerun_outputs/mortality_tage.csv           (per sample)
        rerun_outputs/mortality_within_study.csv   (condition effects + immortalisation)

Usage: 18_mortality_clock.py <rerun_dir> <model_dir>
"""
import sys, warnings, itertools
import joblib, numpy as np, pandas as pd
from scipy.stats import mannwhitneyu
warnings.filterwarnings("ignore")

MORT = "EN_Mortality_Multispecies_Multitissue_scaleddiff.pkl"
SEED, NPERM = 1, 20000
CONDS = {"Contact_inhibited CQ": "CICQ", "Serum_starved CQ": "SSCQ",
         "Replicative CS": "RS", "Stress-induced CS": "SIPS",
         "Oncogene-induced CS": "OIS"}


def _patch(imp):
    if not hasattr(imp, "_fill_dtype"):
        imp._fill_dtype = imp.statistics_.dtype if hasattr(imp, "statistics_") else np.float64


def predict(m, path):
    feats = list(map(str, m.feature_names_in_))
    e = pd.read_csv(path)
    sid = e["sample_id"].values
    e = e.drop(columns=["sample_id"]); e.columns = e.columns.map(str)
    for g in [g for g in feats if g not in e.columns]:
        e[g] = np.nan
    X = e.loc[:, feats]
    Z = m.named_steps["scaler"].transform(m.named_steps["imputation"].transform(X))
    coef = m.named_steps["estimator"].coef_
    return sid, m.named_steps["estimator"].intercept_ + Z @ coef, Z * coef[None, :]


def main(rerun_dir, model_dir):
    PT = f"{rerun_dir}/partial_tage"
    rng = np.random.default_rng(SEED)
    m = joblib.load(f"{model_dir}/{MORT}"); _patch(m.named_steps["imputation"])

    meta = pd.read_csv(f"{rerun_dir}/sample_metadata_RERUN.csv")
    ann = pd.read_csv(f"{rerun_dir}/immortalisation_annotation_corrected.csv").set_index("external_id")
    groups = pd.read_csv(f"{PT}/meta_groups.csv")

    sid, pred, C = predict(m, f"{PT}/meta_scaled_diff.csv")
    d = pd.DataFrame({"external_id": sid, "mortality_tAge": pred})
    d["study"] = meta.set_index("external_id").loc[d.external_id, "study"].values
    d["group"] = groups.set_index("sample_id").loc[d.external_id, "group"].values
    d["condition"] = d.group.map(CONDS).fillna("Proliferating")
    d["immortalised"] = ann.loc[d.external_id, "immortalised"].values
    d["cell_line"] = ann.loc[d.external_id, "cell_line_resolved"].values
    d.to_csv(f"{rerun_dir}/mortality_tage.csv", index=False)
    print(f"mortality tAge computed for {len(d)} samples "
          f"(median {d.mortality_tAge.median():.3f}, range "
          f"{d.mortality_tAge.min():.3f} to {d.mortality_tAge.max():.3f})")
    print(f"  Proliferating median {d.loc[d.condition=='Proliferating','mortality_tAge'].median():.3f}"
          f"  (intercept {m.named_steps['estimator'].intercept_:.3f})")

    rows = []
    # ---- within-study condition effects, same estimator as script 13 --------
    for cond in ["CICQ", "SSCQ", "RS", "SIPS", "OIS"]:
        s = d[d.condition.isin([cond, "Proliferating"])]
        blocks = []
        for st, g in s.groupby("study"):
            x = g.mortality_tAge[g.condition == cond].values
            y = g.mortality_tAge[g.condition == "Proliferating"].values
            if len(x) and len(y):
                blocks.append((x, y, len(x) * len(y) / (len(x) + len(y))))
        W = sum(b[2] for b in blocks)
        obs = sum(b[2] * (b[0].mean() - b[1].mean()) for b in blocks) / W
        perm = np.empty(NPERM)
        for i in range(NPERM):
            acc = 0.0
            for x, y, w in blocks:
                pool = np.concatenate([x, y]); rng.shuffle(pool)
                acc += w * (pool[:len(x)].mean() - pool[len(x):].mean())
            perm[i] = acc / W
        p = (1 + (np.abs(perm) >= abs(obs)).sum()) / (NPERM + 1)
        npos = sum((b[0].mean() - b[1].mean()) > 0 for b in blocks)
        rows.append(dict(test="condition_within_study", condition=cond,
                         n_studies=len(blocks), diff_within_study=obs,
                         diff_pooled=d.loc[d.condition==cond,"mortality_tAge"].median()
                                     - d.loc[d.condition=="Proliferating","mortality_tAge"].median(),
                         studies_positive=npos, p_perm=p))
    res = pd.DataFrame(rows)
    # raw p_perm is the reported statistic; see meta_analysis/13
    res["p_perm_adj_DEPRECATED"] = _bh(res.p_perm.values)
    print("\n=== within-study condition effects, MORTALITY clock ===")
    print(res.round(4).to_string(index=False))

    # ---- immortalisation, the comparison the paper made -------------------
    irows = []
    for cond in ["Proliferating", "CICQ", "SSCQ", "OIS"]:
        s = d[d.condition == cond]
        y = s.mortality_tAge[s.immortalised == "yes"]; n = s.mortality_tAge[s.immortalised == "no"]
        if len(y) < 2 or len(n) < 2: continue
        irows.append(dict(test="immortalised_vs_primary", stratum=cond,
                          n_imm=len(y), n_prim=len(n),
                          median_imm=y.median(), median_prim=n.median(),
                          diff=y.median() - n.median(),
                          p=mannwhitneyu(y, n).pvalue))
    imm = pd.DataFrame(irows); imm["p_adj"] = _bh(imm.p.values)
    print("\n=== immortalised vs primary, MORTALITY clock (study-confounded) ===")
    print(imm.round(4).to_string(index=False))

    # WI-38 specifically: the strain the paper used
    w = d[(d.cell_line == "WI-38")]
    print(f"\n=== WI-38, the strain used in the paper's hTERT experiment ===")
    print(w.groupby("condition").mortality_tAge.agg(["size", "median"]).round(3).to_string())

    # ---- temporal arm: does the mortality clock reproduce the trajectories? ----
    trows = []
    for ct in ["Fibroblast", "Keratinocyte", "Melanocyte"]:
        for tp in ["4_days", "10_days", "20_days"]:
            stem = f"{ct}_{tp}"
            s2, p2, _ = predict(m, f"{PT}/{stem}_scaled_diff.csv")
            g = pd.read_csv(f"{PT}/{stem}_groups.csv").set_index("sample_id")["group"]
            lab = np.array([g.get(x) for x in s2])
            a, b = p2[lab == tp], p2[lab == "none"]
            trows.append(dict(test="temporal_vs_baseline", cell_type=ct, timepoint=tp,
                              n_t=len(a), n_c=len(b),
                              diff=np.median(a) - np.median(b),
                              p=mannwhitneyu(a, b).pvalue))
    tmp = pd.DataFrame(trows); tmp["p_adj"] = _bh(tmp.p.values)
    print("\n=== temporal arm, MORTALITY clock (vs each cell type's own baseline) ===")
    print(tmp.round(4).to_string(index=False))
    for ct in ["Fibroblast", "Keratinocyte", "Melanocyte"]:
        v = tmp.loc[tmp.cell_type == ct, "diff"].values
        shape = ("monotonic rise" if v[2] > v[1] > v[0] else
                 "monotonic decline" if v[2] < v[1] < v[0] else
                 "early peak then partial reversal")
        print(f"  {ct:<13} {np.round(v,3)}  -> {shape}")

    pd.concat([res, imm, tmp]).to_csv(f"{rerun_dir}/mortality_within_study.csv", index=False)
    print(f"\nSaved -> {rerun_dir}/mortality_tage.csv and mortality_within_study.csv")


def _bh(p):
    p = np.asarray(p, float); n = len(p); o = np.argsort(p); adj = np.empty(n)
    adj[o] = np.minimum.accumulate((p[o] * n / (np.arange(n) + 1))[::-1])[::-1]
    return np.clip(adj, 0, 1)


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
