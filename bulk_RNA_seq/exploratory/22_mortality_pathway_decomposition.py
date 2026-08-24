#!/usr/bin/env python
"""22_mortality_pathway_decomposition.py

Runs the gene-set decomposition on the MORTALITY clock for BOTH arms, so the
pathway-level sections can be reported on that clock alone.

WHY MORTALITY ONLY FOR SET-LEVEL WORK. The choice is about decomposability, not
about preferring mortality as an endpoint. A set's contribution is a sum of
coef_i x z_i over its genes, so the decomposition is only informative if the
model spreads weight across the genes of a set:

               features   non-zero   l1_ratio   consequence for a 150-gene set
  chronoage      10,487      1,839      0.20    carried by 25-57 genes; elastic
                                                net picks semi-arbitrarily among
                                                co-expressed genes
  mortality      10,487     10,487      0.00    every gene contributes; ridge
                                                shrinks correlated genes together

Empirically this is what separates them: against weight-matched random sets the
chronological clock resolved 2 of 85 interpretable tests with a median |z| of
0.22-0.26, the mortality clock 6 of 85 with a median |z| of 1.01. The whole-
transcriptome sections use both clocks; the set-level sections use mortality.

THE COST, stated plainly: only a scaled-difference mortality model exists, so
set-level results cannot be corroborated across normalisations the way the
whole-transcriptome results are. That is the trade for decomposability.

WHAT THIS COMPUTES
  per-sample set scores for every group in both arms, so group differences use the
  same Wilcoxon tests as the chronological pipeline;
  within-study combination for the meta-analysis (per meta_analysis/13), and
  vs-own-baseline for the temporal arm, which is a single study;
  the weight-share yardstick with a matched null for the temporal groups (the
  meta-analysis yardstick is already produced by exploratory/20 --mortality).

Output: rerun_outputs/mortality_partial_tage_ALL.csv
        rerun_outputs/pathway_specificity_yardstick_mortality_temporal.csv

Usage: 22_mortality_pathway_decomposition.py <rerun_dir> <model_dir> [n_null]
"""
import sys, warnings
import joblib, numpy as np, pandas as pd
from scipy.stats import mannwhitneyu
warnings.filterwarnings("ignore")

MORT = "EN_Mortality_Multispecies_Multitissue_scaleddiff.pkl"
SEED = 1
CONDS = {"Contact_inhibited CQ": "CICQ", "Serum_starved CQ": "SSCQ",
         "Replicative CS": "RS", "Stress-induced CS": "SIPS",
         "Oncogene-induced CS": "OIS"}
CELLS = ["Fibroblast", "Keratinocyte", "Melanocyte"]
TPS = ["4_days", "10_days", "20_days"]


def _patch(imp):
    if not hasattr(imp, "_fill_dtype"):
        imp._fill_dtype = imp.statistics_.dtype if hasattr(imp, "statistics_") else np.float64


def contribs(m, path, feats):
    e = pd.read_csv(path); sid = e["sample_id"].values
    e = e.drop(columns=["sample_id"]); e.columns = e.columns.map(str)
    for g in [g for g in feats if g not in e.columns]:
        e[g] = np.nan
    Z = m.named_steps["scaler"].transform(m.named_steps["imputation"].transform(e.loc[:, feats]))
    return sid, Z * m.named_steps["estimator"].coef_[None, :]


def bh(p):
    p = np.asarray(p, float); n = len(p); o = np.argsort(p); a = np.empty(n)
    a[o] = np.minimum.accumulate((p[o] * n / (np.arange(n) + 1))[::-1])[::-1]
    return np.clip(a, 0, 1)


def main(rerun_dir, model_dir, n_null=20000):
    PT = f"{rerun_dir}/partial_tage"
    rng = np.random.default_rng(SEED)
    m = joblib.load(f"{model_dir}/{MORT}"); _patch(m.named_steps["imputation"])
    feats = list(map(str, m.feature_names_in_))
    idx = {g: i for i, g in enumerate(feats)}
    coef = m.named_steps["estimator"].coef_; aco = np.abs(coef); wtot = aco.sum()
    pw = pd.read_csv(f"{PT}/hallmark_pathway_mouse_ids.csv")
    sets = {n: np.array([idx[x] for x in s.mouse_gene_id.astype(str) if x in idx])
            for n, s in pw.groupby("pathway")}
    sets = {k: v for k, v in sets.items() if len(v)}
    meta_md = pd.read_csv(f"{rerun_dir}/sample_metadata_RERUN.csv").set_index("external_id")

    # coefficient-decile strata for the matched null (no zero stratum: ridge)
    q = np.quantile(aco, np.linspace(0, 1, 11)[1:-1])
    strat = np.searchsorted(q, aco)
    by_strat = {s: np.where(strat == s)[0] for s in np.unique(strat)}

    rows, yrows = [], []

    # ---------- meta-analysis: within study ----------
    sid, C = contribs(m, f"{PT}/meta_scaled_diff.csv", feats)
    grp = pd.read_csv(f"{PT}/meta_groups.csv").set_index("sample_id")["group"]
    g = np.array([grp.get(s) for s in sid])
    study = np.array([meta_md.loc[s, "study"] for s in sid])
    S = np.column_stack([C[:, ix].sum(axis=1) for ix in sets.values()])
    names = list(sets)
    for cond_value, lab in CONDS.items():
        blocks = []
        for st in np.unique(study):
            ti = np.where((study == st) & (g == cond_value))[0]
            ci = np.where((study == st) & (g == "Proliferating"))[0]
            if len(ti) and len(ci):
                blocks.append((ti, ci, len(ti) * len(ci) / (len(ti) + len(ci))))
        W = sum(b[2] for b in blocks)
        for k, nm in enumerate(names):
            eff = sum(b[2] * (S[b[0], k].mean() - S[b[1], k].mean()) for b in blocks) / W
            x = S[g == cond_value, k]; y = S[g == "Proliferating", k]
            rows.append(dict(analysis="meta_analysis", label=lab, pathway=nm,
                             n_studies=len(blocks), contrib_diff_within_study=eff,
                             wilcox_p=mannwhitneyu(x, y).pvalue))

    # ---------- temporal: per cell type x timepoint, and pooled ----------
    for ct in CELLS:
        for tp in TPS + ["pooled"]:
            stem = f"{ct}_{tp}" if tp != "pooled" else ct
            sid2, C2 = contribs(m, f"{PT}/{stem}_scaled_diff.csv", feats)
            gg = pd.read_csv(f"{PT}/{stem}_groups.csv").set_index("sample_id")["group"]
            lab2 = np.array([gg.get(s) for s in sid2])
            test_lab = tp if tp != "pooled" else "irradiated"
            a = lab2 == test_lab; b = lab2 == "none"
            S2 = np.column_stack([C2[:, ix].sum(axis=1) for ix in sets.values()])
            d = C2[a].mean(axis=0) - C2[b].mean(axis=0)      # per-gene difference
            total = d.sum()
            for k, nm in enumerate(names):
                obs = d[sets[nm]].sum()
                rows.append(dict(
                    analysis="temporal_pooled" if tp == "pooled" else "temporal_bytimepoint",
                    label=stem, cell_type=ct, timepoint=None if tp == "pooled" else tp,
                    pathway=nm, contrib_diff_within_study=obs,
                    wilcox_p=mannwhitneyu(S2[a, k], S2[b, k]).pvalue))
                if tp == "pooled":
                    continue
                cols = sets[nm]
                share = aco[cols].sum() / wtot
                counts = {s: int((strat[cols] == s).sum()) for s in np.unique(strat[cols])}
                draws = np.zeros(n_null)
                for st_id, kk in counts.items():
                    draws += rng.choice(d[by_strat[st_id]], size=(n_null, kk),
                                        replace=True).sum(axis=1)
                p = (1 + np.sum(np.abs(draws - draws.mean()) >= abs(obs - draws.mean()))) / (n_null + 1)
                yrows.append(dict(analysis="temporal_bytimepoint", label=stem,
                                  cell_type=ct, timepoint=tp, pathway=nm,
                                  observed=obs, expected_from_weight=share * total,
                                  total_shift=total, z=(obs - draws.mean()) / draws.std(),
                                  p_emp=p))
            print(f"  {stem} done (total shift {total:+.4f})")

    out = pd.DataFrame(rows)
    out["p_adj"] = out.groupby("analysis").wilcox_p.transform(lambda x: bh(x.values))
    out.to_csv(f"{rerun_dir}/mortality_partial_tage_ALL.csv", index=False)

    y = pd.DataFrame(yrows)
    rep = pd.read_csv(f"{rerun_dir}/pathway_representation_mortality.csv")
    interp = set(rep.loc[rep.tier == "INTERPRETABLE", "pathway"])
    y["interpretable"] = y.pathway.isin(interp)
    # NOT the reported statistic - see header. Kept for provenance only.
    y["p_emp_adj_DEPRECATED"] = np.nan
    sel = y.interpretable
    y.loc[sel, "p_emp_adj_DEPRECATED"] = bh(y.loc[sel, "p_emp"].values)
    y["p_floor"] = 1.0 / n_null
    print(f"\n  temporal yardstick: {len(y)} comparisons, floor {1/n_null:.1e}")
    for thr in (0.05, 0.01, 0.001):
        print(f"    raw p < {thr:<6}: {(y.p_emp < thr).sum():>3}"
              f"   expected by chance {len(y)*thr:.0f}")
    rec = y[y.p_emp < 0.05].pathway.value_counts()
    print(f"    sets passing in 4+ of the 9 groups: "
          f"{ {k.replace('HALLMARK ',''): int(v) for k, v in rec[rec >= 4].items()} }")
    y.to_csv(f"{rerun_dir}/pathway_specificity_yardstick_mortality_temporal.csv", index=False)
    print(f"\nSaved -> mortality_partial_tage_ALL.csv ({len(out)} rows)")
    print(f"Saved -> pathway_specificity_yardstick_mortality_temporal.csv "
          f"(BH family {int(sel.sum())} interpretable tests)")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], int(sys.argv[3]) if len(sys.argv) > 3 else 20000)
