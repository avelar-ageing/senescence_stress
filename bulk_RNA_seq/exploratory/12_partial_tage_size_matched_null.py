#!/usr/bin/env python
"""11_partial_tage_size_matched_null.py

Size-matched null for the partial-tAge decomposition.

Motivation: a pathway's partial-tAge effect size partly tracks how many of the
clock's features it contains (Spearman rho = 0.33 between signed Cohen's d and
n_genes_in_model across the 50 Hallmark sets, meta-analysis, yugene model).
Larger sets sum more coefficient x z terms and so can reach larger |d| for
reasons that have nothing to do with the biology of the set. This script asks,
for each pathway: does its observed effect exceed what a RANDOM gene set of the
same size, drawn from the clock's own features, achieves on the same samples?

Method. The contribution matrix contrib = Z * coef is computed exactly as in
04_partial_tage_decompose.py (same imputer/scaler/coefficients, same
whole-transcriptome preprocessed input). For a pathway of n clock features we
draw N random size-n feature sets without replacement from the model's 10,487
features, sum each to a per-sample score, and compute Cohen's d between the
test and control groups. The empirical two-sided p is the fraction of null
draws with |d_null| >= |d_observed|.

Reporting convention: raw empirical p, the p-floor (1/(N+1)), and the effect
size are reported; NO multiple-testing correction is applied to Monte-Carlo
null p-values, which are a pre-registered robustness check rather than a
discovery family.
"""
import sys
import warnings
import joblib
import numpy as np
import pandas as pd
from sklearn.exceptions import InconsistentVersionWarning

warnings.filterwarnings("ignore", category=InconsistentVersionWarning)
warnings.filterwarnings("ignore", category=UserWarning)

SPECIES_ADJ = {"human": 122.5, "mouse": 48, "rat": 50.4, "monkey": 39}
N_DRAWS = 1000
SEED = 20260816


def _patch_simple_imputer(imputer):
    if not hasattr(imputer, "_fill_dtype"):
        imputer._fill_dtype = (imputer.statistics_.dtype
                               if hasattr(imputer, "statistics_") else np.float64)


def cohens_d(x, y):
    nx, ny = len(x), len(y)
    pooled = np.sqrt(((nx - 1) * np.var(x, ddof=1) + (ny - 1) * np.var(y, ddof=1)) / (nx + ny - 2))
    return (np.mean(x) - np.mean(y)) / pooled if pooled > 0 else np.nan


def run(model_path, expr_csv, pathway_csv, groups_csv, test_label, control_label,
        label, analysis, out_rows, species="human"):
    model = joblib.load(model_path)
    for _, step in model.steps:
        if type(step).__name__ == "SimpleImputer":
            _patch_simple_imputer(step)
    clock_genes = list(model.feature_names_in_)

    expr = pd.read_csv(expr_csv)
    sample_ids = expr["sample_id"].values
    expr = expr.drop(columns=["sample_id"])
    expr.columns = expr.columns.map(str)
    for g in [g for g in clock_genes if g not in expr.columns]:
        expr[g] = np.nan
    X = expr.loc[:, clock_genes]

    Z = model.named_steps["scaler"].transform(model.named_steps["imputation"].transform(X))
    coef = model.named_steps["estimator"].coef_
    contrib = Z * coef[np.newaxis, :] * SPECIES_ADJ.get(species, 1.0)

    groups = pd.read_csv(groups_csv).set_index("sample_id")["group"]
    grp = np.array([groups.get(s, None) for s in sample_ids])
    is_test, is_ctrl = grp == test_label, grp == control_label
    if is_test.sum() == 0 or is_ctrl.sum() == 0:
        print(f"  [skip] {label}: no samples for {test_label}/{control_label}", file=sys.stderr)
        return

    pathways = pd.read_csv(pathway_csv)
    gene_to_col = {g: i for i, g in enumerate(clock_genes)}
    rng = np.random.default_rng(SEED)
    n_feat = len(clock_genes)

    for pw, sub in pathways.groupby("pathway"):
        cols = [gene_to_col[g] for g in sub["mouse_gene_id"].astype(str) if g in gene_to_col]
        n = len(cols)
        if n == 0:
            continue
        obs_scores = contrib[:, cols].sum(axis=1)
        d_obs = cohens_d(obs_scores[is_test], obs_scores[is_ctrl])

        null_d = np.empty(N_DRAWS)
        for k in range(N_DRAWS):
            pick = rng.choice(n_feat, size=n, replace=False)
            s = contrib[:, pick].sum(axis=1)
            null_d[k] = cohens_d(s[is_test], s[is_ctrl])
        null_d = null_d[np.isfinite(null_d)]

        p_emp = (np.sum(np.abs(null_d) >= abs(d_obs)) + 1) / (len(null_d) + 1)
        out_rows.append(dict(
            analysis=analysis, label=label, pathway=pw, n_genes_in_model=n,
            cohens_d=d_obs, null_mean_abs_d=float(np.mean(np.abs(null_d))),
            null_q95_abs_d=float(np.quantile(np.abs(null_d), 0.95)),
            p_empirical=p_emp, p_floor=1.0 / (len(null_d) + 1), n_draws=len(null_d),
        ))


if __name__ == "__main__":
    rerun_dir, model_dir = sys.argv[1], sys.argv[2]
    # optional: alternative grouping (e.g. the paper's modules) + output tag
    mapping_name = sys.argv[3] if len(sys.argv) > 3 else "hallmark_pathway_mouse_ids.csv"
    out_tag = sys.argv[4] if len(sys.argv) > 4 else "partial_tage_size_matched_null"
    PT = f"{rerun_dir}/partial_tage"
    pathway_csv = f"{PT}/{mapping_name}"
    models = {"scaled": f"{model_dir}/EN_Chronoage_Multispecies_Multitissue_scaleddiff.pkl",
              "yugene": f"{model_dir}/EN_Chronoage_Multispecies_Multitissue_yugenediff.pkl"}
    expr_suffix = {"scaled": "scaled_diff", "yugene": "yugene_diff"}

    meta_conditions = {"Contact_inhibited CQ": "CICQ", "Serum_starved CQ": "SSCQ",
                       "Replicative CS": "RS", "Stress-induced CS": "SIPS",
                       "Oncogene-induced CS": "OIS"}

    all_rows = []
    for mdl, mpath in models.items():
        rows = []
        for cond_value, short in meta_conditions.items():
            print(f"[{mdl}] meta: {short}", file=sys.stderr)
            run(mpath, f"{PT}/meta_{expr_suffix[mdl]}.csv", pathway_csv, f"{PT}/meta_groups.csv",
                cond_value, "Proliferating", short, "meta_analysis", rows)
        for ct in ["Fibroblast", "Keratinocyte", "Melanocyte"]:
            print(f"[{mdl}] temporal pooled: {ct}", file=sys.stderr)
            run(mpath, f"{PT}/{ct}_{expr_suffix[mdl]}.csv", pathway_csv, f"{PT}/{ct}_groups.csv",
                "irradiated", "none", ct, "temporal_pooled", rows)
        # per-timepoint: 3 cell types x 3 timepoints, 6v6 vs that cell type's own baseline.
        # Needed because the per-timepoint results are what the temporal write-up leans on;
        # the pooled null does not cover them.
        for ct in ["Fibroblast", "Keratinocyte", "Melanocyte"]:
            for tp in ["4_days", "10_days", "20_days"]:
                grp = f"{ct}_{tp}"
                print(f"[{mdl}] temporal bytimepoint: {grp}", file=sys.stderr)
                run(mpath, f"{PT}/{grp}_{expr_suffix[mdl]}.csv", pathway_csv,
                    f"{PT}/{grp}_groups.csv", tp, "none", grp, "temporal_bytimepoint", rows)
        for r in rows:
            r["model"] = mdl
        all_rows.extend(rows)

    out = pd.DataFrame(all_rows)
    out_path = f"{rerun_dir}/{out_tag}.csv"
    out.to_csv(out_path, index=False)
    print(f"\nSaved {out_path} ({len(out)} rows)", file=sys.stderr)
    sig = out[out.p_empirical < 0.05]
    print(f"Pathways exceeding the size-matched null (p_emp<0.05): {len(sig)} of {len(out)}",
          file=sys.stderr)
