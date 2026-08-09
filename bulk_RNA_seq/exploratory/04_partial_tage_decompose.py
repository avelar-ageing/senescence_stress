#!/usr/bin/env python
"""04_partial_tage_decompose.py

Exact linear decomposition of the fitted ElasticNet tAge model's prediction
by MSigDB Hallmark pathway (the "partial tAge" method described in
Tyshkovskiy/Gladyshev et al. 2026, Nature -- "partial tAge differences
predicted using only genes from the respective module").

ElasticNet is linear: prediction = intercept + sum(coef_i * z_i), where z_i
is the imputed+standardized feature value (exactly what the pipeline's
imputer+scaler steps produce). A pathway's exact contribution to that sum is
just sum(coef_i * z_i) over the genes in that pathway -- computed once from
the correctly, fully-preprocessed data exported by 01/02, with no
re-imputation and no re-normalization of a restricted gene subset.

Also reconstructs the FULL prediction (sum over ALL features + intercept,
species-adjusted) as an exact sanity check against a direct model.predict()
call -- see PARTIAL_TAGE_METHOD.md for the verification (max diff = 0.0).
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


def _patch_simple_imputer(imputer):
    """Same compat patch tage_predict.py applies (sklearn <1.3.0 model
    loaded into a >=1.3.0 environment is missing _fill_dtype)."""
    if not hasattr(imputer, "_fill_dtype"):
        if hasattr(imputer, "statistics_"):
            imputer._fill_dtype = imputer.statistics_.dtype
        else:
            imputer._fill_dtype = np.float64


def decompose(model_path, expr_csv, pathway_csv, out_csv, species="human"):
    model = joblib.load(model_path)
    for _, step in model.steps:
        if type(step).__name__ == "SimpleImputer":
            _patch_simple_imputer(step)
    clock_genes = list(model.feature_names_in_)

    expr = pd.read_csv(expr_csv)
    sample_ids = expr["sample_id"].values
    expr = expr.drop(columns=["sample_id"])
    expr.columns = expr.columns.map(str)

    # Align exactly to the model's expected feature order (same as tage_predict.py).
    missing = [g for g in clock_genes if g not in expr.columns]
    for g in missing:
        expr[g] = np.nan
    X = expr.loc[:, clock_genes]

    imputer = model.named_steps["imputation"]
    scaler = model.named_steps["scaler"]
    estimator = model.named_steps["estimator"]

    X_imp = imputer.transform(X)
    Z = scaler.transform(X_imp)  # samples x genes, standardized -- same space as coef_

    coef = estimator.coef_
    intercept = estimator.intercept_

    contrib = Z * coef[np.newaxis, :]  # samples x genes, per-feature contribution

    # Sanity check: full reconstruction vs a direct model.predict() call.
    full_pred_direct = model.predict(X) * SPECIES_ADJ.get(species, 1.0)
    full_pred_reconstructed = (contrib.sum(axis=1) + intercept) * SPECIES_ADJ.get(species, 1.0)
    max_diff = np.max(np.abs(full_pred_direct - full_pred_reconstructed))
    print(f"[sanity check] max |direct predict() - reconstructed sum| = {max_diff:.10f} "
          f"(should be ~0)", file=sys.stderr)

    pathways = pd.read_csv(pathway_csv)
    gene_to_col = {g: i for i, g in enumerate(clock_genes)}

    results = {"sample_id": sample_ids, "full_tAge_direct": full_pred_direct,
               "full_tAge_reconstructed": full_pred_reconstructed}
    n_genes_used = {}
    for pw, grp in pathways.groupby("pathway"):
        cols = [gene_to_col[g] for g in grp["mouse_gene_id"].astype(str) if g in gene_to_col]
        n_genes_used[pw] = len(cols)
        if len(cols) == 0:
            results[pw] = np.zeros(len(sample_ids))
        else:
            results[pw] = contrib[:, cols].sum(axis=1) * SPECIES_ADJ.get(species, 1.0)

    out_df = pd.DataFrame(results)
    out_df.to_csv(out_csv, index=False)

    n_genes_df = pd.DataFrame({"pathway": list(n_genes_used.keys()), "n_genes_in_model": list(n_genes_used.values())})
    n_genes_df.to_csv(out_csv.replace(".csv", "_ngenes.csv"), index=False)
    print(f"Saved {out_csv} ({len(sample_ids)} samples x {len(pathways['pathway'].unique())} pathways)",
          file=sys.stderr)


if __name__ == "__main__":
    model_path, expr_csv, pathway_csv, out_csv, species = sys.argv[1:6]
    decompose(model_path, expr_csv, pathway_csv, out_csv, species)
