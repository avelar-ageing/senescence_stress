#!/usr/bin/env python
"""14_pathway_representation.py

How well is each gene set actually REPRESENTED in the tAge clock?

A partial-tAge score is a sum of coef_i * z_i over a gene set. If most of that
sum comes from one or two genes, "pathway X moves tAge" is really "gene Y moves
tAge". This script quantifies that per gene set, so pathway-level claims can be
filtered to the sets that can support them.

Metrics, computed on the meta-analysis matrix (n=230, the largest):
  n_clock      genes of the set present among the clock's 10,487 features
  n_nonzero    of those, how many have a non-zero elastic-net coefficient
               (the rest contribute exactly nothing, ever)
  eff_n        effective number of contributing genes = 1 / sum(share^2),
               where share_i = mean|coef_i * z_i| / sum over the set.
               Inverse-Simpson: eff_n = 10 means the set behaves like ~10
               equally-weighted genes; eff_n = 2 means two genes carry it.
  top1, top5   share of the set's total |contribution| from its largest
               1 and 5 genes

Suggested tiers (used in DISCREPANCY_REPORT/PARTIAL_TAGE_VALIDITY.md):
  WELL     eff_n >= 14 on BOTH models
  MODERATE eff_n >= 9  on both
  POOR     otherwise -- do not make pathway-level claims from these

Output: rerun_outputs/pathway_representation.csv
"""
import sys
import warnings
import joblib
import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")

MODELS = {"scaled": "EN_Chronoage_Multispecies_Multitissue_scaleddiff.pkl",
          "yugene": "EN_Chronoage_Multispecies_Multitissue_yugenediff.pkl"}
EXPR = {"scaled": "meta_scaled_diff.csv", "yugene": "meta_yugene_diff.csv"}


def _patch(imp):
    if not hasattr(imp, "_fill_dtype"):
        imp._fill_dtype = imp.statistics_.dtype if hasattr(imp, "statistics_") else np.float64


def representation(rerun_dir, model_dir, mapping="hallmark_pathway_mouse_ids.csv"):
    PT = f"{rerun_dir}/partial_tage"
    pw = pd.read_csv(f"{PT}/{mapping}")
    frames = {}
    for mdl, mfile in MODELS.items():
        m = joblib.load(f"{model_dir}/{mfile}")
        _patch(m.named_steps["imputation"])
        feats = list(map(str, m.feature_names_in_))
        idx = {g: i for i, g in enumerate(feats)}
        coef = m.named_steps["estimator"].coef_

        e = pd.read_csv(f"{PT}/{EXPR[mdl]}").drop(columns=["sample_id"])
        e.columns = e.columns.map(str)
        for g in [g for g in feats if g not in e.columns]:
            e[g] = np.nan
        X = e.loc[:, feats]
        Z = m.named_steps["scaler"].transform(m.named_steps["imputation"].transform(X))
        # mean absolute per-gene contribution across samples
        mc = np.abs(Z * coef[np.newaxis, :]).mean(axis=0)

        rows = []
        for name, sub in pw.groupby("pathway"):
            cols = [idx[g] for g in sub["mouse_gene_id"].astype(str) if g in idx]
            if not cols:
                continue
            v = mc[cols]
            tot = v.sum()
            share = v / tot if tot > 0 else v
            srt = np.sort(share)[::-1]
            rows.append(dict(
                pathway=name, n_clock=len(cols),
                n_nonzero=int(sum(coef[c] != 0 for c in cols)),
                eff_n=float(1.0 / np.sum(share ** 2)) if tot > 0 else 0.0,
                top1=float(srt[0]) if tot > 0 else np.nan,
                top5=float(srt[:5].sum()) if tot > 0 else np.nan,
            ))
        frames[mdl] = pd.DataFrame(rows).set_index("pathway")

    d = frames["scaled"].join(frames["yugene"], lsuffix="_scaled", rsuffix="_yugene")
    d["eff_n_min"] = d[["eff_n_scaled", "eff_n_yugene"]].min(axis=1)
    d["tier"] = np.where(d.eff_n_min >= 14, "WELL",
                np.where(d.eff_n_min >= 9, "MODERATE", "POOR"))
    return d.sort_values("eff_n_min", ascending=False)


if __name__ == "__main__":
    rerun_dir, model_dir = sys.argv[1], sys.argv[2]
    mapping = sys.argv[3] if len(sys.argv) > 3 else "hallmark_pathway_mouse_ids.csv"
    d = representation(rerun_dir, model_dir, mapping)
    out = f"{rerun_dir}/pathway_representation.csv"
    d.to_csv(out)
    print(d[["n_clock_scaled", "n_nonzero_scaled", "n_nonzero_yugene",
             "eff_n_scaled", "eff_n_yugene", "eff_n_min", "tier"]].to_string(
             float_format=lambda x: f"{x:7.1f}"))
    print(f"\ntiers: {dict(d.tier.value_counts())}")
    print(f"Saved -> {out}")
