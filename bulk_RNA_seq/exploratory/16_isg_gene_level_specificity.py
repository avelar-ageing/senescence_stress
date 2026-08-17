#!/usr/bin/env python
"""16_isg_gene_level_specificity.py

Gene-level specificity test for the melanocyte interferon result.

WHY. `HALLMARK INTERFERON ALPHA RESPONSE` fails the representation filter
(eff_n 5.9; `Isg15` alone carries 38-51% of the melanocyte effect), so it cannot
support a pathway-level claim -- but the effect itself is real and reproducible
(label-permutation p = 0.0005 pooled, 0.0011-0.0043 per timepoint, both models).
The defensible claim is therefore gene-level. That claim only holds if the genes
carrying it are doing something DIFFERENT in melanocytes than in the other cell
types and the meta-analysis conditions. This script tests exactly that.

WHAT IT COMPUTES. For every interferon-alpha clock gene with a non-zero
coefficient, in every group comparison and both models:
  contrib_diff  mean(coef*z | test) - mean(coef*z | control) -- the gene's
                contribution to the tAge difference, in the model's own units
  cohens_d      standardised effect on that gene's per-sample contribution
  wilcox_p      two-sided Wilcoxon on the per-sample contributions
BH correction is applied within each analysis x model stratum, matching the
convention used for the pathway-level tests.

Output: rerun_outputs/isg_gene_level_specificity.csv
"""
import sys
import warnings

import joblib
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu

warnings.filterwarnings("ignore")

SPECIES_ADJ = 122.5
PATHWAY = "HALLMARK INTERFERON ALPHA RESPONSE"
MODELS = {"scaled": ("EN_Chronoage_Multispecies_Multitissue_scaleddiff.pkl", "scaled_diff"),
          "yugene": ("EN_Chronoage_Multispecies_Multitissue_yugenediff.pkl", "yugene_diff")}
META = {"Contact_inhibited CQ": "CICQ", "Serum_starved CQ": "SSCQ",
        "Replicative CS": "RS", "Stress-induced CS": "SIPS",
        "Oncogene-induced CS": "OIS"}
CTS = ["Fibroblast", "Keratinocyte", "Melanocyte"]
TPS = ["4_days", "10_days", "20_days"]


def _patch(imp):
    if not hasattr(imp, "_fill_dtype"):
        imp._fill_dtype = imp.statistics_.dtype if hasattr(imp, "statistics_") else np.float64


def cohens_d(x, y):
    nx, ny = len(x), len(y)
    p = np.sqrt(((nx - 1) * np.var(x, ddof=1) + (ny - 1) * np.var(y, ddof=1)) / (nx + ny - 2))
    return (np.mean(x) - np.mean(y)) / p if p > 0 else np.nan


def main(rerun_dir, model_dir):
    PT = f"{rerun_dir}/partial_tage"
    pw = pd.read_csv(f"{PT}/hallmark_pathway_mouse_ids.csv")
    ifn_ids = set(pw.loc[pw.pathway == PATHWAY, "mouse_gene_id"].astype(str))
    gt = pd.read_csv("/home/ro/R/x86_64-pc-linux-gnu-library/4.6/tAge/extdata/metadata/Gene_table_mouse.csv")
    sym = dict(zip(gt.Entrez.astype(str), gt["Gene.Symbol"]))

    jobs = [("meta", "meta", c, "Proliferating", s, "meta_analysis")
            for c, s in META.items()]
    jobs += [(ct, ct, "irradiated", "none", ct, "temporal_pooled") for ct in CTS]
    jobs += [(f"{ct}_{tp}", f"{ct}_{tp}", tp, "none", f"{ct}_{tp}", "temporal_bytimepoint")
             for ct in CTS for tp in TPS]

    rows = []
    for mdl, (mfile, suffix) in MODELS.items():
        m = joblib.load(f"{model_dir}/{mfile}")
        _patch(m.named_steps["imputation"])
        feats = list(map(str, m.feature_names_in_))
        idx = {g: i for i, g in enumerate(feats)}
        coef = m.named_steps["estimator"].coef_
        # interferon clock genes that can contribute at all
        genes = [g for g in ifn_ids if g in idx and coef[idx[g]] != 0]

        for prefix, grpfile, test_lab, ctrl_lab, label, analysis in jobs:
            e = pd.read_csv(f"{PT}/{prefix}_{suffix}.csv")
            sid = e["sample_id"].values
            e = e.drop(columns=["sample_id"]); e.columns = e.columns.map(str)
            for g in [g for g in feats if g not in e.columns]:
                e[g] = np.nan
            X = e.loc[:, feats]
            Z = m.named_steps["scaler"].transform(m.named_steps["imputation"].transform(X))
            contrib = Z * coef[np.newaxis, :] * SPECIES_ADJ

            grp = pd.read_csv(f"{PT}/{grpfile}_groups.csv").set_index("sample_id")["group"]
            g_arr = np.array([grp.get(s, None) for s in sid])
            a, b = g_arr == test_lab, g_arr == ctrl_lab
            if a.sum() < 2 or b.sum() < 2:
                continue
            for gid in genes:
                v = contrib[:, idx[gid]]
                try:
                    p = mannwhitneyu(v[a], v[b]).pvalue
                except ValueError:
                    p = np.nan
                rows.append(dict(analysis=analysis, label=label, model=mdl,
                                 entrez=gid, gene=sym.get(gid, gid),
                                 coefficient=coef[idx[gid]],
                                 contrib_diff=v[a].mean() - v[b].mean(),
                                 cohens_d=cohens_d(v[a], v[b]), wilcox_p=p,
                                 n_test=int(a.sum()), n_control=int(b.sum())))

    out = pd.DataFrame(rows)
    def bh(p):
        """Benjamini-Hochberg, same as R p.adjust(method='BH')."""
        p = np.asarray(p, dtype=float)
        n = len(p)
        o = np.argsort(p)
        ranked = p[o] * n / (np.arange(n) + 1)
        ranked = np.minimum.accumulate(ranked[::-1])[::-1]
        out_ = np.empty(n); out_[o] = np.minimum(ranked, 1.0)
        return out_

    out["p_adj"] = np.nan
    for (an, mdl), sub in out.groupby(["analysis", "model"]):
        ok = sub.wilcox_p.notna()
        out.loc[sub.index[ok], "p_adj"] = bh(sub.wilcox_p[ok].values)
    path = f"{rerun_dir}/isg_gene_level_specificity.csv"
    out.to_csv(path, index=False)
    print(f"Saved {path} ({len(out)} rows, {out.gene.nunique()} genes)", file=sys.stderr)
    return out


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
