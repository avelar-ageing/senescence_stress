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

GATE (top5_max <= 0.65). A set is interpreted at set level only if its five
largest-contributing genes carry no more than 65% of its total contribution
under BOTH models. 17 of 50 sets qualify. This replaced an earlier eff_n >= 9
threshold, which selected an almost identical list but could not be justified as
a number and was non-monotone in the quantity of interest: it retained MTORC1
(top5 60%) while excluding KRAS SIGNALING UP (top5 52%). Filtering directly on
top5 removes that artefact, is stateable in one clause, and is a stipulation
about what may be called a set-level effect rather than an empirical threshold.
eff_n is still reported alongside (rho 0.95 with -top5).

Output: rerun_outputs/pathway_representation.csv
"""
import sys
import warnings
import joblib
import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")

TOP5_CUT = 0.65   # five largest genes may carry at most this share, both models
MODELS = {"scaled": "EN_Chronoage_Multispecies_Multitissue_scaleddiff.pkl",
          "yugene": "EN_Chronoage_Multispecies_Multitissue_yugenediff.pkl"}
EXPR = {"scaled": "meta_scaled_diff.csv", "yugene": "meta_yugene_diff.csv"}
# THE GATE IS CLOCK-SPECIFIC. Gene domination is a property of the model's
# sparsity, so the gate has to be recomputed for whichever clock a section
# reports. On the chronological clocks a typical set is carried by 20-22
# non-zero genes and its five largest carry a median 65-73% of the
# contribution, so the gate excludes most sets. The mortality clock is dense
# ridge: a typical set is carried by 125 genes, median top-5 share is 25%, and
# 49 of 50 sets pass - only HALLMARK PANCREAS BETA CELLS fails. Using the
# chronological gate for mortality results would discard 32 sets that are well
# represented in the model actually being used, and would also set the wrong
# BH family. Run with --mortality for the mortality gate.
MORTALITY_MODELS = {"mortality": "EN_Mortality_Multispecies_Multitissue_scaleddiff.pkl"}
MORTALITY_EXPR = {"mortality": "meta_scaled_diff.csv"}


def _patch(imp):
    if not hasattr(imp, "_fill_dtype"):
        imp._fill_dtype = imp.statistics_.dtype if hasattr(imp, "statistics_") else np.float64


def representation(rerun_dir, model_dir, mapping="hallmark_pathway_mouse_ids.csv",
                   which="chronoage"):
    PT = f"{rerun_dir}/partial_tage"
    pw = pd.read_csv(f"{PT}/{mapping}")
    models = MORTALITY_MODELS if which == "mortality" else MODELS
    exprs = MORTALITY_EXPR if which == "mortality" else EXPR
    frames = {}
    for mdl, mfile in models.items():
        m = joblib.load(f"{model_dir}/{mfile}")
        _patch(m.named_steps["imputation"])
        feats = list(map(str, m.feature_names_in_))
        idx = {g: i for i, g in enumerate(feats)}
        coef = m.named_steps["estimator"].coef_

        e = pd.read_csv(f"{PT}/{exprs[mdl]}").drop(columns=["sample_id"])
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

    if which == "mortality":
        d = frames["mortality"].add_suffix("_mortality")
        d["top5_max"] = d.top5_mortality
        d["top1_max"] = d.top1_mortality
        d["eff_n_min"] = d.eff_n_mortality
        d["interpretable"] = d.top5_max <= TOP5_CUT
        d["tier"] = np.where(d.interpretable, "INTERPRETABLE", "GENE-DOMINATED")
        return d.sort_values("top5_max")
    d = frames["scaled"].join(frames["yugene"], lsuffix="_scaled", rsuffix="_yugene")
    d["eff_n_min"] = d[["eff_n_scaled", "eff_n_yugene"]].min(axis=1)
    # PRIMARY CRITERION (top5_max): worst-case share of the set's total contribution
    # carried by its five largest genes, across the two models. Used as the gate
    # because it is directly interpretable -- "five genes carry X% of this set's
    # signal" -- and because it is monotone in the quantity of interest, which
    # eff_n is not (eff_n retained MTORC1 at top5 60% while excluding KRAS
    # SIGNALING UP at top5 52%). eff_n is retained as a reported companion
    # statistic; the two agree at rho 0.95.
    d["top5_max"] = d[["top5_scaled", "top5_yugene"]].max(axis=1)
    d["top1_max"] = d[["top1_scaled", "top1_yugene"]].max(axis=1)
    d["interpretable"] = d.top5_max <= TOP5_CUT
    d["tier"] = np.where(d.interpretable, "INTERPRETABLE", "GENE-DOMINATED")
    return d.sort_values("top5_max")


if __name__ == "__main__":
    rerun_dir, model_dir = sys.argv[1], sys.argv[2]
    which = "mortality" if "--mortality" in sys.argv else "chronoage"
    args = [a for a in sys.argv[3:] if not a.startswith("--")]
    mapping = args[0] if args else "hallmark_pathway_mouse_ids.csv"
    d = representation(rerun_dir, model_dir, mapping, which)
    out = (f"{rerun_dir}/pathway_representation_mortality.csv" if which == "mortality"
           else f"{rerun_dir}/pathway_representation.csv")
    d.to_csv(out)
    cols = ([c for c in ["n_clock_mortality", "n_nonzero_mortality"] if c in d]
            if which == "mortality"
            else ["n_clock_scaled", "n_nonzero_scaled", "n_nonzero_yugene"])
    print(d[cols + ["top1_max", "top5_max", "eff_n_min", "tier"]].to_string(
          float_format=lambda x: f"{x:7.2f}"))
    print(f"\ntiers (top5 <= {TOP5_CUT:.0%}): {dict(d.tier.value_counts())}")
    print(f"Saved -> {out}")
