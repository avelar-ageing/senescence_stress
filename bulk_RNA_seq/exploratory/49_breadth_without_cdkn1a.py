"""49_breadth_without_cdkn1a.py

WHY. 2.2.4 names the sets whose contribution shifts significantly, in the same direction,
in all three cell types (irradiated vs untreated, each cell type pooled over timepoints;
Wilcoxon, BH over the 150 set x cell-type tests; exploratory/22 -> 21). exploratory/47
shows that Cdkn1a is the largest single contributor to several of these sets. This
script repeats the same test with Cdkn1a removed from every set, to show which of those
cross-cell-type results depend on it.

WHAT. For each of the 50 sets and 3 cell types: the set's per-sample contribution with
and without Cdkn1a, the irradiated-minus-untreated difference, the Wilcoxon p and the BH
q over the 150 tests (with Cdkn1a: reproduces mortality_partial_tage_ALL.csv, checked).
Then, with and without Cdkn1a: significant (q < 0.05) in all three cell types, and in the
same direction in all three.

USAGE: 49_breadth_without_cdkn1a.py <rerun_dir> <model_dir>
OUT:   <rerun_dir>/breadth_without_cdkn1a.csv   one row per set x cell type, plus
       columns all3_sig_same_direction_{with,without}_cdkn1a per set
"""
import os
import subprocess
import sys
import warnings

import joblib
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu

warnings.filterwarnings("ignore")
MORT = "EN_Mortality_Multispecies_Multitissue_scaleddiff.pkl"
CELLS = ["Fibroblast", "Keratinocyte", "Melanocyte"]


def tage_gene_table():
    p = os.environ.get("TAGE_GENE_TABLE")
    if not p:
        p = subprocess.run(["Rscript", "-e", "cat(system.file('extdata/metadata/Gene_table_mouse.csv', package='tAge'))"],
                           capture_output=True, text=True, check=True).stdout.strip()
    if not p or not os.path.exists(p):
        sys.exit("Gene_table_mouse.csv not found: install the tAge R package or set TAGE_GENE_TABLE")
    return p


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
    feats = list(map(str, m.feature_names_in_)); idx = {g: i for i, g in enumerate(feats)}
    coef = m.named_steps["estimator"].coef_
    gt = pd.read_csv(tage_gene_table())
    cdk = idx[str(gt.loc[gt["Gene.Symbol"] == "Cdkn1a", "Entrez"].iloc[0])]
    pw = pd.read_csv(f"{PT}/hallmark_pathway_mouse_ids.csv")
    sets = {n: np.array([idx[x] for x in s.mouse_gene_id.astype(str) if x in idx]) for n, s in pw.groupby("pathway")}
    sets = {k: v for k, v in sets.items() if len(v)}

    rows = []
    for ct in CELLS:
        e = pd.read_csv(f"{PT}/{ct}_scaled_diff.csv"); sid = e.sample_id.values
        e = e.drop(columns=["sample_id"]); e.columns = e.columns.map(str)
        for g in [g for g in feats if g not in e.columns]:
            e[g] = np.nan
        C = m.named_steps["scaler"].transform(imp.transform(e.loc[:, feats])) * coef[None, :]
        gg = pd.read_csv(f"{PT}/{ct}_groups.csv").set_index("sample_id")["group"]
        lab = np.array([gg.get(s) for s in sid]); a = lab == "irradiated"; b = lab == "none"
        for nm, cols in sets.items():
            r = dict(cell_type=ct, pathway=nm, cdkn1a_in_set=bool(cdk in cols))
            for tag, cc in (("with", cols), ("without", cols[cols != cdk])):
                s = C[:, cc].sum(axis=1)
                r[f"diff_{tag}_cdkn1a"] = s[a].mean() - s[b].mean()
                r[f"p_{tag}_cdkn1a"] = mannwhitneyu(s[a], s[b]).pvalue
            rows.append(r)
    d = pd.DataFrame(rows)
    for tag in ("with", "without"):
        d[f"q_{tag}_cdkn1a"] = bh(d[f"p_{tag}_cdkn1a"].values)

    pub = pd.read_csv(f"{rerun_dir}/mortality_partial_tage_ALL.csv")
    pub = pub[pub.analysis == "temporal_pooled"].set_index(["label", "pathway"]).p_adj
    chk = np.array([pub.loc[(r.cell_type, r.pathway)] for r in d.itertuples()])
    assert np.allclose(chk, d.q_with_cdkn1a.values, rtol=1e-9, atol=1e-12), "with-Cdkn1a q does not reproduce exploratory/22"

    for tag in ("with", "without"):
        g = d.groupby("pathway").apply(lambda x: (x[f"q_{tag}_cdkn1a"] < 0.05).all()
                                       and np.sign(x[f"diff_{tag}_cdkn1a"]).nunique() == 1)
        d[f"all3_sig_same_direction_{tag}_cdkn1a"] = d.pathway.map(g)
    d.to_csv(f"{rerun_dir}/breadth_without_cdkn1a.csv", index=False)

    for tag in ("with", "without"):
        sig3 = d.groupby("pathway").apply(lambda x: (x[f"q_{tag}_cdkn1a"] < 0.05).all())
        same = d.groupby("pathway")[f"all3_sig_same_direction_{tag}_cdkn1a"].first()
        print(f"{tag:7s} Cdkn1a: significant in all 3 cell types: {sorted(s.replace('HALLMARK ', '') for s in sig3[sig3].index)}")
        print(f"          of these, same direction in all 3: {sorted(s.replace('HALLMARK ', '') for s in same[same].index)}")
    print("Saved -> breadth_without_cdkn1a.csv")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
