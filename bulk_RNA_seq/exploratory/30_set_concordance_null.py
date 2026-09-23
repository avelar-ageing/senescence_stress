#!/usr/bin/env python3
"""30_set_concordance_null.py

Is a gene set's cancellation structure MORE consistent across conditions than a
random group of genes of the same size would be?

29_set_cancellation_structure.py measured, for each Hallmark set and each pair of
conditions, the Spearman correlation between the two conditions' per-gene
contribution vectors. Because a gene's clock coefficient is fixed, that asks
whether the same genes moved the same way in both conditions.

A raw correlation cannot be read on its own. Any two arrest conditions share a
great deal of biology, so ANY group of genes will show some concordance between
them - that shared baseline is a property of the comparison, not of the set. The
null here supplies it: for each set, draw random gene groups of the same size from
the clock's own features and take their between-condition correlation. The
empirical p is the fraction of draws reaching the observed correlation, so a small
p means the set is more internally consistent across conditions than its size
alone explains.

CONVENTIONS. 20,000 draws and the clock's feature pool, inherited from
20_pathway_specificity_yardstick.py. Raw empirical p and its floor are reported;
BH is applied within each dataset (arrest conditions, time course) over that
dataset's set x condition-pair tests, matching how the yardstick is corrected.
Draws are bucketed by (condition pair, set size) because the null depends only on
those two things, so sets of equal size share a null.

Output: rerun_outputs/set_concordance_null.csv
"""
import sys, warnings, itertools
import joblib, numpy as np, pandas as pd
from sklearn.exceptions import InconsistentVersionWarning
warnings.filterwarnings("ignore", category=InconsistentVersionWarning)

N_DRAWS = 20000
SEED = 1
DATASETS = [("arrest_conditions", ["meta"], "Proliferating"),
            ("time_course", [f"{ct}_{tp}" for ct in ("Fibroblast", "Keratinocyte", "Melanocyte")
                             for tp in ("4_days", "10_days", "20_days")], "none")]


def bh(p):
    p = np.asarray(p, float); n = len(p); o = np.argsort(p); a = np.empty(n)
    a[o] = np.minimum.accumulate((p[o] * n / (np.arange(n) + 1))[::-1])[::-1]
    return np.clip(a, 0, 1)


def rowwise_spearman(A, B):
    """Spearman per row of two (N, k) matrices, via Pearson on within-row ranks."""
    ra = np.argsort(np.argsort(A, axis=1), axis=1).astype(np.float64)
    rb = np.argsort(np.argsort(B, axis=1), axis=1).astype(np.float64)
    ra -= ra.mean(axis=1, keepdims=True); rb -= rb.mean(axis=1, keepdims=True)
    num = (ra * rb).sum(axis=1)
    den = np.sqrt((ra ** 2).sum(axis=1) * (rb ** 2).sum(axis=1))
    return np.divide(num, den, out=np.zeros_like(num), where=den > 0)


def main(rerun_dir, model_dir):
    PT = f"{rerun_dir}/partial_tage"
    m = joblib.load(f"{model_dir}/EN_Mortality_Multispecies_Multitissue_scaleddiff.pkl")
    imp = m.named_steps["imputation"]
    if not hasattr(imp, "_fill_dtype"):
        imp._fill_dtype = imp.statistics_.dtype
    feats = list(map(str, m.feature_names_in_))
    coef = m.named_steps["estimator"].coef_

    pw = pd.read_csv(f"{PT}/hallmark_pathway_mouse_ids.csv")
    pw["mouse_gene_id"] = pw["mouse_gene_id"].astype(str)
    fi = {g: i for i, g in enumerate(feats)}
    sets = {p: np.array([fi[g] for g in d["mouse_gene_id"] if g in fi]) for p, d in pw.groupby("pathway")}
    sets = {p: ix for p, ix in sets.items() if len(ix) >= 5}

    rows = []
    for dset, stems, control in DATASETS:
        delta = {}
        for stem in stems:
            expr = pd.read_csv(f"{PT}/{stem}_scaled_diff.csv")
            sid = expr["sample_id"].values
            expr = expr.drop(columns=["sample_id"]); expr.columns = expr.columns.map(str)
            for g in [g for g in feats if g not in expr.columns]:
                expr[g] = np.nan
            Z = m.named_steps["scaler"].transform(imp.transform(expr.loc[:, feats]))
            contrib = Z * coef[np.newaxis, :]
            grp = pd.read_csv(f"{PT}/{stem}_groups.csv").set_index("sample_id").loc[sid, "group"].values
            ctrl = contrib[grp == control].mean(axis=0)
            for c in [x for x in pd.unique(grp) if x != control]:
                key = c if stem == "meta" else f"{stem.split('_')[0]}_{c}"
                delta[key] = contrib[grp == c].mean(axis=0) - ctrl

        pairs = list(itertools.combinations(sorted(delta), 2))
        sizes = sorted({len(ix) for ix in sets.values()})
        print(f"{dset}: {len(pairs)} condition pairs x {len(sets)} sets, {len(sizes)} distinct sizes")
        rng = np.random.default_rng(SEED)
        for a, b in pairs:
            va, vb = delta[a], delta[b]
            nulls = {}
            for k in sizes:
                idx = rng.choice(len(feats), size=(N_DRAWS, k), replace=True)
                nulls[k] = rowwise_spearman(va[idx], vb[idx])
            for p, ix in sets.items():
                k = len(ix)
                obs = rowwise_spearman(va[ix][None, :], vb[ix][None, :])[0]
                nul = nulls[k]
                rows.append(dict(dataset=dset, pathway=p, condition_a=a, condition_b=b,
                                 n_genes=k, rho=obs,
                                 null_median=float(np.median(nul)),
                                 p_emp=(1 + (nul >= obs).sum()) / (N_DRAWS + 1),
                                 p_floor=1.0 / (N_DRAWS + 1)))
    D = pd.DataFrame(rows)
    for dset, g in D.groupby("dataset"):
        D.loc[g.index, "q_emp"] = bh(g["p_emp"].values)
        print(f"  {dset}: {len(g)} tests, {(D.loc[g.index,'q_emp'] < 0.05).sum()} at 5% FDR, "
              f"{(D.loc[g.index,'q_emp'] < 0.10).sum()} at 10%")
    D.to_csv(f"{rerun_dir}/set_concordance_null.csv", index=False)
    print(f"Saved -> set_concordance_null.csv ({len(D)} rows)")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
