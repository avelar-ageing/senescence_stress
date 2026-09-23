#!/usr/bin/env python3
"""29_set_cancellation_structure.py

TWO QUESTIONS about the near-cancellation inside gene sets, both asked on the
per-gene contributions the existing decomposition computes but never saves.

Recall the quantity. For gene g the contribution to a condition's mortality-score
difference is coef_g x (mean z of the condition - mean z of proliferating
controls). coef_g is a property of the CLOCK and is identical in every condition;
only the expression term changes. A set's reported contribution is the sum over
its genes, and 2.1.5.3 shows that sum is a small residue of much larger opposing
movements.

QUESTION 1 - IS EACH SIDE LARGER THAN CHANCE?
The published test asks whether a set's NET exceeds what a random set of the same
size achieves. A set whose genes move a lot in both directions has a small net and
is therefore indistinguishable from a set where nothing moved. Here the summed
positive contributions and the summed negative contributions are each tested
against the same weight-matched null, so "moves a lot in both directions" becomes
a statement that can be made rather than an invisible case.

QUESTION 2 - DO SETS CANCEL VIA THE SAME GENES IN EVERY CONDITION?
Because coef_g is fixed, a gene switches from pushing the score up to pushing it
down only if its expression change reverses. So asking whether two conditions
cancel "the same way" is asking whether the same genes moved the same way, with
the clock's coefficients as the weighting. Two measures per set per condition
pair: the Spearman correlation of the per-gene contribution vectors, and the
Jaccard overlap of the up-pushing gene sets. A set that cancels through the same
genes everywhere is showing a property of the clock's coefficient structure; one
that cancels through different genes is showing condition-specific biology.

BOTH DATASETS. Run over the five arrest conditions (each against the 91
proliferating controls) and the nine cell-type x timepoint groups of the
irradiation time course (each against its own cell type's untreated samples).

MULTIPLE TESTING. These 500 and 900 side tests are a discovery scan, not a
pre-registered robustness check, so unlike the net test they carry BH. ONE FAMILY
PER DATASET, matching how the gene-set yardstick is corrected: the arrest
conditions and the time course are separate analyses reported in separate
sections. Within a dataset the up and down tests share a family, because they are
the same question asked in two directions over the same sets. The raw empirical p
and its floor are kept alongside so a value sitting at 1/(N+1) is still visible.

Output: rerun_outputs/set_cancellation_sides.csv      (question 1)
        rerun_outputs/set_cancellation_concordance.csv (question 2)
        rerun_outputs/set_gene_contributions.csv       (the per-gene layer)
"""
import sys, warnings, itertools
import joblib, numpy as np, pandas as pd
from scipy.stats import spearmanr
from sklearn.exceptions import InconsistentVersionWarning
warnings.filterwarnings("ignore", category=InconsistentVersionWarning)

N_DRAWS = 20000
SEED = 1
CONTROL = "Proliferating"


def bh(p):
    p = np.asarray(p, float); n = len(p); o = np.argsort(p); a = np.empty(n)
    a[o] = np.minimum.accumulate((p[o] * n / (np.arange(n) + 1))[::-1])[::-1]
    return np.clip(a, 0, 1)


# (dataset label, expression stem, control label). The time course keeps each
# cell type separate because its baseline is that cell type's own untreated cells.
DATASETS = [("arrest_conditions", ["meta"], "Proliferating")] + \
           [("time_course", [f"{ct}_{tp}" for ct in ("Fibroblast", "Keratinocyte", "Melanocyte")
                             for tp in ("4_days", "10_days", "20_days")], "none")]


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
    fidx = {g: i for i, g in enumerate(feats)}
    sets = {p: np.array([fidx[g] for g in d["mouse_gene_id"] if g in fidx])
            for p, d in pw.groupby("pathway")}
    sets = {p: ix for p, ix in sets.items() if len(ix) >= 5}

    gene_rows, side_rows, conc_rows = [], [], []
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
        print(f"{dset}: {len(delta)} groups x {len(sets)} sets")

        rng = np.random.default_rng(SEED)
        for c, d in delta.items():
            for p, ix in sets.items():
                for i in ix:
                    gene_rows.append((dset, c, p, feats[i], coef[i], d[i]))
                k = len(ix)
                v = d[ix]
                obs_up, obs_dn, obs_net = v[v > 0].sum(), v[v < 0].sum(), v.sum()
                dd = d[rng.choice(len(feats), size=(N_DRAWS, k), replace=True)]
                nu = np.where(dd > 0, dd, 0).sum(axis=1)
                nd = np.where(dd < 0, dd, 0).sum(axis=1)
                nn = dd.sum(axis=1)
                side_rows.append(dict(
                    dataset=dset, condition=c, pathway=p, n_genes=k,
                    up=obs_up, down=obs_dn, net=obs_net,
                    surviving_pct=100*abs(obs_net)/(obs_up - obs_dn) if obs_up != obs_dn else np.nan,
                    p_up=(1 + (nu >= obs_up).sum()) / (N_DRAWS + 1),
                    p_down=(1 + (nd <= obs_dn).sum()) / (N_DRAWS + 1),
                    p_net=(1 + (np.abs(nn) >= abs(obs_net)).sum()) / (N_DRAWS + 1),
                    p_floor=1.0 / (N_DRAWS + 1)))
        for p, ix in sets.items():
            for a, b in itertools.combinations(sorted(delta), 2):
                va, vb = delta[a][ix], delta[b][ix]
                ua, ub = set(np.where(va > 0)[0]), set(np.where(vb > 0)[0])
                conc_rows.append(dict(
                    dataset=dset, pathway=p, condition_a=a, condition_b=b, n_genes=len(ix),
                    rho_contribution=spearmanr(va, vb).statistic if len(ix) > 2 else np.nan,
                    jaccard_up_genes=len(ua & ub) / len(ua | ub) if (ua | ub) else np.nan,
                    frac_up_a=len(ua)/len(ix), frac_up_b=len(ub)/len(ix)))

    S = pd.DataFrame(side_rows)
    # BH within each dataset, over that dataset's up AND down tests together
    for dset, g in S.groupby("dataset"):
        both = np.concatenate([g["p_up"].values, g["p_down"].values])
        q = bh(both); half = len(g)
        S.loc[g.index, "q_up"] = q[:half]
        S.loc[g.index, "q_down"] = q[half:]
        print(f"  {dset}: BH over {len(both)} side tests "
              f"({(q < 0.05).sum()} at 5%, {(q < 0.10).sum()} at 10%)")
    S.to_csv(f"{rerun_dir}/set_cancellation_sides.csv", index=False)
    pd.DataFrame(conc_rows).to_csv(f"{rerun_dir}/set_cancellation_concordance.csv", index=False)
    pd.DataFrame(gene_rows, columns=["dataset","condition","pathway","gene","coef","contribution"]) \
      .to_csv(f"{rerun_dir}/set_gene_contributions.csv", index=False)
    print(f"Saved -> set_cancellation_sides.csv ({len(S)}), "
          f"set_cancellation_concordance.csv ({len(conc_rows)}), set_gene_contributions.csv")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
