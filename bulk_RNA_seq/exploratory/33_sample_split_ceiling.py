#!/usr/bin/env python3
"""33_sample_split_ceiling.py

Does the per-gene contribution reproduce in the TIME COURSE, and is the temporal
concordance therefore worth reporting at all?

WHY THIS EXISTS SEPARATELY FROM 32. 32_set_reproducibility_ceiling.py splits
STUDIES, so its ceiling contains batch noise and is the right bound for a
cross-study comparison. The time course is a single study (ERP021140), so that
split does not exist there. The only replicate structure available is the six
samples per cell type x timepoint, split three against three, with the untreated
samples split too so the two halves share nothing.

WHAT THIS IS AND IS NOT. A sample split measures SAMPLING noise only. It omits
the batch and protocol variation that a study split includes, so it is an
optimistic bound and its value is NOT comparable with 32's. It is computed here
for BOTH datasets on the same footing precisely so the two are comparable with
each other, which is the only comparison it supports.

WHAT IT DECIDES. If the time course reproduces about as well as the arrest
conditions under an identical sample split, its weaker between-condition
concordance (E2F 0.37 against 0.83) reflects the comparison and not the
measurement, and the temporal concordance can be reported. If it reproduces
badly, the temporal concordance should be dropped rather than qualified.

Splits are enumerated exhaustively: C(6,3)/2 = 10 per group. Every group in both
datasets has exactly six treated samples, and controls are split to match.

Output: rerun_outputs/set_sample_split_ceiling.csv
"""
import sys, warnings, itertools
import joblib, numpy as np, pandas as pd
from sklearn.exceptions import InconsistentVersionWarning
warnings.filterwarnings("ignore", category=InconsistentVersionWarning)

SEED = 1
MIN_N = 6   # a 3 v 3 split needs six of each


def spear(a, b):
    ra = np.argsort(np.argsort(a)).astype(float)
    rb = np.argsort(np.argsort(b)).astype(float)
    ra -= ra.mean(); rb -= rb.mean()
    d = np.sqrt((ra ** 2).sum() * (rb ** 2).sum())
    return float((ra * rb).sum() / d) if d > 0 else np.nan


def halves(n):
    """all distinct 2-way splits of n items, complement de-duplicated"""
    out, seen = [], set()
    for c in itertools.combinations(range(n), n // 2):
        comp = frozenset(set(range(n)) - set(c))
        if comp in seen:
            continue
        seen.add(frozenset(c)); out.append(list(c))
    return out


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

    STEMS = [("arrest_conditions", ["meta"], "Proliferating"),
             ("time_course", [f"{ct}_{tp}" for ct in ("Fibroblast", "Keratinocyte", "Melanocyte")
                              for tp in ("4_days", "10_days", "20_days")], "none")]
    rng = np.random.default_rng(SEED)
    rows = []
    for dset, stems, control in STEMS:
        acc = {p: [] for p in sets}
        ngroups = 0
        for stem in stems:
            expr = pd.read_csv(f"{PT}/{stem}_scaled_diff.csv")
            sid = expr["sample_id"].values
            expr = expr.drop(columns=["sample_id"]); expr.columns = expr.columns.map(str)
            for g in [g for g in feats if g not in expr.columns]:
                expr[g] = np.nan
            Z = m.named_steps["scaler"].transform(imp.transform(expr.loc[:, feats]))
            contrib = Z * coef[np.newaxis, :]
            grp = pd.read_csv(f"{PT}/{stem}_groups.csv").set_index("sample_id").loc[sid, "group"].values
            ci = np.where(grp == control)[0]
            for c in [x for x in pd.unique(grp) if x != control]:
                ti = np.where(grp == c)[0]
                if len(ti) < MIN_N or len(ci) < MIN_N:
                    continue
                ngroups += 1
                # controls are subsampled to six then split, so both halves are
                # independent of each other in treated AND untreated samples
                cs = rng.permutation(ci)[:MIN_N]
                for hT in halves(MIN_N):
                    for hC in halves(MIN_N):
                        tA, tB = ti[hT], ti[[i for i in range(MIN_N) if i not in hT]]
                        cA, cB = cs[hC], cs[[i for i in range(MIN_N) if i not in hC]]
                        A = contrib[tA].mean(axis=0) - contrib[cA].mean(axis=0)
                        Bv = contrib[tB].mean(axis=0) - contrib[cB].mean(axis=0)
                        for p, ix in sets.items():
                            acc[p].append(spear(A[ix], Bv[ix]))
        for p, v in acc.items():
            rows.append(dict(dataset=dset, pathway=p, n_genes=len(sets[p]),
                             n_groups=ngroups, n_estimates=len(v),
                             ceiling_median=float(np.median(v)),
                             lo=float(np.percentile(v, 5)), hi=float(np.percentile(v, 95))))
        print(f"  {dset}: {ngroups} groups with >= {MIN_N} treated and {MIN_N} control, "
              f"{len(acc[list(sets)[0]])} estimates per set")
    D = pd.DataFrame(rows)
    D.to_csv(f"{rerun_dir}/set_sample_split_ceiling.csv", index=False)
    print(f"\nSaved -> set_sample_split_ceiling.csv ({len(D)} rows)\n")
    piv = D.pivot(index="pathway", columns="dataset", values="ceiling_median")
    print(f"  median sample-split ceiling over the 50 sets:")
    for d in piv.columns:
        print(f"    {d:<20} {piv[d].median():.2f}   (range {piv[d].min():.2f} to {piv[d].max():.2f})")
    print(f"\n  {'set':<32}{'arrest':>9}{'time course':>14}")
    for p in piv.sort_values(piv.columns[0], ascending=False).index[:6]:
        print(f"  {p.replace('HALLMARK ','')[:31]:<32}{piv.loc[p].iloc[0]:>9.2f}{piv.loc[p].iloc[1]:>14.2f}")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
