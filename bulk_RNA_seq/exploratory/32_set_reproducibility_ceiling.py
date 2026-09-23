#!/usr/bin/env python3
"""32_set_reproducibility_ceiling.py

Which gene sets can we say anything about at GENE level, and for which is the
per-gene attribution too noisy to interpret?

THE PROBLEM. 30_set_concordance_null.py found that only E2F TARGETS and G2M
CHECKPOINT have per-gene contribution patterns that repeat across conditions
(rho 0.83 and 0.74 against a 0.13 null). P53 PATHWAY sits at 0.31, significant in
1 of 10 pairs. Two opposite readings give that same number:

  (a) the same set is engaged in every condition but through DIFFERENT genes -
      which is what a convergent stress network would predict, the pathway being
      a hub that distinct stressors enter through distinct effectors; or
  (b) the per-gene measurement is not reproducible at all, in which case the
      between-condition correlation is measuring noise and neither reading holds.

THE TEST. A reproducibility ceiling. Split each condition's STUDIES into two
halves, estimate the per-gene contributions independently in each, and correlate
them within each set. That number is how well the measurement reproduces when the
biology is held constant and only the samples change, so it bounds what any
between-condition correlation can mean:

  ceiling low                      -> gene-level attribution is noise; exclude the
                                      set from any gene-level claim
  ceiling high, between-cond high  -> stereotyped: the same genes carry the set
                                      in every condition
  ceiling high, between-cond low   -> condition-specific effectors: the set
                                      recurs, the genes carrying it do not

ESTIMATOR. Each half uses the WITHIN-STUDY difference - the mean over its studies
of (treated mean - that study's own control mean) - so the two halves share no
control samples and batch noise sits inside the ceiling rather than beside it.
This differs from 04_partial_tage_decompose.py, which references all 91
proliferating samples; the between-condition correlation is recomputed here on
the same within-study footing so the two numbers are comparable.

Only studies carrying both treated and proliferating samples can contribute
(CICQ 6/7, SSCQ 6/6, SIPS 4/5, OIS 15/16, RS 6/6).

Output: rerun_outputs/set_reproducibility_ceiling.csv
"""
import sys, warnings, itertools
import joblib, numpy as np, pandas as pd
from sklearn.exceptions import InconsistentVersionWarning
warnings.filterwarnings("ignore", category=InconsistentVersionWarning)

# EXHAUSTIVE, not sampled. Every condition has few enough usable studies that all
# two-way splits can be enumerated: 3 for SIPS (4 studies), 10 each for CICQ, SSCQ
# and RS (6), and 6,435 for OIS (15 studies split 7 v 8, so no split is its own
# complement). Drawing 200 random partitions, as the first
# version did, re-drew the same 3 or 10 splits repeatedly and reported a "median
# over 200 partitions" that implied resolution it did not have. Enumerating gives
# the exact median and lets the split count and the spread be reported alongside,
# so a ceiling resting on three splits is visibly weaker than one resting on
# thousands.
MAX_ENUM = 20000
LAB = {"Contact_inhibited CQ": "CICQ", "Serum_starved CQ": "SSCQ", "Stress-induced CS": "SIPS",
       "Oncogene-induced CS": "OIS", "Replicative CS": "RS"}


def spear(a, b):
    ra = np.argsort(np.argsort(a)).astype(float)
    rb = np.argsort(np.argsort(b)).astype(float)
    ra -= ra.mean(); rb -= rb.mean()
    d = np.sqrt((ra ** 2).sum() * (rb ** 2).sum())
    return float((ra * rb).sum() / d) if d > 0 else np.nan


def main(rerun_dir, model_dir):
    PT = f"{rerun_dir}/partial_tage"
    m = joblib.load(f"{model_dir}/EN_Mortality_Multispecies_Multitissue_scaleddiff.pkl")
    imp = m.named_steps["imputation"]
    if not hasattr(imp, "_fill_dtype"):
        imp._fill_dtype = imp.statistics_.dtype
    feats = list(map(str, m.feature_names_in_))
    coef = m.named_steps["estimator"].coef_

    expr = pd.read_csv(f"{PT}/meta_scaled_diff.csv")
    sid = expr["sample_id"].values
    expr = expr.drop(columns=["sample_id"]); expr.columns = expr.columns.map(str)
    for g in [g for g in feats if g not in expr.columns]:
        expr[g] = np.nan
    Z = m.named_steps["scaler"].transform(imp.transform(expr.loc[:, feats]))
    contrib = Z * coef[np.newaxis, :]

    ann = pd.read_csv(f"{rerun_dir}/immortalisation_annotation_corrected.csv").set_index("external_id")
    sub = ann.reindex(sid)["cell_substate"].values
    std = ann.reindex(sid)["study"].values

    pw = pd.read_csv(f"{PT}/hallmark_pathway_mouse_ids.csv")
    pw["mouse_gene_id"] = pw["mouse_gene_id"].astype(str)
    fi = {g: i for i, g in enumerate(feats)}
    sets = {p: np.array([fi[g] for g in d["mouse_gene_id"] if g in fi]) for p, d in pw.groupby("pathway")}
    sets = {p: ix for p, ix in sets.items() if len(ix) >= 5}

    # studies that carry both the condition and its own proliferating controls
    usable = {}
    for full, lab in LAB.items():
        ss = []
        for s in pd.unique(std[sub == full]):
            if ((std == s) & (sub == "Proliferating")).sum() >= 1:
                ss.append(s)
        usable[lab] = ss

    def within_study(lab, studies):
        """mean over the given studies of (treated mean - that study's control mean)"""
        full = [k for k, v in LAB.items() if v == lab][0]
        d = [contrib[(std == s) & (sub == full)].mean(axis=0)
             - contrib[(std == s) & (sub == "Proliferating")].mean(axis=0) for s in studies]
        return np.mean(d, axis=0) if d else None

    ceil = {lab: {p: [] for p in sets} for lab in LAB.values()}
    n_splits = {}
    for lab, ss in usable.items():
        if len(ss) < 4:
            print(f"  {lab}: only {len(ss)} usable studies - ceiling not estimable")
            n_splits[lab] = 0
            continue
        half = len(ss) // 2
        # all two-way splits, de-duplicated (a split and its complement are one)
        splits, seen = [], set()
        for comb_ in itertools.combinations(range(len(ss)), half):
            key = frozenset(comb_)
            comp = frozenset(set(range(len(ss))) - key)
            if comp in seen:
                continue
            seen.add(key); splits.append(comb_)
            if len(splits) >= MAX_ENUM:
                break
        n_splits[lab] = len(splits)
        for comb_ in splits:
            A = within_study(lab, [ss[i] for i in comb_])
            Bv = within_study(lab, [ss[i] for i in range(len(ss)) if i not in comb_])
            if A is None or Bv is None:
                continue
            for p, ix in sets.items():
                ceil[lab][p].append(spear(A[ix], Bv[ix]))
        print(f"  {lab}: {len(ss)} usable studies, {len(splits)} distinct splits (all enumerated)")

    # between-condition correlation, on the SAME within-study estimator
    full_d = {lab: within_study(lab, ss) for lab, ss in usable.items() if len(ss) >= 1}
    rows = []
    for p, ix in sets.items():
        btw = [spear(full_d[a][ix], full_d[b][ix])
               for a, b in itertools.combinations(sorted(full_d), 2)]
        per_cond = {lab: (float(np.median(v)) if v else np.nan) for lab, v in
                    ((l, ceil[l][p]) for l in ceil)}
        vals = [v for v in per_cond.values() if not np.isnan(v)]
        # spread of the ceiling ACROSS SPLITS within the best-resolved condition,
        # so the reader can see how firm the ceiling is
        best = max(n_splits, key=lambda k: n_splits[k])
        sp = ceil[best][p]
        rows.append(dict(pathway=p, n_genes=len(ix),
                         ceiling_median=float(np.median(vals)) if vals else np.nan,
                         ceiling_min=float(np.min(vals)) if vals else np.nan,
                         ceiling_max=float(np.max(vals)) if vals else np.nan,
                         ceiling_spread_lo=float(np.percentile(sp, 5)) if sp else np.nan,
                         ceiling_spread_hi=float(np.percentile(sp, 95)) if sp else np.nan,
                         spread_from=best, n_splits_used=n_splits[best],
                         between_median=float(np.median(btw)),
                         **{f"ceiling_{k}": v for k, v in per_cond.items()}))
    D = pd.DataFrame(rows)
    # gap: how far the between-condition correlation falls below what the data can support
    D["gap"] = D["ceiling_median"] - D["between_median"]
    D["verdict"] = np.where(D["ceiling_median"] < 0.30, "noise-limited",
                    np.where(D["between_median"] >= 0.6 * D["ceiling_median"], "stereotyped",
                             "condition-specific"))
    D = D.sort_values("ceiling_median", ascending=False)
    D.to_csv(f"{rerun_dir}/set_reproducibility_ceiling.csv", index=False)
    print(f"\nSaved -> set_reproducibility_ceiling.csv ({len(D)} sets)")
    print(D["verdict"].value_counts().to_string())
    print("\n  {:<34}{:>7}{:>10}{:>10}{:>7}  {}".format(
        "set", "n", "ceiling", "between", "gap", "verdict"))
    for _, r in D.iterrows():
        print("  {:<34}{:>7}{:>10.2f}{:>10.2f}{:>7.2f}  {}".format(
            r["pathway"].replace("HALLMARK ", "")[:33], int(r["n_genes"]),
            r["ceiling_median"], r["between_median"], r["gap"], r["verdict"]))


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
