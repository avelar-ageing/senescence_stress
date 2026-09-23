#!/usr/bin/env python3
"""36_core_structure_tests.py

Three questions about the "fixed cell-cycle core" that 32/34 turned up.

Q1  IS THE CORE THE SAME IN QUIESCENCE AND SENESCENCE? The core was defined as the
    sets whose per-gene contribution pattern repeats across ALL conditions, which
    pooled quiescence and senescence together. If the core is genuinely fixed, a
    set's concordance should be no higher between two quiescence conditions, or
    between two senescences, than it is across the divide. Four within-class pairs
    (CICQ-SSCQ; SIPS-OIS, SIPS-RS, OIS-RS) against six across it, so a rank test
    has floor 2/C(10,4) = 0.0095 and is testable.

Q2  IS THE CORE'S MAGNITUDE STABLE WHILE THE REST VARIES? The fixed-core reading
    predicts that the cell-cycle sets contribute a similar amount in every
    condition while the remainder swings. Tested as the coefficient of variation
    of a set's absolute contribution across conditions, six core sets against the
    other forty-four.

Q3  THE SAME TWO QUESTIONS IN THE TIME COURSE, where "class" is cell type: nine
    within-cell-type timepoint pairs against twenty-seven across cell types. This
    also answers whether the weak temporal composition result is a sample-size
    problem. It is not obviously one: 33/34 showed the temporal per-gene
    measurement reproduces as well as the arrest conditions' (0.84 against 0.87),
    so a low between-group concordance there cannot be blamed on six samples per
    group without first checking it against that ceiling.

CORE SETS are taken from 32's output as the six with the smallest ceiling-minus-
between gap, not chosen by eye.

Raw p, no BH: three pre-specified comparisons per dataset, not a discovery scan.
Floors reported.

Output: rerun_outputs/set_core_structure_tests.csv
"""
import sys, warnings, itertools, collections
import joblib, numpy as np, pandas as pd
from scipy.stats import mannwhitneyu
from math import comb
from sklearn.exceptions import InconsistentVersionWarning
warnings.filterwarnings("ignore", category=InconsistentVersionWarning)

CQ = {"CICQ", "SSCQ"}
CS = {"SIPS", "OIS", "RS"}
LAB = {"Contact_inhibited CQ": "CICQ", "Serum_starved CQ": "SSCQ", "Stress-induced CS": "SIPS",
       "Oncogene-induced CS": "OIS", "Replicative CS": "RS"}


def spear(a, b):
    ra = np.argsort(np.argsort(a)).astype(float)
    rb = np.argsort(np.argsort(b)).astype(float)
    ra -= ra.mean(); rb -= rb.mean()
    d = np.sqrt((ra ** 2).sum() * (rb ** 2).sum())
    return float((ra * rb).sum() / d) if d > 0 else np.nan


def floor_mw(a, b):
    return 2 / comb(a + b, min(a, b))


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

    ceil = pd.read_csv(f"{rerun_dir}/set_reproducibility_ceiling.csv")
    core = set(ceil.nsmallest(6, "gap")["pathway"])
    print("core sets (smallest ceiling-minus-between gap in 32):")
    for c in sorted(core):
        print(f"  {c.replace('HALLMARK ','')}")

    def load(stem):
        expr = pd.read_csv(f"{PT}/{stem}_scaled_diff.csv")
        sid = expr["sample_id"].values
        expr = expr.drop(columns=["sample_id"]); expr.columns = expr.columns.map(str)
        for g in [g for g in feats if g not in expr.columns]:
            expr[g] = np.nan
        Z = m.named_steps["scaler"].transform(imp.transform(expr.loc[:, feats]))
        grp = pd.read_csv(f"{PT}/{stem}_groups.csv").set_index("sample_id").loc[sid, "group"].values
        return sid, Z * coef[np.newaxis, :], grp

    # ---- arrest conditions: within-study delta per condition -------------
    sid, contrib, grp = load("meta")
    ann = pd.read_csv(f"{rerun_dir}/immortalisation_annotation_corrected.csv").set_index("external_id")
    std = ann.reindex(sid)["study"].values
    delta = {}
    for full, lab in LAB.items():
        ds = []
        for s in pd.unique(std[grp == full]):
            if ((std == s) & (grp == "Proliferating")).sum() >= 1:
                ds.append(contrib[(std == s) & (grp == full)].mean(axis=0)
                          - contrib[(std == s) & (grp == "Proliferating")].mean(axis=0))
        delta[lab] = np.mean(ds, axis=0)

    # ---- time course ------------------------------------------------------
    tdelta = {}
    for ct in ("Fibroblast", "Keratinocyte", "Melanocyte"):
        for tp in ("4_days", "10_days", "20_days"):
            _, c2, g2 = load(f"{ct}_{tp}")
            tdelta[f"{ct}|{tp}"] = c2[g2 == tp].mean(axis=0) - c2[g2 == "none"].mean(axis=0)

    rows = []

    def q1(name, d, cls_of):
        wi, bw = [], []
        per_set = {}
        for p, ix in sets.items():
            w, b = [], []
            for a, bb in itertools.combinations(sorted(d), 2):
                r = spear(d[a][ix], d[bb][ix])
                (w if cls_of(a) == cls_of(bb) else b).append(r)
            per_set[p] = (float(np.median(w)), float(np.median(b)), len(w), len(b))
            wi.append(np.median(w)); bw.append(np.median(b))
        nw, nb = per_set[list(sets)[0]][2], per_set[list(sets)[0]][3]
        print(f"\n=== Q1 {name}: same genes within a class as across it? "
              f"({nw} within-class pairs, {nb} across) ===\n")
        print(f"  {'set':<34}{'within':>9}{'across':>9}{'diff':>8}   class")
        for p in sorted(sets, key=lambda z: -(per_set[z][0] - per_set[z][1]))[:5] + \
                 sorted(sets, key=lambda z: (per_set[z][0] - per_set[z][1]))[:2]:
            w, b, _, _ = per_set[p]
            print(f"  {p.replace('HALLMARK ','')[:33]:<34}{w:>9.2f}{b:>9.2f}{w-b:>8.2f}"
                  f"   {'CORE' if p in core else ''}")
            rows.append(dict(test=f"q1_{name}", pathway=p, within_class=w, across_class=b,
                             diff=w - b, is_core=p in core, n_within=nw, n_across=nb))
        for p in sets:
            w, b, _, _ = per_set[p]
            if not any(r["pathway"] == p and r["test"] == f"q1_{name}" for r in rows):
                rows.append(dict(test=f"q1_{name}", pathway=p, within_class=w, across_class=b,
                                 diff=w - b, is_core=p in core, n_within=nw, n_across=nb))
        cd = [per_set[p][0] - per_set[p][1] for p in sets if p in core]
        od = [per_set[p][0] - per_set[p][1] for p in sets if p not in core]
        u = mannwhitneyu(cd, od)
        print(f"\n  core sets' within-minus-across gap: median {np.median(cd):+.2f}")
        print(f"  other sets':                        median {np.median(od):+.2f}")
        print(f"  Mann-Whitney p = {u.pvalue:.4f} (floor {floor_mw(len(cd), len(od)):.2e})")
        return per_set

    q1("arrest", delta, lambda k: "CQ" if k in CQ else "CS")
    q1("temporal", tdelta, lambda k: k.split("|")[0])

    # ---- Q2: is the core's magnitude stable across conditions? -----------
    for name, d in (("arrest", delta), ("temporal", tdelta)):
        cv = {}
        for p, ix in sets.items():
            vals = [abs(d[k][ix].sum()) for k in d]
            cv[p] = float(np.std(vals) / np.mean(vals)) if np.mean(vals) else np.nan
        cc = [cv[p] for p in sets if p in core]; oo = [cv[p] for p in sets if p not in core]
        u = mannwhitneyu(cc, oo)
        print(f"\n=== Q2 {name}: variability of a set's contribution across groups ===\n")
        print(f"  core sets   median CV {np.median(cc):.2f}")
        print(f"  other sets  median CV {np.median(oo):.2f}")
        print(f"  Mann-Whitney p = {u.pvalue:.4f} (floor {floor_mw(len(cc), len(oo)):.2e})")
        for p in sets:
            rows.append(dict(test=f"q2_{name}", pathway=p, cv_across_groups=cv[p], is_core=p in core))

    pd.DataFrame(rows).to_csv(f"{rerun_dir}/set_core_structure_tests.csv", index=False)
    print(f"\nSaved -> set_core_structure_tests.csv")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
