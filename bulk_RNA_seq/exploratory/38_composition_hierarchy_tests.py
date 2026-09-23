#!/usr/bin/env python3
"""38_composition_hierarchy_tests.py

Three tests the composition results needed and did not have.

T1  IS THE HIERARCHY REAL? 36 reported that a set's per-gene contributions are
    more alike within an arrest class (quiescence, senescence) than across it
    (median 0.48 against 0.06), and more alike within a cell type than across
    (0.74 against 0.15), but tested neither. Paired by set, n = 50, Wilcoxon
    signed-rank.

T2  IS CELL TYPE THE RIGHT PARTITION, or would any grouping of the nine
    time-course groups into three threes do as well? There are exactly 280 such
    partitions, so the cell-type one can be ranked against all of them and the
    floor is 1/280 = 0.0036. This is the control the divergence analysis in 2.2.4
    never had - and unlike the metric used there, which random gene sets beat,
    this compares the real partition against alternative partitions of the same
    data, so gene-level response cannot supply the answer on its own.

    The arrest-condition analogue is NOT run: five conditions split two-versus-
    three gives only 10 partitions, so the best attainable p is 0.10 and the
    comparison is not testable. Where the true split ranks among the ten is
    reported descriptively instead.

T3  DOES COMPOSITION DRIFT WITH TIME? Within a cell type, adjacent timepoint
    pairs (4-10, 10-20) against the distant one (4-20): 6 against 3, rank-test
    floor 0.024. If composition moves progressively, adjacent pairs are more
    alike; if the response reorganises once and holds, they are not.

Raw p throughout - three pre-specified tests, not a scan. Floors reported.

Output: rerun_outputs/set_composition_hierarchy.csv
"""
import sys, warnings, itertools
import joblib, numpy as np, pandas as pd
from scipy.stats import wilcoxon, mannwhitneyu
from math import comb
from sklearn.exceptions import InconsistentVersionWarning
warnings.filterwarnings("ignore", category=InconsistentVersionWarning)

CQ = {"CICQ", "SSCQ"}
LAB = {"Contact_inhibited CQ": "CICQ", "Serum_starved CQ": "SSCQ", "Stress-induced CS": "SIPS",
       "Oncogene-induced CS": "OIS", "Replicative CS": "RS"}
CTS = ("Fibroblast", "Keratinocyte", "Melanocyte")
TPS = ("4_days", "10_days", "20_days")


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

    pw = pd.read_csv(f"{PT}/hallmark_pathway_mouse_ids.csv")
    pw["mouse_gene_id"] = pw["mouse_gene_id"].astype(str)
    fi = {g: i for i, g in enumerate(feats)}
    sets = {p: np.array([fi[g] for g in d["mouse_gene_id"] if g in fi]) for p, d in pw.groupby("pathway")}
    sets = {p: ix for p, ix in sets.items() if len(ix) >= 5}

    def load(stem):
        e = pd.read_csv(f"{PT}/{stem}_scaled_diff.csv")
        s = e["sample_id"].values
        e = e.drop(columns=["sample_id"]); e.columns = e.columns.map(str)
        for g in [g for g in feats if g not in e.columns]:
            e[g] = np.nan
        Z = m.named_steps["scaler"].transform(imp.transform(e.loc[:, feats]))
        g_ = pd.read_csv(f"{PT}/{stem}_groups.csv").set_index("sample_id").loc[s, "group"].values
        return s, Z * coef[np.newaxis, :], g_

    # arrest: within-study delta per condition
    sid, contrib, grp = load("meta")
    ann = pd.read_csv(f"{rerun_dir}/immortalisation_annotation_corrected.csv").set_index("external_id")
    std = ann.reindex(sid)["study"].values
    A = {}
    for full, lab in LAB.items():
        ds = [contrib[(std == s) & (grp == full)].mean(axis=0)
              - contrib[(std == s) & (grp == "Proliferating")].mean(axis=0)
              for s in pd.unique(std[grp == full])
              if ((std == s) & (grp == "Proliferating")).sum() >= 1]
        A[lab] = np.mean(ds, axis=0)
    # time course
    T = {}
    for ct in CTS:
        for tp in TPS:
            _, c2, g2 = load(f"{ct}_{tp}")
            T[(ct, tp)] = c2[g2 == tp].mean(axis=0) - c2[g2 == "none"].mean(axis=0)

    rows = []

    # ---- T1 -------------------------------------------------------------
    print("=== T1  is the hierarchy real? paired by set, n = 50 ===\n")
    for name, d, cls in (("arrest class (CQ vs CS)", A, lambda k: "CQ" if k in CQ else "CS"),
                         ("cell type", T, lambda k: k[0])):
        wi, bw = [], []
        for p, ix in sets.items():
            w = [spear(d[a][ix], d[b][ix]) for a, b in itertools.combinations(sorted(d), 2) if cls(a) == cls(b)]
            b_ = [spear(d[a][ix], d[b][ix]) for a, b in itertools.combinations(sorted(d), 2) if cls(a) != cls(b)]
            wi.append(np.median(w)); bw.append(np.median(b_))
            rows.append(dict(test="T1", grouping=name, pathway=p,
                             within=float(np.median(w)), across=float(np.median(b_))))
        st = wilcoxon(wi, bw)
        print(f"  {name:<26} within {np.median(wi):+.2f}  across {np.median(bw):+.2f}  "
              f"n = {len(wi)} sets  p = {st.pvalue:.2e}")

    # ---- T2 -------------------------------------------------------------
    keys = sorted(T)
    # the 36 pairwise correlations per set are the same for every partition, so
    # compute them once; each of the 280 partitions then only re-aggregates them
    pairs = list(itertools.combinations(range(9), 2))
    RHO = np.array([[spear(T[keys[a]][ix], T[keys[b]][ix]) for a, b in pairs]
                    for ix in sets.values()])          # sets x pairs

    def gap_for(partition):
        """median over sets of (within-group concordance - across-group)"""
        lookup = {k: gi for gi, grp_ in enumerate(partition) for k in grp_}
        same = np.array([lookup[keys[a]] == lookup[keys[b]] for a, b in pairs])
        return float(np.median(np.median(RHO[:, same], axis=1)
                               - np.median(RHO[:, ~same], axis=1)))

    parts, seen = [], set()
    for c1 in itertools.combinations(range(9), 3):
        rest = [i for i in range(9) if i not in c1]
        for c2 in itertools.combinations(rest, 3):
            c3 = tuple(i for i in rest if i not in c2)
            key = frozenset([frozenset(c1), frozenset(c2), frozenset(c3)])
            if key in seen:
                continue
            seen.add(key)
            parts.append([[keys[i] for i in c] for c in (c1, c2, c3)])
    true_part = [[(ct, tp) for tp in TPS] for ct in CTS]
    obs = gap_for(true_part)
    null = np.array([gap_for(pp) for pp in parts])
    p2 = (1 + (null >= obs).sum()) / (len(parts) + 1)
    print(f"\n=== T2  is cell type the right partition of the nine groups? ===\n")
    print(f"  {len(parts)} distinct 3x3x3 partitions enumerated; floor {1/(len(parts)+1):.4f}")
    print(f"  cell-type partition gap {obs:+.2f}; best alternative {null.max():+.2f}; "
          f"median alternative {np.median(null):+.2f}")
    print(f"  rank of the true partition: {int((null >= obs).sum())+1} of {len(parts)+1}   p = {p2:.4f}")
    rows.append(dict(test="T2", grouping="cell type vs 280 alternatives", observed=obs,
                     null_median=float(np.median(null)), null_max=float(null.max()),
                     p_emp=p2, n_partitions=len(parts)))

    # arrest analogue: describable only
    ak = sorted(A)
    aparts = [c for c in itertools.combinations(ak, 2)]
    def agap(pair):
        out = []
        for p, ix in sets.items():
            w, b = [], []
            for a, bb in itertools.combinations(ak, 2):
                r = spear(A[a][ix], A[bb][ix])
                same = (a in pair) == (bb in pair)
                (w if same else b).append(r)
            out.append(np.median(w) - np.median(b))
        return float(np.median(out))
    ag = {p: agap(p) for p in aparts}
    tr = tuple(sorted(CQ))
    rank = sorted(ag.values(), reverse=True).index(ag[tr]) + 1
    print(f"\n  arrest analogue (NOT testable, only {len(aparts)} partitions, floor "
          f"{1/(len(aparts)+1):.2f}):")
    print(f"    CQ/CS gap {ag[tr]:+.2f} ranks {rank} of {len(aparts)}; "
          f"best alternative {max(ag.values()):+.2f} ({max(ag, key=ag.get)})")

    # ---- T3 -------------------------------------------------------------
    adj, dist = [], []
    for p, ix in sets.items():
        for ct in CTS:
            adj.append(spear(T[(ct, "4_days")][ix], T[(ct, "10_days")][ix]))
            adj.append(spear(T[(ct, "10_days")][ix], T[(ct, "20_days")][ix]))
            dist.append(spear(T[(ct, "4_days")][ix], T[(ct, "20_days")][ix]))
    u = mannwhitneyu(adj, dist)
    print(f"\n=== T3  does composition drift with time? ===\n")
    print(f"  adjacent timepoint pairs  median {np.median(adj):+.2f}  (n = {len(adj)})")
    print(f"  distant pair (4 v 20 d)   median {np.median(dist):+.2f}  (n = {len(dist)})")
    print(f"  Mann-Whitney p = {u.pvalue:.4f}  (floor {2/comb(6,3):.3f} at one set; "
          f"pooled over 50 sets here)")
    rows.append(dict(test="T3", grouping="adjacent vs distant timepoints",
                     observed=float(np.median(adj)), null_median=float(np.median(dist)),
                     p_emp=float(u.pvalue), n_partitions=len(adj)))

    pd.DataFrame(rows).to_csv(f"{rerun_dir}/set_composition_hierarchy.csv", index=False)
    print(f"\nSaved -> set_composition_hierarchy.csv")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
