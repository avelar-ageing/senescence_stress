#!/usr/bin/env python3
"""34_ceiling_comparison.py

How reproducible is the per-gene contribution, and does that differ between the
arrest conditions and the irradiation time course once the comparison is fair?

THREE SETTINGS, all a 3-versus-3 sample split with the controls split too, so the
two halves of every estimate share no samples. Splits are enumerated exhaustively
(10 treated x 10 control per cell).

  A  arrest, pooled           BOTH sides pooled across studies: the six treated drawn
                              from all samples of that condition, the six controls from
                              all 91 proliferating samples. Hence 5 cells (one per
                              condition) against B's 6 (study x condition). The output
                              key still reads "arrest_pooled_controls" for backward
                              compatibility, but the TREATED side is pooled too - this
                              is not a controls-only difference, and describing it that
                              way in the manuscript was an error corrected 2026-09-18.
                              This is the footing 04_partial_tage_decompose.py uses.
  B  arrest, within-study      treated and controls from the SAME study. Only six
                               study x condition cells carry six of each:
                               SSCQ in ERP021140 and SRP065206, SIPS in ERP021140
                               and SRP062872, OIS in SRP062872 and SRP066947.
  C  time course, within-study every group is one study by construction
                               (ERP021140), nine cell type x timepoint cells.

WHY ALL THREE. B and C are the only pair that can be compared: same split, same
estimator, controls from one study in both. A is included because it is what the
existing decomposition rests on, and the gap between A and B measures what
pooling controls across studies costs. ERP021140 appears in both B and C, so part
of the comparison is within a single study, platform and laboratory.

WHAT THIS IS NOT. A sample split omits batch variation, so none of these numbers
is comparable with 32_set_reproducibility_ceiling.py's cross-study ceiling. That
one bounds a cross-study claim; these bound a within-study one.

STATISTICS. The ceilings are descriptive per set, so the test is across the 50
sets, paired by set: Wilcoxon signed-rank on the 50 differences. n = 50 pairs
gives a floor far below 0.05, so the comparison is testable; the floor is
reported. Raw p, no BH - three pre-specified comparisons, not a discovery scan.

Output: rerun_outputs/set_ceiling_comparison.csv
        rerun_outputs/set_ceiling_paired_tests.csv
"""
import sys, warnings, itertools, collections
import joblib, numpy as np, pandas as pd
from scipy.stats import wilcoxon
from sklearn.exceptions import InconsistentVersionWarning
warnings.filterwarnings("ignore", category=InconsistentVersionWarning)

SEED = 1
N = 6
LAB = {"Contact_inhibited CQ": "CICQ", "Serum_starved CQ": "SSCQ", "Stress-induced CS": "SIPS",
       "Oncogene-induced CS": "OIS", "Replicative CS": "RS"}


def spear(a, b):
    ra = np.argsort(np.argsort(a)).astype(float)
    rb = np.argsort(np.argsort(b)).astype(float)
    ra -= ra.mean(); rb -= rb.mean()
    d = np.sqrt((ra ** 2).sum() * (rb ** 2).sum())
    return float((ra * rb).sum() / d) if d > 0 else np.nan


def halves(n):
    out, seen = [], set()
    for c in itertools.combinations(range(n), n // 2):
        if frozenset(set(range(n)) - set(c)) in seen:
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
    HT, HC = halves(N), halves(N)

    def load(stem):
        expr = pd.read_csv(f"{PT}/{stem}_scaled_diff.csv")
        sid = expr["sample_id"].values
        expr = expr.drop(columns=["sample_id"]); expr.columns = expr.columns.map(str)
        for g in [g for g in feats if g not in expr.columns]:
            expr[g] = np.nan
        Z = m.named_steps["scaler"].transform(imp.transform(expr.loc[:, feats]))
        grp = pd.read_csv(f"{PT}/{stem}_groups.csv").set_index("sample_id").loc[sid, "group"].values
        return sid, Z * coef[np.newaxis, :], grp

    rng = np.random.default_rng(SEED)

    def accumulate(contrib, ti, ci, acc):
        cs = rng.permutation(ci)[:N]
        ts = rng.permutation(ti)[:N]
        for hT in HT:
            for hC in HC:
                A = contrib[ts[hT]].mean(axis=0) - contrib[cs[hC]].mean(axis=0)
                Bv = (contrib[ts[[i for i in range(N) if i not in hT]]].mean(axis=0)
                      - contrib[cs[[i for i in range(N) if i not in hC]]].mean(axis=0))
                for p, ix in sets.items():
                    acc[p].append(spear(A[ix], Bv[ix]))

    settings = {}

    # ---- A and B: arrest conditions --------------------------------------
    sid, contrib, grp = load("meta")
    ann = pd.read_csv(f"{rerun_dir}/immortalisation_annotation_corrected.csv").set_index("external_id")
    std = ann.reindex(sid)["study"].values
    accA = {p: [] for p in sets}
    accB = {p: [] for p in sets}
    cells_A = cells_B = 0
    ctrl_all = np.where(grp == "Proliferating")[0]
    for full, lab in LAB.items():
        ti = np.where(grp == full)[0]
        if len(ti) >= N and len(ctrl_all) >= N:
            accumulate(contrib, ti, ctrl_all, accA); cells_A += 1
        for s in pd.unique(std[grp == full]):
            t = np.where((grp == full) & (std == s))[0]
            c = np.where((grp == "Proliferating") & (std == s))[0]
            if len(t) >= N and len(c) >= N:
                accumulate(contrib, t, c, accB); cells_B += 1
    settings["arrest_pooled_controls"] = (accA, cells_A)
    settings["arrest_within_study"] = (accB, cells_B)

    # ---- C: time course, within study by construction ---------------------
    accC = {p: [] for p in sets}; cells_C = 0
    for ct in ("Fibroblast", "Keratinocyte", "Melanocyte"):
        for tp in ("4_days", "10_days", "20_days"):
            sid2, contrib2, grp2 = load(f"{ct}_{tp}")
            ti = np.where(grp2 == tp)[0]; ci = np.where(grp2 == "none")[0]
            if len(ti) >= N and len(ci) >= N:
                accumulate(contrib2, ti, ci, accC); cells_C += 1
    settings["time_course_within_study"] = (accC, cells_C)

    rows = []
    for name, (acc, cells) in settings.items():
        for p, v in acc.items():
            rows.append(dict(setting=name, pathway=p, n_genes=len(sets[p]), n_cells=cells,
                             n_estimates=len(v), ceiling_median=float(np.median(v)),
                             lo=float(np.percentile(v, 5)), hi=float(np.percentile(v, 95))))
        print(f"  {name:<28} {cells} cells, {len(acc[list(sets)[0]])} estimates per set")
    D = pd.DataFrame(rows)
    D.to_csv(f"{rerun_dir}/set_ceiling_comparison.csv", index=False)

    piv = D.pivot(index="pathway", columns="setting", values="ceiling_median")
    print(f"\n  median over the 50 sets:")
    for c in piv.columns:
        print(f"    {c:<28} {piv[c].median():.2f}  (range {piv[c].min():.2f} to {piv[c].max():.2f})")
    # the paired test was printed but never written; 2.1.5.3 quotes its p, so it is
    # persisted alongside the per-set ceilings
    wrows = []
    print(f"\n  paired Wilcoxon across the 50 sets (floor 2/2^50, so testable):")
    for a, b in itertools.combinations(piv.columns, 2):
        w = wilcoxon(piv[a], piv[b])
        print(f"    {a[:26]:<28} vs {b[:26]:<28} median diff {(piv[a]-piv[b]).median():+.2f}  p = {w.pvalue:.2e}")
        wrows.append(dict(setting_a=a, setting_b=b, n_sets=int(len(piv)),
                          median_a=float(piv[a].median()), median_b=float(piv[b].median()),
                          median_paired_diff=float((piv[a] - piv[b]).median()),
                          wilcoxon_statistic=float(w.statistic), p_value=float(w.pvalue)))


    pd.DataFrame(wrows).to_csv(f"{rerun_dir}/set_ceiling_paired_tests.csv", index=False)
    print(f"\nSaved -> set_ceiling_comparison.csv, set_ceiling_paired_tests.csv")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
