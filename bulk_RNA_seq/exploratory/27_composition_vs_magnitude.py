#!/usr/bin/env python
"""27_composition_vs_magnitude.py

Tests the claim in 2.1.5.2 that how much a condition raises the mortality score tells
us little about which gene sets produce the rise.

THE CLAIM, AS A STATISTIC. For each of the 10 condition pairs, take the difference in
the size of their within-study shift, and the Spearman correlation between their
50-set contribution profiles. If magnitude predicted composition, pairs far apart in
size would be less alike in composition, i.e. those two quantities would correlate
negatively across the 10 pairs. The statistic is that correlation.

THE NULL. The 10 pairs are not independent (each condition appears in four of them),
so a nominal p on a correlation over 10 points is wrong. Instead the assignment of
profiles to conditions is permuted: the five magnitudes stay put and the five profile
vectors are shuffled between them, recomputing the statistic each time. All 5! = 120
assignments are enumerated, so the p is exact and its floor is 1/120 = 0.0083.

WHAT A NULL RESULT WOULD AND WOULD NOT SHOW. A weak, non-extreme statistic is
consistent with magnitude and composition being unrelated, but failing to reject is
not proof of independence - with five conditions the test has little power. The
observed statistic is therefore reported as an effect size alongside the p, and the
text should claim no more than the pair supports.

SECOND TEST, on whether the profile correlations mean anything at all. A rho between
two conditions is computed over 50 Hallmark sets that share genes, so the sets are not
independent and no nominal p applies. The set-shuffle null asks the answerable
question: reassign genes to sets at random, preserving set sizes, recompute both
conditions' 50 contributions from the same per-gene contributions, and recompute rho.
This gives the profile correlation expected from gene-level agreement alone, with no
real pathway structure. Observed rho above that null means Hallmark grouping
concentrates the agreement; at or below it means the correlation is carried by
gene-level similarity that any grouping would reproduce.

Usage: 27_composition_vs_magnitude.py <rerun_dir> <model_dir> [n_shuffle]
Output: <rerun_dir>/composition_vs_magnitude.csv
"""
import itertools
import json
import os
import sys
import warnings

import joblib
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

warnings.filterwarnings("ignore")
SEED = 1
CONDS = ["CICQ", "SSCQ", "RS", "SIPS", "OIS"]
MORT = "EN_Mortality_Multispecies_Multitissue_scaleddiff.pkl"


def main(rerun_dir, model_dir, n_shuffle=10000):
    rng = np.random.default_rng(SEED)

    # ---- profiles and magnitudes ------------------------------------------
    d = pd.read_csv(f"{rerun_dir}/mortality_partial_tage_ALL.csv")
    d = d[d.analysis == "meta_analysis"]
    P = d.pivot_table(index="pathway", columns="label",
                      values="contrib_diff_within_study")[CONDS]
    mw = pd.read_csv(f"{rerun_dir}/mortality_within_study.csv")
    mw = mw[mw.test == "condition_within_study"].set_index("condition")
    mag = mw.loc[CONDS, "diff_within_study"]
    print(f"{P.shape[0]} gene sets x {P.shape[1]} conditions")
    print("magnitudes:", {c: round(float(mag[c]), 3) for c in CONDS})

    pairs = list(itertools.combinations(range(5), 2))

    def stat(order):
        """Spearman between pairwise magnitude gap and pairwise profile rho."""
        gaps, rhos = [], []
        for i, j in pairs:
            gaps.append(abs(mag.iloc[i] - mag.iloc[j]))
            rhos.append(spearmanr(P.iloc[:, order[i]], P.iloc[:, order[j]]).statistic)
        return spearmanr(gaps, rhos).statistic, gaps, rhos

    obs, gaps, rhos = stat(list(range(5)))
    print("\n=== TEST 1: does magnitude predict composition? ===")
    for (i, j), g, r in zip(pairs, gaps, rhos):
        fold = max(mag.iloc[i], mag.iloc[j]) / min(mag.iloc[i], mag.iloc[j])
        print(f"  {CONDS[i]:>4} vs {CONDS[j]:<5} gap {g:5.3f} ({fold:4.1f}x)   profile rho {r:5.2f}")
    null = np.array([stat(list(o))[0] for o in itertools.permutations(range(5))])
    p_two = (np.abs(null) >= abs(obs)).mean()
    print(f"\n  observed Spearman(gap, rho) = {obs:+.3f}")
    print(f"  exact p over all {len(null)} assignments = {p_two:.4f} "
          f"(floor {1/len(null):.4f})")
    print("  negative would mean bigger magnitude gaps go with less similar composition")

    rows = [dict(test="magnitude_predicts_composition", statistic="spearman_gap_vs_rho",
                 observed=obs, p_exact=p_two, n_permutations=len(null),
                 p_floor=1 / len(null))]

    # ---- TEST 2: set-shuffle null for the two headline pairs --------------
    print("\n=== TEST 2: is the profile correlation more than gene-level agreement? ===")
    m = joblib.load(os.path.join(model_dir, MORT))
    imp = m.named_steps["imputation"]           # sklearn version shim, as in script 18
    if not hasattr(imp, "_fill_dtype"):
        imp._fill_dtype = imp.statistics_.dtype if hasattr(imp, "statistics_") else np.float64
    feats = list(map(str, m.feature_names_in_))
    pw = pd.read_csv(f"{rerun_dir}/partial_tage/hallmark_pathway_mouse_ids.csv")
    sets = {n: [str(x) for x in g.mouse_gene_id] for n, g in pw.groupby("pathway")}
    idx = {f: i for i, f in enumerate(feats)}
    setidx = {k: np.array([idx[g] for g in v if g in idx]) for k, v in sets.items()}
    setidx = {k: v for k, v in setidx.items() if len(v) > 0}
    print(f"  {len(setidx)} sets mapped onto the clock's {len(feats)} features")

    # per-gene contributions per condition, from the same matrices the
    # decomposition uses
    PT = f"{rerun_dir}/partial_tage"
    e = pd.read_csv(f"{PT}/meta_scaled_diff.csv")
    sid = e["sample_id"].values
    e = e.drop(columns=["sample_id"]); e.columns = e.columns.map(str)
    for g in [g for g in feats if g not in e.columns]:
        e[g] = np.nan
    Z = m.named_steps["scaler"].transform(
        m.named_steps["imputation"].transform(e.loc[:, feats]))
    C = Z * m.named_steps["estimator"].coef_[None, :]
    grp = pd.read_csv(f"{PT}/meta_groups.csv").set_index("sample_id")["group"]
    lab = np.array([grp.get(s) for s in sid])
    SUB = {"Contact_inhibited CQ": "CICQ", "Serum_starved CQ": "SSCQ",
           "Replicative CS": "RS", "Stress-induced CS": "SIPS",
           "Oncogene-induced CS": "OIS"}
    cond = np.array([SUB.get(x, "Proliferating") for x in lab])
    meta_md = pd.read_csv(f"{rerun_dir}/sample_metadata_RERUN.csv").set_index("external_id")
    study = np.array([meta_md.loc[x, "study"] for x in sid])

    # PER-GENE within-study contribution difference, the SAME estimator as
    # meta_analysis/13 and exploratory/22: per study, mean(condition) minus
    # mean(control), combined with weight n_t*n_c/(n_t+n_c). Using pooled means
    # here instead would give rho values that are not the ones the text reports.
    def gene_delta_within_study(c):
        num = np.zeros(C.shape[1]); den = 0.0
        for st in np.unique(study):
            m1 = (study == st) & (cond == c)
            m0 = (study == st) & (cond == "Proliferating")
            if m1.sum() and m0.sum():
                w = m1.sum() * m0.sum() / (m1.sum() + m0.sum())
                num += w * (C[m1].mean(axis=0) - C[m0].mean(axis=0)); den += w
        return num / den
    gene_delta = {c: gene_delta_within_study(c) for c in CONDS}

    # check the reconstruction against the reported profile correlations
    chk = {c: np.array([gene_delta[c][v].sum() for v in setidx.values()]) for c in CONDS}
    names = list(setidx)
    rep = P.reindex(names)
    for a, b in [("CICQ", "SSCQ"), ("CICQ", "OIS")]:
        r_new = spearmanr(chk[a], chk[b]).statistic
        r_rep = spearmanr(rep[a], rep[b]).statistic
        print(f"  reconstruction check {a}-{b}: rho {r_new:+.3f} vs "
              f"{r_rep:+.3f} from the decomposition CSV")

    sizes = [len(v) for v in setidx.values()]
    for a, b in [("CICQ", "SSCQ"), ("CICQ", "OIS")]:
        ga, gb = gene_delta[a], gene_delta[b]
        pa, pb = chk[a], chk[b]
        obs_r = spearmanr(pa, pb).statistic
        nullr = np.empty(n_shuffle)
        n_feat = len(feats)
        for k in range(n_shuffle):
            perm = rng.permutation(n_feat)
            cut = np.cumsum(sizes)[:-1]
            groups = np.split(perm[:sum(sizes)], cut)
            nullr[k] = spearmanr([ga[g].sum() for g in groups],
                                 [gb[g].sum() for g in groups]).statistic
        p = (1 + (nullr >= obs_r).sum()) / (n_shuffle + 1)
        print(f"  {a} vs {b}: observed rho {obs_r:+.3f}   random-set null "
              f"{nullr.mean():+.3f} (sd {nullr.std():.3f})   p(obs > null) = {p:.4f}")
        rows.append(dict(test="set_shuffle_profile_rho", statistic=f"{a}_vs_{b}",
                         observed=obs_r, null_mean=float(nullr.mean()),
                         null_sd=float(nullr.std()), p_exact=p,
                         n_permutations=n_shuffle, p_floor=1 / (n_shuffle + 1)))

    # ---- TEST 3: the same null on the temporal arm ------------------------
    # 2.2.4 reports that profiles from one cell type at different timepoints
    # correlate at rho = 0.80 against 0.35 between cell types, with a p from
    # shuffling cell-type labels. That null asks whether the grouping is real; it
    # does not ask whether the correlations reflect PATHWAY structure. This does:
    # reassign genes to sets at random, preserving set sizes, and recompute the
    # same within-minus-between gap.
    print("\n=== TEST 3: temporal within- vs between-cell-type gap, set-shuffle null ===")
    CT = ["Fibroblast", "Keratinocyte", "Melanocyte"]
    TP = ["4_days", "10_days", "20_days"]
    # the by-timepoint export is one matrix per cell type PER timepoint
    # (exploratory/02), each holding that timepoint against its own untreated
    # baseline. The pooled {ct}_ files label groups irradiated/none, not by day.
    gdelta = {}
    for ct in CT:
        for tp in TP:
            stem = f"{PT}/{ct}_{tp}"
            et = pd.read_csv(f"{stem}_scaled_diff.csv")
            sid2 = et["sample_id"].values
            et = et.drop(columns=["sample_id"]); et.columns = et.columns.map(str)
            for g in [g for g in feats if g not in et.columns]:
                et[g] = np.nan
            Z2 = m.named_steps["scaler"].transform(
                m.named_steps["imputation"].transform(et.loc[:, feats]))
            C2 = Z2 * m.named_steps["estimator"].coef_[None, :]
            g2 = pd.read_csv(f"{stem}_groups.csv").set_index("sample_id")["group"]
            lab2 = np.array([g2.get(x) for x in sid2])
            gdelta[(ct, tp)] = C2[lab2 == tp].mean(axis=0) - C2[lab2 == "none"].mean(axis=0)
    keys = list(gdelta)
    print(f"  {len(keys)} cell-type/timepoint groups")

    # vectorised: build all 9 profiles, rank them once, and take the whole
    # correlation matrix in one call. Looping scipy.spearmanr over 36 pairs per
    # shuffle is ~40x slower and made 5,000 draws intractable.
    from scipy.stats import rankdata
    same_ct = np.array([[a[0] == b[0] for b in keys] for a in keys])
    triu = np.triu(np.ones((len(keys), len(keys)), bool), 1)
    wi_mask = triu & same_ct
    bw_mask = triu & ~same_ct

    def gap_from(groups_idx):
        prof = np.array([[gdelta[k][g].sum() for g in groups_idx] for k in keys])
        R = np.corrcoef(np.apply_along_axis(rankdata, 1, prof))
        return R[wi_mask].mean() - R[bw_mask].mean(), R[wi_mask].mean(), R[bw_mask].mean()

    real_groups = [v for v in setidx.values()]
    obs_gap, obs_wi, obs_bw = gap_from(real_groups)
    print(f"  observed: within {obs_wi:.3f}, between {obs_bw:.3f}, gap {obs_gap:+.3f}")
    nullg = np.empty(n_shuffle); n_feat2 = len(feats)
    cut = np.cumsum(sizes)[:-1]
    for k in range(n_shuffle):
        perm = rng.permutation(n_feat2)
        nullg[k] = gap_from(np.split(perm[:sum(sizes)], cut))[0]
    p3 = (1 + (nullg >= obs_gap).sum()) / (n_shuffle + 1)
    print(f"  random-set null gap: {nullg.mean():+.3f} (sd {nullg.std():.3f})")
    print(f"  p(observed gap > random-set gap) = {p3:.4f}")
    rows.append(dict(test="temporal_set_shuffle_gap", statistic="within_minus_between_celltype",
                     observed=obs_gap, null_mean=float(nullg.mean()),
                     null_sd=float(nullg.std()), p_exact=p3,
                     n_permutations=n_shuffle, p_floor=1 / (n_shuffle + 1)))

    # ---- TESTS 4 and 5: the other composition claims in 2.2.4 --------------
    # Same principle as test 3. Any statistic computed from per-set contributions
    # can be reproduced by random sets if it really reflects gene-level structure,
    # so each of these claims needs the set-shuffle null before it can be called
    # a pathway result.
    print("\n=== TEST 4: are the leading axes of variation between cell types? ===")

    def pc_between_frac(groups_idx):
        prof = np.array([[gdelta[k][g].sum() for g in groups_idx] for k in keys])
        X = (prof - prof.mean(0)) / (prof.std(0) + 1e-12)
        U, S, Vt = np.linalg.svd(X - X.mean(0), full_matrices=False)
        sc = U[:, :2] * S[:2]
        ct = np.array([k[0] for k in keys])
        fr = []
        for j in range(2):
            y = sc[:, j]; gm = y.mean()
            ssb = sum(((y[ct == c].mean() - gm) ** 2) * (ct == c).sum() for c in np.unique(ct))
            fr.append(ssb / max(((y - gm) ** 2).sum(), 1e-12))
        return float(np.mean(fr))

    obs4 = pc_between_frac(real_groups)
    null4 = np.empty(n_shuffle)
    for k in range(n_shuffle):
        perm = rng.permutation(len(feats))
        null4[k] = pc_between_frac(np.split(perm[:sum(sizes)], cut))
    p4 = (1 + (null4 >= obs4).sum()) / (n_shuffle + 1)
    print(f"  observed between-cell-type variance fraction of PC1-2: {obs4:.3f}")
    print(f"  random-set null: {null4.mean():.3f} (sd {null4.std():.3f})   p = {p4:.4f}")
    rows.append(dict(test="temporal_set_shuffle_pc_between_frac",
                     statistic="mean_between_celltype_ss_frac_PC1_PC2",
                     observed=obs4, null_mean=float(null4.mean()),
                     null_sd=float(null4.std()), p_exact=p4,
                     n_permutations=n_shuffle, p_floor=1 / (n_shuffle + 1)))

    # TEST 5 uses the POOLED temporal contributions - irradiated against untreated
    # per cell type, one value per cell type - because that is what the 31-of-50
    # and melanocyte-20 figures in the text are computed from. The by-timepoint
    # averages used above are a different statistic (27 and 19) and testing the
    # wrong one is how the first version of this script went astray.
    print("\n=== TEST 5/6: sign changes between cell types, and their concentration ===")
    gpool = {}
    for ct in CT:
        et = pd.read_csv(f"{PT}/{ct}_scaled_diff.csv")
        sid3 = et["sample_id"].values
        et = et.drop(columns=["sample_id"]); et.columns = et.columns.map(str)
        for g in [g for g in feats if g not in et.columns]:
            et[g] = np.nan
        Z3 = m.named_steps["scaler"].transform(
            m.named_steps["imputation"].transform(et.loc[:, feats]))
        C3 = Z3 * m.named_steps["estimator"].coef_[None, :]
        g3 = pd.read_csv(f"{PT}/{ct}_groups.csv").set_index("sample_id")["group"]
        lab3 = np.array([g3.get(x) for x in sid3])
        gpool[ct] = C3[lab3 == "irradiated"].mean(axis=0) - C3[lab3 == "none"].mean(axis=0)

    def flips_and_dissent(groups_idx):
        M = np.array([[gpool[c][g].sum() for g in groups_idx] for c in CT])
        flip = (M > 0).any(axis=0) & (M < 0).any(axis=0)
        dis = np.zeros(3, int)
        for j in np.where(flip)[0]:
            sg = np.sign(M[:, j])
            for i in range(3):
                if sg[i] != 0 and all(sg[k] == -sg[i] for k in range(3) if k != i):
                    dis[i] += 1
        return int(flip.sum()), dis

    obs5, obs_dis = flips_and_dissent(real_groups)
    print(f"  observed: {obs5} of {len(sizes)} sets change sign; lone dissent "
          f"{dict(zip(CT, obs_dis))}")
    null5 = np.empty(n_shuffle); null_max = np.empty(n_shuffle); null_mel = np.empty(n_shuffle)
    for k in range(n_shuffle):
        perm = rng.permutation(len(feats))
        gi = np.split(perm[:sum(sizes)], cut)
        fl, dd = flips_and_dissent(gi)
        null5[k] = fl; null_max[k] = dd.max(); null_mel[k] = dd[CT.index("Melanocyte")]
    p5 = (1 + (null5 >= obs5).sum()) / (n_shuffle + 1)
    mel = int(obs_dis[CT.index("Melanocyte")])
    p6 = (1 + (null_mel >= mel).sum()) / (n_shuffle + 1)
    p6b = (1 + (null_max >= mel).sum()) / (n_shuffle + 1)
    print(f"  number of sign changes: null {null5.mean():.1f} (sd {null5.std():.1f})"
          f"   p(obs > null) = {p5:.4f}")
    print(f"  melanocyte lone-dissent count {mel}: null {null_mel.mean():.1f} "
          f"(sd {null_mel.std():.1f})   p = {p6:.4f}")
    print(f"    same, against the LARGEST dissent count of any cell type under the null "
          f"(guards against picking melanocytes because they are the extreme): p = {p6b:.4f}")
    for nm, ob, nl, pv in [("n_sets_changing_sign", obs5, null5, p5),
                           ("melanocyte_lone_dissent", mel, null_mel, p6),
                           ("melanocyte_vs_max_dissent", mel, null_max, p6b)]:
        rows.append(dict(test="temporal_set_shuffle_signflips_pooled", statistic=nm,
                         observed=ob, null_mean=float(nl.mean()), null_sd=float(nl.std()),
                         p_exact=pv, n_permutations=n_shuffle, p_floor=1 / (n_shuffle + 1)))

    out = f"{rerun_dir}/composition_vs_magnitude.csv"
    pd.DataFrame(rows).to_csv(out, index=False)
    print(f"\nSaved -> {out}")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2],
         int(sys.argv[3]) if len(sys.argv) > 3 else 10000)
