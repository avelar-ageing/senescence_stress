#!/usr/bin/env python
"""23_decomposition_diagnostics.py

Diagnostics for the gene-set decomposition, replacing claims that had been asserted
in the text with tests. Covers BOTH arms.

FOUR QUESTIONS, each previously asserted without evidence:

  A  CANCELLATION. Is the net tAge difference a coherent shift, or a residual of
     opposing gene contributions? Reports the summed positive and negative per-gene
     contributions, their ratio to the net, the fraction of genes moving in the net
     direction, and how many genes account for half the net.

  B  DOES THE TRANSCRIPTOME ACTUALLY CHANGE? Reports the differentially expressed
     gene count against the tested background, which is the claim that IS supported
     and should not be confused with A.

  C  WHAT DRIVES THE COUNT OF SIGNIFICANT SETS? Correlation cannot answer this in
     the cross-sectional arm, because sample size and effect size are themselves
     correlated across conditions. Subsampling to a common n can, so every
     condition is subsampled to the smallest available size (11 v 11, where the
     smallest attainable Mann-Whitney p is 2.8e-6 and the test is therefore usable).

     The temporal arm CANNOT be tested this way. It is 6 v 6 throughout, and
     subsampling to 3 v 3 puts the smallest attainable two-sided p at
     2/C(6,3) = 0.100, so after correction nothing can reach significance whatever
     the data - an earlier version of this script did exactly that and returned
     nine zeros that were arithmetic rather than biology. Instead the temporal arm
     is handled by the fact that n is CONSTANT there: any variation in counts
     between groups cannot be a power effect, so the counts are compared directly
     against effect size (D).

  D  EFFECT SIZE, independent of n. Median |AUC - 0.5| across sets per group, so
     counts and effect sizes can be compared directly.

All computation is on the MORTALITY clock, matching the set-level sections, and on
the within-study contrasts for the meta-analysis.

ONE BH FAMILY, matching the rest of the arm. The convention in exploratory/05, 18 and
22 is to correct WITHIN AN ANALYSIS: 250 tests for the five arrest conditions (50 sets
x 5) and 450 for the nine temporal groups. This script uses the same family, including
inside the subsampling loop - for each draw, the resampled condition's 50 p-values are
combined with the other four conditions' unchanged p-values to form the full 250, BH is
applied to that, and the count is taken among the resampled condition's sets. The other
conditions supply the rest of the family exactly as they do in the real analysis, so the
subsampled and observed counts are on the same footing and comparable with the figures
quoted in the Results.

Output: rerun_outputs/decomposition_diagnostics.csv
        rerun_outputs/decomposition_subsampling.csv

NUMBER OF DRAWS. B = 10,000, chosen so the Monte Carlo error is negligible against
the quantity being estimated rather than by habit. The estimand is a median count of
significant sets, whose Monte Carlo standard error is about 1.253 x SD / sqrt(B);
with a between-draw SD of roughly 5 sets that is +-0.06 sets at B = 10,000 against
+-0.4 at B = 200. The achieved MCSE is reported per group in the output so the
choice can be checked. An earlier version used B = 200 for no reason other than that
the Mann-Whitney call had not been vectorised; scipy computes all 50 sets in one
call across an axis, which makes B = 10,000 cost about three minutes.

Usage: 23_decomposition_diagnostics.py <rerun_dir> <model_dir> [n_draws]
"""
import sys, warnings
import joblib, numpy as np, pandas as pd
from scipy.stats import mannwhitneyu
warnings.filterwarnings("ignore")

MORT = "EN_Mortality_Multispecies_Multitissue_scaleddiff.pkl"
SEED = 1
CONDS = {"Contact_inhibited CQ": "CICQ", "Serum_starved CQ": "SSCQ",
         "Replicative CS": "RS", "Stress-induced CS": "SIPS",
         "Oncogene-induced CS": "OIS"}
CELLS = ["Fibroblast", "Keratinocyte", "Melanocyte"]
TPS = ["4_days", "10_days", "20_days"]


def bh(p):
    p = np.asarray(p, float); n = len(p); o = np.argsort(p); a = np.empty(n)
    a[o] = np.minimum.accumulate((p[o] * n / (np.arange(n) + 1))[::-1])[::-1]
    return np.clip(a, 0, 1)


def load(model_dir):
    m = joblib.load(f"{model_dir}/{MORT}")
    imp = m.named_steps["imputation"]
    if not hasattr(imp, "_fill_dtype"):
        imp._fill_dtype = imp.statistics_.dtype
    return m


def contribs(m, path):
    feats = list(map(str, m.feature_names_in_))
    e = pd.read_csv(path); sid = e["sample_id"].values
    e = e.drop(columns=["sample_id"]); e.columns = e.columns.map(str)
    present = {c for c in e.columns if e[c].notna().any()}
    for g in [g for g in feats if g not in e.columns]:
        e[g] = np.nan
    Z = m.named_steps["scaler"].transform(m.named_steps["imputation"].transform(e.loc[:, feats]))
    obs = np.array([f in present for f in feats])
    return sid, Z * m.named_steps["estimator"].coef_[None, :], obs, feats


def cancellation(d):
    """d = per-gene contribution difference, measured genes only."""
    pos, neg, net = d[d > 0].sum(), d[d < 0].sum(), d.sum()
    sgn = np.sign(net)
    srt = np.sort(d * sgn)[::-1]
    n50 = int(np.searchsorted(np.cumsum(srt) / abs(net), 1.0)) + 1
    return dict(sum_positive=pos, sum_negative=neg, net=net,
                cancellation_ratio=(pos - neg) / abs(net),
                frac_genes_in_net_direction=float((np.sign(d) == sgn).mean()),
                n_genes_for_half_the_net=n50, n_genes=len(d))


def main(rerun_dir, model_dir, n_draws=10000):
    PT = f"{rerun_dir}/partial_tage"
    rng = np.random.default_rng(SEED)
    m = load(model_dir)
    pw = pd.read_csv(f"{PT}/hallmark_pathway_mouse_ids.csv")
    md = pd.read_csv(f"{rerun_dir}/sample_metadata_RERUN.csv").set_index("external_id")
    rows, sub = [], []

    # ---------------- cross-sectional ----------------
    sid, C, obs, feats = contribs(m, f"{PT}/meta_scaled_diff.csv")
    idx = {g: i for i, g in enumerate(feats)}
    sets = {n: np.array([idx[x] for x in s.mouse_gene_id.astype(str) if x in idx])
            for n, s in pw.groupby("pathway")}
    sets = {k: v for k, v in sets.items() if len(v)}
    S = np.column_stack([C[:, ix].sum(1) for ix in sets.values()])
    grp = pd.read_csv(f"{PT}/meta_groups.csv").set_index("sample_id")["group"]
    g = np.array([grp.get(s) for s in sid])
    study = np.array([md.loc[s, "study"] for s in sid])
    ctrl = np.where(g == "Proliferating")[0]

    def pvals(ia, ib, M=None):
        M = S if M is None else M
        return mannwhitneyu(M[ia], M[ib], axis=0).pvalue

    def count_sig_family(target_p, other_p):
        """BH over the whole analysis family, count significant among the target."""
        allp = np.concatenate([target_p] + other_p)
        adj = bh(allp)
        return int((adj[:len(target_p)] < 0.05).sum())

    def auc_dev(ia, ib, M=None):
        M = S if M is None else M
        u = mannwhitneyu(M[ia], M[ib], axis=0).statistic / (len(ia) * len(ib))
        return float(np.median(np.abs(u - 0.5)))

    nmin = min(int((g == cv).sum()) for cv in CONDS)
    print(f"cross-sectional: subsampling every condition to {nmin} vs {nmin}")
    for cv, lab in CONDS.items():
        # A: cancellation, on the within-study weighted per-gene difference
        acc = np.zeros(C.shape[1]); w = 0.0
        for s in np.unique(study):
            ti = np.where((study == s) & (g == cv))[0]
            ci = np.where((study == s) & (g == "Proliferating"))[0]
            if len(ti) and len(ci):
                ww = len(ti) * len(ci) / (len(ti) + len(ci))
                acc += ww * (C[ti].mean(0) - C[ci].mean(0)); w += ww
        r = cancellation((acc / w)[obs])
        ia = np.where(g == cv)[0]
        # observed p-values for every condition, so the analysis-wide family can be built
        obs_p = {l2: pvals(np.where(g == c2)[0], ctrl) for c2, l2 in CONDS.items()}
        others = [v for k2, v in obs_p.items() if k2 != lab]
        r.update(arm="cross_sectional", group=lab, n_test=len(ia), n_control=len(ctrl),
                 n_sig_sets=count_sig_family(obs_p[lab], others),
                 median_abs_auc_dev=auc_dev(ia, ctrl))
        rows.append(r)
        # C: subsample to the smallest condition
        cs = np.array([count_sig_family(
                           pvals(rng.choice(ia, nmin, replace=False),
                                 rng.choice(ctrl, nmin, replace=False)), others)
                       for _ in range(n_draws)])
        sub.append(dict(arm="cross_sectional", group=lab, n_matched=nmin, n_draws=n_draws,
                        n_sig_observed=r["n_sig_sets"], n_sig_matched_median=float(np.median(cs)),
                        p10=float(np.percentile(cs, 10)), p90=float(np.percentile(cs, 90)),
                        sd=float(cs.std()), mcse_median=float(1.253 * cs.std() / np.sqrt(n_draws))))
        print(f"  {lab:<5} observed {r['n_sig_sets']:>2} -> matched {np.median(cs):>4.0f}"
              f"   (MCSE {1.253*cs.std()/np.sqrt(n_draws):.2f})"
              f"   cancellation {r['cancellation_ratio']:>5.0f}x")

    # ---------------- temporal ----------------
    # p-value floor guard: refuse to subsample where the floor exceeds 0.05
    from math import comb
    floor_3v3 = 2 / comb(6, 3)
    print(f"\ntemporal: n is 6v6 throughout. Subsampling to 3v3 is NOT attempted, since"
          f" its p-floor is {floor_3v3:.3f} and no result could reach significance."
          f" Counts are compared against effect size instead, n being constant.")
    for ct in CELLS:
        for tp in TPS:
            stem = f"{ct}_{tp}"
            sid2, C2, obs2, feats2 = contribs(m, f"{PT}/{stem}_scaled_diff.csv")
            idx2 = {gg: i for i, gg in enumerate(feats2)}
            sets2 = {n: np.array([idx2[x] for x in s.mouse_gene_id.astype(str) if x in idx2])
                     for n, s in pw.groupby("pathway")}
            sets2 = {k: v for k, v in sets2.items() if len(v)}
            S2 = np.column_stack([C2[:, ix].sum(1) for ix in sets2.values()])
            gg = pd.read_csv(f"{PT}/{stem}_groups.csv").set_index("sample_id")["group"]
            lab2 = np.array([gg.get(s) for s in sid2])
            a = np.where(lab2 == tp)[0]; b = np.where(lab2 == "none")[0]
            r = cancellation((C2[a].mean(0) - C2[b].mean(0))[obs2])
            def cs_t(ia, ib):
                return int((bh(mannwhitneyu(S2[ia], S2[ib], axis=0).pvalue) < 0.05).sum())
            def auc_t(ia, ib):
                u = mannwhitneyu(S2[ia], S2[ib], axis=0).statistic / (len(ia) * len(ib))
                return float(np.median(np.abs(u - 0.5)))
            r.update(arm="temporal", group=stem, n_test=len(a), n_control=len(b),
                     median_abs_auc_dev=auc_t(a, b), _p=mannwhitneyu(S2[a], S2[b], axis=0).pvalue)
            rows.append(r)
            print(f"  {stem:<22} effect-size row done;"
                  f"   effect |AUC-0.5| {r['median_abs_auc_dev']:.3f}"
                  f"   cancellation {r['cancellation_ratio']:>5.0f}x"
                  f"   net {r['net']:+.3f}")

    out = pd.DataFrame(rows)
    # temporal counts on the 450-test family: all nine groups corrected together
    tmask = out.arm == "temporal"
    if tmask.any():
        tp = np.concatenate(list(out.loc[tmask, "_p"]))
        tadj = bh(tp); k = len(out.loc[tmask, "_p"].iloc[0])
        out.loc[tmask, "n_sig_sets"] = [
            int((tadj[i * k:(i + 1) * k] < 0.05).sum()) for i in range(int(tmask.sum()))]
    out = out.drop(columns=["_p"], errors="ignore")
    out["n_sig_sets"] = out.n_sig_sets.astype(int)
    # B: differential expression, the claim that IS supported
    try:
        deg = pd.read_csv(f"{rerun_dir}/deg_count_RERUN.csv").groupby("group_1").n.sum()
        bg = len(pd.read_csv(f"{rerun_dir}/../../Final/SI_tables/enrichment_background.csv"))
        print(f"\ndifferentially expressed genes, of {bg} tested:")
        for k, v in deg.items():
            print(f"  {k:<24} {v:>5} ({100*v/bg:.1f}%)")
        out.attrs["deg_background"] = bg
    except Exception as e:
        print(f"\nDEG counts not read: {e}")
    # D as a test, not a description: with n constant in the temporal arm, does the
    # count track effect size? Reported here so the claim is not made in prose alone.
    from scipy.stats import spearmanr
    t = out[out.arm == "temporal"]
    c = out[out.arm == "cross_sectional"]
    print("\n== what the counts track ==")
    print(f"  temporal (n constant at 6v6): count vs effect size rho = "
          f"{spearmanr(t.n_sig_sets, t.median_abs_auc_dev).statistic:+.2f}, "
          f"count vs net shift rho = {spearmanr(t.n_sig_sets, t.net).statistic:+.2f}")
    print(f"  cross-sectional (n varies)  : count vs n rho = "
          f"{spearmanr(c.n_sig_sets, c.n_test).statistic:+.2f}, "
          f"count vs effect size rho = "
          f"{spearmanr(c.n_sig_sets, c.median_abs_auc_dev).statistic:+.2f}"
          f"  (n and effect size themselves correlate rho = "
          f"{spearmanr(c.n_test, c.median_abs_auc_dev).statistic:+.2f}, which is why"
          f" the subsampling above is needed)")
    out.to_csv(f"{rerun_dir}/decomposition_diagnostics.csv", index=False)
    pd.DataFrame(sub).to_csv(f"{rerun_dir}/decomposition_subsampling.csv", index=False)
    print(f"\nSaved -> decomposition_diagnostics.csv, decomposition_subsampling.csv")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], int(sys.argv[3]) if len(sys.argv) > 3 else 10000)
