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

TWO BH FAMILIES, deliberately, and they must not be confused. The project convention
elsewhere (exploratory/05, 18, 22) corrects WITHIN AN ANALYSIS - 250 tests for the
five arrest conditions, 450 for the nine temporal groups - and that is the family the
Results sections quote. The subsampling test here cannot use it, because under
subsampling the other groups' p-values are not recomputed, so there is no
analysis-wide family to correct within; it therefore corrects within each group over
50 tests, applied identically to the full and subsampled data so the comparison is
internally consistent. Both counts are written out: n_sig_sets_analysiswide is the
figure to quote, n_sig_sets_pergroup the one the subsampling uses. They differ by at
most 5 sets and give the same correlations (rho = 0.94 versus 0.92 against effect
size in the temporal arm), so nothing rests on the choice, but the Results text must
use one consistently.

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

    def count_sig(ia, ib, M=None):
        M = S if M is None else M
        return int((bh(mannwhitneyu(M[ia], M[ib], axis=0).pvalue) < 0.05).sum())

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
        r.update(arm="cross_sectional", group=lab, n_test=len(ia), n_control=len(ctrl),
                 n_sig_sets_full=count_sig(ia, ctrl), median_abs_auc_dev=auc_dev(ia, ctrl))
        rows.append(r)
        # C: subsample to the smallest condition
        cs = np.array([count_sig(rng.choice(ia, nmin, replace=False),
                                 rng.choice(ctrl, nmin, replace=False))
                       for _ in range(n_draws)])
        sub.append(dict(arm="cross_sectional", group=lab, n_matched=nmin, n_draws=n_draws,
                        n_sig_full=r["n_sig_sets_full"], n_sig_matched_median=float(np.median(cs)),
                        p10=float(np.percentile(cs, 10)), p90=float(np.percentile(cs, 90)),
                        sd=float(cs.std()), mcse_median=float(1.253 * cs.std() / np.sqrt(n_draws))))
        print(f"  {lab:<5} full {r['n_sig_sets_full']:>2} -> matched {np.median(cs):>4.0f}"
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
                     n_sig_sets_full=cs_t(a, b), median_abs_auc_dev=auc_t(a, b))
            rows.append(r)
            print(f"  {stem:<22} sets {r['n_sig_sets_full']:>2}"
                  f"   effect |AUC-0.5| {r['median_abs_auc_dev']:.3f}"
                  f"   cancellation {r['cancellation_ratio']:>5.0f}x"
                  f"   net {r['net']:+.3f}")

    out = pd.DataFrame(rows).rename(columns={"n_sig_sets_full": "n_sig_sets_pergroup"})
    # add the analysis-wide count, which is the project convention and the figure to quote
    try:
        allsets = pd.read_csv(f"{rerun_dir}/mortality_partial_tage_ALL.csv")
        conv = {}
        for a, gg in allsets.groupby("analysis"):
            for lb, g2 in gg.groupby("label"):
                conv[lb] = int((g2.p_adj < 0.05).sum())
        out["n_sig_sets_analysiswide"] = out.group.map(conv)
    except Exception as e:
        print(f"analysis-wide counts not added: {e}")
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
          f"{spearmanr(t.n_sig_sets_analysiswide, t.median_abs_auc_dev).statistic:+.2f}, "
          f"count vs net shift rho = {spearmanr(t.n_sig_sets_analysiswide, t.net).statistic:+.2f}")
    print(f"  cross-sectional (n varies)  : count vs n rho = "
          f"{spearmanr(c.n_sig_sets_analysiswide, c.n_test).statistic:+.2f}, "
          f"count vs effect size rho = "
          f"{spearmanr(c.n_sig_sets_analysiswide, c.median_abs_auc_dev).statistic:+.2f}"
          f"  (n and effect size themselves correlate rho = "
          f"{spearmanr(c.n_test, c.median_abs_auc_dev).statistic:+.2f}, which is why"
          f" the subsampling above is needed)")
    out.to_csv(f"{rerun_dir}/decomposition_diagnostics.csv", index=False)
    pd.DataFrame(sub).to_csv(f"{rerun_dir}/decomposition_subsampling.csv", index=False)
    print(f"\nSaved -> decomposition_diagnostics.csv, decomposition_subsampling.csv")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], int(sys.argv[3]) if len(sys.argv) > 3 else 10000)
