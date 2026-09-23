#!/usr/bin/env python3
"""42_coverage_threshold_free.py

Replaces the rank-cut coverage tests with two quantities that need no cutoff.

WHY THE EARLIER VERSIONS WERE BIASED
  half-net lists   "the genes needed to reach half a condition's net" has a length set by
                   how small the residue is, not by how concentrated the movement is. SSCQ
                   cancels 71-fold and gets 3 genes; RS cancels 12-fold and gets 20. The
                   apparent quiescence/senescence coverage difference (Fisher p = 0.020)
                   came from that and vanishes on a fixed-size list (p = 0.534).
  top-N lists      N is arbitrary, the union's size depends on how many groups a dataset
                   has, and ranking by |contribution| selects on |coefficient|. Annotated
                   clock genes carry 1.05x the median |coefficient| of unannotated ones
                   (Mann-Whitney p = 0.0074) and are 35.2% of the top |coefficient| decile
                   against 29.3% overall, so a top-N list reads as better annotated before
                   any biology enters.

WHAT IS COMPUTED INSTEAD, per condition or timepoint group
  n_half      genes needed to account for half the TOTAL MOVEMENT, sum|c|, taken largest
              first. The legible concentration statistic, and unlike the half-net count of
              earlier drafts its denominator does not shrink with cancellation: 767 to 971
              genes for half the movement, against 3 to 20 for half the net.
  n_eff       (sum|c|)^2 / sum(c^2) over all measured clock genes: the number of
              equal-sized contributions that would give the same spread. No cutoff, and
              comparable across groups. Reported alongside n_half for the figure.
  ann_pct     share of the total |contribution| carried by genes in >= 1 Hallmark set, i.e.
              the part a set-level decomposition can see. Every gene enters, weighted by how
              much it moves. unann_pct is its complement and is kept for reference.
  null        unann_pct with the annotation labels permuted WITHIN |coefficient| decile,
              20,000 draws, which holds the coefficient confound fixed. Raw empirical p
              with its floor, two-sided, no BH: one pre-specified test per group.

Contributions are the same per-gene within-study differences as 2.1.5.1 (arrest) and the
timepoint-against-baseline differences of 2.2.3 (time course).

Also writes a per-gene table so the same question can be asked of annotation sources other
than Hallmark: 43_gene_annotation_coverage.R joins GO biological-process terms onto it, since
absence from the 50 Hallmark sets is not absence of known function and the two should not be
conflated.

Output: rerun_outputs/coverage_threshold_free.csv
        rerun_outputs/gene_contribution_table.csv
        rerun_outputs/clock_unaccounted_shares.csv   (panel (a) of the figure)
"""
import sys, warnings
import joblib, numpy as np, pandas as pd
from sklearn.exceptions import InconsistentVersionWarning
warnings.filterwarnings("ignore", category=InconsistentVersionWarning)

NDRAW = 20000
CHUNK = 2500
SEED = 1
LAB = {"Contact_inhibited CQ": "CICQ", "Serum_starved CQ": "SSCQ", "Stress-induced CS": "SIPS",
       "Oncogene-induced CS": "OIS", "Replicative CS": "RS"}
CTS = ("Fibroblast", "Keratinocyte", "Melanocyte")
TPS = ("4_days", "10_days", "20_days")


def stratified_null(w, is_ann, strata, rng, ndraw=NDRAW, chunk=CHUNK):
    """unannotated share of sum|w| with labels permuted inside each stratum"""
    tot = w.sum()
    ann_sum = np.zeros(ndraw)
    for st in np.unique(strata):
        ix = np.where(strata == st)[0]
        k = int(is_ann[ix].sum())
        if k == 0:
            continue
        ws = w[ix]
        for a in range(0, ndraw, chunk):
            b = min(a + chunk, ndraw)
            pick = np.argsort(rng.random((b - a, len(ix))), axis=1)[:, :k]
            ann_sum[a:b] += ws[pick].sum(axis=1)
    return 100 * (tot - ann_sum) / tot


def main(rerun_dir, model_dir):
    PT = f"{rerun_dir}/partial_tage"
    m = joblib.load(f"{model_dir}/EN_Mortality_Multispecies_Multitissue_scaleddiff.pkl")
    imp = m.named_steps["imputation"]
    if not hasattr(imp, "_fill_dtype"):
        imp._fill_dtype = imp.statistics_.dtype
    feats = np.array(list(map(str, m.feature_names_in_)))
    coef = m.named_steps["estimator"].coef_
    hall = set(pd.read_csv(f"{PT}/hallmark_pathway_mouse_ids.csv")["mouse_gene_id"].astype(str))
    is_ann_all = np.array([f in hall for f in feats])

    def load(stem):
        e = pd.read_csv(f"{PT}/{stem}_scaled_diff.csv")
        s = e["sample_id"].values
        e = e.drop(columns=["sample_id"]); e.columns = e.columns.map(str)
        present = {c for c in e.columns if e[c].notna().any()}
        for g in [g for g in feats if g not in e.columns]:
            e[g] = np.nan
        Z = m.named_steps["scaler"].transform(imp.transform(e.loc[:, feats]))
        ob = np.array([f in present for f in feats])
        g_ = pd.read_csv(f"{PT}/{stem}_groups.csv").set_index("sample_id").loc[s, "group"].values
        return s, Z * coef[np.newaxis, :], g_, ob

    rng = np.random.default_rng(SEED)
    rows = []
    per_gene, MASK = {}, {}

    def run(dset, group, v, ob):
        per_gene[(dset, group)] = v
        MASK[(dset, group)] = ob
        w = np.abs(v)
        tot = w.sum()
        A = is_ann_all[ob]
        ac = np.abs(coef)[ob]
        strata = np.digitize(ac, np.quantile(ac, np.linspace(0, 1, 11)[1:-1]))
        neff = tot ** 2 / (v ** 2).sum()
        cum = np.cumsum(np.sort(w)[::-1]) / tot
        nhalf = int(np.searchsorted(cum, 0.5)) + 1
        top10 = 100 * cum[9]
        o = 100 * w[~A].sum() / tot
        nul = stratified_null(w, A, strata, rng)
        z = (o - nul.mean()) / nul.std()
        p = (1 + 2 * min((nul <= o).sum(), (nul >= o).sum())) / (NDRAW + 1)
        # the annotated share is what the section reports; same test, sign flipped
        a_o, a_nul = 100 - o, 100 - nul
        a_z = (a_o - a_nul.mean()) / a_nul.std()
        print(f"  {group:<21}{len(v):>7}{nhalf:>8}{neff:>8.0f}{a_o:>9.1f}{a_nul.mean():>9.1f}"
              f"{a_nul.std():>8.2f}{a_z:>+7.1f}{p:>10.5f}")
        rows.append(dict(dataset=dset, group=group, n_measured=len(v),
                         n_half_movement=nhalf, top10_pct=top10, n_eff=neff,
                         ann_pct=a_o, ann_null_mean=a_nul.mean(), z_ann=a_z,
                         net=v.sum(), total_abs_contribution=tot, unann_pct=o,
                         null_mean=nul.mean(), null_sd=nul.std(), z=z, p_emp=p,
                         p_floor=1 / (NDRAW + 1),
                         pct_genes_unann=100 * (~A).mean(),
                         pct_weight_unann=100 * ac[~A].sum() / ac.sum()))

    hdr = (f"\n  {'group':<21}{'genes':>7}{'n_half':>8}{'n_eff':>8}{'Hall %':>9}{'null %':>9}"
           f"{'null sd':>8}{'z':>7}{'p':>10}")

    # ---- arrest conditions ------------------------------------------------
    sid, contrib, _, ob = load("meta")
    ann = pd.read_csv(f"{rerun_dir}/immortalisation_annotation_corrected.csv").set_index("external_id")
    grp = ann.reindex(sid)["cell_substate"].values
    std = ann.reindex(sid)["study"].values
    print("=== arrest conditions: within-study difference from own-study controls ===" + hdr)
    for full, lab in LAB.items():
        ds, ws = [], []
        for s in pd.unique(std[grp == full]):
            t = contrib[(std == s) & (grp == full)]
            c = contrib[(std == s) & (grp == "Proliferating")]
            if len(t) < 1 or len(c) < 1:
                continue
            wt = len(t) * len(c) / (len(t) + len(c))
            ds.append(wt * (t.mean(axis=0) - c.mean(axis=0))); ws.append(wt)
        run("arrest", lab, (np.sum(ds, axis=0) / sum(ws))[ob], ob)

    # ---- irradiation time course ------------------------------------------
    print("\n=== irradiation time course: each timepoint against its own baseline ===" + hdr)
    for ct in CTS:
        for tp in TPS:
            _, c2, g2, ob2 = load(f"{ct}_{tp}")
            run("time_course", f"{ct}_{tp}",
                (c2[g2 == tp].mean(axis=0) - c2[g2 == "none"].mean(axis=0))[ob2], ob2)

    D = pd.DataFrame(rows)
    D.to_csv(f"{rerun_dir}/coverage_threshold_free.csv", index=False)

    # per-gene table: one row per clock feature, |contribution| per group
    G = pd.DataFrame({"gene": feats, "coef": coef, "in_hallmark": is_ann_all})
    for (dset, group), v in per_gene.items():
        col = np.full(len(feats), np.nan)
        col[MASK[(dset, group)]] = np.abs(v)
        G[f"abs_contrib_{group}"] = col
    G["measured_in_arrest"] = MASK[("arrest", "CICQ")]
    G.to_csv(f"{rerun_dir}/gene_contribution_table.csv", index=False)

    # panel (a): what the 50 sets account for, measured three ways on the same model.
    # Gene count and weight are properties of the trained clock; movement is the mean
    # over the 14 groups. Written here rather than in the figure script so the figure
    # plots numbers it does not itself compute.
    ac = np.abs(coef)
    pd.DataFrame([
        dict(measure="by gene count", pct_in_hallmark=100 * is_ann_all.mean(),
             pct_not_in_hallmark=100 * (~is_ann_all).mean(),
             n_in_hallmark=int(is_ann_all.sum()), n_not_in_hallmark=int((~is_ann_all).sum())),
        dict(measure="by clock weight", pct_in_hallmark=100 * ac[is_ann_all].sum() / ac.sum(),
             pct_not_in_hallmark=100 * ac[~is_ann_all].sum() / ac.sum(),
             n_in_hallmark=np.nan, n_not_in_hallmark=np.nan),
        dict(measure="by tAge movement", pct_in_hallmark=D.ann_pct.mean(),
             pct_not_in_hallmark=100 - D.ann_pct.mean(),
             n_in_hallmark=np.nan, n_not_in_hallmark=np.nan),
    ]).to_csv(f"{rerun_dir}/clock_unaccounted_shares.csv", index=False)
    print("Saved -> clock_unaccounted_shares.csv")
    print(f"Saved -> gene_contribution_table.csv  ({len(G)} clock features)")
    print(f"\n  reference: unannotated share of measured genes "
          f"{D.pct_genes_unann.min():.1f}-{D.pct_genes_unann.max():.1f}%, "
          f"of the clock's |weight| {D.pct_weight_unann.min():.1f}-{D.pct_weight_unann.max():.1f}%")
    print(f"\nSaved -> coverage_threshold_free.csv")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
