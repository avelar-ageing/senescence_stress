"""47_set_driver_genes.py

WHICH GENES MAKE A GENE SET STAND OUT, for every set that beats its matched null.

WHY. exploratory/20 (arrest conditions, --mortality) and exploratory/22 (time course)
test each gene set against 20,000 random sets of clock genes with the same number of
genes and the same spread of clock weights. 40 of 250 and 74 of 450 comparisons beat
that null at p < 0.05. Within those sets, about half the genes push the same way as
the set, as in the sets that do not beat it (exploratory/28,
withinset_cancellation_by_null.csv). So the excess comes from how far some genes move.
This script finds those genes.

WHAT IT COMPUTES, per set x group that beats its null (published p_emp < 0.05):

  excess of a gene = its contribution  -  the average contribution of all clock genes in
                     the same clock-weight tenth (the stratum the null draws from)
  The set's excess over its matched null is the sum of its genes' excesses, so each
  gene's SHARE OF THE SET'S EXCESS is exact and needs no threshold. Shares sum to 1;
  a negative share means the gene works against the set's excess.

  same_weight_rank_pct: the percentage of clock genes of the same weight that move less
                     far in the direction of the set's excess. 99.9 = this gene moves
                     further than 999 of 1,000 genes of equal clock weight.

  p_without_gene:    the set re-tested against its matched null with this one gene
                     removed. The null is the same 20,000 draws with one draw from the
                     gene's weight stratum dropped, so it is the matched null of the
                     reduced set. A set whose p rises to >= 0.05 depends on that gene.

  n_removed_to_lose: remove the set's genes one at a time, largest excess first, until
                     the set no longer beats its null (p >= 0.05). 1 = one gene carries
                     the set; a large number = the excess is spread.

  n_genes_half_excess: fewest genes, largest first, whose shares reach half the excess.

NULL. The draws are replayed with the same seed, order and calls as exploratory/20
(--mortality) and exploratory/22, and the script stops unless every replayed p equals
the published p_emp exactly, so these results refer to the same null as the text.
Empirical p = (1 + k)/(B + 1), B = 20,000, floor 5e-5, two-sided as in 20/22. No BH:
house convention for simulation nulls. Clock: mortality only, as in the set-level text.

A gene belongs to up to several Hallmark sets. Its excess in a group is the same in
every set that contains it (it depends on the gene and its weight stratum only), so a
gene that carries one set carries equally every other set it belongs to; the figure
(exploratory/48) marks such genes.

USAGE: 47_set_driver_genes.py <rerun_dir> <model_dir> [n_null]
OUT:   <rerun_dir>/set_driver_genes.csv     one row per gene x set x group (sets that
                                             beat their null only)
       <rerun_dir>/set_driver_summary.csv   one row per set x group
"""
import os
import subprocess
import sys
import warnings

import joblib
import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")

MORT = "EN_Mortality_Multispecies_Multitissue_scaleddiff.pkl"
SEED = 1
CONDS = {"Contact_inhibited CQ": "CICQ", "Serum_starved CQ": "SSCQ",
         "Replicative CS": "RS", "Stress-induced CS": "SIPS",
         "Oncogene-induced CS": "OIS"}
CELLS = ["Fibroblast", "Keratinocyte", "Melanocyte"]
TPS = ["4_days", "10_days", "20_days"]
ALPHA = 0.05


def _patch(imp):
    if not hasattr(imp, "_fill_dtype"):
        imp._fill_dtype = imp.statistics_.dtype if hasattr(imp, "statistics_") else np.float64


def tage_gene_table():
    p = os.environ.get("TAGE_GENE_TABLE")
    if not p:
        p = subprocess.run(["Rscript", "-e", "cat(system.file('extdata/metadata/Gene_table_mouse.csv', package='tAge'))"],
                           capture_output=True, text=True, check=True).stdout.strip()
    if not p or not os.path.exists(p):
        sys.exit("Gene_table_mouse.csv not found: install the tAge R package or set TAGE_GENE_TABLE")
    return p


def p_two_sided(null, obs):
    """exploratory/20 and 22: (1 + #|null - mean| >= |obs - mean|) / (B + 1)."""
    c = null.mean()
    return (1 + np.sum(np.abs(null - c) >= abs(obs - c))) / (len(null) + 1)


def analyse(arm, group, pathway, cols, d, strat, by_strat, mats, p_pub, sym, feats):
    """mats: {stratum: (B x k_s) draw matrix} exactly as drawn by 20/22."""
    draws = np.zeros(next(iter(mats.values())).shape[0])
    for M in mats.values():
        draws += M.sum(axis=1)
    obs = d[cols].sum()
    p_full = p_two_sided(draws, obs)
    # compare the counts k in p = (1 + k)/(B + 1): exact, unlike the CSV's printed float
    B1 = len(draws) + 1
    assert round(p_full * B1) == round(p_pub * B1), \
        f"{arm} {group} {pathway}: replayed p {p_full} != published {p_pub}"

    mu = {s: d[by_strat[s]].mean() for s in mats}
    s_of = strat[cols]
    exc = d[cols] - np.array([mu[s] for s in s_of])
    total_exc = exc.sum()
    sgn = np.sign(total_exc) if total_exc != 0 else 1.0
    share = exc / total_exc

    # rank among same-weight clock genes, in the direction of the set's excess
    rank = np.array([100 * np.mean(sgn * d[by_strat[s]] < sgn * v) for s, v in zip(s_of, d[cols])])

    # leave one gene out: drop one draw from the gene's stratum
    p_wo = np.empty(len(cols))
    for s in np.unique(s_of):
        null_minus = draws - mats[s][:, 0]
        for i in np.where(s_of == s)[0]:
            p_wo[i] = p_two_sided(null_minus, obs - d[cols[i]])

    # remove genes one at a time, largest excess (in the set's direction) first
    order = np.argsort(-sgn * exc)
    removed = {s: 0 for s in mats}
    null_k, obs_k, n_lose, p_after = draws.copy(), obs, np.nan, np.nan
    for step, i in enumerate(order, start=1):
        s = s_of[i]
        null_k = null_k - mats[s][:, removed[s]]
        removed[s] += 1
        obs_k = obs_k - d[cols[i]]
        pk = p_two_sided(null_k, obs_k)
        if pk >= ALPHA:
            n_lose, p_after = step, pk
            break

    cum = np.cumsum(share[order])
    n_half = int(np.searchsorted(cum, 0.5) + 1) if cum[-1] >= 0.5 else np.nan
    pos = np.clip(share, 0, None)
    n_eff = pos.sum() ** 2 / (pos ** 2).sum() if (pos ** 2).sum() > 0 else np.nan

    genes = pd.DataFrame(dict(
        arm=arm, group=group, pathway=pathway, gene_id=[feats[c] for c in cols],
        gene=[sym.get(feats[c], feats[c]) for c in cols],
        contribution=d[cols], expected_for_weight=[mu[s] for s in s_of], excess=exc,
        share_of_set_excess=share, same_weight_rank_pct=rank,
        p_without_gene=p_wo, set_depends_on_gene=p_wo >= ALPHA))
    genes["rank_in_set"] = (-sgn * genes.excess).rank(method="first").astype(int)
    top = genes.sort_values("rank_in_set").iloc[0]
    summary = dict(
        arm=arm, group=group, pathway=pathway, n_genes=len(cols), observed=obs,
        null_mean=draws.mean(), excess=total_exc, p_emp=p_full,
        direction="up" if sgn > 0 else "down",
        top_gene=top.gene, top_gene_share=top.share_of_set_excess,
        top_gene_same_weight_rank_pct=top.same_weight_rank_pct,
        p_without_top_gene=top.p_without_gene,
        n_genes_set_depends_on=int(genes.set_depends_on_gene.sum()),
        n_removed_to_lose=n_lose, p_after_removal=p_after,
        n_genes_half_excess=n_half, n_eff_positive=n_eff)
    return genes, summary


def main(rerun_dir, model_dir, n_null=20000):
    PT = f"{rerun_dir}/partial_tage"
    m = joblib.load(f"{model_dir}/{MORT}"); _patch(m.named_steps["imputation"])
    feats = list(map(str, m.feature_names_in_))
    idx = {g: i for i, g in enumerate(feats)}
    coef = m.named_steps["estimator"].coef_; aco = np.abs(coef)
    gt = pd.read_csv(tage_gene_table())
    sym = dict(zip(gt.Entrez.astype(str), gt["Gene.Symbol"]))
    pw = pd.read_csv(f"{PT}/hallmark_pathway_mouse_ids.csv")
    genes_out, summ_out = [], []

    # ---------- arrest conditions: replay exploratory/20 --mortality -------------------
    ys = pd.read_csv(f"{rerun_dir}/pathway_specificity_yardstick_mortality.csv")
    pub = {(r.label, r.pathway): r.p_emp for r in ys.itertuples()}
    rng = np.random.default_rng(SEED)
    e = pd.read_csv(f"{PT}/meta_scaled_diff.csv"); sid = e["sample_id"].values
    e = e.drop(columns=["sample_id"]); e.columns = e.columns.map(str)
    for g in [g for g in feats if g not in e.columns]:
        e[g] = np.nan
    Z = m.named_steps["scaler"].transform(m.named_steps["imputation"].transform(e.loc[:, feats]))
    C = Z * coef[np.newaxis, :] * 1.0
    groups = pd.read_csv(f"{PT}/meta_groups.csv")
    groups = groups.set_index(pd.Index(range(len(groups))))
    assert list(groups.sample_id) == list(sid)
    meta = pd.read_csv(f"{rerun_dir}/sample_metadata_RERUN.csv")
    groups["study"] = groups.sample_id.map(dict(zip(meta.external_id, meta.study)))
    nz = aco > 0
    strat = np.zeros(len(feats), dtype=int)
    q = np.quantile(aco[nz], np.linspace(0, 1, 11)[1:-1])
    strat[nz] = 1 + np.searchsorted(q, aco[nz])
    by_strat = {s: np.where(strat == s)[0] for s in np.unique(strat)}
    sets = {n: [idx[x] for x in s.mouse_gene_id.astype(str) if x in idx] for n, s in pw.groupby("pathway")}
    for cond_value, lab in CONDS.items():
        acc = np.zeros(C.shape[1]); wsum = 0.0
        for st, g in groups.groupby("study"):
            ti = g.index[g.group == cond_value].to_numpy(); ci = g.index[g.group == "Proliferating"].to_numpy()
            if not len(ti) or not len(ci):
                continue
            w = len(ti) * len(ci) / (len(ti) + len(ci))
            acc += w * (C[ti].mean(axis=0) - C[ci].mean(axis=0)); wsum += w
        d = acc / wsum
        for name, cols in sets.items():
            if not cols:
                continue
            cols = np.asarray(cols)
            keep = pub[(lab, name)] < ALPHA
            mats = {}
            for st_id in np.unique(strat[cols]):
                k = int((strat[cols] == st_id).sum())
                M = rng.choice(d[by_strat[st_id]], size=(n_null, k), replace=True)
                if keep:
                    mats[st_id] = M
            if keep:
                gdf, s = analyse("arrest", lab, name, cols, d, strat, by_strat, mats, pub[(lab, name)], sym, feats)
                genes_out.append(gdf); summ_out.append(s)
        print(f"  arrest {lab} done")

    # ---------- time course: replay exploratory/22 ------------------------------------
    yt = pd.read_csv(f"{rerun_dir}/pathway_specificity_yardstick_mortality_temporal.csv")
    pubt = {(r.label, r.pathway): r.p_emp for r in yt.itertuples()}
    rng = np.random.default_rng(SEED)
    sets = {n: np.array([idx[x] for x in s.mouse_gene_id.astype(str) if x in idx]) for n, s in pw.groupby("pathway")}
    sets = {k: v for k, v in sets.items() if len(v)}
    q = np.quantile(aco, np.linspace(0, 1, 11)[1:-1])
    strat = np.searchsorted(q, aco)
    by_strat = {s: np.where(strat == s)[0] for s in np.unique(strat)}
    for ct in CELLS:
        for tp in TPS:
            stem = f"{ct}_{tp}"
            e = pd.read_csv(f"{PT}/{stem}_scaled_diff.csv"); sid2 = e["sample_id"].values
            e = e.drop(columns=["sample_id"]); e.columns = e.columns.map(str)
            for g in [g for g in feats if g not in e.columns]:
                e[g] = np.nan
            C2 = m.named_steps["scaler"].transform(m.named_steps["imputation"].transform(e.loc[:, feats])) * coef[None, :]
            gg = pd.read_csv(f"{PT}/{stem}_groups.csv").set_index("sample_id")["group"]
            lab2 = np.array([gg.get(s) for s in sid2])
            d = C2[lab2 == tp].mean(axis=0) - C2[lab2 == "none"].mean(axis=0)
            for nm, cols in sets.items():
                keep = pubt[(stem, nm)] < ALPHA
                mats = {}
                for st_id in np.unique(strat[cols]):
                    kk = int((strat[cols] == st_id).sum())
                    M = rng.choice(d[by_strat[st_id]], size=(n_null, kk), replace=True)
                    if keep:
                        mats[st_id] = M
                if keep:
                    gdf, s = analyse("time_course", stem, nm, cols, d, strat, by_strat, mats, pubt[(stem, nm)], sym, feats)
                    genes_out.append(gdf); summ_out.append(s)
            print(f"  time course {stem} done")

    G = pd.concat(genes_out, ignore_index=True)
    S = pd.DataFrame(summ_out)
    G.to_csv(f"{rerun_dir}/set_driver_genes.csv", index=False)
    S.to_csv(f"{rerun_dir}/set_driver_summary.csv", index=False)

    print(f"\nsets that beat their null: {len(S)} ({(S.arm == 'arrest').sum()} arrest, "
          f"{(S.arm == 'time_course').sum()} time course); replayed p equal to published in all")
    for arm, s in S.groupby("arm"):
        print(f"\n== {arm} ==")
        print(f"  set no longer beats its null without its single largest gene: "
              f"{int((s.p_without_top_gene >= ALPHA).sum())} of {len(s)}")
        print("  genes removed (largest first) before the set stops beating its null:",
              s.n_removed_to_lose.describe()[["min", "25%", "50%", "75%", "max"]].round(1).to_dict())
        print("  share of the excess carried by the largest gene: median "
              f"{100 * s.top_gene_share.median():.0f}% (range {100 * s.top_gene_share.min():.0f}-{100 * s.top_gene_share.max():.0f}%)")
        print("  genes needed for half the excess: median", s.n_genes_half_excess.median())
    print(f"\nSaved -> set_driver_genes.csv ({len(G)} rows), set_driver_summary.csv ({len(S)} rows)")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], int(sys.argv[3]) if len(sys.argv) > 3 else 20000)
