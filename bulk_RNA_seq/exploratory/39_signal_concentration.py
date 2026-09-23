#!/usr/bin/env python3
"""39_signal_concentration.py

Each condition's mortality-clock shift is the small residue of much larger opposing
sums, and half of that residue is reached within 3 to 20 genes of about 9,000 (the
n_genes_for_the_whole_net column of decomposition_diagnostics.csv counts to the whole
net, not half, so its 4-54 range is a different statistic). 37 showed the 545 cell-cycle
features contribute about 2% and can be deleted without changing any condition's
result. So the signal is neither proliferation shutdown nor spread across the
clock: it concentrates into a few dozen genes.

THE QUESTION. Are they the SAME few genes in every condition? If one small group
carries the shift in all five arrest conditions, that is the arrest signature and
it is nameable. If the genes differ, the concentration is condition-specific and
only the set-level description generalises.

WHAT IS MEASURED
  overlap        Jaccard between conditions of the top-k contributing genes, at a
                 fixed k = 50 and at each condition's own half-net k.
  recurrence     how many genes appear in the top-k of all five conditions, and
                 of at least three.
  null           the same statistics for random gene groups of matched size drawn
                 from the clock's own features, 20,000 draws. Two conditions'
                 top-50 lists overlap somewhat by chance alone given they are
                 drawn from the same 9,000 features, so the null supplies that
                 baseline. Raw empirical p with its floor.

Genes are ranked by SIGNED contribution to the within-study difference, the same
estimator as 2.1.5.1. BOTH DIRECTIONS are analysed: the top k are the genes pushing
the score up hardest, the bottom k those pulling it down hardest. An earlier version
reported only the upward half, which made the convergence look like a property of
upward movement; it is not. The two lists cannot overlap, so the down side carries
its own recurrence counts and its own null.

Both datasets are run. Mouse Entrez ids are reported as ids; no symbol mapping is
available locally (org.Mm.eg.db is not installed), so naming is left to the
caller.

CLASS SIGNATURES. The same question inside an arrest class: genes in the top-k of every
senescence condition, of both quiescence conditions, and the overlap between those two
lists. Nulls are matched random groups drawn from the clock's own features, NDRAW
draws like every other null here. These were computed ad hoc for an earlier draft and never written
to a file; they are persisted here so 2.1.5.3 has a table behind every number in it.

Output: rerun_outputs/signal_concentration_directions.csv  (recurrence, up and down)
        rerun_outputs/signal_concentration.csv
        rerun_outputs/signal_concentration_top_genes.csv
        rerun_outputs/signal_concentration_classes.csv
"""
import sys, warnings, itertools
import joblib, numpy as np, pandas as pd
from sklearn.exceptions import InconsistentVersionWarning
warnings.filterwarnings("ignore", category=InconsistentVersionWarning)

NDRAW = 20000   # EVERY null in this script. Earlier versions hardcoded 2,000 for the
                # recurrence counts and 5,000 for the class signatures while this
                # constant said 20,000, so neighbouring p-values in 2.1.5.4 sat on three
                # different floors (5e-4, 2e-4, 5e-5) and read as different strengths of
                # evidence. None is BH-adjusted, so 20,000 here is for consistency with
                # 37 and 42, not for BH reachability (see 22's header for that).
SEED = 1
TOPK = 50
LAB = {"Contact_inhibited CQ": "CICQ", "Serum_starved CQ": "SSCQ", "Stress-induced CS": "SIPS",
       "Oncogene-induced CS": "OIS", "Replicative CS": "RS"}
CTS = ("Fibroblast", "Keratinocyte", "Melanocyte")
TPS = ("4_days", "10_days", "20_days")


def jac(a, b):
    a, b = set(a), set(b)
    return len(a & b) / len(a | b) if (a | b) else np.nan


def half_net_k(v):
    """genes needed, taken largest first, to reach half the net"""
    s = np.sort(v)[::-1]
    net = v.sum()
    if net <= 0:
        return len(v)
    c = np.cumsum(s)
    idx = np.searchsorted(c, net / 2.0)
    return int(min(idx + 1, len(v)))


def main(rerun_dir, model_dir):
    PT = f"{rerun_dir}/partial_tage"
    m = joblib.load(f"{model_dir}/EN_Mortality_Multispecies_Multitissue_scaleddiff.pkl")
    imp = m.named_steps["imputation"]
    if not hasattr(imp, "_fill_dtype"):
        imp._fill_dtype = imp.statistics_.dtype
    feats = np.array(list(map(str, m.feature_names_in_)))
    coef = m.named_steps["estimator"].coef_

    def load(stem):
        e = pd.read_csv(f"{PT}/{stem}_scaled_diff.csv")
        s = e["sample_id"].values
        e = e.drop(columns=["sample_id"]); e.columns = e.columns.map(str)
        for g in [g for g in feats if g not in e.columns]:
            e[g] = np.nan
        Z = m.named_steps["scaler"].transform(imp.transform(e.loc[:, feats]))
        g_ = pd.read_csv(f"{PT}/{stem}_groups.csv").set_index("sample_id").loc[s, "group"].values
        return s, Z * coef[np.newaxis, :], g_

    # ---- per-gene within-study contribution, both datasets ---------------
    sid, contrib, grp = load("meta")
    ann = pd.read_csv(f"{rerun_dir}/immortalisation_annotation_corrected.csv").set_index("external_id")
    std = ann.reindex(sid)["study"].values
    D = {}
    for full, lab in LAB.items():
        ds, ws = [], []
        for s in pd.unique(std[grp == full]):
            t = contrib[(std == s) & (grp == full)]
            c = contrib[(std == s) & (grp == "Proliferating")]
            if len(t) < 1 or len(c) < 1:
                continue
            w = len(t) * len(c) / (len(t) + len(c))
            ds.append(w * (t.mean(axis=0) - c.mean(axis=0))); ws.append(w)
        D[("arrest", lab)] = np.sum(ds, axis=0) / sum(ws)
    for ct in CTS:
        for tp in TPS:
            _, c2, g2 = load(f"{ct}_{tp}")
            D[("time_course", f"{ct}_{tp}")] = c2[g2 == tp].mean(axis=0) - c2[g2 == "none"].mean(axis=0)

    rng = np.random.default_rng(SEED)
    rows, gene_rows, dir_rows = [], [], []
    for dset in ("arrest", "time_course"):
        keys = [k for k in D if k[0] == dset]
        tops, hks = {}, {}
        for k in keys:
            v = D[k]
            order = np.argsort(v)[::-1]
            tops[k] = order[:TOPK]
            hk = half_net_k(v)
            hks[k] = order[:hk]
            for r, gi in enumerate(order[:TOPK], 1):
                gene_rows.append(dict(dataset=dset, group=k[1], rank=r, gene=feats[gi],
                                      contribution=v[gi], coef=coef[gi],
                                      in_half_net=r <= hk))
        print(f"\n=== {dset}: {len(keys)} groups ===")
        print(f"  genes needed for half the net: "
              f"{ {k[1]: len(hks[k]) for k in keys} }")

        # ---- recurrence on BOTH sides, same k, same draws --------------------
        bots = {k: np.argsort(D[k])[:TOPK] for k in keys}
        for dirn, lists in (("up", tops), ("down", bots)):
            cnt_d = {}
            for k in keys:
                for gi in lists[k]:
                    cnt_d[gi] = cnt_d.get(gi, 0) + 1
            allg_d = [gi for gi, c in cnt_d.items() if c == len(keys)]
            most_d = [gi for gi, c in cnt_d.items() if c >= max(3, len(keys) // 2)]
            nd_all, nd_most = [], []
            for _ in range(NDRAW):
                dr = [set(rng.choice(len(feats), TOPK, replace=False)) for _ in keys]
                nd_all.append(len(set.intersection(*dr)))
                kc = {}
                for d_ in dr:
                    for gi in d_:
                        kc[gi] = kc.get(gi, 0) + 1
                nd_most.append(sum(1 for v in kc.values() if v >= max(3, len(keys) // 2)))
            nd_all = np.array(nd_all); nd_most = np.array(nd_most)
            pa = (1 + (nd_all >= len(allg_d)).sum()) / (NDRAW + 1)
            pm = (1 + (nd_most >= len(most_d)).sum()) / (NDRAW + 1)
            print(f"  [{dirn:<4}] in all {len(keys)}: {len(allg_d):>3} (null {nd_all.mean():.3f}, "
                  f"p = {pa:.5f})   in >=3: {len(most_d):>3} (null {nd_most.mean():.3f}, p = {pm:.5f})")
            if allg_d:
                print(f"         genes: {', '.join(feats[g] for g in sorted(allg_d))}")
            dir_rows.append(dict(dataset=dset, direction=dirn, topk=TOPK, n_groups=len(keys),
                                 n_in_all=len(allg_d), null_in_all_mean=float(nd_all.mean()),
                                 p_in_all=pa, n_in_half=len(most_d),
                                 null_in_half_mean=float(nd_most.mean()), p_in_half=pm,
                                 p_floor=1 / (NDRAW + 1), n_draws=NDRAW,
                                 genes_in_all=";".join(feats[g] for g in sorted(allg_d))))

        # observed overlap
        obs = [jac(tops[a], tops[b]) for a, b in itertools.combinations(keys, 2)]
        # null: random groups of size TOPK from the clock's features
        nul = np.array([jac(rng.choice(len(feats), TOPK, replace=False),
                            rng.choice(len(feats), TOPK, replace=False)) for _ in range(NDRAW)])
        p = (1 + (nul >= np.median(obs)).sum()) / (NDRAW + 1)
        print(f"  top-{TOPK} Jaccard between groups: median {np.median(obs):.3f}  "
              f"(range {min(obs):.3f}-{max(obs):.3f})")
        print(f"  random groups of {TOPK}:            median {np.median(nul):.3f}   p = {p:.5f}")

        # recurrence
        cnt = {}
        for k in keys:
            for gi in tops[k]:
                cnt[gi] = cnt.get(gi, 0) + 1
        allg = [gi for gi, c in cnt.items() if c == len(keys)]
        most = [gi for gi, c in cnt.items() if c >= max(3, len(keys) // 2)]
        # null for BOTH recurrence counts in one loop, so they share a draw count and floor
        nul_all, nul_most = [], []
        for _ in range(NDRAW):
            draws = [set(rng.choice(len(feats), TOPK, replace=False)) for _ in keys]
            nul_all.append(len(set.intersection(*draws)))
            kc = {}
            for d in draws:
                for gi in d:
                    kc[gi] = kc.get(gi, 0) + 1
            nul_most.append(sum(1 for v in kc.values() if v >= max(3, len(keys) // 2)))
        nul_all = np.array(nul_all); nul_most = np.array(nul_most)
        p_all = (1 + (nul_all >= len(allg)).sum()) / (NDRAW + 1)
        p_most = (1 + (nul_most >= len(most)).sum()) / (NDRAW + 1)
        print(f"  genes in at least half the groups: {len(most)}   "
              f"(random expectation {nul_most.mean():.3f}, max {nul_most.max()}, p = {p_most:.5f})")
        print(f"  genes in the top-{TOPK} of ALL {len(keys)} groups: {len(allg)}   "
              f"(random expectation {nul_all.mean():.2f}, p = {p_all:.5f})")
        if allg:
            print(f"    recurrent genes: {', '.join(feats[g] for g in allg)}")
        rows.append(dict(dataset=dset, n_groups=len(keys), topk=TOPK,
                         jaccard_median=float(np.median(obs)),
                         jaccard_min=float(min(obs)), jaccard_max=float(max(obs)),
                         null_jaccard_median=float(np.median(nul)), p_jaccard=p,
                         n_in_all=len(allg), n_in_half=len(most),
                         null_in_all_mean=float(nul_all.mean()), p_in_all=p_all,
                         null_in_half_mean=float(nul_most.mean()),
                         null_in_half_max=int(nul_most.max()), p_in_half=p_most,
                         genes_in_all=";".join(feats[g] for g in allg),
                         p_floor=1 / (NDRAW + 1)))

    # ---- class-level signatures, arrest conditions only ------------------
    NCLASS = NDRAW   # one draw count for every null in this script
    CLASSES = {"senescence (SIPS, OIS, RS)": ["Stress-induced CS", "Oncogene-induced CS", "Replicative CS"],
               "quiescence (CICQ, SSCQ)": ["Contact_inhibited CQ", "Serum_starved CQ"]}
    crows, shared = [], {}
    print("\n=== class signatures, top-%d of every condition in the class ===" % TOPK)
    for cname, members in CLASSES.items():
        idx = [np.argsort(D[("arrest", LAB[m])])[::-1][:TOPK] for m in members]
        inter = set(idx[0])
        for x in idx[1:]:
            inter &= set(x)
        shared[cname] = inter
        nul = np.array([len(set.intersection(*[set(rng.choice(len(feats), TOPK, replace=False))
                                               for _ in members])) for _ in range(NCLASS)])
        pv = (1 + (nul >= len(inter)).sum()) / (NCLASS + 1)
        print(f"  {cname:<28}{len(inter):>3} genes shared   null mean {nul.mean():.2f}   p = {pv:.5f}")
        crows.append(dict(test="genes shared across the class", subset=cname,
                          n_conditions=len(members), topk=TOPK, observed=len(inter),
                          null_mean=float(nul.mean()), p_emp=pv, p_floor=1 / (NCLASS + 1),
                          n_draws=NCLASS, genes=";".join(feats[g] for g in sorted(inter))))
    a, b = (shared[k] for k in CLASSES)
    j = len(a & b) / len(a | b) if (a | b) else np.nan
    njac = []
    for _ in range(NCLASS):
        ra = set.intersection(*[set(rng.choice(len(feats), TOPK, replace=False)) for _ in range(3)])
        rb = set.intersection(*[set(rng.choice(len(feats), TOPK, replace=False)) for _ in range(2)])
        njac.append(len(ra & rb) / len(ra | rb) if (ra | rb) else 0.0)
    njac = np.array(njac)
    pj = (1 + (njac >= j).sum()) / (NCLASS + 1)
    print(f"  overlap between the two class lists: Jaccard {j:.3f}   null mean {njac.mean():.4f}"
          f"   p = {pj:.5f}   (shared genes: {len(a & b)})")
    crows.append(dict(test="overlap between the class lists", subset="senescence vs quiescence",
                      n_conditions=5, topk=TOPK, observed=j, null_mean=float(njac.mean()),
                      p_emp=pj, p_floor=1 / (NCLASS + 1), n_draws=NCLASS,
                      genes=";".join(feats[g] for g in sorted(a & b))))
    # how the >=3-condition genes divide by class: "at least three of five" pools genes
    # that are senescence-only, quiescence-heavy, or shared by everything
    a_top = {LAB[m]: set(np.argsort(D[("arrest", LAB[m])])[::-1][:TOPK]) for m in LAB}
    CSL, CQL = ["SIPS", "OIS", "RS"], ["CICQ", "SSCQ"]
    cnt2 = {}
    for lab, st in a_top.items():
        for gi in st:
            cnt2[gi] = cnt2.get(gi, 0) + 1
    print("\n=== composition of the genes recurring in at least three conditions ===")
    for gi in sorted([g for g, c in cnt2.items() if c >= 3]):
        pass
    comp = {}
    for gi in [g for g, c in cnt2.items() if c >= 3]:
        k = (sum(gi in a_top[x] for x in CSL), sum(gi in a_top[x] for x in CQL))
        comp[k] = comp.get(k, 0) + 1
    for (ncs, ncq), n_ in sorted(comp.items(), key=lambda t: (-t[0][0], -t[0][1])):
        print(f"  {ncs} senescence + {ncq} quiescence conditions: {n_} genes")
        crows.append(dict(test="composition of genes in >=3 conditions",
                          subset=f"{ncs} senescence + {ncq} quiescence", n_conditions=ncs + ncq,
                          topk=TOPK, observed=n_, null_mean=np.nan, p_emp=np.nan,
                          p_floor=np.nan, n_draws=np.nan, genes=""))
    pd.DataFrame(crows).to_csv(f"{rerun_dir}/signal_concentration_classes.csv", index=False)

    pd.DataFrame(dir_rows).to_csv(f"{rerun_dir}/signal_concentration_directions.csv", index=False)
    pd.DataFrame(rows).to_csv(f"{rerun_dir}/signal_concentration.csv", index=False)
    pd.DataFrame(gene_rows).to_csv(f"{rerun_dir}/signal_concentration_top_genes.csv", index=False)
    print(f"\nSaved -> signal_concentration.csv, signal_concentration_top_genes.csv")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
