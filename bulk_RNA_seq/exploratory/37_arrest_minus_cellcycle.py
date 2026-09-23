#!/usr/bin/env python3
"""37_arrest_minus_cellcycle.py

Is the tAge rise in arrested cells reducible to the shutdown of proliferation?

Every condition in 2.1.5 is a growth arrest, and 32/36 showed the cell-cycle sets
are the ones carried by a fixed group of genes in every one of them. So the
question the section never asks is whether the whole effect is cell-cycle exit
read by a clock. If it is, "arrested cells are transcriptionally aged" means
something much narrower than it appears to.

THE TEST. Recompute each condition's difference from its own study's proliferating
controls using only the clock features that are NOT in a cell-cycle set, and test
it the way 2.1.5.1 tests the whole clock: the condition label permuted within each
study 20,000 times, so no sample is ever compared with a control from another
study. If the five conditions stay elevated on the reduced clock, the signal is not
just arrest. If they collapse, it is.

WHAT COUNTS AS CELL-CYCLE. The four Hallmark sets that behave as one block in
32/36 - E2F TARGETS, G2M CHECKPOINT, MITOTIC SPINDLE, MYC TARGETS - plus DNA
REPAIR, which sits with them on every measure. Their union is removed and the
number of features dropped is reported, because a reduction could otherwise do its
work simply by shrinking the clock. Two controls guard against that:

  size control  the same NUMBER of features removed at random, 20,000 draws, so
                the observed drop can be compared with an arbitrary block's.
  reverse       only the cell-cycle features kept. If arrest is the whole story
                this should carry the full effect.

IMPLEMENTATION NOTE. The within-study estimator is linear in the per-sample score,
so for fixed labels it is a dot product with a weight vector: +w_s/W for a treated
sample, -w_s/W for a control, zero elsewhere, with w_s = n_t n_c /(n_t + n_c) as
in 2.1.5.1. Two consequences are used here. A permutation needs only a fresh
weight vector and one dot product rather than a loop over studies. And because the
score is a sum over genes, any gene subset's difference is the SUM OF ITS GENES'
individual differences, so the 20,000-draw size control costs one gather per draw
instead of recomputing the clock. Both give identical numbers to the naive form.

Raw permutation p with its floor, no BH, matching 13_condition_within_study.R.

Output: rerun_outputs/arrest_minus_cellcycle.csv
"""
import sys, warnings
import joblib, numpy as np, pandas as pd
from sklearn.exceptions import InconsistentVersionWarning
warnings.filterwarnings("ignore", category=InconsistentVersionWarning)

NPERM = 20000
NSIZE = 20000
SEED = 1
CELLCYCLE = ["HALLMARK E2F TARGETS", "HALLMARK G2M CHECKPOINT", "HALLMARK MITOTIC SPINDLE",
             "HALLMARK MYC TARGETS", "HALLMARK DNA REPAIR"]
LAB = {"Contact_inhibited CQ": "CICQ", "Serum_starved CQ": "SSCQ", "Stress-induced CS": "SIPS",
       "Oncogene-induced CS": "OIS", "Replicative CS": "RS"}
CONTROL = "Proliferating"


def weights(is_t, is_c, study_idx, nstudies):
    """weight vector w such that w . score == the within-study mean difference"""
    w = np.zeros(len(is_t))
    den = 0.0
    for s in range(nstudies):
        ti = np.where(is_t & (study_idx == s))[0]
        ci = np.where(is_c & (study_idx == s))[0]
        if len(ti) == 0 or len(ci) == 0:
            continue
        ws = len(ti) * len(ci) / (len(ti) + len(ci))
        w[ti] += ws / len(ti); w[ci] -= ws / len(ci)
        den += ws
    return w / den if den else w


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
    grp = ann.reindex(sid)["cell_substate"].values
    studies = pd.unique(ann.reindex(sid)["study"].values)
    smap = {s: i for i, s in enumerate(studies)}
    sidx = np.array([smap[s] for s in ann.reindex(sid)["study"].values])

    pw = pd.read_csv(f"{PT}/hallmark_pathway_mouse_ids.csv")
    pw["mouse_gene_id"] = pw["mouse_gene_id"].astype(str)
    fi = {g: i for i, g in enumerate(feats)}
    cc = np.array(sorted({fi[g] for g in pw.loc[pw.pathway.isin(CELLCYCLE), "mouse_gene_id"] if g in fi}))
    keep = np.setdiff1d(np.arange(len(feats)), cc)
    print(f"clock features {len(feats)} | in a cell-cycle set {len(cc)} | remaining {len(keep)}")

    rng = np.random.default_rng(SEED)
    rows = []
    hdr = f"\n{'cond':<7}{'full diff':>11}{'p':>9}{'no-CC diff':>13}{'p':>9}{'CC-only diff':>15}{'p':>9}{'retained':>10}"
    print(hdr)
    for full, lab in LAB.items():
        is_t = grp == full
        is_c = grp == CONTROL
        w = weights(is_t, is_c, sidx, len(studies))
        # per-gene within-study difference: any subset's diff is the sum of these
        per_gene = w @ contrib
        o_full, o_red, o_cc = per_gene.sum(), per_gene[keep].sum(), per_gene[cc].sum()

        # permutation: relabel within study, rebuild w, one dot product per draw
        pool = np.where(is_t | is_c)[0]
        lab_pool = is_t[pool].copy()
        by_study = [np.where(sidx[pool] == s)[0] for s in range(len(studies))]
        cf = cr = cq = 0
        sf = contrib.sum(axis=1); sr = contrib[:, keep].sum(axis=1); sq = contrib[:, cc].sum(axis=1)
        for _ in range(NPERM):
            perm = lab_pool.copy()
            for ii in by_study:
                if len(ii) > 1:
                    perm[ii] = rng.permutation(lab_pool[ii])
            t2 = np.zeros(len(grp), bool); c2 = np.zeros(len(grp), bool)
            t2[pool[perm]] = True; c2[pool[~perm]] = True
            w2 = weights(t2, c2, sidx, len(studies))
            cf += abs(w2 @ sf) >= abs(o_full)
            cr += abs(w2 @ sr) >= abs(o_red)
            cq += abs(w2 @ sq) >= abs(o_cc)
        pf, pr, pq = [(1 + c) / (NPERM + 1) for c in (cf, cr, cq)]
        print(f"{lab:<7}{o_full:>+11.3f}{pf:>9.5f}{o_red:>+13.3f}{pr:>9.5f}"
              f"{o_cc:>+15.3f}{pq:>9.5f}{o_red/o_full:>10.2f}")

        # size control, using the per-gene decomposition: drop |cc| genes at random
        drops = np.array([rng.choice(len(feats), size=len(cc), replace=False) for _ in range(NSIZE)])
        null_red = o_full - per_gene[drops].sum(axis=1)
        p_size = (1 + (np.abs(null_red - o_full) >= abs(o_red - o_full)).sum()) / (NSIZE + 1)
        rows.append(dict(condition=lab, n_features=len(feats), n_cellcycle=len(cc),
                         diff_full=o_full, p_full=pf,
                         diff_no_cellcycle=o_red, p_no_cellcycle=pr,
                         diff_cellcycle_only=o_cc, p_cellcycle_only=pq,
                         frac_retained=o_red / o_full,
                         random_removal_median=float(np.median(null_red)),
                         p_size_control=p_size, p_floor=1 / (NPERM + 1)))

    D = pd.DataFrame(rows)
    print(f"\nsize control: removing {len(cc)} features at random, {NSIZE} draws\n")
    print(f"  {'cond':<7}{'observed drop':>15}{'random drop (median)':>23}{'p':>10}")
    for _, r in D.iterrows():
        print(f"  {r['condition']:<7}{r['diff_full']-r['diff_no_cellcycle']:>+15.3f}"
              f"{r['diff_full']-r['random_removal_median']:>+23.3f}{r['p_size_control']:>10.5f}")
    D.to_csv(f"{rerun_dir}/arrest_minus_cellcycle.csv", index=False)
    print(f"\nSaved -> arrest_minus_cellcycle.csv")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
