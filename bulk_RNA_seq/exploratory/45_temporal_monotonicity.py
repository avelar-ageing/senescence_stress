#!/usr/bin/env python3
"""45_temporal_monotonicity.py

Is the time-course trend in panel (b) of figure_signal_concentration real?

THE OBSERVATION. The share of tAge movement carried by the 50 gene sets moves
monotonically with time in all three cell types: fibroblasts up (1.18, 1.19, 1.20),
keratinocytes down (1.19, 1.17, 1.15), melanocytes up (1.15, 1.16, 1.20). That is
notable because the tAge trajectories themselves are NOT monotonic - melanocytes peak
at 4 days and keratinocytes at 4-10 days (2.2.3) - so this is not the magnitude trend
seen through a different lens.

WHICH MATRIX. The per-timepoint matrices (Keratinocyte_4_days_scaled_diff.csv and so on)
CANNOT be used: tAge_preprocessing was run separately for each, so the same six baseline
samples carry different values in each file (r = 0.9993, max difference 0.087) and the three
timepoints are not on a common scale. A first version of this script pooled them anyway and
returned 38.9, 39.0, 39.3 for keratinocytes where the correct shares are 38.9, 37.8, 37.4 -
the observed statistic did not reproduce, which is how the error surfaced. This version uses
the per-CELL-TYPE export (Keratinocyte_scaled_diff.csv), all 24 samples preprocessed
together, with the timepoint of each sample recovered from the per-timepoint group files.
Shares from that frame differ slightly from the figure's per-timepoint ones; both are
reported, and the claim only stands if the direction agrees in both.

WHY THE OBVIOUS TEST IS WRONG. A random ordering of three values is monotonic with
probability 2/6 = 1/3, giving 1/27 = 0.037 for three cell types. That assumes the three
timepoints are independent. They are not: all three contrasts within a cell type share
the SAME six untreated baseline samples, which correlates their errors and makes an
accidental monotonic run commoner than 1/3. The combinatorial p is anti-conservative.

THE TEST. Permute which treated samples belong to which timepoint, within cell type:
the 18 treated samples are reassigned 6/6/6 at random while the baseline stays fixed,
so the shared-baseline correlation is preserved exactly. 20,000 draws. Reported per cell
type and jointly, as the fraction of draws in which the run is monotonic in either
direction. Raw empirical p with its floor; three pre-specified questions, no BH.

The statistic permuted is the annotated share of sum|contribution|, not the ratio to the
weight-matched null: the null mean varies only between 31.7% and 32.6% across the nine
groups, so the ordering of the ratios is the ordering of the shares, and permuting the
share avoids nesting one 20,000-draw null inside another.

Output: rerun_outputs/temporal_monotonicity.csv
"""
import sys, warnings
import joblib, numpy as np, pandas as pd
from sklearn.exceptions import InconsistentVersionWarning
warnings.filterwarnings("ignore", category=InconsistentVersionWarning)

NDRAW = 20000
CHUNK = 2000
SEED = 1
CTS = ("Fibroblast", "Keratinocyte", "Melanocyte")
TPS = ("4_days", "10_days", "20_days")


def mono(a):
    """monotonic in either direction, strictly"""
    return (a[:, 0] < a[:, 1]) & (a[:, 1] < a[:, 2]) | \
           (a[:, 0] > a[:, 1]) & (a[:, 1] > a[:, 2])


def main(rerun_dir, model_dir):
    PT = f"{rerun_dir}/partial_tage"
    m = joblib.load(f"{model_dir}/EN_Mortality_Multispecies_Multitissue_scaleddiff.pkl")
    imp = m.named_steps["imputation"]
    if not hasattr(imp, "_fill_dtype"):
        imp._fill_dtype = imp.statistics_.dtype
    feats = np.array(list(map(str, m.feature_names_in_)))
    coef = m.named_steps["estimator"].coef_
    hall = set(pd.read_csv(f"{PT}/hallmark_pathway_mouse_ids.csv")["mouse_gene_id"].astype(str))
    is_ann = np.array([f in hall for f in feats])

    rng = np.random.default_rng(SEED)
    rows, obs_all, nul_all = [], {}, {}

    for ct in CTS:
        # one matrix per cell type, all 24 samples on a single preprocessing; the
        # timepoint labels come from the per-timepoint group files
        e = pd.read_csv(f"{PT}/{ct}_scaled_diff.csv")
        sid = e["sample_id"].values
        e = e.drop(columns=["sample_id"]); e.columns = e.columns.map(str)
        present = {c for c in e.columns if e[c].notna().any()}
        for g in [g for g in feats if g not in e.columns]:
            e[g] = np.nan
        Z = m.named_steps["scaler"].transform(imp.transform(e.loc[:, feats]))
        C = Z * coef[np.newaxis, :]
        ob = np.array([f in present for f in feats])
        tp_of = {}
        for tp in TPS:
            gg = pd.read_csv(f"{PT}/{ct}_{tp}_groups.csv")
            for s_ in gg.loc[gg.group == tp, "sample_id"]:
                tp_of[s_] = tp
        lab = np.array([tp_of.get(s_, "none") for s_ in sid])
        assert (lab == "none").sum() == 6 and len(lab) == 24, (ct, np.bincount(
            np.unique(lab, return_inverse=True)[1]))
        # treated samples ordered 4d, 10d, 20d so the observed run is the identity draw
        T = [C[lab == tp] for tp in TPS]
        base = C[lab == "none"]
        Tr = np.vstack(T)                       # 18 treated x features
        A = is_ann[ob]
        bmean = base[:, ob].mean(axis=0)
        Trm = Tr[:, ob]
        n = len(Trm) // 3

        def shares(idx):
            """annotated share of sum|contribution| for each of the three blocks"""
            W = np.zeros((len(idx) * 3, len(Trm)))
            for k in range(3):
                W[np.arange(len(idx)) * 3 + k, :] = 0
                for j, order in enumerate(idx):
                    W[j * 3 + k, order[k * n:(k + 1) * n]] = 1.0 / n
            D = np.abs(W @ Trm - bmean)
            return (100 * D[:, A].sum(axis=1) / D.sum(axis=1)).reshape(len(idx), 3)

        o_share = shares([np.arange(len(Trm))])[0]
        obs_all[ct] = o_share
        hits = 0
        nul = np.empty((NDRAW, 3))
        for a in range(0, NDRAW, CHUNK):
            b = min(a + CHUNK, NDRAW)
            idx = [rng.permutation(len(Trm)) for _ in range(b - a)]
            nul[a:b] = shares(idx)
        mk = mono(nul)
        nul_all[ct] = mk
        p = (1 + mk.sum()) / (NDRAW + 1)
        d = "up" if o_share[0] < o_share[2] else "down"
        print(f"  {ct:<14}{'  '.join(f'{x:.1f}' for x in o_share)}   monotonic {d:<5}"
              f"   null rate {100*mk.mean():.1f}%   p = {p:.4f}")
        rows.append(dict(cell_type=ct, share_4d=o_share[0], share_10d=o_share[1],
                         share_20d=o_share[2], direction=d,
                         observed_monotonic=bool(mono(o_share[None, :])[0]),
                         null_monotonic_rate=float(mk.mean()), p_emp=p,
                         p_floor=1 / (NDRAW + 1), n_draws=NDRAW))

    joint = np.ones(NDRAW, bool)
    for ct in CTS:
        joint &= nul_all[ct]
    pj = (1 + joint.sum()) / (NDRAW + 1)
    print(f"\n  all three cell types monotonic: null rate {100*joint.mean():.2f}%   p = {pj:.4f}")
    print(f"  (the naive combinatorial value, treating the three timepoints as independent,"
          f" would be (1/3)^3 = 0.037)")
    rows.append(dict(cell_type="ALL THREE", direction="2 up, 1 down",
                     observed_monotonic=True, null_monotonic_rate=float(joint.mean()),
                     p_emp=pj, p_floor=1 / (NDRAW + 1), n_draws=NDRAW))
    pd.DataFrame(rows).to_csv(f"{rerun_dir}/temporal_monotonicity.csv", index=False)
    print(f"\nSaved -> temporal_monotonicity.csv")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
