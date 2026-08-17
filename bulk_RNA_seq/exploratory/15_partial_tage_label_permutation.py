#!/usr/bin/env python
"""15_partial_tage_label_permutation.py

Label-permutation null for the partial-tAge decomposition.

WHAT IT ASKS. Hold the gene set fixed and permute the group labels: does this
set's partial score separate THESE groups more than chance? This is the question
the Results sections actually assert, and unlike the size-matched null
(12_..., excluded -- see DISCREPANCY_REPORT/PARTIAL_TAGE_VALIDITY.md Test 1) it
requires no choice of background gene pool, so there is no
"is-the-pool-neutral" objection.

WHAT IT DOES NOT DO -- read this before using it as a filter. It does NOT
control for gene-set size or for how concentrated a set's signal is, because the
set is held fixed. It therefore does not replace what the size-matched null was
attempting; it only removes that null's unfairness. And because the Wilcoxon
rank-sum test is itself an exact permutation test (on ranks rather than means),
this test shares the Wilcoxon's null hypothesis and is expected to agree with it
closely. The script reports that agreement explicitly so the redundancy is
visible rather than assumed. Its value, if any, is as a distribution-free check
on Cohen's d that makes no asymptotic assumption.

Resolution limit: for a 6v6 comparison only C(12,6) = 924 distinct label splits
exist, so the smallest attainable two-sided p is ~1/924 regardless of how many
permutations are drawn. The p-floor is reported per row.

Reads the already-computed per-sample partial scores, so no model loading and no
re-decomposition. Raw empirical p is reported with the floor and effect size; NO
multiple-testing correction is applied to Monte-Carlo null p-values, per project
convention (a pre-registered robustness check is not a discovery family).

Output: rerun_outputs/partial_tage_label_permutation.csv
"""
import sys
import warnings
from math import comb

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")

N_PERM = 2000
SEED = 20260817
META_CONDITIONS = {"Contact_inhibited CQ": "CICQ", "Serum_starved CQ": "SSCQ",
                   "Replicative CS": "RS", "Stress-induced CS": "SIPS",
                   "Oncogene-induced CS": "OIS"}
CELL_TYPES = ["Fibroblast", "Keratinocyte", "Melanocyte"]
TIMEPOINTS = ["4_days", "10_days", "20_days"]


def cohens_d_matrix(M, a, b):
    """Cohen's d per column of M, between row-masks a and b."""
    na, nb = a.sum(), b.sum()
    ma, mb = M[a].mean(axis=0), M[b].mean(axis=0)
    va, vb = M[a].var(axis=0, ddof=1), M[b].var(axis=0, ddof=1)
    pooled = np.sqrt(((na - 1) * va + (nb - 1) * vb) / (na + nb - 2))
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.where(pooled > 0, (ma - mb) / pooled, np.nan)


def run(scores_csv, groups_csv, test_label, control_label, label, analysis,
        model, rows, rng):
    scores = pd.read_csv(scores_csv)
    # pathway names are column headers and may contain spaces; pandas preserves them
    drop = ["sample_id", "full_tAge_direct", "full_tAge_reconstructed"]
    pathways = [c for c in scores.columns if c not in drop]
    groups = pd.read_csv(groups_csv).set_index("sample_id")["group"]
    g = np.array([groups.get(s, None) for s in scores["sample_id"]])

    keep = np.isin(g, [test_label, control_label])
    M = scores.loc[keep, pathways].to_numpy(dtype=float)
    gk = g[keep]
    a = gk == test_label
    b = gk == control_label
    if a.sum() < 2 or b.sum() < 2:
        print(f"  [skip] {label}/{model}: n too small", file=sys.stderr)
        return

    d_obs = cohens_d_matrix(M, a, b)

    n_tot, n_a = len(gk), int(a.sum())
    n_distinct = comb(n_tot, n_a)
    n_perm = min(N_PERM, max(200, n_distinct))

    count = np.zeros(len(pathways), dtype=int)
    idx = np.arange(n_tot)
    for _ in range(n_perm):
        perm = rng.permutation(idx)
        pa = np.zeros(n_tot, dtype=bool)
        pa[perm[:n_a]] = True
        d_p = cohens_d_matrix(M, pa, ~pa)
        count += np.abs(d_p) >= np.abs(d_obs)

    p_emp = (count + 1) / (n_perm + 1)
    p_floor = 2.0 / n_distinct
    for k, pw in enumerate(pathways):
        rows.append(dict(analysis=analysis, label=label, model=model, pathway=pw,
                         cohens_d=d_obs[k], p_permutation=p_emp[k],
                         p_floor=p_floor, n_perm=n_perm,
                         n_distinct_splits=n_distinct,
                         n_test=int(a.sum()), n_control=int(b.sum())))


if __name__ == "__main__":
    rerun_dir = sys.argv[1]
    PT = f"{rerun_dir}/partial_tage"
    rng = np.random.default_rng(SEED)
    rows = []
    for model in ["scaled", "yugene"]:
        for cond, short in META_CONDITIONS.items():
            print(f"[{model}] meta: {short}", file=sys.stderr)
            run(f"{PT}/meta_partial_scores_{model}.csv", f"{PT}/meta_groups.csv",
                cond, "Proliferating", short, "meta_analysis", model, rows, rng)
        for ct in CELL_TYPES:
            print(f"[{model}] pooled: {ct}", file=sys.stderr)
            run(f"{PT}/{ct}_partial_scores_{model}.csv", f"{PT}/{ct}_groups.csv",
                "irradiated", "none", ct, "temporal_pooled", model, rows, rng)
        for ct in CELL_TYPES:
            for tp in TIMEPOINTS:
                grp = f"{ct}_{tp}"
                print(f"[{model}] bytimepoint: {grp}", file=sys.stderr)
                run(f"{PT}/{grp}_partial_scores_{model}.csv", f"{PT}/{grp}_groups.csv",
                    tp, "none", grp, "temporal_bytimepoint", model, rows, rng)

    out = pd.DataFrame(rows)
    path = f"{rerun_dir}/partial_tage_label_permutation.csv"
    out.to_csv(path, index=False)
    print(f"\nSaved {path} ({len(out)} rows)", file=sys.stderr)
    print(f"p_permutation < 0.05: {(out.p_permutation < 0.05).sum()}/{len(out)} "
          f"({100 * (out.p_permutation < 0.05).mean():.1f}%)", file=sys.stderr)

    # Redundancy check against the Wilcoxon p already in partial_tage_ALL.csv
    try:
        w = pd.read_csv(f"{rerun_dir}/partial_tage_ALL.csv")
        m = out.merge(w[["analysis", "label", "model", "pathway", "wilcox_p", "p_adj"]],
                      on=["analysis", "label", "model", "pathway"], how="inner")
        if len(m):
            rho = m[["p_permutation", "wilcox_p"]].corr(method="spearman").iloc[0, 1]
            agree = ((m.p_permutation < 0.05) == (m.wilcox_p < 0.05)).mean()
            print(f"\nREDUNDANCY CHECK vs uncorrected Wilcoxon p (n={len(m)}):",
                  file=sys.stderr)
            print(f"  Spearman rho(p_permutation, wilcox_p) = {rho:.3f}", file=sys.stderr)
            print(f"  same call at 0.05 (uncorrected): {100 * agree:.1f}%", file=sys.stderr)
            agree_adj = ((m.p_permutation < 0.05) == (m.p_adj < 0.05)).mean()
            print(f"  same call vs BH-adjusted p_adj:  {100 * agree_adj:.1f}%", file=sys.stderr)
    except FileNotFoundError:
        pass
