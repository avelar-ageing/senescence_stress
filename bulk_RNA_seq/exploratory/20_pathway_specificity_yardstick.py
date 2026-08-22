#!/usr/bin/env python
"""20_pathway_specificity_yardstick.py

Gives each gene set a YARDSTICK, so its contribution can be read as more or less
than expected rather than merely non-zero.

THE PROBLEM. A set's partial-tAge contribution is the slice of a condition's tAge
shift carried by that set's genes. When the whole transcriptome shifts - OIS is
+27 units - nearly every set's genes shift with it, so "is this set's contribution
non-zero?" is almost always yes: 26-43 of 50 sets reach FDR < 0.05 in a given
condition. That is one global shift detected forty times, not forty findings.

THE YARDSTICK. Contribution is a sum of coef_i * z_i over the set, so what a set
should contribute if it were unremarkable is fixed by how much of the clock's
coefficient weight it holds:

    expected_s = (sum |coef| over set s / sum |coef| over all features) * total shift

The ratio observed/expected is then interpretable directly: 1.0 means the set
carries exactly its share of the shift, 2.0 twice its share. This is deterministic
and needs no null.

SIGNIFICANCE, and why the earlier null was rejected but this one is not. Random
sets are drawn from the clock's own features, matched to the real set on gene
count AND on the distribution of |coef| (stratified by coefficient decile, with
zero-coefficient genes forming their own stratum, since those contribute exactly
nothing and a set's effective size is its non-zero count). Matching on weight
rather than count alone matters because contribution scales with |coef|.

The size-matched null in 12_partial_tage_size_matched_null.py was excluded because
its pool is genes already selected for age-association across mammals, so it is
not centred on "no effect" and non-survival could not be read as absence of
signal. That objection applies to the question it was asked - "does this set have
an effect" - and NOT to the question here, which is "does this set carry more of
the shift than an arbitrary equally-weighted slice of the same clock". For a
question about disproportion, other clock genes are the correct comparator, and
there is no alternative: a non-clock gene has coefficient exactly 0 and cannot
contribute. So the same construction is valid here and invalid there.

DESIGN. Everything is computed on the WITHIN-STUDY contrasts of
exploratory/18: per gene, the precision-weighted mean of per-study
(arrested - control) differences. Any set's contribution is then a sum over that
per-gene vector, which makes the real sets, the matched random sets and the whole
transcriptome all the same operation.

PERMUTATION COUNT AND THE BH FAMILY. An empirical p cannot fall below
1/(B+1), so B must be large enough that BH can reach 0.05. The discovery family is
the 17 INTERPRETABLE sets x 5 conditions = 85 tests per model - the sets whose
contributions are interpreted at set level at all - so the smallest adjusted value
needs p_emp near 0.05/85 = 5.9e-4 and B = 20,000 leaves ample margin. The other 33
sets are computed and reported for the record but excluded from the family and left
unadjusted, because they are not interpreted. An earlier run used B = 500 with BH
over all 250 tests per model, which made significance arithmetically impossible
(floor 2.0e-3 against a required 2.0e-4) and returned zero hits for that reason
alone rather than for any property of the data.

Matched sets are drawn WITH replacement within each stratum, which allows the null
to be vectorised. For the pool and draw sizes here this inflates the null variance
by under 2% relative to sampling without replacement (finite-population correction
1 - (k-1)/(P-1)), so it is very slightly conservative.

Reported per set x condition x model: observed contribution, expected from weight
share, ratio, z against the matched null, empirical two-sided p, and BH within
model across the interpretable family.

Usage: 20_pathway_specificity_yardstick.py <rerun_dir> <model_dir> [n_null]
Output: <rerun_dir>/pathway_specificity_yardstick.csv
"""
import sys
import warnings

import joblib
import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")

ADJ = 122.5
SEED = 1
MODELS = {"scaled": ("EN_Chronoage_Multispecies_Multitissue_scaleddiff.pkl",
                     "meta_scaled_diff.csv"),
          "yugene": ("EN_Chronoage_Multispecies_Multitissue_yugenediff.pkl",
                     "meta_yugene_diff.csv")}
# The mortality clock is a separate instrument, not a third normalisation: pure
# ridge (l1_ratio 0) so all 10,487 features carry weight, against 1,839 for the
# sparse chronological clock. Run with --mortality. No species adjustment is
# applied to it (see meta_analysis/18_mortality_clock.py).
MORTALITY = {"mortality": ("EN_Mortality_Multispecies_Multitissue_scaleddiff.pkl",
                           "meta_scaled_diff.csv")}
CONDS = {"Contact_inhibited CQ": "CICQ", "Serum_starved CQ": "SSCQ",
         "Replicative CS": "RS", "Stress-induced CS": "SIPS",
         "Oncogene-induced CS": "OIS"}


def _patch(imp):
    if not hasattr(imp, "_fill_dtype"):
        imp._fill_dtype = imp.statistics_.dtype if hasattr(imp, "statistics_") else np.float64


def per_gene_within_study(C, groups, cond):
    """Precision-weighted mean of per-study (cond - Proliferating) per-gene diffs."""
    acc = np.zeros(C.shape[1])
    wsum = 0.0
    for st, g in groups.groupby("study"):
        ti = g.index[g.group == cond].to_numpy()
        ci = g.index[g.group == "Proliferating"].to_numpy()
        if not len(ti) or not len(ci):
            continue
        w = len(ti) * len(ci) / (len(ti) + len(ci))
        acc += w * (C[ti].mean(axis=0) - C[ci].mean(axis=0))
        wsum += w
    return acc / wsum


def main(rerun_dir, model_dir, n_null=2000, which="chronoage"):
    models = MORTALITY if which == "mortality" else MODELS
    adj = 1.0 if which == "mortality" else ADJ
    suffix = "_mortality" if which == "mortality" else ""
    PT = f"{rerun_dir}/partial_tage"
    rng = np.random.default_rng(SEED)
    pw = pd.read_csv(f"{PT}/hallmark_pathway_mouse_ids.csv")
    meta = pd.read_csv(f"{rerun_dir}/sample_metadata_RERUN.csv")
    study_of = dict(zip(meta.external_id, meta.study))

    rows = []
    for mdl, (pkl, expr) in models.items():
        m = joblib.load(f"{model_dir}/{pkl}")
        _patch(m.named_steps["imputation"])
        feats = list(map(str, m.feature_names_in_))
        idx = {g: i for i, g in enumerate(feats)}
        coef = m.named_steps["estimator"].coef_
        aco = np.abs(coef)
        wtot = aco.sum()

        e = pd.read_csv(f"{PT}/{expr}")
        sid = e["sample_id"].values
        e = e.drop(columns=["sample_id"])
        e.columns = e.columns.map(str)
        for g in [g for g in feats if g not in e.columns]:
            e[g] = np.nan
        X = e.loc[:, feats]
        Z = m.named_steps["scaler"].transform(m.named_steps["imputation"].transform(X))
        C = Z * coef[np.newaxis, :] * adj

        groups = pd.read_csv(f"{PT}/meta_groups.csv")
        groups = groups.set_index(pd.Index(range(len(groups))))
        assert list(groups.sample_id) == list(sid), "row order mismatch"
        groups["study"] = groups.sample_id.map(study_of)

        # strata for matching: zero-coefficient genes, then deciles of |coef|
        nz = aco > 0
        strat = np.zeros(len(feats), dtype=int)
        q = np.quantile(aco[nz], np.linspace(0, 1, 11)[1:-1])
        strat[nz] = 1 + np.searchsorted(q, aco[nz])
        by_strat = {s: np.where(strat == s)[0] for s in np.unique(strat)}

        sets = {n: [idx[x] for x in s.mouse_gene_id.astype(str) if x in idx]
                for n, s in pw.groupby("pathway")}

        for cond_value, lab in CONDS.items():
            d = per_gene_within_study(C, groups, cond_value)
            total = d.sum()
            for name, cols in sets.items():
                if not cols:
                    continue
                cols = np.asarray(cols)
                obs = d[cols].sum()
                share = aco[cols].sum() / wtot
                exp = share * total
                # matched random sets: same count drawn from each |coef| stratum
                counts = {s: int((strat[cols] == s).sum()) for s in np.unique(strat[cols])}
                # vectorised: sum of k draws from each stratum's own d values
                draws = np.zeros(n_null)
                for st_id, k in counts.items():
                    vals = d[by_strat[st_id]]
                    draws += rng.choice(vals, size=(n_null, k), replace=True).sum(axis=1)
                # two-sided empirical p for the observed being extreme vs matched sets
                p = (1 + np.sum(np.abs(draws - draws.mean()) >= abs(obs - draws.mean()))) \
                    / (n_null + 1)
                rows.append(dict(
                    analysis="meta_analysis", label=lab, model=mdl, pathway=name,
                    n_clock=len(cols), n_nonzero=int(nz[cols].sum()),
                    weight_share=share, total_shift=total,
                    observed=obs, expected_from_weight=exp,
                    ratio=obs / exp if exp != 0 else np.nan,
                    null_mean=draws.mean(), null_sd=draws.std(),
                    z=(obs - draws.mean()) / draws.std() if draws.std() > 0 else np.nan,
                    p_emp=p))
            print(f"  {mdl} {lab} done (total shift {total:+.2f})")

    out = pd.DataFrame(rows)
    # RATIO RELIABILITY. observed/expected divides by (weight share x total shift),
    # so when a condition's whole-transcriptome shift is near zero the ratio
    # explodes and changes sign meaninglessly - CICQ on scaled_diff (+5.9 units)
    # gives ratios up to 52 and SSCQ on scaled_diff (-3.4) down to -63. Those are
    # arithmetic, not signal. z against the matched null is the stable statistic
    # and should be used instead; ratio is flagged unreliable below 10 units.
    out["ratio_reliable"] = out.total_shift.abs() >= 10
    # BH within model, over the interpretable family only (see header)
    # the gate is clock-specific (see exploratory/14): use the mortality gate for
    # the mortality clock, or the BH family is built on the wrong set list
    rep = pd.read_csv(f"{rerun_dir}/pathway_representation"
                      f"{'_mortality' if which == 'mortality' else ''}.csv")
    interp = set(rep.loc[rep.tier == "INTERPRETABLE", "pathway"])
    out["interpretable"] = out.pathway.isin(interp)
    out["p_emp_adj"] = np.nan
    for mdl in out.model.unique():
        sel = (out.model == mdl) & out.interpretable
        out.loc[sel, "p_emp_adj"] = _bh(out.loc[sel, "p_emp"].values)
    print(f"\nBH family per model: {int(out[out.model == out.model.iloc[0]].interpretable.sum())} "
          f"interpretable tests; {int((~out.interpretable).sum() / 2)} sets reported unadjusted")
    out.to_csv(f"{rerun_dir}/pathway_specificity_yardstick{suffix}.csv", index=False)
    print(f"\nSaved -> {rerun_dir}/pathway_specificity_yardstick{suffix}.csv")


def _bh(p):
    p = np.asarray(p, dtype=float)
    n = len(p)
    o = np.argsort(p)
    adj = np.empty(n)
    adj[o] = np.minimum.accumulate((p[o] * n / (np.arange(n) + 1))[::-1])[::-1]
    return np.clip(adj, 0, 1)


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2],
         int(sys.argv[3]) if len(sys.argv) > 3 else 2000,
         "mortality" if "--mortality" in sys.argv else "chronoage")
