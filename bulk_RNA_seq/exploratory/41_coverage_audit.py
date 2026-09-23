#!/usr/bin/env python3
"""41_coverage_audit.py

Audits the one number in 2.1.5.3 that no script derived.

"Only 3,071 of the clock's 10,487 features sit in any Hallmark set" was a hardcoded
literal in 40_figure_signal_concentration.R (lines 57, 88). This recomputes it from the
model and the pathway map, and settles the comparison the text now makes.

WHAT IS COMPUTED
  clock coverage      how many of the mortality clock's features are in >= 1 Hallmark
                      set, over all 10,487 features and over the subset actually
                      MEASURED in each dataset. The measured subset is the right
                      denominator for any statement about the genes that carry the
                      shift, because a gene absent from the expression matrix gets the
                      imputed mean in both groups and so contributes exactly zero to
                      every difference; it can never enter a top-contributor list.
  top-gene coverage   the same for the genes carrying half of each condition's net,
                      pooled across conditions within a dataset.
  binomial            P(X <= observed unannotated) under the measured-clock coverage.
                      The point is DIRECTION: the earlier draft of 2.1.5.3 stated the
                      coverage gap as biting harder on the shift-carrying genes than on
                      the clock at large. If the top genes sit at or below the clock's
                      own unannotated share, that is one limitation appearing twice, not
                      a second and worse one.

Output: rerun_outputs/coverage_audit.csv
"""
import sys, warnings
import joblib, numpy as np, pandas as pd
from scipy.stats import binom, fisher_exact
from sklearn.exceptions import InconsistentVersionWarning
warnings.filterwarnings("ignore", category=InconsistentVersionWarning)


def main(rerun_dir, model_dir):
    PT = f"{rerun_dir}/partial_tage"
    m = joblib.load(f"{model_dir}/EN_Mortality_Multispecies_Multitissue_scaleddiff.pkl")
    feats = set(map(str, m.feature_names_in_))

    pw = pd.read_csv(f"{PT}/hallmark_pathway_mouse_ids.csv")
    hall = set(pw["mouse_gene_id"].astype(str))
    ann_clock = feats & hall

    print(f"clock features                       {len(feats):>6}")
    print(f"unique mouse ids in the pathway map  {len(hall):>6}")
    print(f"clock features in >= 1 Hallmark set  {len(ann_clock):>6}"
          f"   <- manuscript says 3,071  [{'MATCH' if len(ann_clock)==3071 else 'MISMATCH'}]")
    print(f"  share in no Hallmark set           {100*(1-len(ann_clock)/len(feats)):>5.1f}%"
          f"   <- manuscript says 71%")

    rows, BASE = [], {}
    rows.append(dict(scope="all clock features", dataset="", subset="", n=len(feats),
                     annotated=len(ann_clock), unannotated=len(feats) - len(ann_clock),
                     pct_unannotated=100 * (1 - len(ann_clock) / len(feats)),
                     baseline_pct_unannotated=np.nan, expected_unannotated=np.nan, p_binom=np.nan))

    top = pd.read_csv(f"{rerun_dir}/signal_concentration_top_genes.csv")
    top["gene"] = top["gene"].astype(str)

    CTS = ("Fibroblast", "Keratinocyte", "Melanocyte")
    TPS = ("4_days", "10_days", "20_days")
    STEM = {"arrest": ["meta"],
            "time_course": [f"{ct}_{tp}" for ct in CTS for tp in TPS]}
    for dset, stems in STEM.items():
        # MEASURED, as 23_decomposition_diagnostics.py defines it: the column exists AND
        # holds at least one non-missing value. Every clock feature is a column in these
        # matrices, but ~1,380 are entirely NaN and get the imputed mean, so they cancel
        # exactly between the two groups of any comparison. Taking the union across the
        # dataset's groups, since a gene measured in any group can enter that group's list.
        meas = set()
        for st in stems:
            e = pd.read_csv(f"{PT}/{st}_scaled_diff.csv")
            e.columns = e.columns.map(str)
            meas |= {c for c in e.columns if c in feats and e[c].notna().any()}
        base = 1 - len(meas & hall) / len(meas)
        BASE[dset] = base
        print(f"\n{dset}")
        print(f"  measured clock features            {len(meas):>6}")
        print(f"  of those, in no Hallmark set       {100*base:>5.1f}%   <- the right baseline")

        # the dataset as a whole, then the arrest classes and each condition on their
        # own, since a pooled union can hide one group carrying the coverage
        subsets = [("all groups pooled", None)]
        if dset == "arrest":
            subsets += [("quiescence (CICQ, SSCQ)", ["CICQ", "SSCQ"]),
                        ("senescence (SIPS, OIS, RS)", ["SIPS", "OIS", "RS"])]
            subsets += [(c, [c]) for c in ("CICQ", "SSCQ", "SIPS", "OIS", "RS")]
        else:
            # cell type is the time course's analogue of the arrest class
            groups9 = [f"{ct}_{tp}" for ct in CTS for tp in TPS]
            subsets += [(ct, [f"{ct}_{tp}" for tp in TPS]) for ct in CTS]
            subsets += [(g, [g]) for g in groups9]
        for sub_name, groups in subsets:
            sel = (top.dataset == dset) & top.in_half_net
            if groups is not None:
                sel &= top.group.isin(groups)
            gs = set(top.loc[sel, "gene"])
            u = len(gs - hall)
            pv = float(binom.cdf(u, len(gs), base))
            print(f"    {sub_name:<27} {len(gs):>3} genes  {u:>3} unannotated "
                  f"({100*u/len(gs):>3.0f}%)  expected {len(gs)*base:>4.1f}  P = {pv:.3f}")
            rows.append(dict(scope="genes carrying half the net", dataset=dset,
                             subset=sub_name, n=len(gs), annotated=len(gs & hall),
                             unannotated=u, pct_unannotated=100 * u / len(gs),
                             baseline_pct_unannotated=100 * base,
                             expected_unannotated=len(gs) * base, p_binom=pv))
        rows.append(dict(scope="measured clock features", dataset=dset, subset="",
                         n=len(meas), annotated=len(meas & hall),
                         unannotated=len(meas) - len(meas & hall),
                         pct_unannotated=100 * base, baseline_pct_unannotated=np.nan,
                         expected_unannotated=np.nan, p_binom=np.nan))

    # ---- the OTHER gene lists in 2.1.5.3: top-50 recurrence, not half-net ----
    # These select differently and are covered differently, which is why the section
    # reports them separately. Same baselines as above.
    print("\n=== top-50 recurrence lists (a different selection from the half-net lists) ===")
    for dset, base in BASE.items():
        keys = sorted(set(top.loc[top.dataset == dset, "group"]))
        if dset == "arrest":
            defs = [("all five conditions", keys),
                    ("both quiescence conditions", ["CICQ", "SSCQ"]),
                    ("all three senescence conditions", ["SIPS", "OIS", "RS"])]
        else:
            defs = [("all nine groups", keys)] + \
                   [(f"all three {ct} timepoints", [f"{ct}_{tp}" for tp in TPS]) for ct in CTS]
        for nm, gs in defs:
            inter = set.intersection(*[set(top.loc[top.group == g, "gene"]) for g in gs])
            if not inter:
                print(f"    {dset:<12} {nm:<34} 0 genes")
                continue
            u = len(inter - hall)
            pv = float(binom.cdf(u, len(inter), base))
            print(f"    {dset:<12} {nm:<34} {len(inter):>3} genes  {u:>3} unannotated "
                  f"({100*u/len(inter):>3.0f}%)  expected {len(inter)*base:>4.1f}  P = {pv:.3f}")
            rows.append(dict(scope="genes in the top 50 of every group", dataset=dset,
                             subset=nm, n=len(inter), annotated=len(inter & hall),
                             unannotated=u, pct_unannotated=100 * u / len(inter),
                             baseline_pct_unannotated=100 * base,
                             expected_unannotated=len(inter) * base, p_binom=pv))

    # is the coverage of the two arrest classes' gene lists different from each other?
    D = pd.DataFrame(rows)
    q = D[(D.subset == "quiescence (CICQ, SSCQ)")].iloc[0]
    c = D[(D.subset == "senescence (SIPS, OIS, RS)")].iloc[0]
    tab = [[int(q.unannotated), int(q.annotated)], [int(c.unannotated), int(c.annotated)]]
    fp = fisher_exact(tab)[1]
    print(f"\nquiescence vs senescence gene lists, unannotated share: "
          f"{q.pct_unannotated:.0f}% ({int(q.n)} genes) against {c.pct_unannotated:.0f}% "
          f"({int(c.n)} genes), Fisher p = {fp:.3f}")

    pd.DataFrame(rows).to_csv(f"{rerun_dir}/coverage_audit.csv", index=False)
    print(f"\nSaved -> coverage_audit.csv")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
