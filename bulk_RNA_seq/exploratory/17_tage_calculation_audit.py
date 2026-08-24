#!/usr/bin/env python
"""17_tage_calculation_audit.py

Audit of how tAge is computed, and whether the reported quantities are
statistically sound. Reproduces every number in
DISCREPANCY_REPORT/TAGE_CALCULATION_AUDIT.md.

Four checks:

  A  COVERAGE. How much of the clock is operating on measured data? Reports
     features observed of 10,487, non-zero-coefficient features observed of
     1,839, and the share of total |coefficient| weight that is measured.
     Run-dependent, because filter_genes applies its threshold within each run.

  B  IMPUTATION IS DIFFERENCE-NEUTRAL. Unmeasured features are filled with the
     training median, so they take the SAME value in every sample of a run and
     cancel from any between-group difference. Verifies their contribution
     difference is exactly zero, and that observed features alone reproduce the
     whole-transcriptome difference.

  C  CONTROL-REFERENCE NOISE IS DIFFERENCE-NEUTRAL. control_subtraction removes
     a per-gene median estimated from the control samples (only 6 in the
     temporal runs). Bootstraps that reference and shows the reported
     differences do not move. Analytically expected: subtracting a per-gene
     constant c shifts test and control means equally, so mean_test - mean_ctrl
     is invariant. Consequently the Wilcoxon test is invariant too (all
     per-sample values shift by the same amount, so ranks are unchanged).

  D  CROSS-RUN COMPARABILITY. The genuine problem. Reports per-run gene
     retention, measured-weight share and input scale, which differ between the
     meta-analysis and temporal runs and make MAGNITUDES non-comparable across
     runs. Correlation-based analyses are unaffected (scale-invariant).

Conclusion: within a run the reported differences are sound and insensitive to
imputation and to control-reference noise. Across runs, magnitudes are not
comparable. Absolute tAge is not interpretable (see the audit doc on units).

Usage: 17_tage_calculation_audit.py <rerun_dir> <model_dir>
"""
import sys
import warnings

import joblib
import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")

ADJ = 122.5
B_BOOT = 200
SEED = 1
RUNS = ["meta", "Fibroblast", "Keratinocyte", "Melanocyte",
        "Fibroblast_4_days", "Fibroblast_20_days", "Melanocyte_20_days"]
GROUPS = {  # run -> (groups_file_stem, test_label, control_label)
    "meta": ("meta", "Stress-induced CS", "Proliferating"),
    "Fibroblast": ("Fibroblast", "irradiated", "none"),
    "Keratinocyte": ("Keratinocyte", "irradiated", "none"),
    "Melanocyte": ("Melanocyte", "irradiated", "none"),
    "Fibroblast_4_days": ("Fibroblast_4_days", "4_days", "none"),
    "Fibroblast_20_days": ("Fibroblast_20_days", "20_days", "none"),
    "Melanocyte_20_days": ("Melanocyte_20_days", "20_days", "none"),
}


def _patch(imp):
    if not hasattr(imp, "_fill_dtype"):
        imp._fill_dtype = imp.statistics_.dtype if hasattr(imp, "statistics_") else np.float64


def load_run(PT, run, feats):
    e = pd.read_csv(f"{PT}/{run}_scaled_diff.csv")
    sid = e["sample_id"].values
    e = e.drop(columns=["sample_id"])
    e.columns = e.columns.map(str)
    present = {c for c in e.columns if e[c].notna().any()}
    for g in [g for g in feats if g not in e.columns]:
        e[g] = np.nan
    return sid, e.loc[:, feats].to_numpy(dtype=float), present


def main(rerun_dir, model_dir):
    PT = f"{rerun_dir}/partial_tage"
    m = joblib.load(f"{model_dir}/EN_Chronoage_Multispecies_Multitissue_scaleddiff.pkl")
    _patch(m.named_steps["imputation"])
    feats = list(map(str, m.feature_names_in_))
    idx = {g: i for i, g in enumerate(feats)}
    coef = m.named_steps["estimator"].coef_
    nz_w = np.abs(coef).sum()

    pw = pd.read_csv(f"{PT}/hallmark_pathway_mouse_ids.csv")
    sets = {n: [idx[x] for x in s.mouse_gene_id.astype(str) if x in idx]
            for n, s in pw.groupby("pathway")}

    def contribs(D):
        Z = m.named_steps["scaler"].transform(
            m.named_steps["imputation"].transform(pd.DataFrame(D, columns=feats)))
        return Z * coef[np.newaxis, :] * ADJ

    print("=" * 78)
    print("A / D  COVERAGE AND CROSS-RUN COMPARABILITY")
    print("=" * 78)
    print(f"model features {len(feats)}, non-zero coefficients {(coef != 0).sum()}")
    print(f"{'run':<20}{'n':>4}{'obs':>7}{'obs%':>7}{'obsNZ':>7}{'weight%':>9}{'med|z|':>9}{'SD z':>8}")
    for run in RUNS:
        sid, D, present = load_run(PT, run, feats)
        obs = np.array([f in present for f in feats])
        w = np.abs(coef[obs & (coef != 0)]).sum()
        fin = D[np.isfinite(D)]
        print(f"{run:<20}{len(sid):>4}{obs.sum():>7}{100*obs.mean():>6.1f}%"
              f"{(obs & (coef!=0)).sum():>7}{100*w/nz_w:>8.1f}%"
              f"{np.median(np.abs(fin)):>9.3f}{np.std(fin):>8.3f}")
    print("\nGene retention, measured-weight share and input scale all differ by run,")
    print("so MAGNITUDES are not comparable across runs. Correlations are.")

    print("\n" + "=" * 78)
    print("B  IS IMPUTATION DIFFERENCE-NEUTRAL?")
    print("=" * 78)
    for run in ["Fibroblast_20_days", "meta"]:
        stem, tl, cl = GROUPS[run]
        sid, D, present = load_run(PT, run, feats)
        grp = pd.read_csv(f"{PT}/{stem}_groups.csv").set_index("sample_id")["group"]
        g = np.array([grp.get(s) for s in sid])
        a, b = g == tl, g == cl
        C = contribs(D)
        d = C[a].mean(axis=0) - C[b].mean(axis=0)
        obs = np.array([f in present for f in feats])
        print(f"{run}: {(~obs).sum()} imputed features")
        print(f"   max |contribution difference| among imputed : {np.abs(d[~obs]).max():.3e}")
        print(f"   observed-only total  {d[obs].sum():+.4f}   whole-transcriptome "
              f"{C.sum(axis=1)[a].mean() - C.sum(axis=1)[b].mean():+.4f}")

    print("\n" + "=" * 78)
    print("C  IS THE CONTROL REFERENCE ESTIMATE DIFFERENCE-NEUTRAL?")
    print("=" * 78)
    rng = np.random.default_rng(SEED)
    for run in ["Fibroblast_20_days", "Melanocyte_20_days"]:
        stem, tl, cl = GROUPS[run]
        sid, D, present = load_run(PT, run, feats)
        grp = pd.read_csv(f"{PT}/{stem}_groups.csv").set_index("sample_id")["group"]
        g = np.array([grp.get(s) for s in sid])
        a, b = g == tl, g == cl
        ci = np.where(b)[0]

        def diffs(Dm):
            C = contribs(Dm)
            dd = C[a].mean(axis=0) - C[b].mean(axis=0)
            return np.array([dd[c].sum() for c in sets.values()]), \
                   C.sum(axis=1)[a].mean() - C.sum(axis=1)[b].mean()

        base, base_w = diffs(D)
        boot, boot_w = [], []
        for _ in range(B_BOOT):
            r = rng.choice(ci, size=len(ci), replace=True)
            c = np.nanmedian(D[r], axis=0)
            s, w = diffs(D - c[np.newaxis, :])
            boot.append(s); boot_w.append(w)
        boot = np.array(boot)
        print(f"{run} (n control = {len(ci)}), B={B_BOOT}")
        print(f"   whole-transcriptome: observed {base_w:+.4f}, bootstrap SD {np.std(boot_w):.3e}")
        print(f"   set contributions  : max bootstrap SD {boot.std(axis=0).max():.3e}, "
              f"max |bias| {np.abs(boot.mean(axis=0) - base).max():.3e}")
    print("\nZero to floating-point precision, as expected analytically: subtracting a")
    print("per-gene constant shifts test and control means equally, so their difference")
    print("is invariant, and per-sample ranks are unchanged so Wilcoxon is invariant too.")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
