#!/usr/bin/env python
"""20_gse175533_export.py

Extracts the WI-38 replicative-senescence time course and its hTERT-immortalised
counterpart from GSE175533 (SRP321317) and writes a single expression matrix plus
sample annotation, for 21_gse175533_tage.R to push through tAge_preprocessing.

WHY THIS DATASET. Script 16 established that our own immortalisation contrast is
not estimable: 0 of 34 studies contain both immortalised and primary cells, so the
+31.1-unit scaled_diff gap between immortalised and primary proliferating controls
is inseparable from a between-study difference (study baselines among proliferating
controls span 87.8 units, script 11). Tyshkovskiy et al. report the opposite-looking
result - hTERT "abolished the increase in tAge" - from a REANALYSIS of this
dataset, not from an experiment of their own. GSE175533 therefore contains the only
one-laboratory, one-strain, immortalised-versus-primary comparison available to us,
and it can settle whether our +31.1 survives when study is held constant.

WHAT IS ACTUALLY IN THE DATASET (verified per sample against each sample's own GEO
characteristics, not inferred from titles):

  parental WI-38  30 bulk samples, PD 20/25/29/33/37/45/46/50/52/53, TP1-TP10, n=3 each
  hTERT           18 bulk samples, PD 46/51/64/73/86/109,            TP1,2,4,5,6,7, n=3 each

Two consequences for the design, both of which constrain what can be asked:

1. PD RANGES BARELY OVERLAP. The hTERT arm STARTS at PD 46, which is essentially
   the parental Hayflick limit here (the parental course ends at PD 53). So a
   PD-matched contrast exists only at PD 46 (hTERT TP1 vs parental PDL46_TP6) and
   near PD 51 (hTERT TP2 vs parental PD 50/52/53) - and at those PDs the parental
   cells are at or approaching replicative senescence, not proliferating.

   NOTE, CORRECTED 2026-08-31: an earlier version of this note read the PD 46
   start as "direct evidence that the hTERT cells were long-established before
   sampling began; it is NOT a fresh transduction." That does not follow. PD is
   counted from the WI-38 origin, so PD 46 says 46 doublings had accumulated in
   the lineage -- it says nothing about when hTERT was introduced, which could
   have been at PD 20 or at PD 46. Nothing in the deposited metadata records the
   transduction timing. What the PD range does bound is the total culture history:
   at most 109 doublings, far short of what BJ-hTERT, Tig3ET or IMR90-hTERT carry
   after decades of passage, so this dataset is not a like-for-like analogue of
   the established lines in our own pool.

2. TP LABELS ARE CALENDAR TIME, NOT PD. hTERT TP1 is PD 46 while parental TP1 is
   PD 20, so the shared TP schedule pairs samples by time in culture while the
   immortalised cells run far ahead in doublings. Pairing by TP and pairing by PD
   are therefore different questions and both are reported downstream.

Because our own claim concerns PROLIFERATING controls, the primary contrast is
hTERT (dividing throughout, PD 46-109) against parental cells while still dividing
(PD 20-37), with the PD-matched late contrasts reported alongside as the
complementary, oppositely-confounded view.

ONE ANNOTATION ERROR, corrected here and recorded: GSM5340976 is titled
hTERT_TP7_C and sits in the hTERT_TPM sheet, but its "cell line" characteristic
reads WI-38 rather than hTERT. It is assigned to the hTERT arm on the strength of
its title and its sheet membership. Taking the characteristic at face value is what
produced an earlier false claim that hTERT and parental cells share six population
doublings; they share one.

INPUT IS TPM, NOT COUNTS. GEO deposits no count matrix for this series (the only
raw archive is a 14.8 GB tar that is mostly ATAC-seq), so the SALMON raw-TPM sheets
are used. TPM_ij = counts_ij / (len_i * L_j) * 1e6, and tAge_preprocessing
z-scores every gene across samples before subtracting the control mean, which
removes any gene-specific constant - the gene-length term included - exactly. The
per-sample factor L_j is what RLE normalisation handles regardless. So the
scaled_diff matrix should be near-invariant to the counts/TPM choice; script 22
tests that rather than assuming it. TPM columns are rescaled to a nominal 20M
library so that the count_threshold=10 detection filter keeps its intended meaning.

Sheets used: hTERT_TPM and RS_TC_raw_TPM (both "RAW TPM values (SALMON)" per the
workbook's own readme). RS_TC_batch_corrected_TPM is NOT used - its readme says
"used only for visualization *Figure 1E".

Usage: 20_gse175533_export.py <xlsx> <out_dir>
Output: <out_dir>/gse175533_expression.csv, gse175533_samples.csv
"""
import os
import sys

import numpy as np
import pandas as pd

from xlsx_strict import rows

LIB = 20e6  # nominal library size the TPM columns are rescaled to

# PD per sample, read from each sample's own GEO characteristics (see docstring).
#
# TWO NAMING DISCREPANCIES between the workbook columns and the GEO titles, both
# resolved in favour of the population doubling, which the two sources agree on:
#   - the workbook drops the TP suffix for four late samples (PDL45, PDL52, PDL53)
#     and numbers PDL50 as TP7 where GEO calls it TP8 (likewise TP9/TP10). The TP
#     index is therefore taken from GEO, keyed on PDL, not from the column name.
#   - RS_PDL28_TP3 is titled PDL28 but its "population doublings" characteristic
#     reads 29. The characteristic is used.
PD = {"hTERT_TP1": 46, "hTERT_TP2": 51, "hTERT_TP4": 64, "hTERT_TP5": 73,
      "hTERT_TP6": 86, "hTERT_TP7": 109,
      "PDL20_TP1": 20, "PDL25_TP2": 25, "PDL28_TP3": 29, "PDL33_TP4": 33,
      "PDL37_TP5": 37, "PDL45": 45, "PDL46_TP6": 46, "PDL50_TP7": 50,
      "PDL52": 52, "PDL53": 53}
# TP index as GEO records it, keyed on (arm, PD). It MUST be keyed on the arm as
# well as the doubling: PD 46 occurs in both arms - hTERT's first timepoint and the
# parental sixth - so a PD-only key silently relabels hTERT TP1 as TP6, which
# collides with hTERT's real TP6 and corrupts any model using time as the x axis.
GEO_TP = {("hTERT", 46): "TP1", ("hTERT", 51): "TP2", ("hTERT", 64): "TP4",
          ("hTERT", 73): "TP5", ("hTERT", 86): "TP6", ("hTERT", 109): "TP7",
          ("parental", 20): "TP1", ("parental", 25): "TP2", ("parental", 29): "TP3",
          ("parental", 33): "TP4", ("parental", 37): "TP5", ("parental", 45): "TP7",
          ("parental", 46): "TP6", ("parental", 50): "TP8", ("parental", 52): "TP9",
          ("parental", 53): "TP10"}
# the parental course reaches its replicative limit at PD ~50; cells at PD <= 37
# are still dividing and are the parental "proliferating" set
PROLIF_MAX_PD = 37


def read_sheet(path, sheet):
    it = rows(path, sheet)
    hdr = next(it)
    cols = [c for c in hdr[1:] if c]
    genes, vals = [], []
    for r in it:
        if not r or r[0] is None:
            continue
        genes.append(str(r[0]))
        v = [np.nan if (i + 1) >= len(r) or r[i + 1] is None else float(r[i + 1])
             for i in range(len(cols))]
        vals.append(v)
    df = pd.DataFrame(vals, index=genes, columns=cols)
    print(f"  {sheet}: {df.shape[0]} genes x {df.shape[1]} samples")
    return df


def main(xlsx, out_dir):
    os.makedirs(out_dir, exist_ok=True)
    print("reading strict-OOXML sheets (openpyxl cannot open this workbook)")
    ht = read_sheet(xlsx, "hTERT_TPM")
    rs = read_sheet(xlsx, "RS_TC_raw_TPM")

    # Duplicated gene symbols must be dropped from each arm BEFORE intersecting:
    # .loc on a de-duplicated intersection still expands them, silently emitting
    # more rows than genes and mis-pairing the two arms.
    for nm, df in (("hTERT", ht), ("parental", rs)):
        d = df.index[df.index.duplicated(keep=False)]
        if len(d):
            print(f"  {nm}: dropping {len(d)} rows with duplicated gene symbols")
    ht = ht[~ht.index.duplicated(keep=False)]
    rs = rs[~rs.index.duplicated(keep=False)]

    shared = ht.index.intersection(rs.index)
    print(f"\ngene overlap after de-duplication: {len(shared)}")
    ht, rs = ht.loc[shared], rs.loc[shared]
    assert ht.index.is_unique and ht.index.equals(rs.index)

    E = pd.concat([ht, rs], axis=1)
    # TPM -> pseudo-counts on a nominal library, so count_threshold=10 means what
    # it means everywhere else in this project
    tot = E.sum(axis=0)
    print(f"\ncolumn sums before rescale: {tot.min():.0f} to {tot.max():.0f} "
          f"(TPM columns should sum to ~1e6)")
    E = E.divide(tot, axis=1) * LIB

    rows_ = []
    for s in E.columns:
        stem = s.rsplit("_", 1)[0]
        arm = "hTERT" if s.startswith("hTERT") else "parental"
        pd_ = PD.get(stem)
        if pd_ is None:
            raise SystemExit(f"no population doubling recorded for {s} (stem {stem})")
        if arm == "hTERT":
            state = "dividing"
        else:
            state = "dividing" if pd_ <= PROLIF_MAX_PD else "late"
        rows_.append(dict(sample_id=s, arm=arm, stem=stem, population_doublings=pd_,
                          timepoint=GEO_TP[(arm, pd_)], state=state,
                          # the control group for tAge_preprocessing: parental cells
                          # while still dividing, i.e. the closest analogue of the
                          # "Proliferating" controls used throughout 2.1.5
                          group=("parental_dividing" if arm == "parental" and state == "dividing"
                                 else f"{arm}_{state}")))
    S = pd.DataFrame(rows_)

    bad = S.groupby(["arm", "timepoint"]).population_doublings.nunique()
    if (bad > 1).any():
        raise SystemExit(f"timepoint label collision:\n{bad[bad > 1]}")
    print(f"\n{S.groupby('arm').timepoint.nunique().to_dict()} distinct timepoints per arm; "
          "no (arm, timepoint) maps to more than one doubling")

    print("\nsample groups:")
    print(S.groupby(["arm", "state"]).agg(
        n=("sample_id", "size"),
        PD=("population_doublings", lambda x: f"{x.min()}-{x.max()}")).to_string())
    print("\nPD-matched pairs available (same PD, both arms):")
    for p in sorted(set(S.population_doublings[S.arm == "hTERT"]) &
                    set(S.population_doublings[S.arm == "parental"])):
        print(f"  PD {p}: hTERT n={sum((S.arm=='hTERT')&(S.population_doublings==p))}, "
              f"parental n={sum((S.arm=='parental')&(S.population_doublings==p))}")

    ep = os.path.join(out_dir, "gse175533_expression.csv")
    E.round(4).to_csv(ep)
    S.to_csv(os.path.join(out_dir, "gse175533_samples.csv"), index=False)
    print(f"\nSaved -> {ep} ({E.shape[0]} genes x {E.shape[1]} samples)")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
