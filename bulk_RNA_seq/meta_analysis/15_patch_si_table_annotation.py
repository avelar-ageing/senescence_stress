#!/usr/bin/env python
"""15_patch_si_table_annotation.py

Corrects the `immortalised` and `cell_line` columns IN PLACE in the supplementary
tables, rather than shipping a second table alongside them.

WHY. The SI tables carry the same two defective columns audited in
meta_analysis/10 and verified per sample in meta_analysis/14:

  `immortalised` - wrong for 26 of 230 samples. hTERT status is a property of the
  cell line and every study here contributes one line, so it must be constant
  within a study; it was not. The errors are all false negatives, in SRP017378
  (whose curated cell_line is literally "BJ hTERT"), SRP066947 (Tig3ET, 18
  samples), SRP123346 and SRP136727.

  `cell_line` - not wrong so much as misleading, in three ways. "Primary" (n=69)
  pools 8 unrelated strains from 8 studies AND includes HCA2-hTert, which is not
  primary at all. "IMR-90" (n=107) includes the 6 IMR90-hTERT samples of
  SRP127037. "TIG-3" is Tig3ET, an hTERT derivative, so the parental name conceals
  the immortalisation. Any reader stratifying on this column gets strata that mix
  immortalised with primary cells.

WHAT IS WRITTEN. `immortalised` is replaced with the GEO-verified value and
`cell_line` with the resolved strain. The original coarse label is preserved as
`cell_line_asubmitted` so nothing is lost, and `immortalisation_evidence` carries
the quotation that establishes each call, so the correction is auditable inside the
same table without adding a new one.

Source of truth: rerun_outputs/immortalisation_annotation_corrected.csv.
Backups written alongside each file as <name>.pre_annotation_fix.csv.

Usage: 15_patch_si_table_annotation.py <si_tables_dir> <rerun_dir>
"""
import os
import shutil
import sys

import pandas as pd

TARGETS = ["sample_metadata.csv", "study_info_all.csv",
           "study_info_cq.csv", "study_info_cs.csv"]


def main(si_dir, rerun_dir):
    ver = pd.read_csv(os.path.join(
        rerun_dir, "immortalisation_annotation_corrected.csv")).set_index("external_id")
    for name in TARGETS:
        path = os.path.join(si_dir, name)
        if not os.path.exists(path):
            print(f"  skip (absent): {name}")
            continue
        d = pd.read_csv(path)
        if "external_id" not in d.columns:
            print(f"  skip (no external_id): {name}")
            continue
        miss = set(d.external_id) - set(ver.index)
        if miss:
            print(f"  WARNING {name}: {len(miss)} ids absent from the verified table")
        shutil.copy(path, path.replace(".csv", ".pre_annotation_fix.csv"))
        v = ver.reindex(d.external_id)
        n_imm = int((d["immortalised"].values != v["immortalised"].values).sum()) \
            if "immortalised" in d.columns else 0
        if "cell_line" in d.columns:
            # keep the submitted label; do not overwrite it silently
            cols = list(d.columns)
            d["cell_line_assubmitted"] = d["cell_line"]
            d["cell_line"] = v["cell_line_resolved"].values
            i = cols.index("cell_line")
            d = d[cols[:i + 1] + ["cell_line_assubmitted"] + cols[i + 1:]]
        n_line = int((d["cell_line"].values != d["cell_line_assubmitted"].values).sum()) \
            if "cell_line_assubmitted" in d.columns else 0
        d["immortalised"] = v["immortalised"].values
        d["immortalisation_evidence"] = v["evidence"].values
        d.to_csv(path, index=False)
        print(f"  {name}: rows={len(d)}  immortalised corrected={n_imm}  "
              f"cell_line refined={n_line}")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
