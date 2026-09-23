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

  `tissue` - wrong for 10 of 230 samples, in two strains, both called "Skin" when
  they are foreskin: HCA2-hTERT (6 samples, SRP089801; Cellosaurus CVCL_E2UW) and
  HDF 12-3 (4 samples, SRP153724; the HDF series' own paper, Mitra 2018, Genome
  Biol). After the fix the only strain left under "Skin" is HDF161, which really is
  adult dermis, and the totals go from Foreskin 73 / Lung 141 / Skin 16 to
  Foreskin 83 / Lung 141 / Skin 6. This column is what the tissue panel of SI
  Figures 1 and 16 is drawn from, so it has to be right at the source and not only
  inside the two figure scripts.

WHAT IS WRITTEN. `immortalised` is replaced with the GEO-verified value, `cell_line`
with the resolved strain, and `tissue` with the verified tissue. The original coarse
labels are preserved as `cell_line_assubmitted` and `tissue_assubmitted` so nothing
is lost, and `immortalisation_evidence` carries the quotation that establishes each
call, so the correction is auditable inside the same table without adding a new one.

IDEMPOTENT. Re-running is safe: a `*_assubmitted` column that already exists is left
alone rather than being overwritten with the already-corrected value, and an existing
backup is not replaced. Both would otherwise destroy the submitted labels on a second
run, which matters because this script has already been run once on these tables.

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
        bk = path.replace(".csv", ".pre_annotation_fix.csv")
        if not os.path.exists(bk):
            shutil.copy(path, bk)
        else:
            print(f"  {name}: keeping existing backup {os.path.basename(bk)}")
        v = ver.reindex(d.external_id)
        n_imm = int((d["immortalised"].values != v["immortalised"].values).sum()) \
            if "immortalised" in d.columns else 0
        def refine(col, source):
            """Overwrite `col` with the verified value, preserving the submitted
            label as <col>_assubmitted. Idempotent: an existing _assubmitted
            column is the true original and is never overwritten."""
            if col not in d.columns:
                return 0, list(d.columns)
            keep = col + "_assubmitted"
            cols = list(d.columns)
            if keep not in cols:
                d[keep] = d[col]
                d[col] = v[source].values
                i = cols.index(col)
                order = cols[:i + 1] + [keep] + cols[i + 1:]
            else:
                d[col] = v[source].values
                order = cols
            n = int((d[col].values != d[keep].values).sum())
            return n, order

        n_line, order = refine("cell_line", "cell_line_resolved")
        d = d[order]
        n_tissue, order = refine("tissue", "tissue_verified")
        d = d[order]
        d["immortalised"] = v["immortalised"].values
        d["immortalisation_evidence"] = v["evidence"].values
        d.to_csv(path, index=False)
        print(f"  {name}: rows={len(d)}  immortalised corrected={n_imm}  "
              f"cell_line refined={n_line}  tissue corrected={n_tissue}")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
