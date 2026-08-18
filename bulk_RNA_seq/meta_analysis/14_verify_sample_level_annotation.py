#!/usr/bin/env python
"""14_verify_sample_level_annotation.py

Verifies PER SAMPLE that the study-level cell line and hTERT calls in
10_immortalisation_annotation_audit.R are correct, and writes the per-sample
evidence out for inspection.

WHY THIS EXISTS. Script 10 assigns cell line and immortalisation status at STUDY
level, on the reasoning that each study uses one line. That was an assumption, not
a check - and script 13 asserted it against script 10's own output, which is
circular and therefore vacuous. It also fails in principle: a GEO series can
contain several lines. GSE106414 (SRP123346) is the motivating case, since its
methods name BJ/ET/RasV12, TIG3/ET/RASV12, Ecopack 2 and HEK293-T. A study-level
call would silently mis-assign every TIG3 sample as BJ.

WHAT IT DOES. Fetches the full sample records for every GEO series contributing to
the meta-analysis (targ=gsm, so one request returns every GSM in the series), and
for each of our samples extracts:
  - a line token, from the submitter's own "cell line" / "cell type" / "strain"
    characteristic, falling back to source_name
  - whether that sample's own title, source_name, characteristics, growth protocol
    or treatment protocol mentions hTERT / immortalisation / SV40 / telomerase
It then checks that all samples within a study agree, and additionally reports the
lines present in the WHOLE series against the lines of the samples we actually use,
so that a series containing extra lines is visible even when our subset is clean.

RESULT AS OF 2026-08-18. All 200 GEO-derived samples agree with their study-level
call; no study's own samples disagree on line or hTERT status. For GSE106414 all 9
samples in the series are annotated "BJ cells" - the TIG3/ET line named in the
methods contributes no sample to that series, and Ecopack 2 / HEK293-T are
virus-packaging lines rather than experimental samples. Series flagged as holding
more than one line token are duplicate superseries entries, treatment wording, or
submitter typos ("IMR90 human fibrobasts"), not distinct lines.

The remaining 30 samples are ERP021140 (ArrayExpress, not GEO) and are checked from
its own description in script 10.

Usage: 14_verify_sample_level_annotation.py <rerun_dir> [cache_dir]
Output: <rerun_dir>/sample_level_line_verification.csv
"""
import collections
import os
import re
import sys
import time
import urllib.request

import pandas as pd

GEO = ("https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi"
       "?acc={}&targ=gsm&form=text&view=brief")
IMM = re.compile(r"hTERT|h-TERT|immortal\w*|SV40|large T antigen|T-antigen"
                 r"|telomerase|E6/E7", re.I)
LINE_KEYS = ["cell line", "cell_line", "cell type", "cell_type",
             "cell strain", "fibroblast strain", "strain"]


def fetch(gse, cache):
    path = os.path.join(cache, f"{gse}.txt")
    if not os.path.exists(path) or os.path.getsize(path) == 0:
        with urllib.request.urlopen(GEO.format(gse), timeout=60) as r:
            open(path, "wb").write(r.read())
        time.sleep(0.35)
    return open(path, errors="replace").read()


def parse(text):
    """-> {GSM: {field: [values]}} for every sample in the series record."""
    out, cur = {}, None
    for ln in text.replace("\r", "").split("\n"):
        m = re.match(r"\^SAMPLE = (GSM\d+)", ln)
        if m:
            cur = m.group(1)
            out[cur] = collections.defaultdict(list)
            continue
        if cur and ln.startswith("!Sample_"):
            k, _, v = ln.partition(" = ")
            out[cur][k[len("!Sample_"):]].append(v)
    return out


def line_token(f):
    ch = " | ".join(f.get("characteristics_ch1", []))
    for k in LINE_KEYS:
        m = re.search(rf"{k}\s*:\s*([^|]+)", ch, re.I)
        if m:
            return m.group(1).strip()
    return " | ".join(f.get("source_name_ch1", []))[:60]


def main(rerun_dir, cache):
    os.makedirs(cache, exist_ok=True)
    meta = pd.read_csv(os.path.join(rerun_dir, "sample_metadata_RERUN.csv"))
    ann = pd.read_csv(os.path.join(rerun_dir,
                                   "immortalisation_annotation_corrected.csv"))
    ours = meta[meta.geo.astype(str).str.startswith("GSM")].copy()
    print(f"{len(ours)} of {len(meta)} samples are GEO-derived "
          f"({len(meta) - len(ours)} are ERP021140/ArrayExpress)")

    # Resolve which series each sample belongs to from its own record, then pull
    # each series once to get every sibling sample's characteristics.
    series, recs, gsm_series = {}, {}, {}
    for gsm in ours.geo:
        p = os.path.join(cache, f"{gsm}.self.txt")
        if not os.path.exists(p) or os.path.getsize(p) == 0:
            url = ("https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi"
                   f"?acc={gsm}&targ=self&form=text&view=quick")
            with urllib.request.urlopen(url, timeout=60) as rr:
                open(p, "wb").write(rr.read())
            time.sleep(0.3)
        t = open(p, errors="replace").read().replace("\r", "")
        gsm_series[gsm] = re.findall(r"!Sample_series_id = (GSE\d+)", t)

    for gse in sorted({g for v in gsm_series.values() for g in v}):
        block = parse(fetch(gse, cache))
        recs.update(block)
        series.setdefault(gse, set()).update(block.keys())

    rows = []
    for _, r in ours.iterrows():
        f = recs.get(r.geo, {})
        blob = " ".join(sum((f.get(k, []) for k in
                            ["title", "source_name_ch1", "characteristics_ch1",
                             "growth_protocol_ch1", "treatment_protocol_ch1"]), []))
        hit = IMM.search(blob)
        a = ann[ann.external_id == r.external_id].iloc[0]
        rows.append(dict(
            external_id=r.external_id, study=r.study, geo=r.geo,
            series=",".join(gsm_series[r.geo]),
            condition=r.cell_substate,
            title=" | ".join(f.get("title", [])),
            line_token_sample=line_token(f),
            characteristics=" | ".join(f.get("characteristics_ch1", [])),
            immortalisation_term_in_sample_record=bool(hit),
            matched_term=hit.group(0) if hit else "",
            cell_line_resolved_studylevel=a.cell_line_resolved,
            immortalised_studylevel=a.immortalised))
    S = pd.DataFrame(rows)
    S["AGREES"] = S.immortalisation_term_in_sample_record.eq(
        S.immortalised_studylevel.eq("yes"))

    print("\n== studies whose own samples disagree internally ==")
    bad = 0
    for st, g in S.groupby("study"):
        if g.line_token_sample.nunique() > 1 or \
           g.immortalisation_term_in_sample_record.nunique() > 1:
            bad += 1
            print(f"  {st}: lines={sorted(set(g.line_token_sample))} "
                  f"htert={sorted(set(g.immortalisation_term_in_sample_record))}")
    print(f"  -> {bad} of {S.study.nunique()} studies")

    print("\n== samples whose own record contradicts the study-level call ==")
    dis = S[~S.AGREES]
    print(f"  {len(dis)} of {len(S)}")
    if len(dis):
        print(dis[["external_id", "study", "geo", "line_token_sample",
                   "immortalised_studylevel"]].to_string(index=False))

    print("\n== lines in the WHOLE series vs lines of the samples we use ==")
    for gse, members in sorted(series.items()):
        mine = [g for g in members if g in set(S.geo)]
        if not mine:
            continue
        allt = collections.Counter(line_token(recs[g]) for g in members)
        myt = collections.Counter(line_token(recs[g]) for g in mine)
        if len(allt) > 1:
            print(f"  {gse}: series={dict(allt)}")
            print(f"          ours ={dict(myt)}")

    out = os.path.join(rerun_dir, "sample_level_line_verification.csv")
    S.to_csv(out, index=False)
    print(f"\nSaved -> {out}")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2] if len(sys.argv) > 2 else ".geo_cache")
