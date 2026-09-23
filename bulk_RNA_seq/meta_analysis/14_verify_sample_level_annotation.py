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
# TISSUE, added 2026-08-31. This script verified cell line and hTERT status only,
# so the tissue column was never checked against the submitters' own records - and
# it is wrong for at least one line: sample_metadata_RERUN.csv calls HDF 12-3
# "Skin", while the paper the samples come from (Mitra et al. 2018, Genome Biol)
# says the whole HDF series is foreskin. Since tissue is one of the four
# predictors of baseline tAge in meta_analysis/11, an unchecked Skin/Foreskin
# split makes that result a partly arbitrary division.
TISSUE_KEYS = ["tissue", "source tissue", "organ", "anatomic site", "body site",
               "cell origin", "origin"]
TISSUE_WORDS = re.compile(r"foreskin|prepuce|lung|dermal|dermis|skin|embryo|fetal|"
                          r"foetal|neonat|newborn", re.I)


def is_soft(text):
    """GEO answers some clients with an HTML reCAPTCHA page (HTTP 200). Such a page
    parses to zero samples, so every sample would be scored as 'no term found'.
    Only text that carries SOFT sample records counts."""
    return "!Sample_" in text and not text.lstrip().lower().startswith(("<!doctype", "<html"))


def cached_get(url, path, pause):
    """Return the SOFT text at url, reading path if it holds a valid record. A cached
    file that is not SOFT (e.g. a stored CAPTCHA page) is ignored and re-fetched; a
    fetch that returns something other than SOFT stops the script and is NOT cached."""
    if os.path.exists(path):
        t = open(path, errors="replace").read()
        if is_soft(t):
            return t
    with urllib.request.urlopen(url, timeout=60) as r:
        raw = r.read()
    t = raw.decode("utf-8", errors="replace")
    if not is_soft(t):
        sys.exit(f"GEO returned a non-SOFT response for {url} (first bytes: {t[:80]!r}). "
                 "This is usually a reCAPTCHA page. Retry later or from another network; "
                 "nothing was written to the cache.")
    open(path, "wb").write(raw)
    time.sleep(pause)
    return t


def fetch(gse, cache):
    return cached_get(GEO.format(gse), os.path.join(cache, f"{gse}.txt"), 0.35)


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


def tissue_token(f):
    """Tissue as the submitter recorded it, plus any provenance words in free text."""
    ch = " | ".join(f.get("characteristics_ch1", []))
    explicit = ""
    for k in TISSUE_KEYS:
        m = re.search(rf"{k}\s*:\s*([^|]+)", ch, re.I)
        if m:
            explicit = m.group(1).strip(); break
    blob = " ".join(sum((f.get(k, []) for k in
                        ["title", "source_name_ch1", "characteristics_ch1",
                         "growth_protocol_ch1", "extract_protocol_ch1"]), []))
    words = sorted({w.lower() for w in TISSUE_WORDS.findall(blob)})
    return explicit, ",".join(words)


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
        url = ("https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi"
               f"?acc={gsm}&targ=self&form=text&view=quick")
        t = cached_get(url, os.path.join(cache, f"{gsm}.self.txt"), 0.3).replace("\r", "")
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
            tissue_stated=tissue_token(f)[0],
            tissue_words_in_record=tissue_token(f)[1],
            tissue_metadata=r.get("tissue", ""),
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

    print("\n== TISSUE: what the metadata says vs what the sample record says ==")
    tt = S[["study", "cell_line_resolved_studylevel", "tissue_metadata",
            "tissue_stated", "tissue_words_in_record"]].drop_duplicates()
    tt = tt.sort_values("cell_line_resolved_studylevel")
    print(tt.to_string(index=False))
    blank = int((tt.tissue_stated.fillna("") == "").sum())
    print(f"\n  sample records stating a tissue explicitly: {len(tt) - blank} of {len(tt)}")
    mism = tt[(tt.tissue_words_in_record != "") &
              tt.apply(lambda r: isinstance(r.tissue_metadata, str) and r.tissue_metadata != ""
                       and r.tissue_metadata.lower() not in r.tissue_words_in_record, axis=1)]
    print(f"  metadata tissue not among the words in the record: {len(mism)}")
    if len(mism):
        print(mism.to_string(index=False))

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
    # NULL-RESULT GUARD (2026-09-01). This script depends on a populated .geo_cache;
    # GEO now serves captchas to curl, so a run with an empty cache finds no
    # characteristics and every immortalisation_term_in_sample_record comes back False.
    # On 2026-08-31 exactly that happened and the null result overwrote a good output
    # file, turning 46 genuine hits into zero and silently breaking the assertion in
    # meta_analysis/13. Refuse to overwrite in that case: a run that finds no evidence
    # at all has failed, and must not be mistaken for evidence of absence.
    _found = int((S["immortalisation_term_in_sample_record"] == True).sum()) \
        if "immortalisation_term_in_sample_record" in S.columns else 0
    if _found == 0 and os.path.exists(out):
        raise SystemExit(
            f"REFUSING TO WRITE {out}: found 0 immortalisation terms across {len(S)} samples, "
            "which means the GEO cache is empty rather than that the terms are absent. "
            "Populate .geo_cache (the 34 SOFT records are on origin/restructure-bulk-rna-seq "
            "under DISCREPANCY_REPORT/evidence/geo_soft/) and re-run. The existing file has "
            "been left untouched.")
    S.to_csv(out, index=False)
    print(f"\nSaved -> {out}")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2] if len(sys.argv) > 2 else ".geo_cache")
