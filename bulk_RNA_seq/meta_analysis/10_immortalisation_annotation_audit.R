# 10_immortalisation_annotation_audit.R
#
# Audits the hand-curated `immortalised` column in sample_metadata_RERUN.csv
# against the source annotation shipped with each study, and writes a corrected
# per-sample annotation.
#
# WHY. 09_immortalisation_confound.R rests entirely on this column. Spot-checking
# it found it internally inconsistent, so it is re-derived here from the raw
# text rather than trusted.
#
# GROUND TRUTH. Immortalisation status is a property of the CELL LINE, and each
# study in this meta-analysis uses one line, so the flag must be constant within
# a study. It is derived from `full_attributes` (the submitter's own
# cell line / cell type / source_name fields, e.g. "IMR90 hTERT",
# "Immortalized primary fibroblasts") and `abstract`, matched at study level for
# immortal|hTERT|h-TERT|SV40|large T|T-antigen|E6/E7.
#
# JUDGEMENT CALL, stated explicitly: an inducible ONCOGENE construct is not
# immortalisation. Studies annotated "IMR90 ER:RAS" (SRP046254, SRP113324,
# SRP113329) or carrying HRASG12V are therefore NOT counted as immortalised
# unless the line is separately described as hTERT/immortalised. This matters
# because the curated column appeared at first glance to be tracking "has a
# vector" rather than "is immortalised"; it is not - its errors are misses.
#
# WHAT THE AUDIT FINDS: 8 false negatives in 3 studies, and no false positives.
# All four studies flagged `yes` are confirmed by explicit hTERT text. The
# curated count of 20 immortalised samples should be 28.
#
# Output: rerun_outputs/immortalisation_annotation_corrected.csv

source("R/config.R")
suppressPackageStartupMessages(library(dplyr))

d <- read.csv(file.path(RERUN_DIR, "sample_metadata_RERUN.csv"))
PAT <- "immortal|hTERT|h-TERT|SV40|large T|T-antigen|E6/E7|E6-E7"

# Derive at STUDY level: a study is immortalised if its own annotation says so.
study_call <- d %>%
  group_by(study) %>%
  summarise(
    n = n(),
    cell_line = paste(sort(unique(cell_line)), collapse = "/"),
    curated = paste(sort(unique(immortalised)), collapse = "/"),
    # evidence taken from the submitter's fields and the study abstract
    evidence = {
      fa <- paste(unique(full_attributes), collapse = " ")
      ab <- paste(unique(abstract), collapse = " ")
      hit <- regmatches(fa, regexpr(paste0(".{0,60}(", PAT, ").{0,60}"), fa,
                                    ignore.case = TRUE, perl = TRUE))
      if (!length(hit)) {
        hit <- regmatches(ab, regexpr(paste0(".{0,60}(", PAT, ").{0,60}"), ab,
                                      ignore.case = TRUE, perl = TRUE))
      }
      if (length(hit)) trimws(gsub("\\s+", " ", hit[1])) else NA_character_
    },
    .groups = "drop"
  ) %>%
  mutate(corrected = ifelse(is.na(evidence), "no", "yes"))

cat("== per-study call ==\n")
print(as.data.frame(study_call %>% select(study, n, cell_line, curated, corrected)),
      row.names = FALSE)

cat("\n== studies where the curated column disagrees with its own source annotation ==\n")
bad <- study_call %>% filter(curated != corrected)
for (i in seq_len(nrow(bad))) {
  cat(sprintf("\n  %s (n=%d, line '%s'): curated '%s' -> corrected '%s'\n     evidence: %s\n",
              bad$study[i], bad$n[i], bad$cell_line[i], bad$curated[i],
              bad$corrected[i], bad$evidence[i]))
}

out <- d %>%
  select(external_id, study, cell_line, tissue, cell_substate,
         immortalised_curated = immortalised) %>%
  left_join(study_call %>% select(study, immortalised = corrected, evidence),
            by = "study") %>%
  mutate(changed = immortalised_curated != immortalised)

cat(sprintf("\n== samples changed: %d of %d ==\n", sum(out$changed), nrow(out)))
print(out %>% filter(changed) %>%
        select(external_id, study, cell_line, cell_substate,
               immortalised_curated, immortalised),
      row.names = FALSE)

cat("\n== counts before / after ==\n")
print(table(curated = out$immortalised_curated, corrected = out$immortalised))

write.csv(out, file.path(RERUN_DIR, "immortalisation_annotation_corrected.csv"),
          row.names = FALSE)
cat(sprintf("\nSaved -> %s\n",
            file.path(RERUN_DIR, "immortalisation_annotation_corrected.csv")))
