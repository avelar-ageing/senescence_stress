# 10_immortalisation_annotation_audit.R
#
# Re-derives two sample attributes that 09/11/12/13 depend on and that the
# curated metadata gets wrong: hTERT-immortalisation status, and cell line at
# strain level.
#
# ============================ WHY, AND WHAT WAS WRONG ========================
#
# (1) THE CURATED `immortalised` COLUMN IS UNRELIABLE. Immortalisation is a
# property of the cell line and each study here uses one line, so the flag must
# be constant within a study. It was not: SRP123346 had 2 of 3 samples flagged
# one way and 1 the other. Checking the rest found 8 samples in 3 studies
# flagged `no` whose own recount3 metadata says otherwise, including SRP017378
# whose cell_line is literally "BJ hTERT".
#
# (2) recount3's `full_attributes` IS NOT COMPLETE, so text-matching the local
# metadata is not sufficient either. An earlier version of this script did that
# and missed SRP066947 (n=18) entirely: its samples are Tig3ET - "Human diploid
# fibroblast Tig3 cells expressing the ecotropic receptor, hTERT, and the
# (pBabe-hygro) empty vector" - but that sentence lives in GEO's
# Sample_treatment_protocol field, which recount3 does not carry. The published
# paper (Genome Res 27:1634, GSE76605) describes the line as "hTERT-immortalized
# Tig3 fibroblasts (Tig3ET)".
#
# Status is therefore taken from the GEO SOFT records (Sample + Series, fetched
# for all 34 studies), not from the local text. Evidence quotes are recorded
# below. Verification covered every study, not only those the local text
# flagged: all 34 sample-level records and all 34 series-level records were
# searched for hTERT|h-TERT|immortal*|SV40|telomerase|E6/E7. 31 of 34 sample
# records carry a growth protocol and every series record carries a summary and
# overall design, so absence of a hit is informative rather than a gap. The
# three studies without a sample growth protocol (SRP040745, SRP045867,
# SRP172671) were resolved at series level.
#
# ERP021140 is ArrayExpress (E-MTAB-5403), not GEO. Its own description names
# primary HCA2 (from O. Pereira-Smith) and BJ (ATCC CRL-2522) with no hTERT, and
# the study includes replicative senescence, which mortal cells are required
# for. Treated as primary.
#
# RESULT: 46 immortalised samples across 8 studies, against 20 in the curated
# column. The 26 remaining studies show no immortalisation language at either
# level.
#
# JUDGEMENT CALL, stated explicitly: an inducible ONCOGENE construct is not
# immortalisation. IMR90 ER:RAS (SRP046254, SRP113324, SRP113329) and HRASG12V
# or BRAFV600E transductions are NOT counted unless the line is separately
# described as hTERT. For these the evidence is absence of a mention rather than
# a positive statement of primary status, which is weaker; it is the best
# available from the deposited records.
#
# (3) THE `cell_line` COLUMN IS TOO COARSE TO STRATIFY ON. "Primary" (n=69)
# spans 8 unrelated strains in 8 studies - HCA2/BJ, foreskin 2DD, HCA2-hTert,
# HDF161, and dermal strains 12-1, 12-3, 10-2, 10-5. "IMR-90" (n=107) spans 15
# studies, one of which (SRP127037) is IMR90-hTERT, so an "IMR-90" stratum mixes
# immortalised and primary cells in both arms. A strain-level `cell_line_resolved`
# is built here. Note it is still nested within study, so it does not remove
# study/batch effects - see script 13.
#
# Output: rerun_outputs/immortalisation_annotation_corrected.csv

source("R/config.R")
suppressPackageStartupMessages(library(dplyr))

d <- read.csv(file.path(RERUN_DIR, "sample_metadata_RERUN.csv"))

# ---- GEO-verified immortalisation, with the quote that establishes it --------
IMMORTAL <- tibble::tribble(
  ~study,      ~line,          ~evidence,
  "SRP017378", "BJ-hTERT",     "Sample_growth_protocol: 'Immortalized human BJ primary fibroblast cells (by hTERT expression)'",
  "SRP066917", "BJ-hTERT",     "Sample_characteristics: 'cell type: hTERT-immortalized BJ cells'",
  "SRP066947", "Tig3ET",       "Sample_treatment_protocol: 'Human diploid fibroblast Tig3 cells expressing the ecotropic receptor, hTERT'; Genome Res 27:1634 calls it 'hTERT-immortalized Tig3 fibroblasts (Tig3ET)'",
  "SRP089801", "HCA2-hTERT",   "Sample_characteristics: 'cell type: HCA2-hTert'",
  "SRP123346", "BJ-hTERT",     "Sample_growth_protocol: 'Immortalized human BJ primary fibroblast cells (by hTERT expression)'",
  "SRP127037", "IMR90-hTERT",  "Sample_characteristics: 'cell line: IMR90 hTERT'",
  "SRP136727", "BJ-hTERT",     "Sample_growth_protocol: 'Immortalized human BJ primary fibroblast cells (by hTERT expression)'",
  "SRP172671", "BJ-hTERT",     "Sample_characteristics: 'cell type: Immortalized BJ-hTERT cells'"
)

# ---- strain-level cell line for the primary studies -------------------------
PRIMARY_LINE <- c(
  ERP021140 = "HCA2/BJ (primary)", SRP017142 = "WI-38",   SRP034163 = "IMR90",
  SRP034541 = "IMR90",   SRP040243 = "IMR90",   SRP040745 = "BJ",
  SRP045867 = "MRC-5",   SRP046254 = "IMR90 ER:RAS", SRP052706 = "foreskin 2DD",
  SRP060598 = "IMR90",   SRP062872 = "IMR90",   SRP064207 = "WI-38",
  SRP065206 = "IMR90",   SRP069768 = "WI-38",   SRP070636 = "IMR90",
  SRP096629 = "HDF161",  SRP098713 = "IMR90",   SRP113324 = "IMR90 ER:RAS",
  SRP113329 = "IMR90 ER:RAS", SRP117883 = "IMR90", SRP121031 = "IMR90",
  SRP136071 = "IMR90",   SRP153205 = "HDF 12-1", SRP153724 = "HDF 12-3",
  SRP154382 = "HDF 10-2", SRP154577 = "HDF 10-5"
)

stopifnot(setequal(unique(d$study), c(IMMORTAL$study, names(PRIMARY_LINE))))

ann <- d %>%
  transmute(external_id, study, tissue, cell_substate,
            cell_line_curated = cell_line,
            immortalised_curated = immortalised) %>%
  mutate(
    immortalised = ifelse(study %in% IMMORTAL$study, "yes", "no"),
    cell_line_resolved = ifelse(
      study %in% IMMORTAL$study,
      IMMORTAL$line[match(study, IMMORTAL$study)],
      unname(PRIMARY_LINE[study])),
    evidence = ifelse(study %in% IMMORTAL$study,
                      IMMORTAL$evidence[match(study, IMMORTAL$study)],
                      "no hTERT/immortal/SV40/telomerase mention in GEO sample or series record"),
    changed = immortalised_curated != immortalised)

cat("== immortalisation: curated vs GEO-verified ==\n")
print(table(curated = ann$immortalised_curated, verified = ann$immortalised))
cat(sprintf("\nsamples changed: %d of %d, in %d studies\n",
            sum(ann$changed), nrow(ann), n_distinct(ann$study[ann$changed])))
print(ann %>% filter(changed) %>% count(study, cell_line_curated, cell_line_resolved,
                                        immortalised_curated, immortalised),
      row.names = FALSE)

cat("\n== immortalised samples per condition (verified) ==\n")
print(table(ann$cell_substate, ann$immortalised))

cat("\n== cell_line: curated label vs resolved strain ==\n")
print(ann %>% count(cell_line_curated, cell_line_resolved, immortalised) %>%
        arrange(cell_line_curated), row.names = FALSE)

write.csv(ann, file.path(RERUN_DIR, "immortalisation_annotation_corrected.csv"),
          row.names = FALSE)
cat(sprintf("\nSaved -> %s\n",
            file.path(RERUN_DIR, "immortalisation_annotation_corrected.csv")))
