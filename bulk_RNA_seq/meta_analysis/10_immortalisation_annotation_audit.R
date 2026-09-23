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

# ---- SRP089801: CLASSIFIED PRIMARY, against its own GEO characteristic -------
# This study is the one case where the submitters contradict themselves and the
# balance of evidence goes against the immortalisation call. Their GEO record
# carries `cell type: HCA2-hTert` on all 6 samples -- the only mention of hTERT
# anywhere in it -- while the same record says `source_name: Human primary
# fibroblasts`, `tissue: primary fibroblasts`, and titles the series "Expression
# level comparison under dividing and quiescent states in human primary
# fibroblasts". Their paper (Mitra et al. 2017, PNAS 114:E4990,
# doi:10.1073/pnas.1710238114) describes only "Human primary fibroblasts ... in
# Eagle's Minimum Essential Medium (ATCC) supplemented with 15% FBS" -- the same
# growth protocol as the GEO record, confirming it is that paper -- and never
# mentions hTERT or telomerase anywhere.
#
# One characteristics token against three statements of "primary" in the same
# record plus the whole paper. Note the curated metadata also called these samples
# immortalised, so this reverses that call as well, and the strain is therefore
# recorded as the parental HCA2 rather than HCA2-hTERT.
#
# HCA2 is neonatal foreskin regardless: ERP021140's own submission states "human
# foreskin fibroblasts HCA2 were obtained from O. Pereira-Smith", so the tissue
# call does not depend on this decision.
#
# EFFECT: the immortalised-vs-primary gap among proliferating controls weakens
# from +31.1 to +27.8 scaled-difference units (p 1.6e-4 -> 3.3e-3) and stays
# non-significant on YuGene. The paper's conclusion -- that this gap is a
# property of the laboratories and lines that use immortalised cells, not of
# immortalisation -- is unaffected either way.

# ---- GEO-verified immortalisation, with the quote that establishes it --------
IMMORTAL <- tibble::tribble(
  ~study,      ~line,          ~evidence,
  "SRP017378", "BJ-hTERT",     "Sample_growth_protocol: 'Immortalized human BJ primary fibroblast cells (by hTERT expression)'",
  "SRP066917", "BJ-hTERT",     "Sample_characteristics: 'cell type: hTERT-immortalized BJ cells'",
  "SRP066947", "Tig3ET",       "Sample_treatment_protocol: 'Human diploid fibroblast Tig3 cells expressing the ecotropic receptor, hTERT'; Genome Res 27:1634 calls it 'hTERT-immortalized Tig3 fibroblasts (Tig3ET)'",
  "SRP123346", "BJ-hTERT",     "Sample_growth_protocol: 'Immortalized human BJ primary fibroblast cells (by hTERT expression)'",
  "SRP127037", "IMR90-hTERT",  "Sample_characteristics: 'cell line: IMR90 hTERT'",
  "SRP136727", "BJ-hTERT",     "Sample_growth_protocol: 'Immortalized human BJ primary fibroblast cells (by hTERT expression)'",
  "SRP172671", "BJ-hTERT",     "Sample_characteristics: 'cell type: Immortalized BJ-hTERT cells'"
)

# ---- three strain calls settled from the papers' Methods (2026-08-31) --------
#
# ERP021140 = "HCA2", not the composite "HCA2/BJ (primary)". The deposition records no
# strain (no SDRF column; one growth protocol shared by all 156 rows, naming MEFs and
# mouse endothelial cells that are not in the accession). The paper's Sample Preparation
# attributes BJ to exactly one experiment -- replicative senescence, "~PD 65 for BJ
# cells" -- and RS is not in E-MTAB-5403 at all, whose only conditions are Proliferation,
# Quiescence and irradiation day 4/10/20. RS in that paper is entirely public (Alspach
# GSE56293 in BJ = our SRP040745; Marthandan GSE64553/GSE63577; Rai GSE53356 = our
# SRP034541), so BJ had no reason to appear in their own deposition. Six samples per
# condition is replicate structure, not two strains: the keratinocytes and melanocytes
# come from one ATCC lot each (PCS-200-010, PCS-200-012) and also have six per condition.
# NOTE: this merges with SRP089801, which is now also HCA2 -- so HCA2 becomes the one
# primary strain in this pool contributed by two independent studies.
#
# SRP070636 = "IMR90 ER:RAS", not plain "IMR90". Capell et al., Genes Dev (MLL1/SASP),
# Methods: "The cells were generated by retrovirally infecting normal IMR90 fibroblasts
# with pLNCX-ER:Ras, and senescence was induced with 4-OHT". Same inducible system as
# SRP046254 / SRP113324 / SRP113329. This also settles the one sample in the pool that
# looked possibly misclassified: GSM2067916 ("RNA WT Control") is ER:Ras WITHOUT 4-OHT,
# so `Proliferating` is correct, and its `genetic_manipulations = HRASG12V` is the
# uninduced construct rather than an active oncogene.
#
# SRP017142 = "WI-38", unchanged, but its three proliferating controls are confirmed
# transduced and selected, not unmodified: Bischof et al. Methods, "Culturing of human
# diploid WI38 fibroblasts (ATCC) and infection by retroviral-mediated gene transfer with
# either pBABE-puro-HRASG12V, pBABE-puro-HA-PIASY, or empty vector as a control ... at a
# physiological oxygen concentration of 3%". The empty vector is that study's own internal
# control, so it is the right comparator; it is simply not an unmodified cell.

# ---- strain-level cell line for the primary studies -------------------------
PRIMARY_LINE <- c(
  ERP021140 = "HCA2",    SRP017142 = "WI-38",   SRP034163 = "IMR90",
  SRP034541 = "IMR90",   SRP040243 = "IMR90",   SRP040745 = "BJ",
  SRP045867 = "MRC-5",   SRP046254 = "IMR90 ER:RAS", SRP052706 = "foreskin 2DD",
  SRP060598 = "IMR90",   SRP062872 = "IMR90",   SRP064207 = "WI-38",
  SRP065206 = "IMR90",   SRP069768 = "WI-38",   SRP070636 = "IMR90 ER:RAS",
  SRP089801 = "HCA2",    SRP096629 = "HDF161",  SRP098713 = "IMR90",
  SRP113324 = "IMR90 ER:RAS",
  SRP113329 = "IMR90 ER:RAS", SRP117883 = "IMR90", SRP121031 = "IMR90",
  SRP136071 = "IMR90",   SRP153205 = "HDF (per-sample)", SRP153724 = "HDF (per-sample)",
  SRP154382 = "HDF (per-sample)", SRP154577 = "HDF (per-sample)"
)

# ---- the four HDF studies carry the strain PER SAMPLE ------------------------
# These four studies each sequenced several strains of the Kaplon/Coller HDF
# series, and the strain is in every sample's own GEO characteristics -- so a
# study-level label is simply wrong here, not merely coarse. The field name is
# not consistent between them, which is how it was missed the first time:
#
#   SRP153205  strain: 12-1 / 12-3                     (2 strains, 4 samples)
#   SRP153724  cell strain: 12-1 / 12-3                (2 strains, 4 samples)
#   SRP154382  strain: 10-2 / 12-1 / 12-2              (3 strains, 9 samples)
#   SRP154577  fibroblast strain: 10-5 / 12-1          (2 strains, 6 samples)
#
# Collapsing them to one strain per study put three different strains under
# "HDF 10-2" and left "HDF 12-2" out of the annotation entirely. Every strain in
# the series is foreskin, so no tissue call moves, but any per-strain figure does.
HDF_STUDIES <- c("SRP153205", "SRP153724", "SRP154382", "SRP154577")

strain_from_attributes <- function(fa) {
  # matches 'strain;;12-1', 'cell strain;;12-1', 'fibroblast strain;;10-5'
  if (is.na(fa)) return(NA_character_)
  m <- regmatches(fa, regexpr("[a-z ]*strain;;[0-9]{2}-[0-9]", fa, ignore.case = TRUE))
  if (length(m) != 1 || !nzchar(m)) return(NA_character_)
  paste("HDF", sub(".*strain;;", "", m))
}

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
    cell_line_resolved = ifelse(
      study %in% HDF_STUDIES,
      vapply(d$full_attributes, strain_from_attributes, character(1), USE.NAMES = FALSE),
      cell_line_resolved),
    evidence = ifelse(study %in% IMMORTAL$study,
                      IMMORTAL$evidence[match(study, IMMORTAL$study)],
                      "no hTERT/immortal/SV40/telomerase mention in GEO sample or series record"),
    changed = immortalised_curated != immortalised)

hdf_missing <- ann$external_id[ann$study %in% HDF_STUDIES & is.na(ann$cell_line_resolved)]
if (length(hdf_missing))
  stop("no per-sample strain parsed for: ", paste(hdf_missing, collapse = ", "))
cat("\n== HDF series: per-sample strain parsed from GEO characteristics ==\n")
print(table(study = ann$study[ann$study %in% HDF_STUDIES],
            strain = ann$cell_line_resolved[ann$study %in% HDF_STUDIES]))

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


# ---- TISSUE, verified per strain against primary sources (added 2026-08-31) ---
# The `tissue` column in sample_metadata_RERUN.csv had never been checked - script
# 14 verifies the cell line and hTERT status per sample, and tissue was not in its
# scope. Two strains were mislabelled "Skin" when their own sources say foreskin.
# Tissue is one of the four predictors of baseline tAge in meta_analysis/11, so an
# unchecked Skin/Foreskin split made part of that result an arbitrary division.
#
# Every entry below was read from a primary record, named in the comment. Strain
# determines tissue, and strain is verified per sample in script 14, so verifying
# the strain-to-tissue map verifies every sample.
TISSUE_VERIFIED <- c(
  # fetal lung
  "IMR90"             = "Lung",      # ATCC CCL-186, lung, 16-week female fetus
  "IMR90 ER:RAS"      = "Lung",      # IMR90 derivative
  "IMR90-hTERT"       = "Lung",      # IMR90 derivative
  "WI-38"             = "Lung",      # ATCC CCL-75, lung, 3 months gestation
  "MRC-5"             = "Lung",      # ATCC CCL-171, lung, 14 weeks gestation
  "Tig3ET"            = "Lung",      # Cellosaurus CVCL_E939 (TIG-3), fetal lung
  # neonatal foreskin
  "BJ"                = "Foreskin",  # ATCC CRL-2522, foreskin, neonate
  "BJ-hTERT"          = "Foreskin",  # BJ derivative
  "HCA2/BJ (primary)" = "Foreskin",  # ENA PRJEB19157 sample titles: "Primary foreskin fibroblasts"
  "HCA2"              = "Foreskin",  # ERP021140's own submission: "human foreskin
                                     # fibroblasts HCA2 were obtained from
                                     # O. Pereira-Smith". SRP089801 is recorded as
                                     # parental HCA2, not HCA2-hTERT -- see the note
                                     # above the IMMORTAL table.               [was "Skin"]
  # NOT IN USE: no study resolves to HCA2-hTERT any more. SRP089801 was the only
  # candidate and is now recorded as parental HCA2 (see the note above the IMMORTAL
  # table). Kept only so that reversing that decision does not reintroduce an
  # unmapped strain and trip the hard stop below.
  "HCA2-hTERT"        = "Foreskin",  # Cellosaurus CVCL_E2UW, site "Foreskin, skin"
                                     # (retained: no study now resolves to it)
  # The HDF series is foreskin on the submitters' own evidence, not only on
  # Mitra 2018 Genome Biol (PMC6203201): SRP153205 states "tissue: Foreskin" for
  # strains 12-1 and 12-3 outright, and SRP154382/SRP154577 give "cell type:
  # human dermal fibroblasts (foreskin)" for 10-2/12-1/12-2/10-5. SRP153724 says
  # only "dermal fibroblasts", but its samples are the same 12-1 and 12-3 that
  # SRP153205 labels Foreskin -- which is why its 4 samples were the ones
  # mislabelled "Skin".
  "HDF 10-2"          = "Foreskin",  # SRP154382: "human dermal fibroblasts (foreskin)"
  "HDF 10-5"          = "Foreskin",  # SRP154577: idem
  "HDF 12-1"          = "Foreskin",  # SRP153205: "tissue: Foreskin"
  "HDF 12-2"          = "Foreskin",  # SRP154382: "human dermal fibroblasts (foreskin)"
  "HDF 12-3"          = "Foreskin",  # SRP153205: "tissue: Foreskin"                [was "Skin"]
  "foreskin 2DD"      = "Foreskin",  # Trost 2015 R Soc Open Sci: "normal human foreskin fibroblasts"
  # adult dermis - the only adult-derived strain here
  "HDF161"            = "Skin")      # Lammermann 2018 (PMC5895844): "skin biopsies of healthy adult donors", dermis

ann$tissue_verified <- TISSUE_VERIFIED[ann$cell_line_resolved]
miss <- unique(ann$cell_line_resolved[is.na(ann$tissue_verified)])
if (length(miss)) stop("no verified tissue for strain(s): ", paste(miss, collapse = ", "))

meta_t <- read.csv(file.path(RERUN_DIR, "sample_metadata_RERUN.csv"))
tcol <- names(meta_t)[tolower(names(meta_t)) == "tissue"][1]
ann$tissue_metadata <- meta_t[[tcol]][match(ann$external_id, meta_t$external_id)]
chg <- ann %>% filter(tissue_metadata != tissue_verified) %>%
  count(cell_line_resolved, tissue_metadata, tissue_verified)
cat("\n== tissue: metadata vs verified ==\n")
if (nrow(chg)) {
  print(as.data.frame(chg), row.names = FALSE)
  cat(sprintf("  %d samples across %d strains relabelled\n",
              sum(chg$n), nrow(chg)))
} else cat("  no disagreements\n")
cat("  verified tissue counts: ")
print(table(ann$tissue_verified))

write.csv(ann, file.path(RERUN_DIR, "immortalisation_annotation_corrected.csv"),
          row.names = FALSE)
cat(sprintf("\nSaved -> %s\n",
            file.path(RERUN_DIR, "immortalisation_annotation_corrected.csv")))
