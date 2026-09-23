# 08_fibroblast_series_confound.R
# ---------------------------------------------------------------------------
# The E-MTAB-5403 fibroblast groups are composites of two sets of cultures.
# From the SDRF (rerun_outputs/E-MTAB-5403.sdrf.txt), source names split C/H/I:
#   baseline = 3 C + 3 H | day 4 = 3 H + 3 I | day 10 = 3 H + 3 I | day 20 = 6 H
# Keratinocytes (K1-K30) and melanocytes (M1-M30) are single series - unaffected.
# ENA runs are one contiguous alphabetically-sorted block (ERR1805188-217) and
# sra.run_published differs by 1 s, so neither indicates a batch; the evidence
# below is expression- and clock-level.
#
# HOUSE RULE: a contrast whose rank-test floor exceeds 0.05 is NOT TESTABLE and
# gets no p-value (cf. BJ 7 v 2 in immortalisation_contrasts.csv). 3 v 3 -> 0.100.
# Series is balanced 3+3 across day 4 and day 10, so centring within timepoint
# gives a legitimate 6 v 6 test of series (floor 0.0022).
# ---------------------------------------------------------------------------
suppressMessages(library(SummarizedExperiment))
# run from bulk_RNA_seq/, like every other script: Rscript temporal_analysis/09_fibroblast_series_confound.R
source('R/config.R')
B <- file.path(PROJECT_DIR, 'bulk_RNA_seq')
fl <- function(a,b) 2/choose(a+b,min(a,b))
wt <- function(x,y) suppressWarnings(wilcox.test(x,y)$p.value)

## ---- clock-level -----------------------------------------------------------
d <- read.csv(file.path(B,'rerun_outputs','tage_temporal_by_celltype.csv'))
m <- read.csv(file.path(B,'rerun_outputs','mortality_temporal_by_celltype.csv'))
d$mort <- m$mortality_tAge[match(d$external_id,m$external_id)]
d$sn <- sub("^E-MTAB-5403:","",d$label); d$ser <- sub("[0-9].*$","",d$sn)
fb <- d[d$cell_type=="Fibroblast",]
CL <- c(yugene_diff_EN_tAge="yugene", scaled_diff_EN_tAge="scaled", mort="mortality")

cat("=== 1. series composition by timepoint (design fact) ===\n"); print(table(fb$ser, fb$time_after_treatment))

cat("\n=== 2. series effect on tAge, timepoint-centred (6 v 6, floor",
    sprintf("%.4f", fl(6,6)), ") ===\n")
s <- fb[fb$time_after_treatment %in% c("4_days","10_days"),]
cat(sprintf("%-10s %9s %9s %9s %9s\n","clock","H","I","diff","p"))
for (v in names(CL)) { r <- s[[v]] - ave(s[[v]], s$time_after_treatment)
  cat(sprintf("%-10s %+9.2f %+9.2f %+9.2f %9.4f\n", CL[[v]], mean(r[s$ser=="H"]),
      mean(r[s$ser=="I"]), mean(r[s$ser=="I"])-mean(r[s$ser=="H"]), wt(r[s$ser=="H"],r[s$ser=="I"]))) }

cat("\n=== 3. composition-matched trajectory, series H only ===\n")
H <- fb[fb$ser=="H",]; base <- H[H$time_after_treatment=="none",]
cat(sprintf("%-22s %-10s %9s %9s %13s\n","contrast","clock","diff","floor","p"))
row <- function(x,y,lab,v){ f <- fl(length(x),length(y))
  cat(sprintf("%-22s %-10s %+9.2f %9.3f %13s\n", lab, CL[[v]], mean(x)-mean(y), f,
      if (f>0.05) "not testable" else sprintf("%.3f", wt(x,y)))) }
for (v in names(CL)) { for (tp in c("4_days","10_days","20_days"))
    row(H[[v]][H$time_after_treatment==tp], base[[v]], paste0(tp," vs baseline"), v)
  row(H[[v]][H$time_after_treatment=="20_days"], H[[v]][H$time_after_treatment=="4_days"], "20d vs 4d", v)
  cat("\n") }

## ---- expression-level ------------------------------------------------------
se <- readRDS(file.path(B,'rerun_outputs','cs_cq_download_raw.rds'))
cd <- as.data.frame(colData(se)); k <- cd$study=="ERP021140"; se <- se[,k]; cd <- cd[k,]
sn <- sub("^E-MTAB-5403:","",cd$sra.sample_title); ser <- sub("[0-9].*$","",sn)
ph <- readRDS(file.path(B,'rerun_outputs','temporal_recount_pheno.rds')); ph$sn <- sub("^E-MTAB-5403:","",ph$label)
ct <- setNames(ph$cell_type,ph$sn)[sn]; tp <- setNames(as.character(ph$time_after_treatment),ph$sn)[sn]
tp[is.na(tp) & !is.na(ct)] <- "baseline"
pca <- function(idx, lab) {
  M <- as.matrix(assay(se,"raw_counts")[, idx, drop=FALSE]); storage.mode(M)<-"double"
  M <- M[rowSums(M>0) >= length(idx)/2, , drop=FALSE]
  z <- log2(sweep(M,2,colSums(M)/1e6,"/")+1); v <- apply(z,1,var); z <- z[v>1e-8,,drop=FALSE]
  v <- v[v>1e-8]; X <- z[order(-v)[1:min(2000,nrow(z))],,drop=FALSE]
  p <- prcomp(t(X), scale.=TRUE); ve <- round(100*p$sdev^2/sum(p$sdev^2),1)
  cat(sprintf("\n=== %s | n=%d | %d genes | PC1-3: %s ===\n", lab, length(idx), nrow(X),
      paste(paste0(ve[1:3],"%"), collapse=" ")))
  data.frame(name=sn[idx], ser=ser[idx], tp=tp[idx], PC1=p$x[,1], PC2=p$x[,2], ve1=ve[1], ve2=ve[2])
}
df <- pca(which(ct=="Fibroblast"), "4. PCA, all fibroblast time-course samples")
ag <- aggregate(cbind(PC1,PC2)~tp+ser, df, mean); ag$n <- aggregate(PC1~tp+ser,df,length)$PC1
print(ag[order(ag$tp,ag$ser),], row.names=FALSE, digits=4)
b <- df[df$tp=="baseline",]
cat(sprintf("\n  baseline C vs H on PC1: %+.1f vs %+.1f (gap %.1f; within-group sd %.1f/%.1f)",
    mean(b$PC1[b$ser=="C"]), mean(b$PC1[b$ser=="H"]),
    abs(diff(tapply(b$PC1,b$ser,mean))), sd(b$PC1[b$ser=="C"]), sd(b$PC1[b$ser=="H"])))
cat("  -> 3 v 3, not testable; the two baseline sets are indistinguishable.\n")

d2 <- pca(which(ct=="Fibroblast" & tp %in% c("4 days","10 days")),
          "5. PCA, fibroblast day4+day10 only (6 H vs 6 I)")
for (i in 1:2) { pc <- paste0("PC",i)
  r <- d2[[pc]] - ave(d2[[pc]], d2$tp); ps <- wt(r[d2$ser=="H"], r[d2$ser=="I"])
  q <- d2[[pc]] - ave(d2[[pc]], d2$ser); pt <- wt(q[d2$tp=="4 days"], q[d2$tp=="10 days"])
  cat(sprintf("  %s (%4.1f%%): series H %+7.2f vs I %+7.2f p=%.4f | timepoint p=%.4f\n",
      pc, if(i==1) d2$ve1[1] else d2$ve2[1], mean(r[d2$ser=="H"]), mean(r[d2$ser=="I"]), ps, pt)) }
cat("\nBoth floors are 0.0022, so p=0.0022 denotes complete separation.\n")
