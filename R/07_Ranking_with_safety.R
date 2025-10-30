library(data.table); library(stringr); library(scales); library(tidyverse)

# Load your current ranking table
dt <- fread("reports/target_ranked_list.tsv", sep = "\t", na.strings = "")

# ---------------- gnomAD constraint (pLI / LOEUF) ----------------
gnomad_path <- "external/gnomad.v4.0.constraint_metrics.tsv"
if (file.exists(gnomad_path)) {
  gnom <- fread(gnomad_path, sep = "\t", na.strings = "")
  # normalize column names
  setnames(gnom, "ene", "gene", skip_absent = TRUE)  # rename ene to gene
  
  # Select only needed columns and handle duplicates
  gnom <- unique(gnom[, .(
    gene = gene,
    pLI = suppressWarnings(as.numeric(lof.pLI)),
    LOEUF = suppressWarnings(as.numeric(lof.oe))  # using lof.oe as LOEUF equivalent
  )], by = "gene")
  
  dt <- merge(dt, gnom, by = "gene", all.x = TRUE)
  
  # Haploinsufficiency risk proxy (adjusted thresholds for lof.oe):
  dt[, haplo_risk := fifelse(!is.na(LOEUF),
                             pmin(1, pmax(0, (0.65 - LOEUF) / 0.65)),  # adjusted threshold for lof.oe
                             pmin(1, pmax(0, pLI)))]
  
  dt <- dt %>% 
    mutate(haplo_risk=replace_na(0))
} else {
  dt[, `:=`(pLI = NA_real_, LOEUF = NA_real_, haplo_risk = 0)]
}

# ---------------- GTEx v10 tissue expression (on-target toxicity) -------------
gtex_path <- "external/GTEx_Analysis_2022-06-06_v10_RNASeQCv2.4.2_gene_median_tpm.gct"
if (file.exists(gtex_path)) {
  # GCT has 2-line header; read from line 3
  gct <- fread(gtex_path, skip = 2)
  # Name = ENSG.version ; Description = gene symbol
  setnames(gct, c("Name", "Description"), c("ensembl_id", "gene"))
  
  # collapse duplicate symbols (use median TPM across transcripts)
  tissue_cols <- setdiff(names(gct), c("ensembl_id", "gene"))
  gct_agg <- gct[, lapply(.SD, median, na.rm = TRUE), by = gene, .SDcols = tissue_cols]
  
  # helper to pull tissue medians safely (updated patterns for v10 naming)
  get_tpm <- function(dt, pat) {
    cols <- grep(pat, names(dt), ignore.case = TRUE, value = TRUE)
    if (length(cols) == 0) return(rep(NA_real_, nrow(dt)))
    rowMeans(as.matrix(dt[, ..cols]), na.rm = TRUE)
  }
  
  # Critical tissues (adjusted patterns for v10 naming)
  gct_agg[, `:=`(
    tpm_brain  = get_tpm(.SD, "^Brain_|Brain.*Cortex|Brain.*Cerebell"),
    tpm_heart  = get_tpm(.SD, "Heart_"),
    tpm_liver  = get_tpm(.SD, "^Liver($|_)"),
    tpm_kidney = get_tpm(.SD, "Kidney_"),
    tpm_lung   = get_tpm(.SD, "^Lung($|_)"),
    tpm_blood  = get_tpm(.SD, "Whole_Blood"),
    tpm_max    = do.call(pmax, c(.SD, list(na.rm=TRUE)))
  ), .SDcols = tissue_cols]
  
  # Ubiquity calculation
  thr <- 3
  expr_mat <- as.matrix(gct_agg[, ..tissue_cols])
  ubiq <- rowSums(expr_mat >= thr, na.rm = TRUE)
  gct_agg[, `:=`(
    ubiq_count = ubiq,
    ubiq_score = scales::rescale(pmin(ubiq, length(tissue_cols)), to = c(0,1))
  )]
  
  # Calculate tissue-specific TPMs
  gct_agg[, `:=`(
    tpm_brain  = get_tpm(.SD, "^Brain_|Brain.*Cortex|Brain.*Cerebell"),
    tpm_heart  = get_tpm(.SD, "Heart_"),
    tpm_liver  = get_tpm(.SD, "^Liver($|_)"),
    tpm_kidney = get_tpm(.SD, "Kidney_"),
    tpm_lung   = get_tpm(.SD, "^Lung($|_)"),
    tpm_blood  = get_tpm(.SD, "Whole_Blood"),
    tpm_max    = do.call(pmax, c(.SD, list(na.rm=TRUE)))
  ), .SDcols = tissue_cols]
  
  # Calculate high expression flags
  hi <- 1.3
  gct_agg[, `:=`(
    brain_hi  = as.integer(tpm_brain  >= hi),
    heart_hi  = as.integer(tpm_heart  >= hi),
    liver_hi  = as.integer(tpm_liver  >= hi),
    kidney_hi = as.integer(tpm_kidney >= hi),
    lung_hi   = as.integer(tpm_lung   >= hi),
    blood_hi  = as.integer(tpm_blood  >= hi)
  )]
  
  # Calculate critical_tox score
  gct_agg[, critical_tox := as.numeric(pmin(1,
                                            0.30*brain_hi + 0.25*heart_hi + 0.15*liver_hi +
                                              0.10*kidney_hi + 0.10*lung_hi + 0.10*blood_hi
  ))]
  
  # Force evaluation of all columns before merge
  gct_agg <- copy(gct_agg)
  
  # Merge with main data table
  dt <- merge(dt, 
              gct_agg[, .(gene, 
                          ubiq_count, ubiq_score, tpm_max,
                          tpm_brain, tpm_heart, tpm_liver, 
                          tpm_kidney, tpm_lung, tpm_blood,
                          brain_hi, heart_hi, liver_hi, 
                          kidney_hi, lung_hi, blood_hi,
                          critical_tox)],
              by = "gene", all.x = TRUE)
  
} else {
  dt[, `:=`(ubiq_count = NA_integer_, ubiq_score = 0, tpm_max = NA_real_,
            tpm_brain = NA_real_, tpm_heart = NA_real_, tpm_liver = NA_real_,
            tpm_kidney = NA_real_, tpm_lung = NA_real_, tpm_blood = NA_real_,
            brain_hi = 0L, heart_hi = 0L, liver_hi = 0L, kidney_hi = 0L, lung_hi = 0L, blood_hi = 0L,
            critical_tox = 0)]
}

# ---------------- Combine into a safety penalty & adjusted priority ----------
# Weights are example defaults—tune them to your risk posture.
dt[, safety_penalty :=
     0.40*haplo_risk      +   # loss-of-function constraint → higher on-target risk
     0.25*critical_tox    +   # high median expression in critical tissues
     0.15*ubiq_score   ]      # strong human trait associations

# Keep penalty in [0,1]
dt[, safety_penalty := pmin(1, pmax(0, safety_penalty))]

# Adjust your composite priority (e.g., subtract a fraction of penalty)
dt[, priority_safety := priority - 1*safety_penalty]   # 0.30 is the penalty weight—tune it

setorder(dt, -priority_safety)
fwrite(dt, "reports/target_ranked_list_with_safety.tsv", sep = "\t")
fwrite(dt[, .(gene, pLI, LOEUF, haplo_risk, ubiq_count, critical_tox,
              safety_penalty, priority, priority_safety)],
       "reports/safety_metrics.tsv", sep = "\t")
