################################################### -
## Title: Preliminary analysis, HMPV co-detection & illness severity
## Author: Ray Pomponio
## Email: pomponiord@upmc.edu
## Project: HMPV Co-detection (NVSN)
## Date Created: 2026-06-30
## Note: Code written with assistance from Claude (Anthropic);
##       all analytical decisions made by the author
################################################### -

library(MASS)
library(data.table)

# ── Data ingest ───────────────────────────────────────────────────────────────

cdc <- fread("Data/Pitt_Anna_HMPV_JUN26.csv")

# ── Design switch ──────────────────────────────────────────────────────────────
# Controls which co-detection/CT handling rule is applied downstream.
#   "A_unrestricted" : all HMPV-positive cases; co-detection by standard lab positivity
#   "B_restricted"   : HMPV CT>30 dropped; co-detection partner CT>30 or inconclusive
#                      results in the CASE BEING DROPPED (not reclassified to mono)
#   "C_reclassify"    : [not yet implemented - placeholder for PI-discussion Analysis C]
#                      HMPV CT<=30 retained; co-detection partner failing CT reclassified
#                      to "absent" rather than dropping the case
DESIGN <- "C_reclassify"
CT.THRESHOLD <- 30

stopifnot(DESIGN %in% c("A_unrestricted", "B_restricted", "C_reclassify"))

# ── Integrity checks ──────────────────────────────────────────────────────────

any(is.na(cdc$Caseid))
range(cdc$scrdate, na.rm=TRUE)
table(cdc$studysite, exclude=NULL)
table(cdc$tmpv, exclude=NULL)
summary(cdc$tmpvCT)

# ── Site restriction: 4 CT-reporting sites ────────────────────────────────────
# Per confirmation: Houston, Pittsburgh, Rochester, Vanderbilt (studysite 5, 8, 2, 1)

CT.SITES <- c(1, 2, 5, 8)
dat <- cdc[studysite %in% CT.SITES]
dat[, d_studysite:=factor(studysite, c(1, 2, 5, 8),
                          c("Vanderbilt", "Rochester", "Houston", "Pittsburgh"))]

# ── Target population: HMPV-positive cases ────────────────────────────────────
# tmpv: 0=Negative, 1=Positive, 2=Inconclusive, 8=Not performed
dat[, d_hmpv_result:=factor(tmpv, c(0, 1, 2, 8),
                            c("Negative", "Positive", "Inconclusive", "Not performed"))]
dat[, d_hmpv_ct:=as.numeric(tmpvCT)]

dat <- dat[d_hmpv_result == "Positive"]

# ── Demographics & covariates (prefix: d_) ────────────────────────────────────

dat[, d_agemonths:=as.numeric(c_agemonths)]
dat[, d_sexch := factor(sexch, levels=c(1, 2), labels=c("Male", "Female"))]
dat[, d_premature:=factor(c_premature, 0:1, c("No", "Yes"))]
dat[, d_anyunderlying:=factor(c_xunderlying, 0:1, c("No", "Yes"))]
dat[, d_scrinsurance := factor(
  scrinsurance,
  levels = c(1, 2, 3, 4),
  labels = c("Public", "Private", "Both", "None/Self-Pay"))]
dat[, d_race_eth := factor(
  c_race_int2,
  levels = 1:7,
  labels = c("White NH", "Black NH", "Hispanic", "NH Pacific Islander",
             "AI/AN NH", "Asian NH", "Multiple/Other NH"))]
dat[, d_ariyear := factor(
  c_ariyear,
  levels = 1:11,
  labels = c("2015-16", "2016-17", "2017-18", "2018-19", "2019-20",
             "2020-21", "2021-22", "2022-23", "2023-24", "2024-25", "2025-26"))]

# ── Pathogen co-detection panel ────────────────────────────────────────────────
# Each pathogen group: lab result variable + CT column(s). For multi-subtype
# pathogens, the CT used downstream is the MINIMUM CT across positive subtypes
# (most conservative - lowest CT corresponds to highest viral load).

PATHOGENS <- list(
  rsv         = list(result="c_rsv_result", ct=c("trsvCT", "trsvACT", "trsvBCT")),
  adenovirus  = list(result="tAdeno",        ct="tAdenoCT"),
  influenza   = list(result="anyflu_result", ct=c("tFluApdmH1CT", "tFluApdmACT", "tFluAH3N2CT",
                                                  "tFluBCT", "tFluBvicCT", "tFluCCT")),
  piv         = list(result="piv14_pos",     ct=c("tpiv1CT", "tpiv2CT", "tpiv3CT", "tpiv4CT")),
  rhino_ent   = list(result="rhent_pos",     ct=c("trhentCT", "trhinoCT", "tenteroCT", "tevd68CT")),
  hcov        = list(result="hcov_pos",      ct=c("tCor229eCT", "tCorhku1CT", "tCorNL63CT", "tCorOC43CT")),
  sarscov2    = list(result="c_sarscov2",    ct="tsarscov2ctrp")
)

# build d_<pathogen>_result (factor) and d_<pathogen>_ct (numeric, min across subtypes)
for (p in names(PATHOGENS)) {
  res.var <- PATHOGENS[[p]]$result
  ct.vars <- PATHOGENS[[p]]$ct
  
  res.col <- paste0("d_", p, "_result")
  ct.col  <- paste0("d_", p, "_ct")
  
  # NOTE: result coding (0/1/2[/8]) is not fully uniform across these variables -
  # verify against the dictionary before trusting labels for sarscov2/hcov/piv/rhent,
  # which only show 0/1/2 (no explicit "not performed" code in the dictionary excerpt)
  dat[, (res.col):=fcase(
    get(res.var) == 1, "Positive",
    get(res.var) == 0, "Negative",
    get(res.var) == 2, "Inconclusive",
    get(res.var) == 8, "Not performed")]
  
  # row-wise min CT across subtype columns, NA if all missing
  dat[, (ct.col):=apply(.SD, 1, function(x) {
    x <- suppressWarnings(as.numeric(x))
    if (all(is.na(x))) NA_real_ else min(x, na.rm=TRUE)
  }), .SDcols=ct.vars]
}

# ── Co-detection classification: Analysis A (lab positivity, unrestricted) ───
# Identifies which (if any) partner pathogen(s) are lab-positive alongside HMPV.
# EXCLUDE cases positive for >1 non-HMPV pathogen, cannot be cleanly assigned
# to a single pairwise reference category. 

RESULT.COLS <- paste0("d_", names(PATHOGENS), "_result")

dat[, d_n_codetect_lab:=rowSums(.SD == "Positive", na.rm=TRUE), .SDcols=RESULT.COLS]

dat <- dat[d_n_codetect_lab <= 1]

dat[, d_codetect_lab:=fcase(
  d_n_codetect_lab == 0, "hmpv-only",
  default=names(PATHOGENS)[apply(.SD, 1, function(x) which(x == "Positive")[1])]),
  .SDcols=RESULT.COLS]
# the line above only resolves correctly for n==1 rows; safe because fcase short-circuits,
# but verify d_codetect_lab is never NA among d_n_codetect_lab==1 rows before proceeding
stopifnot(!any(dat[d_n_codetect_lab == 1, is.na(d_codetect_lab)]))

dat[, d_codetect_lab:=factor(d_codetect_lab, c("hmpv-only", names(PATHOGENS)))]
dat[, d_codetect_lab:=droplevels(d_codetect_lab)]

# ── Cohort assembly by design ──────────────────────────────────────────────────
# d_hospitalized derived here on dat so it is available on both dat and prelim
# (prelim is a copy/subset of dat; generateTables.R uses it on dat for the overall column)
dat[, d_hospitalized:=fifelse(c_finalstatus == 1, "Hospitalized", "Not Hospitalized")]

if (DESIGN == "A_unrestricted") {
  
  prelim <- copy(dat)
  
} else if (DESIGN %in% c("B_restricted", "C_reclassify")) {
  
  # step 1: drop if HMPV CT missing or fails threshold (same for B and C)
  keep <- dat[!is.na(d_hmpv_ct) & d_hmpv_ct <= CT.THRESHOLD]
  
  # step 2: for single co-detections, evaluate partner CT and result
  keep[, d_partner_ct:=as.numeric(NA)]
  keep[, d_partner_inconclusive:=FALSE]
  for (i in which(keep[, d_n_codetect_lab == 1])) {
    pth <- as.character(keep$d_codetect_lab[i])
    keep[i, d_partner_ct        :=get(paste0("d_", pth, "_ct"))]
    keep[i, d_partner_inconclusive:=(get(paste0("d_", pth, "_result")) == "Inconclusive")]
  }
  
  # flag: partner fails CT screen (CT missing, >threshold, or inconclusive result)
  keep[, d_partner_fails_ct:=fcase(
    d_n_codetect_lab == 1 & (is.na(d_partner_ct) | d_partner_ct > CT.THRESHOLD), TRUE,
    d_n_codetect_lab == 1 & d_partner_inconclusive,                              TRUE,
    default=FALSE)]
  
  if (DESIGN == "B_restricted") {
    
    # B: drop the case if partner fails CT screen
    prelim <- keep[d_partner_fails_ct == FALSE]
    prelim[, d_codetect:=d_codetect_lab]
    
  } else if (DESIGN == "C_reclassify") {
    
    # C: retain the case; reclassify co-detection as absent (hmpv-only) if partner
    # fails CT screen. HMPV monoinfection group now includes two sub-types:
    # (a) true lab monoinfections and (b) lab co-detections reclassified after CT screen.
    # d_reclassified flags the latter for sensitivity/audit use.
    prelim <- copy(keep)
    prelim[, d_reclassified:=d_partner_fails_ct]
    prelim[d_partner_fails_ct == TRUE, d_codetect:="hmpv-only"]
    prelim[d_partner_fails_ct == FALSE, d_codetect:=d_codetect_lab]
    
    cat(sprintf(
      "Design C: %d co-detections reclassified to hmpv-only after CT screen\n",
      sum(prelim$d_reclassified)))
  }
}

# categorize co-detection status
prelim[, d_codetect:=factor(d_codetect_lab, c("hmpv-only", names(PATHOGENS)))]
prelim[, d_codetect:=droplevels(d_codetect)]

cat(sprintf("Design: %s | N retained: %d (of %d HMPV-positive at CT sites)\n",
            DESIGN, nrow(prelim), nrow(dat)))
table(prelim$d_codetect, exclude=NULL)

# ── Outcome: 6-level ordinal illness severity ────────────────────────────────
# c_finalstatus: 1=Inpatient, 2=ED, 3=Outpatient, 5=Urgent Care
# d_hospitalized is derived on dat above; available here via prelim as a copy/subset
prelim[, d_severity:=fcase(
  c_died      == 1L, 6L,
  c_intubated == 1L, 5L,
  inptEcmo    == 1L, 5L,
  inptICU     == 1L, 4L,
  d_hospitalized == "Hospitalized" & (c_suppoxy == 1L | c_blowby == 1L |
                                        c_hfnc == 1L | c_cpap == 1L), 3L,
  d_hospitalized == "Hospitalized", 2L,
  default=1L)]

prelim[, d_severity:=factor(d_severity,
                            levels=1:6,
                            labels=c("Discharge", "Hospitalized", "Hospitalized + O2",
                                     "ICU", "Ventilated / ECMO", "Died"),
                            ordered=TRUE)]

table(prelim$d_severity, exclude=NULL)

# ── Missingness check on key analytic variables (flag >10% missing) ──────────

KEY.VARS <- c("d_agemonths", "d_premature", "d_anyunderlying", "d_severity",
              "d_codetect", "d_studysite")
miss.pct <- sapply(KEY.VARS, function(v) mean(is.na(prelim[[v]])) * 100)
miss.flag <- miss.pct[miss.pct > 10]
if (length(miss.flag) > 0) {
  cat("Variables exceeding 10% missing:\n")
  print(round(miss.flag, 1))
} else {
  cat("No key variables exceed 10% missing.\n")
}