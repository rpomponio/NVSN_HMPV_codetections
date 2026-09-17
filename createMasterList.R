################################################### -
## Title: Create a master list
## Author: Ray Pomponio
## Email: pomponiord@upmc.edu
## Project: HMPV Co-detection (NVSN)
## Date Created: 2026-07-08
################################################### -

library(data.table)

# ── Data ingest ──────────────────────────────────────────────────────────── -

cdc <- fread("Data/Pitt_Anna_HMPV_SEP26.csv")
dat <- copy(cdc)

dat[, d_studysite:=factor(studysite, c(8, 1:6),
                        c("Pittsburgh", "Nashville", "Rochester", "Cincinnati",
                          "Seattle", "Houston", "Kansas City"))]

# ── Target population: HMPV-positive cases, CT-restricted ──────────────── -
# tmpv: 0=Negative, 1=Positive, 2=Inconclusive, 8=Not performed
# HMPV CT restriction (REMOVE.HIGH.CT / CT.THRESHOLD) is applied here, on the
# shared `dat` object, so it is now IN EFFECT FOR ALL THREE DESIGNS (A/B/C) -
# not just the CT-restricted designs. This is a change from the original
# three-way comparison, where Design A was meant to be fully unrestricted;
# Design A is now unrestricted only with respect to the PARTNER pathogen's
# CT, not HMPV's own CT. See Table 1 caption and Figure 1 in downstream
# scripts, both of which reflect this shared restriction.
dat[, d_hmpv_result := factor(
  tmpv, c(0, 1, 2, 8),
  c("Negative", "Positive", "Inconclusive", "Not performed"))]
dat[, d_hmpv_ct := as.numeric(tmpvCT)]

dat <- dat[d_hmpv_result == "Positive"]

# ── Demographics & covariates (prefix: d_) ──────────────────────────────── -

dat[, d_agemonths := as.numeric(c_agemonths)]
dat[, d_sexch := factor(sexch, levels=c(1, 2), labels=c("Male", "Female"))]
dat[, d_premature := factor(c_premature, 0:1, c("No", "Yes"))]
dat[, d_anyunderlying := factor(c_xunderlying, 0:1, c("No", "Yes"))]
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
             "2020-21", "2021-22", "2022-23", "2023-24", "2024-25",
             "2025-26"))]

# ── Pathogen co-detection panel ──────────────────────────────────────────── -
# Each pathogen group: lab result variable + CT column(s). For multi-subtype
# pathogens, the CT used downstream is the MINIMUM CT across positive
# subtypes (most conservative - lowest CT corresponds to highest viral
# load).

PATHOGENS <- list(
  rsv = list(
    result = "c_rsv_result",
    ct     = c("trsvCT", "trsvACT", "trsvBCT")),
  adenovirus = list(
    result = "tAdeno",
    ct     = "tAdenoCT"),
  influenza = list(
    result = "anyflu_result",
    ct     = c("tFluACT", "tFluApdmH1CT", "tFluApdmACT", "tFluAH3N2CT",
               "tFluBCT", "tFluBvicCT", "tFluCCT",
               "sFluApdmH1CT", "sFluApdmACT", "sFluAH3N2CT",
               "sFluBvicCT", "sFluCCT")),
  piv = list(
    result = "piv14_pos",
    ct     = c("tpiv1CT", "tpiv2CT", "tpiv3CT", "tpiv4CT")),
  rhino_ent = list(
    result = "rhent_pos",
    ct     = c("trhentCT", "trhinoCT", "tenteroCT", "tevd68CT")),
  hcov = list(
    result = "hcov_pos",
    ct     = c("tCor229eCT", "tCorhku1CT", "tCorNL63CT", "tCorOC43CT")),
  sarscov2 = list(
    result = "c_sarscov2",
    ct     = "tsarscov2ctrp", "tsarscov2p1ct", "tsarscov2p2ct")
)

# build d_<pathogen>_result (factor) and d_<pathogen>_ct (numeric, min
# across subtypes)
for (p in names(PATHOGENS)) {
  res.var <- PATHOGENS[[p]]$result
  ct.vars <- PATHOGENS[[p]]$ct
  
  res.col <- paste0("d_", p, "_result")
  ct.col  <- paste0("d_", p, "_ct")
  
  # NOTE: result coding (0/1/2[/8]) is not fully uniform across these
  # variables - verify against the dictionary before trusting labels for
  # sarscov2/hcov/piv/rhent, which only show 0/1/2 (no explicit "not
  # performed" code in the dictionary excerpt)
  dat[, (res.col) := fcase(
    get(res.var) == 1, "Positive",
    get(res.var) == 0, "Negative",
    get(res.var) == 2, "Inconclusive",
    default="Not performed")]
  
  # row-wise min CT across subtype columns, NA if all missing
  dat[, (ct.col) := apply(.SD, 1, function(x) {
    x <- suppressWarnings(as.numeric(x))
    if (all(is.na(x))) NA_real_ else min(x, na.rm=TRUE)
  }), .SDcols=ct.vars]
  
  # remove inconclusives
  dat <- dat[get(res.col) != "Inconclusive"]
}

# ── Co-detection classification: Analysis A (lab positivity) ───────────── -
# Identifies which (if any) partner pathogen(s) are lab-positive alongside
# HMPV. EXCLUDE cases positive for >1 non-HMPV pathogen, cannot be cleanly
# assigned to a single pairwise reference category.

RESULT.COLS <- paste0("d_", names(PATHOGENS), "_result")

dat[, d_n_codetect_lab := rowSums(.SD == "Positive", na.rm=TRUE),
    .SDcols=RESULT.COLS]

dat[, d_codetect_lab := fcase(
  d_n_codetect_lab == 0, "hmpv-only",
  default=names(PATHOGENS)[
    apply(.SD, 1, function(x) which(x == "Positive")[1])]),
  .SDcols=RESULT.COLS]
dat[, d_codetect_lab := factor(
  d_codetect_lab, c("hmpv-only", names(PATHOGENS)))]
dat[, d_codetect_lab := droplevels(d_codetect_lab)]

# ── Cohort assembly ──────────────────────────────────────────────────────── -
# d_hospitalized derived here on dat so it is available on both dat and
# prelim (prelim is a copy/subset of dat; generateTables.R uses dat for the
# Table 1 overall column, which requires d_hospitalized on the pre-partner-CT
# cohort - note this cohort is already HMPV CT-restricted, see note above)
dat[, d_hospitalized := fifelse(
  c_finalstatus == 1, "Hospitalized", "Not Hospitalized")]

sel.cols <- c("Caseid", "scrdate", "respvisitdt",
              colnames(dat)[grepl("^d_", colnames(dat))])
fwrite(
  dat[, ..sel.cols],
  file.path("Archive/all_mpvPositive_analyticalSample_withCodetectionLabels.csv"))

# ── Line List Agreement ──────────────────────────────────────────────────── -

llist <- readxl::read_xlsx("Data/PITT HMPV specimen line list 9.25.24.xlsx",
                           col_types=c("text", "text", "text", "text", "date",
                                       "date", "numeric"))

# Which samples in the line list are not in the master list?
! llist$Caseid %in% dat$Caseid

# What about the source list?
all(! llist$Caseid %in% cdc$Caseid == ! llist$Caseid %in% dat$Caseid)
