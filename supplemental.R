################################################### -
## Title: Supplemental analysis with RSV-infected cases
## Author: Ray Pomponio
## Email: pomponiord@upmc.edu
## Project: HMPV Co-detection (NVSN)
## Date Created: 2026-07-06
################################################### -

library(MASS)
library(data.table)
library(mice)
library(flextable)
library(gtsummary)

# ── Data ingest ──────────────────────────────────────────────────────────── -

cdc <- fread("Data/Pitt_Anna_HMPV_JUL26.csv")

# ── Site restriction: 4 CT-reporting sites ──────────────────────────────── -
# Per confirmation: Houston, Pittsburgh, Rochester, Vanderbilt
# (studysite 5, 8, 2, 1)

CT.THRESHOLD <- 30
CT.SITES <- c(1, 2, 5, 8)
dat <- cdc[studysite %in% CT.SITES]
dat[, d_studysite := factor(
  studysite, c(1, 2, 5, 8),
  c("Vanderbilt", "Rochester", "Houston", "Pittsburgh"))]

# row-wise min CT across RSV columns, NA if all missing
dat[, d_rsv_ct := apply(.SD, 1, function(x) {
  x <- suppressWarnings(as.numeric(x))
  if (all(is.na(x))) NA_real_ else min(x, na.rm=TRUE)
}), .SDcols=c("trsvCT", "trsvACT", "trsvBCT")]
dat[, d_hmpv_ct := as.numeric(tmpvCT)]

# restrict to samples less than threshold on either pathogen type
dat <- dat[d_hmpv_ct <= CT.THRESHOLD | d_rsv_ct <= CT.THRESHOLD]

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

# ── Co-detection classification: Analysis A (lab positivity) ───────────── -

RESULT.COLS <- c("tAdeno", "anyflu_result", "piv14_pos", "rhent_pos",
                 "hcov_pos", "c_sarscov2")

# exclude if positive on one or more other pathogens
dat[, d_n_codetect_lab := rowSums(.SD == 1, na.rm=TRUE), .SDcols=RESULT.COLS]
dat <- dat[d_n_codetect_lab < 1]

dat[, d_codetect:=fcase(
  d_hmpv_ct <= CT.THRESHOLD & d_rsv_ct <= CT.THRESHOLD, "HMPV/RSV",
  d_hmpv_ct <= CT.THRESHOLD, "HMPV-only",
  d_rsv_ct <= CT.THRESHOLD, "RSV-only")]
dat[, d_codetect:=factor(d_codetect, c("RSV-only", "HMPV/RSV", "HMPV-only"))]

# ── Cohort assembly ──────────────────────────────────────────────────────── -
# derive outcome: 6-level ordinal illness severity
# c_finalstatus: 1=Inpatient, 2=ED, 3=Outpatient, 5=Urgent Care
# d_hospitalized is on dat before this call
dat[, d_hospitalized := fifelse(
  c_finalstatus == 1, "Hospitalized", "Not Hospitalized")]

dat[, d_severity := fcase(
  c_died      == 1L, 6L,
  c_intubated == 1L, 5L,
  inptEcmo    == 1L, 5L,
  inptICU     == 1L, 4L,
  d_hospitalized == "Hospitalized" & (c_suppoxy == 1L | c_blowby == 1L |
                                        c_hfnc == 1L | c_cpap == 1L), 3L,
  d_hospitalized == "Hospitalized", 2L,
  default=1L)]
dat[, d_severity := factor(
  d_severity, levels=1:6,
  labels=c("Discharge", "Hospitalized", "Hospitalized + O2",
           "ICU", "Ventilated / ECMO", "Died"),
  ordered=TRUE)]

table(dat$d_codetect, exclude=NULL)
table(dat$d_severity, exclude=NULL)

# ── Variable lists ───────────────────────────────────────────────────────── -

DEMO.VARS <- c(
  "d_agemonths", "d_hmpv_ct", "d_sexch", "d_race_eth", "d_scrinsurance",
  "d_premature", "d_anyunderlying", "d_hospitalized", "d_studysite", "d_ariyear")

DEMO.LABELS <- list(
  d_agemonths     = "Age, months",
  d_hmpv_ct       = "HMPV CT value",
  d_sexch         = "Sex",
  d_race_eth      = "Race/ethnicity",
  d_scrinsurance  = "Insurance type",
  d_premature     = "Premature birth (<37 weeks)",
  d_anyunderlying = "Any underlying condition",
  d_hospitalized  = "Hospitalized",
  d_studysite     = "Site",
  d_ariyear       = "Study year")

# ── SUPPLEMENTARY TABLE 1: Demographic and clinical characteristics ──────── -

sup1 <- tbl_summary(
  dat,
  by=d_codetect,
  label=DEMO.LABELS,
  include=all_of(DEMO.VARS),
  type=list(
    d_premature     ~ "dichotomous",
    d_anyunderlying ~ "dichotomous",
    d_hospitalized  ~ "dichotomous"),
  value=list(
    d_premature     ~ "Yes",
    d_anyunderlying ~ "Yes",
    d_hospitalized  ~ "Hospitalized"),
  statistic=list(
    all_continuous() ~ "{median} ({p25}, {p75})")) |>
  add_overall(last=FALSE)

sup1 |> as_flex_table() |> save_as_docx(path="Output/sup_table1.docx")
sup1

# ── MI + ordinal regression ──────────────────────────────────────────────── -
# cases with unknown exposure (d_codetect) or outcome (d_severity) are
# excluded via complete-case filter before imputation; imputation methods
# are also set explicitly to "" for those two variables as a
# belt-and-suspenders safeguard. remaining covariates (d_agemonths,
# d_premature, d_anyunderlying) are imputed using standard mice defaults
# for their variable types.
# NOTE: verify imputation diagnostics (trace plots, density overlay) before
# treating pooled results as final. m=5 is adequate for exploration;
# increase to m>=20 for publication-ready inference.

MODEL.VARS <- c("d_severity", "d_codetect", "d_agemonths",
                "d_premature", "d_anyunderlying", "d_studysite")

IMP.METHODS <- c(
  d_severity      = "",        # outcome:  complete-case only, not imputed
  d_codetect      = "",        # exposure: complete-case only, not imputed
  d_agemonths     = "pmm",
  d_premature     = "logreg",
  d_anyunderlying = "logreg",
  d_studysite     = "")        # fully observed enrollment variable, not
# imputed

run.mi.polr <- function(prelim, m=5, seed=42) {
  mod.dat <- prelim[, .SD, .SDcols=MODEL.VARS]
  mod.dat <- mod.dat[!is.na(d_codetect) & !is.na(d_severity)]
  mids <- mice(mod.dat, m=m, seed=seed, printFlag=FALSE, method=IMP.METHODS)
  with(mids, MASS::polr(
    d_severity ~ d_codetect + d_agemonths + d_premature + d_anyunderlying +
      d_studysite,
    Hess=TRUE))
}

# returns a list of mira objects
fit <- run.mi.polr(dat)

# ── SUPPLEMENTAL TABLE 2 ─────────────────────────────────────────────────── -

sup2 <- tbl_regression(
  fit,
  exponentiate=TRUE,
  include="d_codetect",
  label=list(d_codetect="Co-detected pathogen")) |>
  modify_header(estimate="**OR (95% CI)**") |>
  add_n(location="level") |>
  bold_p(t=0.05)

sup2 |> as_flex_table() |> save_as_docx(path="Output/sup_table2.docx")
sup2