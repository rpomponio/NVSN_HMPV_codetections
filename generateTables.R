################################################### -
## Title: Tables for HMPV co-detection & illness severity analysis
## Author: Ray Pomponio
## Email: pomponiord@upmc.edu
## Project: HMPV Co-detection (NVSN)
## Date Created: 2026-06-30
## Note: Code written with assistance from Claude (Anthropic);
##       all analytical decisions made by the author
################################################### -

library(gtsummary)

# runs data ingest, cohort assembly, and all d_ derivations;
# produces: dat (pre-CT, HMPV-positive), prelim (analytic cohort), DESIGN, CT.THRESHOLD
# d_hospitalized and d_codetect_lab (factor) are derived on dat inside prelim.R
source("prelim.R")

# ── CT inclusion flag ─────────────────────────────────────────────────────────
# TRUE  = case is in prelim (proceeds to downstream analysis under active DESIGN)
# FALSE = case was excluded for any CT-related reason
# Under Design A no cases are excluded; flag is still constructed for consistency
dat[, d_ct_included:=factor(
  Caseid %in% prelim$Caseid,
  c(TRUE, FALSE),
  c("Included", "Excluded (CT)"))]

N.EXCLUDED <- sum(dat$d_ct_included == "Excluded (CT)")

# ── Variable lists ────────────────────────────────────────────────────────────

DEMO.VARS <- c(
  "d_agemonths", "d_hmpv_ct", "d_sexch", "d_race_eth", "d_scrinsurance",
  "d_premature", "d_anyunderlying", "d_hospitalized",
  "d_codetect_lab", "d_studysite", "d_ariyear")

DEMO.LABELS <- list(
  d_agemonths     = "Age, months",
  d_hmpv_ct       = "HMPV CT value",
  d_sexch         = "Sex",
  d_race_eth      = "Race/ethnicity",
  d_scrinsurance  = "Insurance type",
  d_premature     = "Premature birth (<37 weeks)",
  d_anyunderlying = "Any underlying condition",
  d_hospitalized  = "Hospitalized",
  d_codetect_lab  = "Co-detected pathogen (lab)",
  d_studysite     = "Site",
  d_ariyear       = "Study year")

# ── Shared helpers ────────────────────────────────────────────────────────────

fmt.stars <- function(tbl) {
  tbl |>
    bold_p() |>
    modify_table_body(
      ~dplyr::mutate(.x,
                     label=ifelse(
                       row_type=="label" & variable %in% .x$variable[!is.na(.x$p.value) & .x$p.value < 0.001],
                       paste0(label, " ***"),
                       ifelse(
                         row_type=="label" & variable %in% .x$variable[!is.na(.x$p.value) & .x$p.value < 0.01],
                         paste0(label, " **"),
                         ifelse(
                           row_type=="label" & variable %in% .x$variable[!is.na(.x$p.value) & .x$p.value < 0.05],
                           paste0(label, " *"),
                           label)))))
}

fmt <- function(tbl, footnote.excl=NULL) {
  tbl <- tbl |>
    modify_footnote_header(
      footnote='Column percentages shown, exclusive of missing ("Unknown") values',
      columns=all_stat_cols()) |>
    modify_indent(columns="label", rows=row_type=="level",   indent=8L) |>
    modify_indent(columns="label", rows=row_type=="missing", indent=8L)
  if (!is.null(footnote.excl)) {
    tbl <- tbl |>
      modify_footnote_header(footnote=footnote.excl, columns=starts_with("stat_"))
  }
  tbl
}

# ── TABLE 1: Demographic and clinical characteristics ─────────────────────────
# Overall column: pre-CT cohort (all HMPV-positive cases at 4 CT sites, single
# co-detection only). Stratified columns (Designs B/C only): Included = proceeds
# to downstream analysis; Excluded (CT) = dropped for any CT-related reason.

if (DESIGN == "A_unrestricted") {
  
  # design A applies no CT exclusions; single overall column only
  tab1 <- tbl_summary(
    dat,
    label=DEMO.LABELS,
    include=all_of(DEMO.VARS),
    type=list(
      d_premature     ~ "dichotomous",
      d_anyunderlying ~ "dichotomous",
      d_hospitalized  ~ "dichotomous"),
    statistic=list(
      all_continuous() ~ "{median} ({p25}, {p75})")) |>
    fmt()
  
} else {
  
  # build design-specific footnote describing what "Excluded (CT)" means
  excl.note <- if (DESIGN == "B_restricted") {
    sprintf(
      paste("\"Excluded (CT)\": HMPV CT >%d or missing; or co-detection partner",
            "CT >%d, missing, or inconclusive. Case dropped entirely (N=%d)."),
      CT.THRESHOLD, CT.THRESHOLD, N.EXCLUDED)
  } else {
    sprintf(
      paste("\"Excluded (CT)\": HMPV CT >%d or missing. Case dropped entirely.",
            "Partner CT failures under Design C result in reclassification to",
            "HMPV monoinfection, not exclusion (N=%d excluded on HMPV CT only)."),
      CT.THRESHOLD, N.EXCLUDED)
  }
  
  tab1 <- tbl_summary(
    dat,
    by=d_ct_included,
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
    add_overall(last=FALSE) |>
    add_p(test=list(
      all_continuous()  ~ "wilcox.test",
      all_dichotomous() ~ "fisher.test",
      all_categorical() ~ "chisq.test",
      d_race_eth     ~ "fisher.test"),
      test.args=list(
        d_race_eth ~ list(simulate.p.value=TRUE, B=10000))) |>
    modify_footnote_header(
      footnote=paste("Wilcoxon rank-sum for continuous variables; Fisher's exact test",
                     "for dichotomous and small-cell categorical variables; chi-square",
                     "for site, insurance type, study year, and co-detected pathogen."),
      columns="p.value") |>
    fmt.stars() |>
    fmt(footnote.excl=excl.note)
}

tab1