################################################### -
## Title: Tables for HMPV co-detection & illness severity analysis
## Author: Ray Pomponio
## Email: pomponiord@upmc.edu
## Project: HMPV Co-detection (NVSN)
## Date Created: 2026-06-30
## Note: Code written with assistance from Claude (Anthropic);
##       all analytical decisions made by the author
################################################### -

library(mice)
library(flextable)
library(gtsummary)

# runs data ingest, cohort assembly, and all d_ derivations;
# produces: dat, prelim, build.prelim(), DESIGN, CT.THRESHOLD, PATHOGENS
# d_hospitalized and d_codetect_lab (factor) are both derived on dat in
# prelim.R. `dat` is already restricted to HMPV-positive, CT <= CT.THRESHOLD
# for ALL designs (A/B/C alike) - see prelim.R for details.
source("prelim.R")

# ── CT inclusion flag ────────────────────────────────────────────────────── -
# TRUE  = case is in prelim (proceeds to downstream analysis under active
#         DESIGN)
# FALSE = case was excluded for a PARTNER-CT-related reason (HMPV CT has
#         already been restricted upstream on `dat`, identically across all
#         designs, so it is not a source of "Excluded (CT)" here)
# Under Design A no cases are excluded; flag is still constructed for
# consistency
dat[, d_ct_included := factor(
  Caseid %in% prelim$Caseid,
  c(TRUE, FALSE),
  c("Included", "Excluded (CT)"))]

N.EXCLUDED <- sum(dat$d_ct_included == "Excluded (CT)")

# ── Variable lists ───────────────────────────────────────────────────────── -

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

# ── Shared helpers ───────────────────────────────────────────────────────── -

fmt.stars <- function(tbl) {
  tbl |>
    bold_p() |>
    modify_table_body(
      ~dplyr::mutate(.x,
                     label=ifelse(
                       row_type=="label" &
                         variable %in% .x$variable[!is.na(.x$p.value) & .x$p.value < 0.001],
                       paste0(label, " ***"),
                       ifelse(
                         row_type=="label" &
                           variable %in% .x$variable[!is.na(.x$p.value) & .x$p.value < 0.01],
                         paste0(label, " **"),
                         ifelse(
                           row_type=="label" &
                             variable %in% .x$variable[
                               !is.na(.x$p.value) & .x$p.value < 0.05],
                           paste0(label, " *"),
                           label)))))
}

fmt <- function(tbl, footnote.excl=NULL) {
  tbl <- tbl |>
    modify_footnote_header(
      footnote=paste('Column percentages shown, exclusive of missing',
                     '("Unknown") values'),
      columns=all_stat_cols()) |>
    modify_indent(columns="label", rows=row_type=="level",   indent=8L) |>
    modify_indent(columns="label", rows=row_type=="missing", indent=8L)
  if (!is.null(footnote.excl)) {
    tbl <- tbl |>
      modify_footnote_header(
        footnote=footnote.excl, columns=starts_with("stat_"))
  }
  tbl
}

# ── TABLE 1: Demographic and clinical characteristics ───────────────────── -
# Cohort for ALL columns (overall and stratified alike) is already restricted
# to HMPV-positive, CT <= CT.THRESHOLD (applied once, upstream, in prelim.R -
# not a per-design step). This differs from the original three-way design:
# Design A is no longer a fully unrestricted comparator, so its caption/
# footnote language below has been updated to describe only the PARTNER
# pathogen's CT handling, not HMPV's own CT.
#   Overall column: HMPV-positive, CT <= CT.THRESHOLD, at 4 CT-reporting
#                    sites, single co-detection only (pre-partner-CT-filter).
#   Stratified columns (Designs B/C only): Included = proceeds to downstream
#                    analysis; Excluded (CT) = dropped for a partner-CT
#                    reason (see caption note below).

TAB1.CAPTION <- sprintf(
  paste("**Table 1.** All cases are HMPV-positive with HMPV CT \u2264 %d;",
        "this restriction is applied identically across every column and",
        "every design (A/B/C)."),
  CT.THRESHOLD)

if (DESIGN == "A_unrestricted") {
  
  # design A applies no further (partner-CT) exclusions beyond the shared
  # HMPV CT <= CT.THRESHOLD restriction; single overall column only
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
    fmt() |>
    modify_caption(TAB1.CAPTION)
  
} else {
  
  # build design-specific footnote describing what "Excluded (CT)" means.
  # NOTE: HMPV CT is intentionally absent from this text - that restriction
  # already applies to every case in `dat`, so it cannot be a reason a case
  # is excluded at this step. Only the PARTNER pathogen's CT differs by
  # design.
  excl.note <- if (DESIGN == "B_restricted") {
    sprintf(
      paste("\"Excluded (CT)\": co-detection partner CT >%d, missing, or",
            "inconclusive. Case dropped entirely (N=%d). (All cases are",
            "already restricted to HMPV CT \u2264%d; see table caption.)"),
      CT.THRESHOLD, N.EXCLUDED, CT.THRESHOLD)
  } else {
    sprintf(
      paste("\"Excluded (CT)\": co-detection partner CT missing or",
            "inconclusive - cases dropped entirely (N=%d). Partner CT >%d",
            "results in reclassification to HMPV monoinfection (retained",
            "in analysis); reclassified cases are not excluded. (All cases",
            "are already restricted to HMPV CT \u2264%d; see table",
            "caption.)"),
      N.EXCLUDED, CT.THRESHOLD, CT.THRESHOLD)
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
      d_race_eth     ~ "fisher.test",
      d_ariyear      ~ "fisher.test",
      d_codetect_lab ~ "fisher.test",
      d_scrinsurance ~ "fisher.test"),
      test.args=list(
        d_race_eth ~ list(simulate.p.value=TRUE, B=10000),
        d_ariyear ~ list(simulate.p.value=TRUE, B=10000),
        d_codetect_lab ~ list(simulate.p.value=TRUE, B=10000),
        d_scrinsurance ~ list(simulate.p.value=TRUE, B=10000))) |>
    modify_footnote_header(
      footnote=paste("Wilcoxon rank-sum for continuous variables; Fisher's",
                     "exact test for dichotomous and small-cell categorical",
                     "variables; chi-square for site."),
      columns="p.value") |>
    fmt.stars() |>
    fmt(footnote.excl=excl.note) |>
    modify_caption(TAB1.CAPTION)
}

tab1 |> as_flex_table() |> save_as_docx(path="Output/table1.docx")
tab1

# ── TABLE 2 setup ────────────────────────────────────────────────────────── -
# build.prelim() is defined in prelim.R (sourced above); calling it here for
# all three designs without re-sourcing or duplicating the assembly logic

# ── Build all three cohorts ──────────────────────────────────────────────── -

DESIGNS <- c("A_unrestricted", "B_restricted", "C_reclassify")

prelim.list <- setNames(lapply(DESIGNS, build.prelim, dat=dat), DESIGNS)

invisible(lapply(names(prelim.list), function(nm)
  cat(sprintf("%-20s N=%d\n", nm, nrow(prelim.list[[nm]])))))

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

# returns a named list of mira objects (one per design)
fit.list <- lapply(prelim.list, run.mi.polr)

# ── TABLE 2: side-by-side pooled ORs by design ───────────────────────────── -
# shows co-detection rows only; full model (including site, age,
# comorbidities) available in fit.list for supplementary reporting

DESIGN.LABELS <- c(
  A_unrestricted = "A: Partner co-detection unrestricted",
  B_restricted   = "B: Partner CT restricted (\u226430)",
  C_reclassify   = "C: Partner CT reclassified (\u226430)")

COMPARE.DESIGNS <- unique(c("A_unrestricted", DESIGN))

make.tbl2 <- function(nm) {
  tbl_regression(
    fit.list[[nm]],
    exponentiate=TRUE,
    include="d_codetect",
    label=list(d_codetect="Co-detected pathogen")) |>
    modify_header(estimate="**OR (95% CI)**") |>
    add_n(location="level") |>
    bold_p(t=0.05)
}

tab2.footnote <- paste(
  "Proportional odds (ordinal logistic) regression; m=5 multiple",
  "imputation, pooled via Rubin's rules. Reference: HMPV monoinfection.",
  "Adjusted for age (months), preterm birth, any underlying condition,",
  "and enrollment site. OR >1 indicates higher odds of more severe",
  "illness. All designs restricted to HMPV CT \u2264", CT.THRESHOLD, ".")

if (length(COMPARE.DESIGNS) == 1) {
  
  tab2 <- make.tbl2("A_unrestricted") |>
    modify_footnote_header(footnote=tab2.footnote, columns="estimate")
  
} else {
  
  tbl2.list <- setNames(lapply(COMPARE.DESIGNS, make.tbl2), COMPARE.DESIGNS)
  
  tab2 <- tbl_merge(
    tbl2.list,
    tab_spanner=paste0("**", DESIGN.LABELS[COMPARE.DESIGNS], "**")) |>
    modify_footnote_header(
      footnote=tab2.footnote, columns=starts_with("estimate"))
}

tab2 |> as_flex_table() |> save_as_docx(path="Output/table2.docx")
tab2

# ── Conclusion-change summary ────────────────────────────────────────────── -
# compares Design A (partner co-detection unrestricted) against the active
# DESIGN (set at top of prelim.R). flags (a) significance flips (p<0.05 in
# one design but not the other) and (b) direction changes (OR crosses 1.0
# between designs). if DESIGN == "A_unrestricted" the comparison is
# trivially identical; a warning is printed rather than a meaningless table.

get.pooled.codetect <- function(fit, design.name) {
  sm <- as.data.table(summary(mice::pool(fit), conf.int=TRUE,
                              exponentiate=FALSE))
  sm <- sm[grepl("^d_codetect", term)]
  sm[, `:=`(
    design = design.name,
    level  = sub("^d_codetect", "", term),
    or     = exp(estimate),
    ci.lo  = exp(estimate - qt(0.975, df) * std.error),
    ci.hi  = exp(estimate + qt(0.975, df) * std.error),
    sig    = p.value < 0.05)]
  sm[, .(design, level, or, ci.lo, ci.hi, p.value, sig)]
}

if (DESIGN == "A_unrestricted") {
  
  warning(paste("DESIGN is 'A_unrestricted'; conclusion-change summary",
                "requires a partner-CT-restricted design for comparison."))
  
} else {
  
  COMPARE.DESIGNS <- c("A_unrestricted", DESIGN)
  
  pooled.all <- rbindlist(
    mapply(get.pooled.codetect,
           fit.list[COMPARE.DESIGNS],
           COMPARE.DESIGNS,
           SIMPLIFY=FALSE))
  
  pooled.wide <- dcast(pooled.all, level ~ design,
                       value.var=c("or", "ci.lo", "ci.hi", "p.value", "sig"))
  
  # column names depend on DESIGN; build references dynamically
  col.a   <- "A_unrestricted"
  col.b   <- DESIGN
  or.a    <- paste0("or_",      col.a);  or.b    <- paste0("or_",      col.b)
  cilo.a  <- paste0("ci.lo_",   col.a);  cilo.b  <- paste0("ci.lo_",   col.b)
  cihi.a  <- paste0("ci.hi_",   col.a);  cihi.b  <- paste0("ci.hi_",   col.b)
  pval.a  <- paste0("p.value_", col.a);  pval.b  <- paste0("p.value_", col.b)
  sig.a   <- paste0("sig_",     col.a);  sig.b   <- paste0("sig_",     col.b)
  
  pooled.wide[, flag_sig_flip := xor(
    get(sig.a) %in% TRUE,
    get(sig.b) %in% TRUE)]
  
  pooled.wide[, flag_dir_change := {
    oa <- get(or.a); ob <- get(or.b)
    !is.na(oa) & !is.na(ob) & ((oa > 1 & ob < 1) | (oa < 1 & ob > 1))
  }]
  
  fmt.cell <- function(or, lo, hi, pv) {
    if (is.na(or)) return(sprintf("%-30s", "(absent from design)"))
    sprintf("%.2f (%.2f-%.2f) [%s]", or, lo, hi,
            ifelse(is.na(pv),  "  NA  ",
                   ifelse(pv < 0.001, "<.001",
                          formatC(pv, digits=3, format="f"))))
  }
  
  label.b <- DESIGN.LABELS[DESIGN]
  cat(sprintf("\n\u2500\u2500 Conclusion changes: A vs %s \u2500\u2500\n",
              label.b))
  cat(sprintf("%-14s  %-30s  %-30s  %s  %s\n",
              "Level", "A: OR (95% CI) [p]",
              paste0(label.b, ": OR (95% CI) [p]"),
              "SigFlip", "DirChg"))
  cat(strrep("-", 100), "\n")
  
  invisible(pooled.wide[, {
    cat(sprintf("%-14s  %-30s  %-30s  %-7s  %s\n",
                level,
                fmt.cell(get(or.a), get(cilo.a), get(cihi.a), get(pval.a)),
                fmt.cell(get(or.b), get(cilo.b), get(cihi.b), get(pval.b)),
                ifelse(flag_sig_flip,   "YES", "\u2014"),
                ifelse(flag_dir_change, "YES", "\u2014")))
  }, by=seq_len(nrow(pooled.wide))])
}