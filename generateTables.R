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
library(gtsummary)

# runs data ingest, cohort assembly, and all d_ derivations;
# produces: dat (pre-CT, HMPV-positive), prelim (analytic cohort), DESIGN, CT.THRESHOLD
source("prelim.R")

# ── Variables needed on dat but only derived on prelim in prelim.R ────────────

# d_hospitalized is constructed from c_finalstatus on prelim inside prelim.R;
# re-derive here on the full pre-CT dat for use in Table 1 overall/excluded rows
dat[, d_hospitalized:=fifelse(c_finalstatus == 1, "Hospitalized", "Not Hospitalized")]

# d_codetect_lab is character on dat; factor for display, retaining all observed levels
dat[, d_codetect_lab:=factor(d_codetect_lab, c("hmpv-only", names(PATHOGENS)))]
dat[, d_codetect_lab:=droplevels(d_codetect_lab)]

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

# ── TABLE 2 setup ─────────────────────────────────────────────────────────────

# build.prelim() encapsulates the cohort assembly logic from prelim.R so that
# all three designs can be built from the same dat object without re-sourcing.
# keep in sync with prelim.R if the assembly or severity derivation logic changes.
build.prelim <- function(dat, design, ct.threshold=CT.THRESHOLD) {
  
  if (design == "A_unrestricted") {
    out <- copy(dat)
    out[, d_codetect:=d_codetect_lab]
    
  } else {
    # step 1: drop if HMPV CT missing or fails threshold (shared between B and C)
    out <- dat[!is.na(d_hmpv_ct) & d_hmpv_ct <= ct.threshold]
    out[, d_partner_ct:=as.numeric(NA)]
    out[, d_partner_inconclusive:=FALSE]
    for (i in which(out[, d_n_codetect_lab == 1])) {
      pth <- as.character(out$d_codetect_lab[i])
      out[i, d_partner_ct        :=get(paste0("d_", pth, "_ct"))]
      out[i, d_partner_inconclusive:=(get(paste0("d_", pth, "_result")) == "Inconclusive")]
    }
    out[, d_partner_fails_ct:=fcase(
      d_n_codetect_lab == 1 & (is.na(d_partner_ct) | d_partner_ct > ct.threshold), TRUE,
      d_n_codetect_lab == 1 & d_partner_inconclusive,                              TRUE,
      default=FALSE)]
    
    if (design == "B_restricted") {
      out <- out[d_partner_fails_ct == FALSE]
      out[, d_codetect:=d_codetect_lab]
      
    } else if (design == "C_reclassify") {
      out[, d_reclassified:=d_partner_fails_ct]
      out[d_partner_fails_ct == TRUE,  d_codetect:="hmpv-only"]
      out[d_partner_fails_ct == FALSE, d_codetect:=d_codetect_lab]
    }
  }
  
  # derive severity; d_hospitalized already on dat (line 21 above), inherited here
  out[, d_severity:=fcase(
    c_died      == 1L, 6L,
    c_intubated == 1L, 5L,
    inptEcmo    == 1L, 5L,
    inptICU     == 1L, 4L,
    d_hospitalized == "Hospitalized" & (c_suppoxy == 1L | c_blowby == 1L |
                                          c_hfnc == 1L | c_cpap == 1L), 3L,
    d_hospitalized == "Hospitalized", 2L,
    default=1L)]
  out[, d_severity:=factor(d_severity, levels=1:6,
                           labels=c("Discharge", "Hospitalized", "Hospitalized + O2",
                                    "ICU", "Ventilated / ECMO", "Died"),
                           ordered=TRUE)]
  
  out[, d_codetect:=factor(d_codetect, c("hmpv-only", names(PATHOGENS)))]
  out[, d_codetect:=droplevels(d_codetect)]
  out
}

# ── Build all three cohorts ───────────────────────────────────────────────────

DESIGNS <- c("A_unrestricted", "B_restricted", "C_reclassify")

prelim.list <- setNames(lapply(DESIGNS, build.prelim, dat=dat), DESIGNS)

invisible(lapply(names(prelim.list), function(nm)
  cat(sprintf("%-20s N=%d\n", nm, nrow(prelim.list[[nm]])))))

# ── MI + ordinal regression ───────────────────────────────────────────────────
# cases with unknown exposure (d_codetect) or outcome (d_severity) are excluded
# via complete-case filter before imputation; imputation methods are also set
# explicitly to "" for those two variables as a belt-and-suspenders safeguard.
# remaining covariates (d_agemonths, d_premature, d_anyunderlying) are imputed
# using standard mice defaults for their variable types.
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
  d_studysite     = "")        # fully observed enrollment variable, not imputed

run.mi.polr <- function(prelim, m=5, seed=42) {
  mod.dat <- prelim[, .SD, .SDcols=MODEL.VARS]
  mod.dat <- mod.dat[!is.na(d_codetect) & !is.na(d_severity)]
  mids <- mice(mod.dat, m=m, seed=seed, printFlag=FALSE, method=IMP.METHODS)
  with(mids, MASS::polr(
    d_severity ~ d_codetect + d_agemonths + d_premature + d_anyunderlying + d_studysite,
    Hess=TRUE))
}

# returns a named list of mira objects (one per design)
fit.list <- lapply(prelim.list, run.mi.polr)

# ── TABLE 2: side-by-side pooled ORs by design ───────────────────────────────
# shows co-detection rows only; full model (including site, age, comorbidities)
# available in fit.list for supplementary reporting

DESIGN.LABELS <- c(
  A_unrestricted = "A: Unrestricted",
  B_restricted   = "B: Restricted (CT \u226430)",
  C_reclassify   = "C: Reclassified (CT \u226430)")

tbl2.list <- lapply(DESIGNS, function(nm) {
  tbl_regression(
    fit.list[[nm]],
    exponentiate=TRUE,
    include="d_codetect",
    label=list(d_codetect="Co-detected pathogen")) |>
    modify_header(estimate="**OR (95% CI)**") |>
    bold_p(t=0.05)
})
names(tbl2.list) <- DESIGNS

tab2 <- tbl_merge(
  tbl2.list,
  tab_spanner=paste0("**", DESIGN.LABELS, "**")) |>
  modify_footnote_header(
    footnote=paste(
      "Proportional odds (ordinal logistic) regression; m=5 multiple imputation,",
      "pooled via Rubin's rules. Reference: HMPV monoinfection.",
      "Adjusted for age (months), preterm birth, any underlying condition,",
      "and enrollment site. OR >1 indicates higher odds of more severe illness."),
    columns=starts_with("estimate"))

tab2

# ── Conclusion-change summary ─────────────────────────────────────────────────
# compares Design A (unrestricted, lab positivity) against the active DESIGN
# (set at top of prelim.R). flags (a) significance flips (p<0.05 in one design
# but not the other) and (b) direction changes (OR crosses 1.0 between designs).
# if DESIGN == "A_unrestricted" the comparison is trivially identical; a warning
# is printed rather than a meaningless table.

get.pooled.codetect <- function(fit, design.name) {
  sm <- as.data.table(summary(mice::pool(fit), conf.int=TRUE, exponentiate=FALSE))
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
  
  warning("DESIGN is 'A_unrestricted'; conclusion-change summary requires a CT-restricted design for comparison.")
  
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
  
  pooled.wide[, flag_sig_flip:=xor(
    get(sig.a) %in% TRUE,
    get(sig.b) %in% TRUE)]
  
  pooled.wide[, flag_dir_change:={
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
  cat(sprintf("\n\u2500\u2500 Conclusion changes: A (unrestricted) vs %s \u2500\u2500\n", label.b))
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