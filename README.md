# HMPV Co-detection & Illness Severity Analysis

**Project:** New Vaccine Surveillance Network (NVSN)  
**PI / Collaborator:** Anna Wang-Erickson  
**Analyst:** Ray Pomponio (rdp57@pitt.edu / pomponiord@upmc.edu)  
**Department:** Med-Pediatrics, University of Pittsburgh

---

## Timeline

* 3/26/2026: Proposal of one-pager (approved)
* Q3 2026: Clinical data anlysis complete
* Q4 2026: Manuscript drafted
* Q2 2027: Manuscript publication

## Overview
 
This pipeline analyzes HMPV co-detection and ordinal illness severity using multi-site
NVSN surveillance data (2015–2026) from four CT-reporting sites: Houston, Pittsburgh,
Rochester, and Vanderbilt. Three analytical designs are implemented via a single design
switch, allowing side-by-side comparison of results under different CT-threshold handling
assumptions.
 
---
 
## Prerequisites
 
### R version
 
**R ≥ 4.1.0** is required. The pipeline uses the native pipe operator (`|>`) introduced
in R 4.1. Verify your version with:
 
```r
R.version.string
```
 
### Required packages
 
Install all dependencies in a single call before running any script:
 
```r
install.packages(c(
  "data.table",   # all data manipulation throughout the pipeline
  "MASS",         # polr() for proportional odds regression
  "mice",         # multiple imputation (Table 2)
  "gtsummary",    # Table 1 and Table 2 formatting
  "ggplot2"       # Figure 1 CONSORT diagram
))
```
 
#### Known version sensitivities
 
| Package | Minimum tested | Notes |
|---|---|---|
| `gtsummary` | 2.0.0 | `modify_footnote_header()`, `modify_indent()`, and `bold_p(t=)` argument require v2.x |
| `mice` | 3.14.0 | Default method selection for ordered factors (`polr`) and binary factors (`logreg`) relies on v3.x behavior |
| `ggplot2` | 3.4.0 | `linewidth=` argument in `geom_segment()` replaces `size=` in earlier versions |
 
---
 
## Directory structure
 
```
project/
├── Data/
│   └── Pitt_Anna_HMPV_JUN26.csv          ← raw NVSN extract (not under version control)
├── Output/                                ← figures written here (create if absent)
├── prelim.R                               ← data ingest, cohort assembly, all derivations
├── generateTables.R                       ← Table 1 and Table 2
├── drawFigures.R                          ← Figure 1 (CONSORT diagram)
└── README.md
```
 
Create the `Output/` directory if it does not exist:
 
```r
dir.create("Output", showWarnings = FALSE)
```
 
---
 
## Configuration
 
**All analytical design choices are controlled by two constants at the top of `prelim.R`.
These are the only lines that should be changed between runs.**
 
```r
DESIGN       <- "B_restricted"   # see options below
CT.THRESHOLD <- 30               # CT value above which a result fails the threshold
```
 
### DESIGN options
 
| Value | Behavior |
|---|---|
| `"A_unrestricted"` | All HMPV-positive cases retained; co-detection defined by standard lab positivity with no CT restriction |
| `"B_restricted"` | HMPV CT > threshold or missing: case **dropped**. Partner CT > threshold, missing, or inconclusive: case **dropped** |
| `"C_reclassify"` | HMPV CT > threshold or missing: case **dropped**. Partner CT missing or inconclusive: case **dropped**. Partner CT > threshold: case **retained**, co-detection reclassified to HMPV-only (`d_reclassified == TRUE`) |
 
> **Note:** `drawFigures.R` will stop with an error if `DESIGN = "A_unrestricted"` since
> Figure 1 is designed to compare Analysis A against a CT-restricted arm.
 
---
 
## Execution order
 
`prelim.R` **must be consistent with the active `DESIGN`** before any downstream script
is run. Both `generateTables.R` and `drawFigures.R` call `source("prelim.R")` internally,
so they automatically use whichever `DESIGN` is set in that file.
 
### Step 1 — Run `prelim.R` (data ingest and QC)
 
Run this script first, independently, to verify the data loads cleanly and check the
console output before generating outputs.
 
```r
source("prelim.R")
```
 
**Console output to verify:**
 
- No `Caseid` NAs (integrity check)
- Site and HMPV result frequency tables look reasonable
- `"Design: ... | N retained: ... (of ... HMPV-positive at CT sites)"` summary line
- `d_codetect` and `d_severity` frequency tables
- Missingness check: any variable exceeding 10% missing is flagged

### Step 2 — Generate tables (`generateTables.R`)
 
```r
source("generateTables.R")
```
 
This script sources `prelim.R` automatically, then produces:
 
- **Table 1** — demographic and clinical characteristics, stratified by CT inclusion
  status (Included vs. Excluded) under the active design. Under `"A_unrestricted"`,
  a single overall column is shown with no stratification.
- **Table 2** — side-by-side pooled proportional odds regression results for Analysis A
  (unrestricted) vs. the active design.
- **Conclusion-change summary** — printed to console; flags co-detection levels where
  significance (p < 0.05) or OR direction changes between Analysis A and the active design.

Both tables are returned as `gtsummary` objects (`tab1`, `tab2`) rendered in the R viewer.
To save tables to file, add the following after each table is produced:
 
```r
gtsummary::as_gt(tab1) |> gt::gtsave("Output/table1.html")
gtsummary::as_gt(tab2) |> gt::gtsave("Output/table2.html")
```
 
> **Note:** Table 2 runs multiple imputation (`mice`, m = 5) for all three designs
> internally (`prelim.list`), which is computationally the most expensive step. Expect
> several minutes of runtime depending on sample size. Increase `m` to ≥ 20 for
> publication-ready inference before submitting.
 
### Step 3 — Generate figures (`drawFigures.R`)
 
```r
source("drawFigures.R")
```
 
This script sources `prelim.R` automatically, then produces:
 
- **Figure 1** — CONSORT-style inclusion waterfall showing the parallel Analysis A
  (unrestricted) and active design (CT-restricted) arms.
The figure is rendered in the R viewer. To save to file, **uncomment** the two `ggsave`
lines near the end of `drawFigures.R`:
 
```r
ggsave("Output/fig1_consort.pdf", fig1, width=7, height=9, units="in")
ggsave("Output/fig1_consort.png", fig1, width=7, height=9, units="in", dpi=300)
```
 
---
 
## Notes on reproducibility
 
- **Seed:** Multiple imputation in `generateTables.R` uses `seed = 42` in `run.mi.polr()`.
  Results will be exactly reproducible across runs on the same R version and `mice`
  version. Upgrading either may produce minor numerical differences.
- **Data file:** `Data/Pitt_Anna_HMPV_JUN26.csv` is not under version control. Contact
  the CDC data manager for the most recent extract.
- **Working directory:** All scripts use relative paths (`"Data/..."`, `"Output/..."`).
  Set the working directory to the project root before sourcing any script:

```r
  setwd("/path/to/project")
```