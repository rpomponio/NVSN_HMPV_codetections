################################################### -
## Title: Figures — HMPV co-detection & illness severity
## Author: Ray Pomponio
## Email: pomponiord@upmc.edu
## Project: HMPV Co-detection (NVSN)
## Date Created: 2026-06-30
## Note: Code written with assistance from Claude (Anthropic);
##       all analytical decisions made by the author
################################################### -

library(mice)
library(ggplot2)

# produces: cdc, dat, prelim, build.prelim(), DESIGN, CT.THRESHOLD, CT.SITES, PATHOGENS
source("prelim.R")

if (DESIGN == "A_unrestricted")
  stop("Figure 1 requires a CT-restricted DESIGN (B or C); set DESIGN in prelim.R")

# ── Build design-specific cohort for figure ────────────────────────────────────

prelim.design <- build.prelim(dat, DESIGN)

DESIGN.LABELS <- c(
  B_restricted = sprintf("Design B: Restricted (CT \u2264%d)", CT.THRESHOLD),
  C_reclassify = sprintf("Design C: Reclassified (CT \u2264%d)", CT.THRESHOLD))

# ── FIGURE 1: CONSORT-style inclusion waterfall ─────────────────────────────────────────────────────
# two parallel arms derived from the pre-CT cohort:
#   left  — Analysis A (unrestricted; no further exclusions)
#   right — active DESIGN (CT-restricted; exclusions or reclassification applied)

# ── N computations ─────────────────────────────────────────────────────────────

N.CT.SITES      <- cdc[studysite %in% CT.SITES, .N]
N.HMPV.POS      <- cdc[studysite %in% CT.SITES & tmpv == 1, .N]
N.EXCL.NOT.HMPV <- N.CT.SITES - N.HMPV.POS
N.EXCL.MULTI    <- N.HMPV.POS - nrow(dat)
N.DAT           <- nrow(dat)
N.AFTER.HMPV.CT <- dat[!is.na(d_hmpv_ct) & d_hmpv_ct <= CT.THRESHOLD, .N]
N.EXCL.HMPV.CT  <- N.DAT - N.AFTER.HMPV.CT
N.DESIGN        <- nrow(prelim.design)

fmt.n <- function(n) formatC(n, format="d", big.mark=",")

# N.DROPPED.PARTNER: cases lost at the partner CT step in both designs
#   Design B: all partner CT failures (missing, inconclusive, >threshold)
#   Design C: only missing/inconclusive (CT >threshold is reclassified, not dropped)
N.DROPPED.PARTNER <- N.AFTER.HMPV.CT - N.DESIGN

if (DESIGN == "C_reclassify")
  N.RECLASSIFIED <- sum(prelim.design$d_reclassified)

# ── Box layout ─────────────────────────────────────────────────────────────────
# coordinate system: x=[0, 11.5], y=[2.0, 11.2] (y=11.2 at top)
# main shared flow:  xc=4.2
# exclusion column:  xc=8.5 (main flow) / xc=9.8 (CT restriction arm)
# Analysis A arm:    xc=1.8
# CT restriction arm: xc=6.8

W.M <- 4.2   # width: main flow boxes
W.E <- 3.2   # width: exclusion boxes
W.A <- 3.0   # width: fork arm boxes
H   <- 0.78  # height: all main/arm boxes
HE  <- 0.72  # height: exclusion boxes

# fill colors
CLR.MAIN <- "#dce8f5"   # blue-grey: shared enrollment flow
CLR.EXCL <- "#f5e6e6"   # red-tinted: exclusion boxes
CLR.A    <- "#e8f5e9"   # green-tinted: Analysis A
CLR.D    <- "#fff3e0"   # orange-tinted: Design B/C

mb <- function(id, xc, yc, w, h, fill, label) {
  data.table(id=id, xc=xc, yc=yc,
             xmin=xc - w/2, xmax=xc + w/2,
             ymin=yc - h/2, ymax=yc + h/2,
             fill=fill, label=label)
}

shared.boxes <- list(
  # shared enrollment funnel (3 boxes, top to bottom)
  mb("b_sites",   4.2, 10.5, W.M, H,  CLR.MAIN,
     sprintf("Enrolled: 4 CT-reporting sites\nN = %s", fmt.n(N.CT.SITES))),
  mb("b_hmpv",    4.2,  8.7, W.M, H,  CLR.MAIN,
     sprintf("HMPV-positive\nN = %s", fmt.n(N.HMPV.POS))),
  mb("b_dat",     4.2,  6.9, W.M, H,  CLR.MAIN,
     sprintf("Pre-CT eligible cohort\nN = %s", fmt.n(N.DAT))),
  
  # exclusion boxes: right of main flow
  mb("e_hmpv",    8.5,  9.6, W.E, HE, CLR.EXCL,
     sprintf("Excluded: not HMPV-positive\nN = %s", fmt.n(N.EXCL.NOT.HMPV))),
  mb("e_multi",   8.5,  7.8, W.E, HE, CLR.EXCL,
     sprintf("Excluded: multiple/inconclusive co-detections\nN = %s", fmt.n(N.EXCL.MULTI))),
  
  # left fork arm: Design A (no further restriction)
  mb("b_arm_a",   1.8,  5.2, W.A, H,  CLR.A,
     sprintf("Design A (Unrestricted)\nN = %s", fmt.n(N.DAT))),
  
  # right fork arm: CT restriction (shared through HMPV CT step)
  mb("b_hmpv_ct", 6.8,  5.2, W.A, H,  CLR.MAIN,
     sprintf("HMPV CT \u2264 %d\nN = %s", CT.THRESHOLD, fmt.n(N.AFTER.HMPV.CT))),
  mb("e_hmpvct",  9.8,  5.9, W.E, HE, CLR.EXCL,
     sprintf("Excluded: HMPV CT >%d or missing\nN = %s", CT.THRESHOLD, fmt.n(N.EXCL.HMPV.CT)))
)

# ── Partner CT step (design-specific) ─────────────────────────────────────────
# Design B: single red exclusion box — all partner CT failures dropped
# Design C: two boxes — red for dropped (CT missing/inconclusive) and
#           yellow for reclassified (CT >threshold; cases are retained)
#           b_design shifts down to accommodate the second box

if (DESIGN == "B_restricted") {
  
  partner.boxes <- list(
    mb("e_partner", 9.8, 3.9, W.E, HE, CLR.EXCL,
       sprintf("Excluded: partner CT >%d or missing\nN = %s",
               CT.THRESHOLD, fmt.n(N.DROPPED.PARTNER))),
    mb("b_design",  6.8, 3.2, W.A, H, CLR.D,
       sprintf("%s\nN = %s", DESIGN.LABELS[DESIGN], fmt.n(N.DESIGN))))
  
  y.junc.partner <- 3.90
  yt.design      <- 3.2 + H / 2
  ylim.bottom    <- 2.2
  
} else {
  
  partner.boxes <- list(
    mb("e_partner_drop",      9.8, 4.3, W.E, HE, CLR.EXCL,
       sprintf("Dropped: partner CT missing\nN = %s",
               fmt.n(N.DROPPED.PARTNER))),
    mb("e_partner_reclassify",9.8, 3.4, W.E, HE, "#fef9e7",
       sprintf("Reclassified: partner CT >%d\nN = %s",
               CT.THRESHOLD, fmt.n(N.RECLASSIFIED))),
    mb("b_design",  6.8, 2.6, W.A, H, CLR.D,
       sprintf("%s\nN = %s", DESIGN.LABELS[DESIGN], fmt.n(N.DESIGN))))
  
  y.junc.drop        <- 4.30
  y.junc.reclassify  <- 3.40
  yt.design          <- 2.6 + H / 2
  ylim.bottom        <- 1.85
  
}

boxes <- rbindlist(c(shared.boxes, partner.boxes))

ARR <- arrow(length=unit(0.18, "cm"), type="closed")

# box edge y-values used below
yb.sites  <- 10.5 - H/2   # = 10.11
yt.hmpv   <-  8.7 + H/2   # =  9.09
yb.hmpv   <-  8.7 - H/2   # =  8.31
yt.dat    <-  6.9 + H/2   # =  7.29
yb.dat    <-  6.9 - H/2   # =  6.51
yt.arm_a  <-  5.2 + H/2   # =  5.59
yt.hmpvct <-  5.2 + H/2   # =  5.59
yb.hmpvct <-  5.2 - H/2   # =  4.81

# junction y-values: midpoints used for exclusion branches
y.junc.e1     <- (yb.sites + yt.hmpv) / 2   # = 9.60
y.junc.e2     <- (yb.hmpv  + yt.dat)  / 2   # = 7.80
y.fork        <- 6.20
y.junc.hmpvct <- 5.90

# left edge of right-arm exclusion boxes (xc=9.8, w=3.2 → left=8.2)
x.excl.right.left <- 9.8 - W.E / 2   # = 8.2

# shared arrowed segments (everything above and including the HMPV CT step)
segs.arr.shared <- rbind(
  data.table(x=4.2, y=yb.sites,      xend=4.2,               yend=yt.hmpv,      col="black"),  # B1→B2
  data.table(x=4.2, y=y.junc.e1,     xend=6.9,               yend=y.junc.e1,    col="grey40"), # →e_hmpv
  data.table(x=4.2, y=yb.hmpv,       xend=4.2,               yend=yt.dat,       col="black"),  # B2→B3
  data.table(x=4.2, y=y.junc.e2,     xend=6.9,               yend=y.junc.e2,    col="grey40"), # →e_multi
  data.table(x=1.8, y=y.fork,        xend=1.8,               yend=yt.arm_a,     col="black"),  # →arm_a
  data.table(x=6.8, y=y.junc.hmpvct, xend=6.8,               yend=yt.hmpvct,    col="black"),  # →hmpv_ct
  data.table(x=6.8, y=y.junc.hmpvct, xend=x.excl.right.left, yend=y.junc.hmpvct,col="grey40") # →e_hmpvct
)

# shared plain segments (no arrowhead)
segs.line.shared <- rbind(
  data.table(x=4.2, y=yb.dat,        xend=4.2, yend=y.fork),        # b_dat→fork
  data.table(x=1.8, y=y.fork,        xend=6.8, yend=y.fork),        # horizontal fork
  data.table(x=6.8, y=y.fork,        xend=6.8, yend=y.junc.hmpvct)  # fork→hmpv_ct junction
)

# ── Partner CT segments (design-specific) ─────────────────────────────────────
if (DESIGN == "B_restricted") {
  
  # single junction: both drop reasons handled together
  partner.segs.arr <- rbind(
    data.table(x=6.8, y=y.junc.partner, xend=6.8,               yend=yt.design,      col="black"),
    data.table(x=6.8, y=y.junc.partner, xend=x.excl.right.left, yend=y.junc.partner, col="grey40"))
  
  partner.segs.line <- data.table(x=6.8, y=yb.hmpvct, xend=6.8, yend=y.junc.partner)
  
} else {
  
  # two junctions: drop junction (CT missing/inconclusive) above reclassify junction (CT >threshold)
  partner.segs.arr <- rbind(
    data.table(x=6.8, y=y.junc.drop,       xend=x.excl.right.left, yend=y.junc.drop,       col="grey40"),
    data.table(x=6.8, y=y.junc.reclassify, xend=x.excl.right.left, yend=y.junc.reclassify, col="grey40"),
    data.table(x=6.8, y=y.junc.reclassify, xend=6.8,               yend=yt.design,         col="black"))
  
  # single line from hmpv_ct bottom through both junctions; branch arrows handle the rest
  partner.segs.line <- data.table(x=6.8, y=yb.hmpvct, xend=6.8, yend=y.junc.reclassify)
  
}

segs.arr  <- rbind(segs.arr.shared,  partner.segs.arr)
segs.line <- rbind(segs.line.shared, partner.segs.line)

# ── Build figure ───────────────────────────────────────────────────────────────

fig1 <- ggplot() +
  # draw segments before boxes so boxes sit on top
  geom_segment(data=segs.line,
               aes(x=x, y=y, xend=xend, yend=yend),
               linewidth=0.45, color="black") +
  geom_segment(data=segs.arr,
               aes(x=x, y=y, xend=xend, yend=yend, color=col),
               arrow=ARR, linewidth=0.45) +
  scale_color_identity() +
  # boxes on top of segments
  geom_rect(data=boxes,
            aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, fill=fill),
            color="black", linewidth=0.35) +
  scale_fill_identity() +
  geom_text(data=boxes,
            aes(x=xc, y=yc, label=label),
            size=2.75, lineheight=1.25) +
  coord_cartesian(xlim=c(0, 11.5), ylim=c(ylim.bottom, 11.2), expand=FALSE) +
  theme_void() +
  theme(
    legend.position  = "none",
    plot.background  = element_rect(fill="white", color=NA),
    plot.margin      = margin(8, 8, 8, 8))

ggsave("Output/fig1_consort.pdf", fig1, width=7, height=9, units="in")
ggsave("Output/fig1_consort.png", fig1, width=7, height=9, units="in", dpi=300)

fig1
# ── FIGURE 2: OR comparison, Analysis A vs active DESIGN ──────────────────────
# each co-detection group shown as an ellipse:
#   center: (OR_A, OR_design) in log scale
#   x semi-axis: half the Analysis A CI width in log space
#   y semi-axis: half the active design CI width in log space
# the 45-degree line indicates perfect agreement between designs.
# ellipses colored by significance pattern (p<0.05) across the two designs.
#
# note: if generateTables.R was sourced earlier in the same session, the fitted
# models (fit.list) and helper functions are reused without re-running MI.
# running this script cold (without a prior generateTables.R source) will
# trigger MI which takes several minutes.



# ── Reuse or rebuild model infrastructure ─────────────────────────────────────

FIG2.DESIGNS <- unique(c("A_unrestricted", DESIGN))

if (!exists("prelim.list"))
  prelim.list <- setNames(lapply(FIG2.DESIGNS, build.prelim, dat=dat), FIG2.DESIGNS)

if (!exists("run.mi.polr")) {
  # mirrors generateTables.R; reproduced here so drawFigures.R is self-contained
  MODEL.VARS  <- c("d_severity", "d_codetect", "d_agemonths",
                   "d_premature", "d_anyunderlying", "d_studysite")
  IMP.METHODS <- c(d_severity="", d_codetect="",
                   d_agemonths="pmm", d_premature="logreg",
                   d_anyunderlying="logreg", d_studysite="")
  run.mi.polr <- function(prelim, m=5, seed=42) {
    mod.dat <- prelim[, .SD, .SDcols=MODEL.VARS]
    mod.dat <- mod.dat[!is.na(d_codetect) & !is.na(d_severity)]
    mids <- mice(mod.dat, m=m, seed=seed, printFlag=FALSE, method=IMP.METHODS)
    with(mids, MASS::polr(
      d_severity ~ d_codetect + d_agemonths + d_premature + d_anyunderlying + d_studysite,
      Hess=TRUE))
  }
}

if (!exists("fit.list")) {
  message("fit.list not found; running MI for Figure 2 (this may take several minutes)...")
  fit.list <- lapply(prelim.list[FIG2.DESIGNS], run.mi.polr)
}

if (!exists("get.pooled.codetect")) {
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
    sm[, .(design, level, or, ci.lo, ci.hi, sig)]
  }
}

# ── Extract and reshape pooled estimates ──────────────────────────────────────

pooled.fig2 <- rbindlist(mapply(
  get.pooled.codetect, fit.list[FIG2.DESIGNS], FIG2.DESIGNS, SIMPLIFY=FALSE))

fig2.wide <- dcast(pooled.fig2, level ~ design, value.var=c("or", "ci.lo", "ci.hi", "sig"))

# rename to design-agnostic short names for figure code below
setnames(fig2.wide,
         old=c(paste0(c("or_", "ci.lo_", "ci.hi_", "sig_"), "A_unrestricted"),
               paste0(c("or_", "ci.lo_", "ci.hi_", "sig_"), DESIGN)),
         new=c("or_A", "ci.lo_A", "ci.hi_A", "sig_A",
               "or_D", "ci.lo_D", "ci.hi_D", "sig_D"))

# significance pattern: compare p<0.05 across both designs
# NA in or_D means the level is absent from the active design; these are
# excluded from the figure since no ellipse can be drawn.
# groups with non-finite or non-positive CI bounds (e.g. near-separation from
# sparse data) are also excluded — they would produce log(-Inf) in the scale.
fig2.wide <- fig2.wide[
  !is.na(or_A) & !is.na(or_D) &
    is.finite(ci.lo_A) & is.finite(ci.hi_A) & ci.lo_A > 0 &
    is.finite(ci.lo_D) & is.finite(ci.hi_D) & ci.lo_D > 0]

if (nrow(fig2.wide) == 0)
  stop("No co-detection groups have finite CI bounds in both designs; cannot draw Figure 2.")

if (nrow(fig2.wide) < nrow(dcast(pooled.fig2, level ~ design)))
  message(sprintf("Figure 2: %d group(s) excluded due to non-finite CI bounds (likely sparse data).",
                  nrow(dcast(pooled.fig2, level ~ design)) - nrow(fig2.wide)))

fig2.wide[, sig_cat:=fcase(
  sig_A &  sig_D,  "Significant in both",
  !sig_A & !sig_D,  "Not significant in either",
  default="Significance differs")]
fig2.wide[, sig_cat:=factor(sig_cat,
                            c("Significant in both", "Significance differs", "Not significant in either"))]

# ── Ellipse polygons ──────────────────────────────────────────────────────────
# generated in log space so the semi-axes correspond to the Wald CI half-widths.
# back-transformed to OR scale for ggplot2 + scale_*_log10().

THETA <- seq(0, 2 * pi, length.out=120)

ellipses <- rbindlist(lapply(seq_len(nrow(fig2.wide)), function(i) {
  r <- fig2.wide[i]
  # semi-axes in log space: half CI width
  cx <- log(r$or_A);   a <- (log(r$ci.hi_A) - log(r$ci.lo_A)) / 2
  cy <- log(r$or_D);   b <- (log(r$ci.hi_D) - log(r$ci.lo_D)) / 2
  data.table(
    level   = r$level,
    sig_cat = r$sig_cat,
    x       = exp(cx + a * cos(THETA)),
    y       = exp(cy + b * sin(THETA)))
}))

# ── Axis range ────────────────────────────────────────────────────────────────
# symmetric across both axes so the agreement line bisects at 45 degrees.
# padded by 15% in log space to give labels and annotation room.

all.bounds <- c(fig2.wide$ci.lo_A, fig2.wide$ci.hi_A,
                fig2.wide$ci.lo_D, fig2.wide$ci.hi_D)
pad      <- exp(diff(log(range(all.bounds, na.rm=TRUE))) * 0.15)
ax.range <- c(min(all.bounds, na.rm=TRUE) / pad,
              max(all.bounds, na.rm=TRUE) * pad)

# axis breaks: a tidy subset of the standard log-OR ladder within ax.range
LOG.BREAKS <- c(0.125, 0.25, 0.5, 1, 2, 4, 8)
ax.breaks  <- LOG.BREAKS[LOG.BREAKS >= ax.range[1] & LOG.BREAKS <= ax.range[2]]
ax.labels  <- ifelse(ax.breaks < 1,
                     formatC(ax.breaks, format="g"),
                     as.character(ax.breaks))

# ── Color palette ─────────────────────────────────────────────────────────────
SIG.COLORS <- c(
  "Significant in both"       = "#c0392b",   # red
  "Significance differs"      = "#e67e22",   # orange
  "Not significant in either" = "#7f8c8d")   # grey

# ── Build figure ───────────────────────────────────────────────────────────────

fig2 <- ggplot() +
  # null-association reference lines (OR=1 on each axis)
  geom_vline(xintercept=1, linewidth=0.35, linetype="dashed", color="grey55") +
  geom_hline(yintercept=1, linewidth=0.35, linetype="dashed", color="grey55") +
  # 45-degree agreement line
  geom_line(data=data.table(x=ax.range, y=ax.range),
            aes(x=x, y=y), linewidth=0.55, color="black") +
  annotate("text",
           x=ax.range[2] / 1.07, y=ax.range[2] / 1.35,
           label="Line of\nagreement", size=2.4, color="black",
           hjust=1, vjust=0, lineheight=1.1) +
  # ellipses
  geom_polygon(data=ellipses,
               aes(x=x, y=y, group=level, fill=sig_cat, color=sig_cat),
               alpha=0.18, linewidth=0.40) +
  # center points
  geom_point(data=fig2.wide,
             aes(x=or_A, y=or_D, color=sig_cat),
             size=2.2, shape=16) +
  # co-detection level labels (offset above center)
  geom_text(data=fig2.wide,
            aes(x=or_A, y=or_D, label=level, color=sig_cat),
            size=2.8, vjust=-1.0, hjust=0.5, fontface="italic") +
  # scales
  scale_x_log10(
    name    = "OR \u2014 Analysis A (Unrestricted)",
    limits  = ax.range,
    breaks  = ax.breaks,
    labels  = ax.labels) +
  scale_y_log10(
    name    = paste0("OR \u2014 ", DESIGN.LABELS[DESIGN]),
    limits  = ax.range,
    breaks  = ax.breaks,
    labels  = ax.labels) +
  scale_fill_manual(values=SIG.COLORS,  name=NULL) +
  scale_color_manual(values=SIG.COLORS, name=NULL) +
  # equal aspect ratio in log space so agreement line renders at true 45 degrees
  coord_equal() +
  theme_bw(base_size=10) +
  theme(
    legend.position  = "bottom",
    legend.key.size  = unit(0.45, "cm"),
    panel.grid.minor = element_blank(),
    plot.background  = element_rect(fill="white", color=NA),
    plot.margin      = margin(6, 6, 6, 6))

ggsave("Output/fig2_or_comparison.pdf", fig2, width=6, height=6.5, units="in")
ggsave("Output/fig2_or_comparison.png", fig2, width=6, height=6.5, units="in", dpi=300)

fig2