################################################### -
## Title: Figures — HMPV co-detection & illness severity
## Author: Ray Pomponio
## Email: pomponiord@upmc.edu
## Project: HMPV Co-detection (NVSN)
## Date Created: 2026-06-30
## Note: Code written with assistance from Claude (Anthropic);
##       all analytical decisions made by the author
################################################### -

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
     sprintf("Excluded: multiple co-detections\nN = %s", fmt.n(N.EXCL.MULTI))),
  
  # left fork arm: Analysis A (no further restriction)
  mb("b_arm_a",   1.8,  5.2, W.A, H,  CLR.A,
     sprintf("Analysis A\n(Unrestricted)\nN = %s", fmt.n(N.DAT))),
  
  # right fork arm: CT restriction (shared through HMPV CT step)
  mb("b_hmpv_ct", 6.8,  5.2, W.A, H,  CLR.MAIN,
     sprintf("HMPV CT \u2264 %d\nN = %s", CT.THRESHOLD, fmt.n(N.AFTER.HMPV.CT))),
  mb("e_hmpvct",  9.8,  5.9, W.E, HE, CLR.EXCL,
     sprintf("Excluded: HMPV CT >%d\nor missing\nN = %s", CT.THRESHOLD, fmt.n(N.EXCL.HMPV.CT)))
)

# ── Partner CT step (design-specific) ─────────────────────────────────────────
# Design B: single red exclusion box — all partner CT failures dropped
# Design C: two boxes — red for dropped (CT missing/inconclusive) and
#           yellow for reclassified (CT >threshold; cases are retained)
#           b_design shifts down to accommodate the second box

if (DESIGN == "B_restricted") {
  
  partner.boxes <- list(
    mb("e_partner", 9.8, 3.9, W.E, HE, CLR.EXCL,
       sprintf("Excluded: partner CT >%d,\nmissing, or inconclusive\nN = %s",
               CT.THRESHOLD, fmt.n(N.DROPPED.PARTNER))),
    mb("b_design",  6.8, 3.2, W.A, H, CLR.D,
       sprintf("%s\nN = %s", DESIGN.LABELS[DESIGN], fmt.n(N.DESIGN))))
  
  y.junc.partner <- 3.90
  yt.design      <- 3.2 + H / 2
  ylim.bottom    <- 2.2
  
} else {
  
  partner.boxes <- list(
    mb("e_partner_drop",      9.8, 4.3, W.E, HE, CLR.EXCL,
       sprintf("Dropped: partner CT missing\nor inconclusive\nN = %s",
               fmt.n(N.DROPPED.PARTNER))),
    mb("e_partner_reclassify",9.8, 3.4, W.E, HE, "#fef9e7",
       sprintf("Reclassified to HMPV-only:\npartner CT >%d\nN = %s",
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