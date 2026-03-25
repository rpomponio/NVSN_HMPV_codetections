################################################### -
## Title: Power Calculation for HMPV co-detection
## Author: Ray Pomponio
## Email: pomponiord@upmc.edu
## PI: Anna Wang-Erickson
## Date Created: 2026-03-25
##
## Ref: Whitehead, J. (1993). Sample size calculations for ordered categorical
##      data. Statistics in Medicine, 12(24), 2257–2271.
################################################### -
library(Hmisc)
library(data.table)

# load data from CDC; perform superficial clean-up
cdc <- fread("Data/Pitt_Anna_HMPV_DEC25.csv")
dat <- cdc[research_tested==1]
dat[is.na(respvisitdt), respvisitdt:=as.IDate(as.Date(scrdate, format="%m/%d/%Y"))]
dat[tmpv>=2, tmpv:=NA]
mpv <- dat[tmpv==1]

# compute marginal cell probabilities of the outcome
mpv[, c_severity:=fcase(
  c_finalstatus==1 & inptICU==1 & (c_intubated==1 | c_died==1), 5,  
  c_finalstatus==1 & inptICU==1, 4,
  c_finalstatus==1 & c_respsup==1 & c_suppoxy==1, 3,
  c_finalstatus==1 & c_respsup==1, 3, ## likely CPAP/BIPAP or Nasal Canula
  c_finalstatus==1 & c_respsup==0, 2,
  c_finalstatus>=2, 1,
  default=0)]
p.marg.sev <- as.numeric(table(mpv$c_severity) / nrow(mpv))

# compute number of subjects in exposure/control groups
n1 <- as.numeric(mpv[c_rsv_result==1, .N])
n2 <- as.numeric(nrow(mpv) - n1)

# conduct power calculation for range of minimally-detectable odds ratios
or.seq <-seq(3, 1, -0.1)
po.res <- popower(p=p.marg.sev, odds.ratio=or.seq, n1=n1, n2=n2)
or.star <- or.seq[which.min(po.res$power[po.res$power > 0.8])]

# plot results
plot(or.seq, po.res$power, type="b", xlab="Effect Size (OR)", ylab="Power")
abline(h=0.8, lty=2, col="red")
abline(v=or.star, lty=2, col="green")
mtext(paste("80% Power Achieved at Approx. Effect Size:", or.star), col="green")
