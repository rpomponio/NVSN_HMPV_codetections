################################################### -
## Title: Simulation to demonstrate selection bias
## Author: Ray Pomponio
## Email: pomponiord@upmc.edu
## Project: HMPV Co-detections
## Date Created: 2026-04-27
################################################### -

library(data.table)
library(ordinal)
set.seed(1818)

# HARDCODED PARAMETERS WHERE RSV IS MORE-VIRULENT AGENT
N.REP <- 1000
RSV.PARAMS <- c("mu"=60, "sigma"=10)
MPV.PARAMS <- c("mu"=40, "sigma"=15)
THRESHOLD <- 50

sim <- replicate(
  N.REP,
  {
    # chance of detection is 20% for both agents, mutually independent
    dat <- data.table(RSV=rbinom(10000, size=1, prob=0.2),
                      MPV=rbinom(10000, size=1, prob=0.2))
    dat[, codetect:=as.integer(RSV * MPV)]
    
    # virulence determined by agent, if co-detected then max virulence used
    dat[RSV==1, virulence:=rnorm(
      sum(RSV==1), mean=RSV.PARAMS["mu"], sd=RSV.PARAMS["sigma"])]
    dat[MPV==1, virulence:=rnorm(
      sum(MPV==1), mean=MPV.PARAMS["mu"], sd=MPV.PARAMS["sigma"])]
    dat[MPV==1 & RSV==1, virulence:=pmax(
      rnorm(sum(MPV==1 & RSV==1), RSV.PARAMS["mu"], RSV.PARAMS["sigma"]),
      rnorm(sum(MPV==1 & RSV==1), MPV.PARAMS["mu"], MPV.PARAMS["sigma"]))]
    
    # severity (ordinal outcome) ultimately determined by virulence
    dat[, severity:=cut(virulence, c(seq(0, 100, 10), 999))]
    
    # induce selection bias by conditioning on virulence above threshold
    enr <- dat[virulence>=THRESHOLD]
    
    # focusing on MPV cases alone, evaluate effect of co-detection
    mpv <- enr[MPV==1]
    fit0 <- clm(severity ~ codetect, data=mpv, link="logit")
    
    # next evaluate effect of co-detection in general population of MPV
    pop <- dat[MPV==1]
    fit1 <- clm(severity ~ codetect, data=pop, link="logit")
    
    # quantify bias based on under-estimation of true co-detection effect
    fit0$beta - fit1$beta
  }
)
hist(sim)