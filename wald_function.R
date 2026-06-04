# ── Custom tidy function: polr() output with Wald p-values ────────────────────
#    tbl_regression(tidy_fun=) expects a function with signature f(x, ...) that
#    returns a data frame in broom::tidy() format. We delegate to the default
#    broom tidier and append p-values derived from a two-sided Wald z-test.
#    This is the asymptotic approximation: z = estimate / std.error,
#    p = 2 * pnorm(-|z|). For polr() objects, broom omits p-values by default
#    because the distribution of the test statistic under the null is only
#    asymptotically normal; the approximation is reasonable at this sample size.
tidy.polr.wald <- function(x, exponentiate=FALSE, conf.level=0.95, ...) {
  tidy.out <- broom::tidy(x,
                          exponentiate = exponentiate,
                          conf.level   = conf.level,
                          ...)
  
  # ── Wald z-statistic and two-sided p-value ────────────────────────────────
  #    broom::tidy.polr() returns the raw log-odds coefficient in `estimate`
  #    regardless of exponentiate=, so we re-extract coef() for the z-statistic
  #    rather than back-transforming from the (possibly exponentiated) estimate
  coef.se <- coef(summary(x))[, c("Value", "Std. Error"), drop=FALSE]
  coef.se <- as.data.frame(coef.se)
  coef.se$term <- rownames(coef.se)
  
  tidy.out <- merge(
    tidy.out,
    coef.se[, c("term", "Value", "Std. Error")],
    by   = "term",
    all.x = TRUE)
  
  tidy.out$p.value <- 2 * pnorm(-abs(tidy.out$Value / tidy.out$`Std. Error`))
  
  # ── Drop the working columns; return clean tidy frame ─────────────────────
  tidy.out$Value      <- NULL
  tidy.out$`Std. Error` <- NULL
  tidy.out
}