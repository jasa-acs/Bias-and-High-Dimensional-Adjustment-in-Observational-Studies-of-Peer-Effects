# Statistical tests comparing overall estimates
library(MASS)
library(plyr)
library(reshape2)
library(dplyr)
library(Hmisc)
library(latticeExtra)
library(directlabels)
library(xtable)

source("ps_strat_prep_data_overall_boot.R")


###
# Seemingly unrelated estimators test

suest.test.math <- function(b1, b0, v1, v0, c12 = NULL) {
  require(MASS)
  if (is.null(c12)) {
    v <- abs(v0 - v1)
  } else {
    v <- v0 + v1 - 2 * c12
  }
  ts <- abs(t((b1 - b0)) %*% ginv(v) %*% (b1 - b0))
  df <- qr(v0 - v1)$rank
  p <- 1 - pchisq(ts, df)
  c(ts = ts, df = df, p = p)
}

pssbo.0 <- pssbo %>%
  filter(k == "var_k_l4") %>%
  select(r, model, rr.rel.error, p.0, p.diff.error)

pssbo.w <- dcast(pssbo.0, r ~ model, value.var = "rr.rel.error")

cor(pssbo.w[, 2:10])

rr.cov <- cov(pssbo.w[-1, 2:11])

pssbo.w.1 <- pssbo.w[1, 2:11]

rr.suest <- foreach(i = 1:nrow(rr.cov), .combine = rbind) %:%
  foreach(j = 1:ncol(rr.cov), .combine = c) %do% {
    r <- suest.test.math(
      pssbo.w.1[, i], pssbo.w.1[, j],
      rr.cov[i, i], rr.cov[j, j],
      rr.cov[i, j]
    )['p']
    r
  }
dimnames(rr.suest) <- list(names(pssbo.w)[2:11], names(pssbo.w)[2:11])

rr.suest

rr.suest.m.1 <- apply(
  rr.suest, 2, function(x) {
    ifelse(is.na(x), "", ifelse(x < 1e-12, "< 1e-12", sprintf("%5.2e", x)))
})
rr.suest.m.1[1:9, 1:9]


rr.suest.m.1[lower.tri(rr.suest.m.1, TRUE)] <- NA
rr.suest.m.1[1:7, 2:8]

# num tests
sum(!is.na(rr.suest.m.1[1:7, 2:8]))
.05 / sum(!is.na(rr.suest.m.1[1:7, 2:8]))

## write table of p-values
print(
  xtable(rr.suest.m.1[1:7, 2:8]),
  include.colnames = TRUE,
  include.rownames = TRUE,
  only.contents = TRUE,
  file = "tables/overall_suest_tests.tex"
)
