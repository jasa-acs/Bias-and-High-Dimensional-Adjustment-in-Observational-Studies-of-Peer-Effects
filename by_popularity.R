# Comparying estimates in buckets formed by domain prior popularity
library(plyr)
library(Hmisc)
library(foreach)
library(reshape2)
library(directlabels)
library(MASS)
library(dplyr)
library(doMC)
registerDoMC(cores = parallel::detectCores())
options(digits = 3)

source("boot_ci.R")
source("helpers.R")
source("ps_strat_prep_data_overall_boot.R")
source("plotting_settings.R")

pops.names <- gsub(
  "_", ".",
  read.table(
    "data/ps_strata_popularity_error_boot_describe.txt",
    header = FALSE)[,1],
  fixed = TRUE)
pops.names <- gsub("shared.url.p.", "p.", pops.names, fixed = TRUE)
pops.classes <- ifelse(pops.names %in% c("model", "k", "pop.measure"), "character", "numeric")

pops <- read.table(
  "data/ps_strata_popularity_error_boot.txt",
  sep='\t', header=FALSE, quote='', comment.char='',
  colClasses = pops.classes, col.names = pops.names,
  na.strings=c('NULL', 'NA', 'NaN', 'Infinity', '-Infinity'),
  stringsAsFactors = FALSE
  )


pops$model <- factor(pops$model, levels = rev(MODELS))
mc <- cbind(model = unique(pops$model), ldply(unique(pops$model), model.to.cats))
pops$model.s <- mc$model.s[pmatch(pops$model, mc$model, duplicates.ok = T)]
pops$model.M <- mc$model.M[pmatch(pops$model, mc$model, duplicates.ok = T)]
pops$model.cat <- factor(mc$model.cat[pmatch(pops$model, mc$model, duplicates.ok = T)],
                         levels = rev(MODEL.CATS))

####
## do suest for each quintile bucket, each measure

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


do.suest <- function(df, variable = "p.diff.error") {
  df.w <- dcast(df, r ~ model, value.var = variable)
  est.cov <- cov(df.w[-1,-1])
  df.w.1 <- df.w[1, -1]

  est.suest <- foreach(i = 1:nrow(est.cov), .combine = rbind) %:%
  foreach(j = 1:i, .combine = rbind) %do% {
    r <- suest.test.math(
      df.w.1[, i], df.w.1[, j],
      est.cov[i, i], est.cov[j, j],
      est.cov[i, j]
    )['p']
    r <- data.frame(
      i = rownames(est.cov)[i],
      j = colnames(est.cov)[j],
      p = r
    )
    r
  }
  est.suest %>% filter(i != j)
}

# do tests
ss.1 <- pops %>%
  group_by(pop.measure, pop.bucket) %>%
  do({
    tests <- do.suest(.)
    tests
  })

save(ss.1, file = "saved_results/suest_by_prior_pop.RData")

ss.1 %>% filter(i == "exp")

ss.1 <- ss.1 %>%
  ungroup() %>%
  mutate(
    pop.measure.nice = factor(
      pop.measure,
      labels = c(
        "NECG pairs with prior sharing",
        "exposed pairs with prior sharing",
        "fraction NECG pairs with prior sharing",
        "fraction exposed pairs with prior sharing",
        "number of prior sharing users"
      )
      )
    )

ss.2 <- ss.1 %>% filter(pop.measure == "num.users.us.ecdft")

pdf("figures/pss_pop_quintiles_p_values_vs_exp.pdf", width = 5.5, height = 4.5)

trellis.par.set(theme = ps.theme.2.m)
xyp1 <- xyplot(
  pmax(p, 1e-12) ~ I(pop.bucket + 1),
  groups = j,
  data = ss.2 %>% filter(i == "exp"),
  ylab = "p-value (vs. experimental estimate)",
  xlab = "Quintile of number of prior domain sharing users",
  xlim = c(.7, 5.5),
  scales = list(y = list(log = 10)),
  lwd = 1.3, lty = 1, pch = 19, cex = .6,
  type = c("o"),
  panel = function(...) {
    panel.abline(h = log10(.05), lty = 2, alpha = .8)
    panel.abline(h = log10(.05 / 9), lty = 3, alpha = .8)
    panel.abline(h = log10(.05 / (9 * 5)), lty = 3, alpha = .3)
    panel.xyplot(...)
  }
  )
direct.label(
  xyp1,
  list(last.qp, dl.trans(y = y + .1, x = x + .1))
  )

dev.off()

pdf("figures/pss_pop_quintiles_p_values_vs_naive.pdf", width = 5.5, height = 4.5)

trellis.par.set(theme = ps.theme.2.m)
xyp1 <- xyplot(
  pmax(p, 1e-12) ~ I(pop.bucket + 1),
  groups = j,
  data = ss.2 %>% filter(i == "naive"),
  ylab = "p-value (vs. naive estimate)",
  xlab = "Quintile of number of prior domain sharing users",
  xlim = c(.55, 5.1),
  scales = list(y = list(log = 10)),
  lwd = 1.3, lty = 1, pch = 19, cex = .6,
  type = c("o"),
  panel = function(...) {
    panel.abline(h = log10(.05), lty = 2, alpha = .8)
    panel.abline(h = log10(.05 / 9), lty = 3, alpha = .8)
    panel.abline(h = log10(.05 / (9 * 5)), lty = 3, alpha = .3)
    panel.xyplot(...)
  }
  )
direct.label(
  xyp1,
  list(first.qp, dl.trans(y = y + .1, x = x - .1))
  )

dev.off()


###
## Plotting summaries by bucket

pops.point <- pops %>%
  filter(r == -1) %>%
  mutate(
    pop.measure.nice = factor(
      pop.measure,
      labels = c(
        "NECG pairs with prior sharing",
        "exposed pairs with prior sharing",
        "fraction NECG pairs with prior sharing",
        "fraction exposed pairs with prior sharing",
        "number of prior sharing users"
      )
      )
    )

pops.point %>%
  filter(pop.measure == 'num.users.us.ecdft') %>%
  group_by(model) %>%
  mutate(frac.1 = n.1 / sum(n.1)) %>%
  select(model, pop.bucket, p.0, p.1, rr, rr.rel.error, frac.1) %>%
  print(n = 100)

weights.by.quintiles <- pops.point %>%
  filter(model == "exp") %>%
  group_by(pop.measure.nice) %>%
  mutate(frac.1 = n.1 / sum(n.1)) %>%
  select(pop.measure, pop.bucket, frac.1)

## match with overall results
pops.point %>%
  filter(pop.measure == 'num.users.us.ecdft') %>%
  group_by(model) %>%
  summarise(
    p.0 = weighted.mean(p.0, n.1),
    n.1 = sum(n.1)
    )


pdf("figures/pss_pop_quintiles_weights.pdf", width = 4.5, height = 3.5)

xyplot(
  frac.1 ~ I(pop.bucket + 1),
  data = weights.by.quintiles %>%
    filter(pop.measure == 'num.users.us.ecdft'),
  ylab = "Fraction of all exposed pairs",
  xlab = "Quintile of number of prior domain sharing users",
  lwd = 1.3, lty = 1, pch = 19,
  cex = .6,
  ylim = c(0, .6),
  scales = list(y = list(at = c(0, .2, .4, .6))),
  type = c("o"),
  col = "black"
  )

dev.off()

pps.1 <- pops.point %>%
  filter(pop.measure == 'num.users.us.ecdft')

pdf("figures/pss_pop_quintiles_p0.pdf", width = 5.5, height = 4.5)

trellis.par.set(theme = ps.theme.2.m)
xyp1 <- xyplot(
  p.0 ~ I(pop.bucket + 1),
  groups = model,
  data = pps.1,
  ylab = p.0.expr,
  xlab = "Quintile of number of prior domain sharing users",
  xlim = c(0.85, 5.57),
  lwd = 1.3, lty = 1, pch = 19,
  cex = .6,
  type = c("o"),
  panel = function(...) {
    panel.xyplot(...)
  }
  )
direct.label(
  xyp1,
  list(last.qp, dl.trans(y = y - .15, x = x + .1))
)

dev.off()

pdf("figures/pss_pop_quintiles_pdiff.pdf", width = 5.5, height = 4.5)

trellis.par.set(theme = ps.theme.2.m)
xyp1 <- xyplot(
  p.diff ~ I(pop.bucket + 1),
  groups = model,
  data = pps.1,
  ylab = "Estimated risk difference",
  xlab = "Quintile of number of prior domain sharing users",
  xlim = c(0.85, 5.57), 
  lwd = 1.3, lty = 1, pch = 19,
  cex = .6,
  type = c("o"),
  panel = function(...) {
    panel.xyplot(...)
  }
  )
direct.label(
  xyp1,
  list(last.qp, dl.trans(y = y - .15, x = x + .1))
)

dev.off()

pdf("figures/pss_pop_quintiles_rr.pdf", width = 5.5, height = 4.5)

trellis.par.set(theme = ps.theme.2.m)
xyp1 <- xyplot(
  rr ~ I(pop.bucket + 1),
  groups = model,
  data = pps.1,
  ylab = "Estimated relative risk",
  xlab = "Quintile of number of prior domain sharing users",
  xlim = c(0.85, 5.57), ylim = c(4, 120),
  scales = list(y = list(
    log = 10,
    at = c(4, 8, 16, 32, 64, 128),
    labels = c(4, 8, 16, 32, 64, 128)
  )),
  lwd = 1.3, lty = 1, pch = 19,
  cex = .6,
  type = c("o"),
  panel = function(...) {
    panel.xyplot(...)
  }
  )
direct.label(
  xyp1,
  list(last.qp, dl.trans(y = y - .15, x = x + .1))
)

dev.off()


pdf("figures/pss_pop_quintiles_rr_rel_error.pdf", width = 5.5, height = 4.5)

trellis.par.set(theme = ps.theme.2.m)
xyp1 <- xyplot(
  100 * rr.rel.error ~ I(pop.bucket + 1),
  groups = model,
  data = pps.1 %>% filter(model != "exp"),
  ylab = "% error in relative risk",
  xlab = "Quintile of number of prior domain sharing users",
  xlim = c(0.55, 5.1), ylim = c(-50, 800),
  lwd = 1.3, lty = 1, pch = 19,
  cex = .6,
  type = c("o"),
  panel = function(...) {
    panel.abline(h = 0, col = ps.theme.2.m$superpose.line$col['exp'], alpha = .9)
    panel.xyplot(...)
  }
  )
direct.label(
  xyp1,
  list(top.bumpup, dl.trans(y = y + .01, x = x - .4))
)

dev.off()
