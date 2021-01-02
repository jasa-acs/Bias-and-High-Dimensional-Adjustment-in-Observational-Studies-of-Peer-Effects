# Overall results, including Figure 2 and Table 2
library(plyr)
library(dplyr)
library(Hmisc)
library(latticeExtra)
library(directlabels)
library(ggplot2)

source("ps_strat_prep_data_overall_boot.R")

source("plotting_settings.R")

pdf("figures/pss_overall_p0_boot.pdf", width = 3.5, height = 2.9)

trellis.par.set(theme = ps.theme.2)
dotplot(
  model ~ p.0,
  lower = pssbo.ci$p.0.ci.clt1,
  upper = pssbo.ci$p.0.ci.clt2,
  horizontal = TRUE,
  groups = factor(model.cat),
  draw.bands = FALSE,
  scales = list(x = list(tick.number = 3, cex = 1)),
  data = pssbo.ci,
  xlab = p.0.expr,
  family = "serif",
  panel = function(x, y, ...) {
    panel.abline(h = y, alpha = .1)
    panel.superpose(x, y, ...)
  },
  panel.groups = function(x, y, lower, upper, group.number, subscripts, ...) {
    panel.segments(lower[subscripts], y, upper[subscripts], y,
                   col = trellis.par.get("superpose.symbol")$col[group.number],
                   alpha = .5, lwd = 1.5)
    panel.xyplot(x, y, ...)
  })

dev.off()


pdf("figures/pss_overall_p0_boot_with_0.pdf", width = 3.5, height = 2.9)

trellis.par.set(theme = ps.theme.2)
dotplot(
  model ~ p.0,
  lower = pssbo.ci$p.0.ci.clt1,
  upper = pssbo.ci$p.0.ci.clt2,
  horizontal = TRUE,
  groups = factor(model.cat),
  draw.bands = FALSE,
  scales = list(x = list(
    at = c(0, 1e-4, 2e-4),
    labels = c("0", "0.0001", "0.0002"),
    cex = .8)),
  data = pssbo.ci,
  #xlab = p.0.expr,
  xlab = "p0",
  xlim = c(0, 2.05e-4),
  family = "serif",
  panel = function(x, y, ...) {
    panel.abline(h = y, alpha = .1)
    panel.superpose(x, y, ...)
  },
  panel.groups = function(x, y, lower, upper, group.number, subscripts, ...) {
    panel.segments(lower[subscripts], y, upper[subscripts], y,
                   col = trellis.par.get("superpose.symbol")$col[group.number],
                   alpha = .5, lwd = 1.5)
    panel.xyplot(x, y, ...)
  })

dev.off()


pdf("figures/pss_overall_rr_boot.pdf", width = 3.5, height = 2.9)

trellis.par.set(theme = ps.theme.2)
dotplot(
  model ~ rr,
  lower = pssbo.ci$rr.ci.clt1,
  upper = pssbo.ci$rr.ci.clt2,
  groups = model.cat,
  data = pssbo.ci,
  xlim = c(4.2, 34),
  xlab = rel.risk.expr,
  panel = function(x, y, ...) {
    panel.abline(h = y, alpha = .1)
    panel.superpose(x, y, ...)
  },
  panel.groups = function(x, y, lower, upper, group.number, subscripts, ...) {
    panel.segments(lower[subscripts], y, upper[subscripts], y,
                   col = trellis.par.get("superpose.symbol")$col[group.number],
                   alpha = .5, lwd = 1.5)
    panel.xyplot(x, y, ...)
  })

dev.off()


pdf("figures/pss_overall_pdiff_boot.pdf", width = 3.5, height = 2.9)

trellis.par.set(theme = ps.theme.2)
dotplot(
  model ~ p.diff,
  lower = pssbo.ci$p.diff.ci.clt1,
  upper = pssbo.ci$p.diff.ci.clt2,
  groups = model.cat,
  xlim = c(.001079, .00129),
  data = pssbo.ci,
  xlab = atet.expr,
  panel = function(x, y, ...) {
    panel.abline(h = y, alpha = .1)
    panel.superpose(x, y, ...)
  },
  panel.groups = function(x, y, lower, upper, group.number, subscripts, ...) {
    panel.segments(lower[subscripts], y, upper[subscripts], y,
                   col = trellis.par.get("superpose.symbol")$col[group.number],
                   alpha = .5, lwd = 1.5)
    panel.xyplot(x, y, ...)
  })

dev.off()

pdf("figures/pss_overall_rr_rel_error_boot.pdf", width = 3.5, height = 2.9)

pc.ne <- subset(pssbo.ci, model != "exp")
trellis.par.set(theme = ps.theme.2)
dotplot(
  model ~ 100 * rr.rel.error,
  lower = 100 * pc.ne$rr.rel.error.ci.clt1,
  upper = 100 * pc.ne$rr.rel.error.ci.clt2,
  groups = model.cat,
  xlim = c(-15, 390),
  data = pssbo.ci,
  scales = list(x = list(at = c(0, 100, 200, 300))),
  xlab = "Percent error in RR",
  panel = function(x, y, ...) {
    panel.abline(h = y, alpha = .1)
    panel.abline(v = 0, alpha = .7,
                 col = trellis.par.get("superpose.symbol")$col['exp'])
    panel.superpose(x, y, ...)
  },
  panel.groups = function(x, y, lower, upper, group.number, subscripts, ...) {
    panel.segments(lower[subscripts], y,
                   upper[subscripts], y,
                   col = trellis.par.get("superpose.symbol")$col[group.number],
                   alpha = .5, lwd = 1.5)
    panel.xyplot(x, y, ...)
  })

dev.off()

pdf("figures/pss_overall_pdiff_rel_error_boot.pdf", width = 3.5, height = 2.9)

pc.ne <- subset(pssbo.ci, model != "exp")
trellis.par.set(theme = ps.theme.2)
dotplot(
  model ~ p.diff.rel.error,
  lower = pc.ne$p.diff.rel.error.ci.clt1,
  upper = pc.ne$p.diff.rel.error.ci.clt2,
  groups = model.cat,
  xlim = c(-.009, .15),
  data = pc.ne,
  xlab = "Relative error in risk difference",
  panel = function(x, y, ...) {
    panel.abline(h = y, alpha = .1)
    panel.abline(v = 0, alpha = .7,
                 col = trellis.par.get("superpose.symbol")$col['exp'])
    panel.superpose(x, y, ...)
  },
  panel.groups = function(x, y, lower, upper, group.number, subscripts, ...) {
    panel.segments(lower[subscripts], y,
                   upper[subscripts], y,
                   col = trellis.par.get("superpose.symbol")$col[group.number],
                   alpha = .5, lwd = 1.5)
    panel.xyplot(x, y, ...)
  })

dev.off()

pdf("figures/pss_overall_pdiff_error_reduction_boot.pdf", width = 3.5, height = 2.9)

trellis.par.set(theme = ps.theme.2)
dotplot(
  model ~ I(100 - p.diff.error.ctn * 100),
  lower = 100 - pssbo.ci$p.diff.error.ctn.ci.clt1 * 100,
  upper = 100 - pssbo.ci$p.diff.error.ctn.ci.clt2 * 100,
  groups = model.cat,
  xlim = c(-5, 105),
  data = pssbo.ci, subset = model != "exp" & model != "naive",
  xlab = percent.red.expr,
  panel = function(x, y, ...) {
    panel.abline(h = y, alpha = .1)
    panel.abline(v = 100, alpha = .5, col = trellis.par.get("superpose.symbol")$col["exp"])
    panel.abline(v = 0, alpha = .5, col = trellis.par.get("superpose.symbol")$col["naive"])
    panel.superpose(x, y, ...)
  },
  panel.groups = function(x, y, lower, upper, group.number, subscripts, ...) {
    panel.segments(lower[subscripts], y,
                   upper[subscripts], y,
                   col = trellis.par.get("superpose.symbol")$col[group.number],
                   alpha = .5, lwd = 1.5)
          panel.xyplot(x, y, ...)
  })

dev.off()

pdf("figures/pss_overall_pdiff_error_rel_to_possible_boot.pdf", width = 3.5, height = 2.9)

trellis.par.set(theme = ps.theme.2)
dotplot(
  model ~ I(100 * p.diff.error.of.max),
  lower = pssbo.ci$p.diff.error.of.max.ci.clt1 * 100,
  upper = pssbo.ci$p.diff.error.of.max.ci.clt2 * 100,
  groups = model.cat,
  xlim = c(-5, 100),
  scales = list(x = list(at = c(0, 20, 40, 60, 80, 100))),
  data = pssbo.ci, subset = model != "exp",
  xlab = percent.poss.expr,
  panel = function(x, y, ...) {
    panel.abline(h = y, alpha = .1)
    panel.abline(v = 0, alpha = .7,
                 col = trellis.par.get("superpose.symbol")$col["exp"])
    panel.superpose(x, y, ...)
  },
  panel.groups = function(x, y, lower, upper, group.number, subscripts, ...) {
    panel.segments(lower[subscripts], y,
                   upper[subscripts], y,
                   col = trellis.par.get("superpose.symbol")$col[group.number],
                   alpha = .5, lwd = 1.5)
    panel.xyplot(x, y, ...)
  })

dev.off()


# change from adding same domain sharing

pdf("figures/pss_overall_pdiff_error_reduction_sds_boot.pdf", width = 4, height = 2)

trellis.par.set(theme = ps.theme.3)
dotplot(
  model.comparison ~ I(-100 * p.diff.error.change.rel),
  groups = c(1, 1, 1, 1),
  lower = -100 * pssbo.sds.ci$p.diff.error.change.rel.ci1,
  upper = -100 * pssbo.sds.ci$p.diff.error.change.rel.ci2,
  data = pssbo.sds.ci,
  col = "black",
  xlim = c(-5, 105),
  xlab = percent.red.sds.expr,
  panel = function(x, y, lower, upper, ...) {
    panel.abline(h = y, alpha = .1)
    panel.segments(lower, y, upper, y, ..., alpha = .5, lwd = 1.5)
    panel.abline(v = 100, alpha = .5, col = ps.theme.2$superpose.symbol$col["exp"])
    panel.abline(v = 0, alpha = .2, col = ps.theme.2$superpose.symbol$col["naive"])
    panel.xyplot(x, y, ...)
  })

dev.off()



####
# Statements in paper

# N

# NECG n
with(pssbo[1,], n.0)

# exp n
with(pssbo[1,], n.1.exp)
with(pssbo[1,], n.0.exp)

with(pssbo[1,], n.1.exp + n.0.exp)

# total n
with(pssbo[1,], n.1 + n.0 + n.1.exp + n.0.exp)

# exp estimates
subset(pssbo.ci, model == "exp")$p.0 * 100
subset(pssbo.ci, model == "exp")$p.1 * 100

subset(pssbo.ci, model == "exp")[,c("p.diff", "p.diff.ci.clt1", "p.diff.ci.clt2")] * 100

subset(pssbo.ci, model == "exp")[,c("rr", "rr.ci.clt1", "rr.ci.clt2")]

# naive
subset(pssbo.ci, model == "naive")[,c("rr", "rr.ci.clt1", "rr.ci.clt2")]

subset(pssbo.ci, model == "naive")[,c("p.diff.rel.error", "p.diff.rel.error.ci.clt1", "p.diff.rel.error.ci.clt2")]

subset(pssbo.ci, model == "naive")[,c("p.diff.error.of.max", "p.diff.error.of.max.ci.clt1", "p.diff.error.of.max.ci.clt2")]

# p.0, RR, p.diff
pssbo.ci %>% select(model, p.0, rr, p.diff)

pssbo.ci %>% select(model, rr, rr.ci.clt1, rr.ci.clt2)

# rel error in delta
pssbo.ci %>% select(model, p.diff.error, p.diff.error.ctn)

# rel error in RR
pssbo.ci %>% select(model, rr.rel.error)

# error reduction in RR, delta
pssbo.ci %>%
  mutate(
    rr.error.red = 1 - rr.rel.error / rr.rel.error[model == "naive"],
    p.diff.error.red = 1 - p.diff.error.ctn
  ) %>%
  select(model, rr.rel.error, rr.error.red, p.diff.error.red)

# of possible error in delta
pssbo.ci %>% select(model, p.diff.error.of.max, p.diff.error.of.max.ci.clt1, p.diff.error.of.max.ci.clt2)

# print a table
library(xtable)

prob.multiplier <- function(p) {
  p * 1e4
}


pssbo.ci.ft <- data.frame(model = pssbo.ci$model)

pssbo.ci.ft$demo <- ifelse(pssbo.ci.ft$model %in% c("exp", "naive"), "", "Yes")
pssbo.ci.ft$general <- ifelse(pssbo.ci.ft$model %in% c("exp", "naive", "D", "Ds"), "", "Yes")
pssbo.ci.ft$sds <- ifelse(grepl("s", pssbo.ci.ft$model), "Yes", "")
pssbo.ci.ft$ods <- ifelse(grepl("M", pssbo.ci.ft$model), "Yes", "")

pssbo.ci.ft$p.0 <- format.num.with.ci(pssbo.ci, "p.0", "%.2f", transform.fnc = prob.multiplier)
pssbo.ci.ft$rr <- format.num.with.ci(pssbo.ci, "rr")
pssbo.ci.ft$p.diff <- format.num.with.ci(pssbo.ci, "p.diff", "%.2f", transform.fnc = prob.multiplier)
#pssbo.ci.ft$p.diff.rel.error <- format.num.with.ci(pssbo.ci, "p.diff.rel.error", "%2.2f",
#                                                   transform.fnc = function(x) {x * 100})

pssbo.ci.ft <- pssbo.ci.ft[nrow(pssbo.ci.ft):1, ]

names(pssbo.ci.ft) <- c(
  "Model",
  "Demographics", "General behavior,\\newline degree",
  "Same domain shares", "Other domain shares",
  "$\\hat{p}^{(0)} \\times 10^4$",
  "$\\widehat{RR}$",
  "$\\hat{\\delta} \\times 10^4$"#,
  #"Percent error in $\\hat{\\delta}$"
  )

print(
  xtable(pssbo.ci.ft),
  include.colnames = TRUE,
  include.rownames = FALSE,
  only.contents = TRUE,
  #hline = c(0, 1, nrow(pssbo.ci.ft)),
  sanitize.text.function = function(x) {x},
  sanitize.colnames.function = function(x) {x},
  file = "tables/overall_ests_by_model.tex"
  )

