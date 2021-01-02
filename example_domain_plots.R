# Illustration with example domain in Figure 1
library(lattice)
library(plyr)
library(dplyr)
source("ps_strat_prep_data_overall_boot.R")
source("plotting_settings.R")

pssds.names <- gsub(
  "_", ".",
  read.table(
    "data/ps_strata_by_domain_means_describe.txt",
    header = FALSE)[,1],
  fixed = TRUE)
pssds.names <- gsub("shared.url.p.", "p.", pssds.names, fixed = TRUE)
pssds.classes <- ifelse(pssds.names %in% c("model", "domain"), "character", "numeric")

pssds <- read.table(
  "data/ex_ps_strata_by_domain_means.txt",
  sep='\t', header=FALSE, quote='', comment.char='',
  colClasses = pssds.classes, col.names = pssds.names,
  na.strings=c('NULL', 'NA', 'NaN', 'Infinity', '-Infinity'),
  stringsAsFactors = FALSE
  )

psspe.names <- gsub(
  "_", ".",
  read.table(
    "data/ps_strata_by_domain_est_error_boot_describe.txt",
    header = FALSE)[,1],
  fixed = TRUE)
psspe.names <- append(psspe.names, "num_users_us_q", after = 10)
psspe.names <- gsub("shared.url.p.", "p.", psspe.names, fixed = TRUE)[-3]
psspe.classes <- ifelse(psspe.names %in% c("model", "domain"), "character", "numeric")

psspe <- read.table(
  "data/ex_ps_strata_by_domain_est_error.txt",
  sep='\t', header=FALSE, quote='', comment.char='',
  colClasses = psspe.classes, col.names = psspe.names,
  na.strings=c('NULL', 'NA', 'NaN', 'Infinity', '-Infinity'),
  stringsAsFactors = FALSE
  )

psdp <- read.table(
  "data/ex_ps_percentiles_by_domain_exploded.txt",
  sep='\t', header=FALSE, quote='', comment.char='',
  na.strings=c('NULL', 'NA', 'NaN', 'Infinity', '-Infinity'),
  stringsAsFactors = FALSE
  )
names(psdp) <- c("domain", "model", "ps", "j")
psdp$domain <- as.character(psdp$domain)
psdp$j <- psdp$j + 1

pssds.1 <- pssds
psspe.1 <- psspe
psdp.1 <- psdp

pssds.1 <- pssds.1[order(pssds.1$model, pssds.1$ps.p),]

pssds.1 <- pssds.1 %>%
  group_by(model, domain) %>%
  mutate(
    freq.1 = n.1 / sum(n.1),
    freq.0 = n.0 / sum(n.0),
    ps.approx = freq.1 / (freq.0 + freq.1)
  )

pssds.2 <- pssds.1
psdp.2 <- psdp.1
psspe.2.AMs <- subset(psspe.1, model == "AMs")

pssds.3 <- pssds.2 %>%
  inner_join(psdp.2 %>% mutate(ps.p = j - 1)) %>%
  arrange(model, j) %>%
  group_by(model) %>%
  mutate(
    cumfreq = cumsum(freq.0 + freq.1) / 2
  )

pssds.3$model <- factor(
  pssds.3$model,
  levels = rev(c("exp", "naive", "D", "A", "M", "AM", "Ds", "As", "Ms", "AMs"))
  )


###
# plots

pdf("figures/pss_ex_percentiles_all_models.pdf", width = 7, height = 5)

trellis.par.set(theme = ps.theme.2.m.noexp)
xyplot(
  100 * cumfreq ~ ps | model,
  groups = model,
  data = pssds.3,
  type = "l",
  ylab = "Stratum (ECDF x 100)",
  xlab = "Estimated propensity score",
  layout = c(4, 2),
  scales = list(y = list(rot = 90),
                x = list(at = c(0, .1, .2))),
  ylim = c(0, 100)
)

dev.off()


prop.1 <- with(psspe.2.AMs, n.1.exp / (n.0.obs + n.1.exp))

pdf("figures/pss_ex_percentiles.pdf", width = 4, height = 3.2)

trellis.par.set(theme = ps.theme.ps)
xyplot(
  j ~ ps,
  data = psdp.2,
  subset = model == "AMs",
  type = "s", cex = .2, pch = 20,
  scales = list(y = list(rot = 90)),
  col = "black",
  ylab = expression(paste("Stratum (", ECDF %*% 100, ")")),
  xlab = "Propensity score",
  panel = function(...) {
    panel.abline(v = prop.1, lty = 2, alpha = .5)
    panel.xyplot(...)
  }
  )

dev.off()


pdf("figures/pss_ex_freq_by_strata.pdf", width = 4, height = 3.2)

trellis.par.set(theme = ps.theme.ps)
xyplot(
  freq.0 + freq.1 ~ ps.p,
  data = pssds.2,
  subset = model == "AMs",
  type = "p", cex = .4,
  scales = list(y = list(rot = 90)),
  auto.key = list(columns = 1, lines = F, points = T,
    text = c("Unexposed", "Exposed"), cex = .75,
    x=0.01, y=0.95, corner=c(0,1)),
  ylab = "Relative frequency",
  xlab = "Stratum"
  )

dev.off()

pdf("figures/pss_ex_implicit_weights_by_strata.pdf", width = 4, height = 3.5)

trellis.par.set(theme = ps.theme.ps)
xyplot(
  I(freq.1 / freq.0) ~ ps.p,
  data = pssds.2,
  subset = model == "AMs",
  type = "p", cex = .4,
  col = "black", pch = 4,
  scales = list(y = list(rot = 90)),
  ylab = "Implicit weights",
  xlab = "Stratum"
  )

dev.off()

pssds.2.comb <- ddply(
  pssds.2, .(domain, model), summarise,
  p.1 = weighted.mean(p.1, n.1),
  p.0.w = weighted.mean(p.0, n.1),
  p.0.uw = weighted.mean(p.0, n.0),
  n.1 = sum(n.1),
  n.0 = sum(n.0)
  )

pssds.2.comb.AMs <- subset(pssds.2.comb, model == "AMs")

pdf("figures/pss_ex_p_sharing_by_strata.pdf", width = 4, height = 3.2)

trellis.par.set(theme = ps.theme.ps)
xyplot(
  ifelse(p.0 == 0, 7e-6, p.0) + p.1 ~ ps.p,
  data = pssds.2,
  subset = model == "AMs",
  type = c("p"),
  cex = .4, lwd = 1.2,
  scales = list(
    y = list(log = 10,
      at = c(7e-6, 1e-4, 1e-3, 1e-2),
      label = c(0, 1e-4, 1e-3, 1e-2),
      rot = 90)),
  #auto.key = list(columns = 2, lines = F, points = T,
  #  text = c("Unexposed", "Exposed"), cex = .8),
  ylab = "Proportion sharing",
  xlab = "Stratum",
  panel = function(..., lwd) {
    panel.xyplot(...)
    panel.abline(h = log10(psspe.2.AMs$p.0.exp), col = model.cat.col["exp"], lwd = lwd)
    panel.abline(h = log10(pssds.2.comb.AMs$p.1), col = exposed.col, lwd = lwd)
    panel.abline(h = log10(pssds.2.comb.AMs$p.0.w), col = model.cat.col["Ms"], lwd = lwd)
    panel.abline(h = log10(pssds.2.comb.AMs$p.0.uw), col = model.cat.col["naive"], lwd = lwd)
  }
  )

dev.off()
