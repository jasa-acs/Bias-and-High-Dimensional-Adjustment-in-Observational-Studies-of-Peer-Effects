library(plyr)
library(Hmisc)
library(foreach)
library(doMC)
library(lattice)
registerDoMC(cores = 3)

source("helpers.R")
source("boot_ci.R")

# prep data
pssbo.names <- gsub(
  "_", ".",
  read.table(
    "data/pssb_overall_describe.txt",
    header = FALSE)[,1],
  fixed = TRUE)
pssbo.names <- gsub("shared.url.p.", "p.", pssbo.names, fixed = TRUE)
pssbo.classes <- ifelse(pssbo.names %in% c("model", "k"), "character", "numeric")
pssbo <- read.table(
  "data/pssb_overall.txt",
  sep='\t', header=FALSE, quote='', comment.char='',
  colClasses = pssbo.classes, col.names = pssbo.names,
  na.strings=c('NULL', 'NA', 'NaN', 'Infinity', '-Infinity'),
  stringsAsFactors = FALSE
  )

# diagnostics
table(pssbo$r, pssbo$model)


pssbo$model <- factor(pssbo$model, levels = rev(c("exp", "naive", "D", "S", "B", "A", "Ds", "Ss", "Bs", "As", "M", "AM", "Ms", "AMs")))
pssbo$model.s <- as.character(pssbo$model)
pssbo$model.s <- ifelse(substr(pssbo$model.s, nchar(pssbo$model.s), nchar(pssbo$model.s)) == 's', TRUE, FALSE)
pssbo$model.M <- as.character(pssbo$model)
pssbo$model.M <- grepl("M", pssbo$model.M)

pssbo$model.cat <- as.character(pssbo$model)
pssbo$model.cat[!pssbo$model.s & !pssbo$model.M & !pssbo$model %in% c("exp", "naive")] <- "basic"
pssbo$model.cat[pssbo$model.s & !pssbo$model.M] <- "s"
pssbo$model.cat[!pssbo$model.s & pssbo$model.M] <- "M"
pssbo$model.cat[pssbo$model.s & pssbo$model.M] <- "Ms"
pssbo$model.cat <- factor(pssbo$model.cat, levels = rev(c("exp", "naive", "basic", "M", "s", "Ms")))

pssbo$rr.rel.error <- with(pssbo, rr / rr[model == "exp"] - 1)

save(pssbo, file = "saved_results/pssbo.RData")
load("saved_results/pssbo.RData")

pssbo.1 <- subset(pssbo, k == "var_k_l4")

pssbo.ci <- ddply(
  pssbo.1, .(k, model, model.s, model.M, model.cat),
  function(df) {
    with(df, 
    c(p.0 = p.0[r == -1],
      p.0.ci = percentile.ci(p.0[r != -1]),
      p.0.ci.clt = clt.ci(p.0[r != -1], t0 = p.0[r == -1]),
      p.1 = p.1.full[r == -1],
      p.1.ci = percentile.ci(p.1.full[r != -1]),
      p.diff = p.diff[r == -1],
      p.diff.ci = percentile.ci(p.diff[r != -1]),
      p.diff.ci.clt = clt.ci(p.diff[r != -1], t0 = p.diff[r == -1]),
      p.diff.error = p.diff.error[r == -1],
      p.diff.error.ci = percentile.ci(p.diff.error[r != -1]),
      p.diff.error.of.max = p.diff.error[r == -1] / p.0.exp[r == -1],
      p.diff.error.of.max.ci = percentile.ci(p.diff.error[r != -1] / p.0.exp[r != -1]),
      p.diff.error.of.max.ci.clt = clt.ci(p.diff.error[r != -1] / p.0.exp[r != -1],
        t0 = p.diff.error[r == -1] / p.0.exp[r == -1]),
      rr = rr[r == -1],
      rr.ci = percentile.ci(rr[r != -1]),
      rr.ci.clt = clt.ci(rr[r != -1], t0 = rr[r == -1]),
      rr.ci.f = fieller.boot.ci(p.1.full[r != -1], p.0[r != -1], p.1.full[r == -1], p.0[r == -1]),
      rr.rel.error = rr.rel.error[r == -1],
      rr.rel.error.ci.clt = clt.ci(rr.rel.error[r != -1], t0 = rr.rel.error[r == -1]),
      rr.rel.error.se = sd(rr.rel.error[r != -1]),
      p.diff.rel.error = p.diff.error[r == -1] / p.diff.exp[r == -1],
      p.diff.rel.error.ci = percentile.ci(p.diff.error[r != -1] / p.diff.exp[r != -1]),
      p.diff.rel.error.ci.clt = clt.ci(p.diff.error[r != -1] / p.diff.exp[r != -1],
                                       t0 = p.diff.error[r == -1] / p.diff.exp[r == -1]),
      p.diff.rel.error.ci.f = fieller.boot.ci(p.diff.error[r != -1], p.diff.exp[r != -1],
                                                p.diff.error[r == -1], p.diff.exp[r == -1]),
      p.diff.error.ctn = abs(p.diff.error[r == -1] / p.diff.obs.error[r == -1]),
      p.diff.error.ctn.ci = percentile.ci(abs(p.diff.error[r != -1] / p.diff.obs.error[r != -1])),
      p.diff.error.ctn.ci.clt = clt.ci(abs(p.diff.error[r != -1] / p.diff.obs.error[r != -1]),
        t0 = abs(p.diff.error[r == -1] / p.diff.obs.error[r == -1]))
      ))
  })

pssbo.ci <- subset(
  pssbo.ci,
  model %in% c("exp", "naive", "D", "A", "M", "AM", "Ds", "As", "Ms", "AMs")
  )
pssbo.ci$model <- factor(
  pssbo.ci$model,
  levels = rev(c("exp", "naive", "D", "A", "M", "AM", "Ds", "As", "Ms", "AMs"))
  )
pssbo.ci$model.cat <- droplevels(pssbo.ci$model.cat)

save(pssbo.ci, file = "saved_results/pssbo.ci.RData")
load("saved_results/pssbo.ci.RData")

write.table(
  format(pssbo.ci, digits=7),
  file = "saved_results/pssbo.ci.tsv",
  quote = FALSE, sep = "\t", row.names = FALSE
)

##
# same domain sharing comparisons

pssbo.sds <- subset(pssbo, model %in% c("D", "Ds", "A", "As", "M", "Ms", "AM", "AMs"))

pssbo.sds$model.comparison <- paste(
  sub("s", "", pssbo.sds$model),
  " vs. ",
  sub("s", "", pssbo.sds$model), "s",
  sep = ""
  )

pssbo.sds.2 <- ddply(
  pssbo.sds, .(model.comparison, r, k),
  function(df) {
    if (nrow(df) == 2) {
      return (c(p.diff.error.0 = df$p.diff.error[!df$model.s],
                p.diff.error.change = (abs(df$p.diff.error[df$model.s]) -
                                       abs(df$p.diff.error[!df$model.s]))
                ))
    }
  })

pssbo.sds <- merge(subset(pssbo.sds, !model.s),
                  pssbo.sds.2,
                  all.x = FALSE, all.y = FALSE)

pssbo.sds.ci <- ddply(pssbo.sds, .(model.comparison, k), function(df) {
  c(p.diff.error.0 = df$p.diff.error.0[df$r == -1],
    p.diff.error.0.ci = percentile.ci(df$p.diff.error.0[df$r != -1]),
    p.diff.error.change = df$p.diff.error.change[df$r == -1],
    p.diff.error.change.ci = percentile.ci(df$p.diff.error.change[df$r != -1]),
    p.diff.error.change.rel = df$p.diff.error.change[df$r == -1] / df$p.diff.error.0[df$r == -1],
    p.diff.error.change.rel.ci = percentile.ci(df$p.diff.error.change[df$r != -1] / df$p.diff.error.0[df$r != -1])
    )
})

save(pssbo.sds.ci, file = "saved_results/pssbo.sds.ci.RData")
load("saved_results/pssbo.sds.ci.RData")

write.table(
  format(pssbo.sds.ci, digits=3),
  file = "saved_results/pssbo.sds.ci.tsv",
  quote = FALSE, sep = "\t", row.names = FALSE
  )

pssbo.sds.ci.all <- pssbo.sds.ci
pssbo.sds.ci <- subset(pssbo.sds.ci, k == "var_k_l4")
