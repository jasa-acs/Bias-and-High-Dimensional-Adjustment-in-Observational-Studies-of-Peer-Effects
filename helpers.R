

# which models to include in analysis, in this order
PS.MODELS <- c("D", "A", "M", "AM", "Ds", "As", "Ms", "AMs")
EXP <- "exp"
NAIVE <- "naive"
MODELS <- c(EXP, NAIVE, PS.MODELS)
MODEL.CATS <- c(EXP, NAIVE, "basic", "M", "s", "Ms")

model.to.cats <- function(m) {
  m <- as.character(m)
  model.s <- grepl("s", m)
  model.M <- grepl("M", m)
  if (m %in% c(EXP, NAIVE)) {
    model.cat <- m
  } else if (model.s & !model.M) {
    model.cat <- "s"
  } else if (!model.s & model.M) {
    model.cat <- "M"
  } else if (model.s & model.M) {
    model.cat <- "Ms"
  } else {
    model.cat <- "basic"
  }
  data.frame(model.s = model.s, model.M = model.M, model.cat = model.cat)
}


format.num.with.ci <- function(
  df, variable,
  one.format = "%2.2f", ci.suffix = ".ci.clt",
  line.break = TRUE,
  transform.fnc = function(x) x
  ) {
  if (line.break) {
    format = sprintf(
      "\\shortstack{ %s \\\\ \\relax [%s, %s] }",
      one.format, one.format, one.format
      )
  } else {
    format = sprintf(
      "%s [%s, %s]",
      one.format, one.format, one.format
      )
  }
  x <- sprintf(
    format,
    transform.fnc(df[,variable]),
    transform.fnc(df[,sprintf("%s%s1", variable, ci.suffix)]),
    transform.fnc(df[,sprintf("%s%s2", variable, ci.suffix)])
    )
  x
}

weighted.mean.ips.se <- function(y, weights, ...) {
  y.bar <- weighted.mean(y, weights, ...)
  weights <- weights / sum(weights)
  y.bar.var.hat <- sum(weights^2 * (y - y.bar)^2, ...)
  sqrt(y.bar.var.hat)
}
