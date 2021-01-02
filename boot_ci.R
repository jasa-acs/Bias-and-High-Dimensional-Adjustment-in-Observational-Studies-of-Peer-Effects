
percentile.ci <- function(t, alpha = .05, p.value = FALSE, h.0 = 0, type = 9) {
  ci <- unname(quantile(
    t,
    c(alpha / 2, 1 - alpha / 2),
    na.rm = TRUE,
    type = type
    ))
  if (p.value) {
    p <- sum(t > h.0) / length(t)
    list(
      ci = ci,
      p = max(2 * min(p, 1 - p), 2 / length(t))
      )
  } else {
    ci
  }
}

bc.ci <- function(t, t0, alpha = .05, names = FALSE, type = 7) {
  t.ecdf <- ecdf(t)
  z0 <- qnorm(t.ecdf(t0))
  zs <- qnorm(c(alpha / 2, 1 - alpha / 2))
  ci <- quantile(
    x = t,
    probs = pnorm(2 * z0 + zs),
    type = type,
    names = names
    )
  ci
}


clt.ci <- function(t, t0 = mean(t), alpha = .05) {
  r <- t0 + c(-1, 1) * qnorm(1 - alpha / 2) * sd(t)
  return(r)
}

log.clt.ci <- function(t, t0 = mean(t), alpha = .05) {
  s <- sd(log(t))
  r <- log(t0) + c(-1, 1) * qnorm(1 - alpha / 2) * s
  return(exp(r))
}

fieller.boot.ci <- function(a, b, a0, b0, alpha = .05) {
  v <- (a0 / b0)^2 * ((var(a) / a0^2) + var(b) / b0^2)
  a0 / b0 + c(-1, 1) * qnorm(1 - alpha / 2) * sqrt(v)
}

bc.ci.smooth <- function(t, t0, alpha = .05,
                         adjust = 1, ngrid = 1e4,
                         from = min(t) - diff(range(t)),
                         to = max(t) + diff(range(t)), ...) {
  loadNamespace("sROC")
  ecdf.l <- sROC::kCDF(
    x = t,
    bw = bw.ucv(t, lower = min(t) / 1000, upper = max(t) / 2),
    adjust = adjust, ngrid = ngrid, from = from, to = to
    )
  t.ecdf <- approxfun(ecdf.l$x, ecdf.l$Fhat,
                      yleft = 0, yright = 1)
  t.inv.ecdf <- approxfun(ecdf.l$Fhat, ecdf.l$x)
  z0 <- qnorm(t.ecdf(t0))
  zs <- qnorm(c(alpha / 2, 1 - alpha / 2))
  ci <- t.inv.ecdf(pnorm(2 * z0 + zs))
  ci
}

xyplot.bootCI <- function(
  formula, data, CIs, poly = FALSE, CI.alpha = .07,
  border = TRUE, ...
  ) {
  xyplot(
    formula, data, CIs = CIs, poly = poly,
    panel = function(...) {
      panel.grid(v=-4, h=-4, lwd = .8, alpha = .02)
      panel.superpose(...)
    },
    panel.groups = function(x, y, CIs = CIs, subscripts = subscripts, col.line = col.line, ...) {
      CIs <- CIs[subscripts,]
      if (poly==TRUE) {
        for (i in 1:(ncol(CIs) / 2)) {
          panel.polygon(
            x = c(sort(x), rev(sort(x))),
            y = c(CIs[order(x),i],
              CIs[order(-x), ncol(CIs) - (i-1)]),
            col = col.line, alpha = CI.alpha,
            border = NULL
            )
             }
      } else {
        for (i in 1:ncol(CIs)) {
          panel.xyplot(sort(x), y=CIs[order(x), i], ...)
        }
      }
      panel.xyplot(sort(x), y[order(x)], col.line = col.line, ...)
    }, ...)
}
