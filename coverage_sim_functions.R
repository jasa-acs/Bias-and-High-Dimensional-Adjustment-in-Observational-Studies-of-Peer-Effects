
# simulate data for sharing of a single URL
gen.one.url.data <- function(A, x.user, p.holdout = 1/16, intercept = -4, beta = 1, n_steps = Inf) {
  n <- nrow(A)
  u.url <- rnorm(1, 0, .5)
  x <- intercept - rowSums(x.user) + rnorm(n) + u.url
  z <- rbinom(n, 1, 1 - p.holdout)

  # potential outcomes
  y0 <- ifelse(x > 0, 1, 0)
  y1 <- ifelse(beta + x > 0, 1, 0)

  # exposure and observed outcome
  y <- y0

  t = 0
  while (TRUE) {
    d <- pmin(as.vector(A %*% y), 1)
    y_new <- ifelse(z == 1 & d == 1, y1, y0)
    y_change = y_new - y
    y = y_new
    
    if (sum(y_change) == 0) break
    t = t + 1
    if (t >= n_steps) break
  }

  data.frame(
    user = 1:n,
    x = x.user,
    u.url = u.url,
    degree = rowSums(A),
    z, d,
    y0 = y0, y1 = y1,
    y = y
    )
}

# estimation using experimental variation
exp.est <- function(df, weights = NULL) {
  if (!is.null(weights)) {
    df <- df[weights > 0, ]
  }
  f <- with(
    df %>% filter(d == 1),
    lm.fit(cbind(1, z), y)
  )
  data.frame(
    p.diff = coef(f)[2],
    rr = sum(coef(f)) / coef(f)[1]
    )
}

# choice of number of strata k based on n
get.num.strata <- function(n, p = 1/2, a = 4) {
  min(n, floor(a * n^p))
}

# estimation using propensity score stratification
pss.est <- function(df, weights = NULL, k = get.num.strata, ...) {
  n <- sum(df$z * df$d)
  if (!is.null(weights)) {
    df <- df[weights > 0, ]
  }
  # exclude experimental data
  df <- df %>% filter(z == 1)

  # sample according to URL prevalence
  dfu = df %>%
    group_by(url) %>%
    summarise(
      n.1 = sum(d),
      n.0 = sum(1 - d),
      ratio = n.0 / n.1,
      x.1.1 = mean(x.1[d == 1]),
      x.1.0 = mean(x.1[d == 0])
    ) %>%
    mutate(
      pinc = min(ratio) / ratio
    )

  df$url_d = ifelse(df$d == 1, 0, df$url)
  df$inc <- randomizr::block_ra(df$url_d, block_prob = pmax(0, pmin(1, c(1, dfu$pinc))))

  df <- df %>% filter(inc == 1)

  # fit model
  if (min(sum(df$d), sum(1 - df$d)) > 1) {
    mm <- sparse.model.matrix(
      ~ x.1 + x.2,
      data = df
    )
    m <- glmnet(
      mm,
      df$d,
      family = "binomial"
    )

    lambda <- max(0.5 / nrow(df), min(m$lambda))
    df$ps <- predict(m, newx = mm, type = "response", s = lambda)[,1]
  } else {
    df$ps = mean(df$d)
  }

  # choose number of strata k
  if (is.function(k)) {
    k <- k(n, ...)
  }
  df$stratum <- ntile(df$ps, k)
  
  tmp <- df %>%
    group_by(stratum) %>%
    mutate(
      n.treated.stratum = sum(d),
      n.control.stratum = sum(1 - d)
    ) %>%
    filter(
      n.treated.stratum > 0
      ) %>%
    group_by(d, stratum) %>%
    summarise(
      y.bar = mean(y),
      n.treated.stratum = n.treated.stratum[1],
      n.control.stratum = n.control.stratum[1]
    ) %>%
    group_by(d) %>%
    summarise(
      y.bar = weighted.mean(y.bar, n.treated.stratum)
    )
  if (nrow(tmp) != 2) {
    warning("No matches.")
    #browser()
    est <- data.frame(p.diff = NA, rr = NA)
  } else {
    est <- tmp %>%
      summarise(
        p.diff = y.bar[d == 1] - y.bar[d == 0],
        rr = y.bar[d == 1] / y.bar[d == 0]
      )
  }
  return(est)
}

inv.logit <- function(x) exp(x) / (1 + exp(x))

latent.space.graph <- function(n, d, alpha = 1, rescaled = TRUE) {
  if (is.function(alpha))
    alpha <- alpha(n, d)
  x <- matrix(rnorm(n * d), ncol = d)
  x <- x[order(x[, 1]), ]
  p <- inv.logit(alpha - as.matrix(dist(x)))
  if (rescaled) {
    p <- p * (10 * log10(n) / n)
  }
  A <- Matrix(1 * (p > matrix(runif(n^2), ncol = n)), sparse = TRUE)
  diag(A) <- 0
  list(A = A, x = x)
}

ed <- foreach(n = 2^(4:10), .combine = rbind) %do% {
  d <- foreach(i = 1:10, .combine = c) %do% {
    latent.space.graph(n, 2, 1)$A %>% rowSums() %>% mean()
  }
  data.frame(
    n = n,
    mean.degree = mean(d)
  ) %>%
    mutate(
      mean.degree.ratio = mean.degree / log10(n)
      )
}

  

do.sim <- function(n.users = 100, n.urls = 100, R = 100, beta = 0, conf.level = .9, with.exp = FALSE) {
  lsg <- latent.space.graph(
    n.users,
    d = 2,
    alpha = 0,
    rescaled = TRUE
  )
  A <- lsg$A
  x.user <- lsg$x
  rm(lsg)

  df <- foreach(i = 1:n.urls, .combine = bind_rows) %do% {
    min.cond <- 0
    while (min.cond == 0) {
      df1 <- gen.one.url.data(A, x.user, beta = beta, p.holdout = with.exp * 1/16)
      min.cond <- min(sum(df1$d), sum(1 - df1$d))
    }
    df1$url <- i
    df1
  }
  
  obs.ests <- pss.est(df)
  

  oracle <- df %>%
    filter(d == 1) %>%
    summarise(
      true.y.bar.0 = mean(y0),
      true.y.bar.1 = mean(y1),
      true.p.diff = true.y.bar.1 - true.y.bar.0,
      true.rr = true.y.bar.1 / true.y.bar.0
    )

  ci.alpha <- 1 - conf.level

  obs.boot.ests <- multiway.boot(
    statistic = pss.est,
    R = R,
    groups = cbind(df$user, df$url),
    df = df,
  ) %>% bind_rows() %>%
    summarise(
      p.diff.se = sd(p.diff),
      rr.se = sd(rr),
      p.diff.ci1 = clt.ci(p.diff, t0 = obs.ests$p.diff, alpha = ci.alpha)[1],
      p.diff.ci2 = clt.ci(p.diff, t0 = obs.ests$p.diff, alpha = ci.alpha)[2],
      rr.ci1 = clt.ci(rr, t0 = obs.ests$rr, alpha = ci.alpha)[1],
      rr.ci2 = clt.ci(rr, t0 = obs.ests$rr, alpha = ci.alpha)[2]
    )

  if (with.exp) {
    exp.ests <- exp.est(df)
  
    exp.boot.ests <- multiway.boot(
      statistic = exp.est,
      R = R,
      groups = cbind(df$user, df$url),
      df = df
    ) %>% bind_rows() %>%
      summarise(
        p.diff.se = sd(p.diff),
        rr.se = sd(rr),
        p.diff.ci1 = clt.ci(p.diff, t0 = exp.ests$p.diff, alpha = ci.alpha)[1],
        p.diff.ci2 = clt.ci(p.diff, t0 = exp.ests$p.diff, alpha = ci.alpha)[2],
        rr.ci1 = clt.ci(rr, t0 = exp.ests$rr, alpha = ci.alpha)[1],
        rr.ci2 = clt.ci(rr, t0 = exp.ests$rr, alpha = ci.alpha)[2]
      )
    ests <- rbind(
      cbind(obs.ests, obs.boot.ests),
      cbind(exp.ests, exp.boot.ests)
    ) %>%
      mutate(
        estimator = c("obs", "exp")
      )
    cbind(ests, oracle)
  } else {
    ests <- cbind(obs.ests, obs.boot.ests) %>%
    mutate(
      estimator = c("obs")
    )
    cbind(ests, oracle)

  }

}

replace.nas <- function(x, replace.with = 0) {
  x[is.na(x)] <- replace.with
  x
}
