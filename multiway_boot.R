library(plyr)
library(foreach)
library(digest)

r.double.or.nothing <- function(n) {
  2 * rbinom(n, 1, .5)
}

make.seed <- function(object, salt = '') {
  as.numeric(paste('0x', substr(digest(paste(object, salt)), 10, 15), sep=''))
}

RNG.with.group.seed <- function(r, grouping.name, groups, RNG) {
  x <- foreach(g = groups, .combine = c) %do% {
    set.seed(make.seed(c(r, grouping.name, g)))
    RNG(1)
  }
  x
}

multiway.boot <- function(
  statistic, R,
  groups = as.matrix(1:N),
  verbose = FALSE,
  RNG = r.double.or.nothing,
  .parallel = FALSE,
  .progress = 'none',
  use.seeds = FALSE,
  require.positive.weights = TRUE,
  ...
  ) {
  N.groupingFactors <- ncol(groups)
  groups.f <- list()
  unique.groups <- list()
  groups.num <- list()
  N.groups <- c()
  for (d in 1:N.groupingFactors) {
    groups.f[[d]] <- as.factor(groups[, d])
    unique.groups[[d]] <- levels(groups.f[[d]])
    groups.num[[d]] <- as.numeric(groups.f[[d]])
    N.groups <- c(N.groups, length(unique.groups[[d]]))
  }

  do.one.replicate <- function(r) {
    # Observation weights are products of weights for each factor
    need.weights <- TRUE
    while (need.weights) {
      w <- foreach(i = 1:N.groupingFactors, .combine = `*`) %do% {
        if (use.seeds) {
          RNG.with.group.seed(r, i, unique.groups[[i]], RNG)[groups.num[[i]]]
        } else {
          RNG(N.groups[i])[groups.num[[i]]]
        }
      }
      need.weights <- sum(w, na.rm = T) == 0 & require.positive.weights
    }

    if (verbose) cat(i, " ")
    statistic(..., weights = w)
  }

  llply(
    1:R, do.one.replicate,
    .parallel = .parallel,
    .progress = .progress
    )
}
