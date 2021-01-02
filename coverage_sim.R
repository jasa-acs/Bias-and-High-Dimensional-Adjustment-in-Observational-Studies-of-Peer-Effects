library(igraph)
library(foreach)
library(plyr)
library(dplyr)
library(tidyr)
library(splines)
library(glmnet)
library(randomizr)
library(doMC)
library(ggplot2)
theme_set(theme_bw())
options(digits = 3)

source("boot_ci.R")
source("multiway_boot.R")
source("coverage_sim_functions.R")

###
# example graph and simulation
latent.space.graph(16, 2, 0)

one.result <- do.sim(n.users = 2^10, n.urls = 2^6, beta = 1)
one.result


###
#  many simulations

# to actually produce the figures, need 100 replicates
n.replicates <- 100

params.df <- expand.grid(
  rep = 1:n.replicates,
  n.users = 2^(10:11),
  n.urls = 2^(5:7),
  beta = c(0, 1)
  )

# registerDoMC(cores = 38)
registerDoSEQ()

results <- foreach(
  i = 1:nrow(params.df),
  .errorhandling = "pass"
) %dopar% {
  params <- params.df[i, ]
  r1 <- do.sim(
    n.users = params$n.users,
    n.urls = params$n.urls,
    beta = params$beta
  )
  r1 <- cbind(r1, params)
  r1
}

saveRDS(results, "saved_results/coverage_sim_results.RDS")

# results <- readRDS("saved_results/coverage_sim_results.RDS")

valid <- results %>% sapply(is.data.frame)
results[!valid]
table(valid)


results.1 <- results[valid] %>% bind_rows()

results.1 <- results.1 %>%
  group_by(estimator, beta) %>%
  mutate(
    p.diff.true.mean = mean(true.p.diff),
    rr.true.mean = mean(true.y.bar.1) / mean(true.y.bar.0),
    p.diff.cov.fp = (true.p.diff > p.diff.ci1 &
                       true.p.diff < p.diff.ci2),
    p.diff.cov.sp = (p.diff.true.mean > p.diff.ci1 &
                       p.diff.true.mean < p.diff.ci2),
    rr.cov.fp = (true.rr > rr.ci1 &
                   true.rr < rr.ci2),
    rr.cov.sp = (rr.true.mean > rr.ci1 &
                   rr.true.mean < rr.ci2)
  ) %>%
  group_by(estimator, beta, n.users, n.urls) %>%
  mutate(
    rep.rank.rr = rank(rr, ties.method = "first"),
    rep.rank.p.diff = rank(p.diff, ties.method = "first")
    )

rs <- results.1 %>%
  group_by(n.users, n.urls, beta, estimator) %>%
  summarise(
    p.diff.true.mean = mean(true.p.diff),
    p.diff.sd = sd(p.diff),
    p.diff.se.mean = mean(p.diff.se),
    p.diff.cov.fp = mean(p.diff.cov.fp),
    p.diff.cov.sp = mean(p.diff.cov.sp),
    rr.true.mean = mean(true.y.bar.1) / mean(true.y.bar.0),
    rr.sd = sd(rr),
    rr.se.mean = mean(rr.se),
    rr.cov.fp = mean(replace.nas(rr.cov.fp)),
    rr.cov.sp = mean(replace.nas(rr.cov.sp))
    )


rs.long <- rs %>%
  gather("variable", "value", p.diff.true.mean:rr.cov.sp)

###
# plot simulation results

results.2 <- results.1 %>% filter(
  estimator == "obs",
  n.users %in% 2^(10:11),
  n.urls %in% 2^(5:7)
)

rs.2 <- rs %>% filter(
  estimator == "obs",
  n.users %in% 2^(10:11),
  n.urls %in% 2^(5:7)
)

pdf("figures/coverage_sim_p_diff_cov_caterpillar.pdf", width = 7, height = 6)

ggplot(
  aes(ymin = p.diff.ci1, ymax = p.diff.ci2,
      y = p.diff,
      x = rep.rank.p.diff,
      group = rep.rank.p.diff,
      color = (p.diff.true.mean > p.diff.ci1 & p.diff.true.mean < p.diff.ci2)
      ),
  data = results.2
) +
  facet_grid(
    n.users ~ n.urls,
    labeller = label_bquote(
      cols = `number of URLs`: .(n.urls),
      rows = `number of users`: .(n.users)
    )
  ) +
  geom_hline(aes(yintercept = p.diff.true.mean), color = "blue") +
  geom_linerange(size = .5) +
  geom_point(size = .3) +
  scale_color_manual(values = c("red", "#444444")) +
  xlab("replication (sorted by estimate)") +
  ylab("estimated risk difference and CI") +
  geom_text(
    data = rs.2,
    aes(label = sprintf("%s%%", 100 * p.diff.cov.sp),
        y = p.diff.true.mean + .02),
    hjust = 0,
    inherit.aes = FALSE,
    x = 2, size = 3
    ) +
  theme(
    strip.background = element_rect(fill="white"),
    legend.position = "none",
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank()
  )

 dev.off()

pdf("figures/coverage_sim_rr_cov_caterpillar.pdf", width = 7, height = 8)

ggplot(
  aes(ymin = rr.ci1, ymax = rr.ci2,
      y = rr,
      x = rep.rank.rr,
      group = rep.rank.rr,
      color = (rr.true.mean > rr.ci1 & rr.true.mean < rr.ci2)
      ),
  data = results.2
) +
  facet_grid(
    n.users + beta ~ n.urls,
    labeller = label_bquote(
      cols = `number of URLs`: .(n.urls),
      rows = `number of users`: .(n.users)
    )
  ) +
  geom_hline(aes(yintercept = rr.true.mean), color = "blue") +
  geom_linerange(size = .5) +
  geom_point(size = .3) +
  scale_color_manual(values = c("red", "#444444")) +
  xlab("replication (sorted by estimate)") +
  ylab("estimated relative risk and CI") +
  coord_cartesian(ylim = c(0, 7.5)) +
  geom_text(
    data = rs.2,
    aes(label = sprintf("%s%% coverage with beta = %s", 100 * rr.cov.sp, beta)),
    hjust = 0,
    inherit.aes = FALSE,
    x = 1, y = 6.5, size = 3
    ) +
  theme(
    strip.background = element_rect(fill="white"),
    legend.position = "none",
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank()
  )


dev.off()


###
# expected degrees in random graph model
ed <- foreach(n = 2^(6:12), .combine = rbind) %do% {
  d <- foreach(i = 1:100, .combine = c) %dopar% {
    latent.space.graph(n, 2, 0)$A %>% rowSums() %>% mean()
  }
  data.frame(
    n = n,
    i = 1:100,
    mean.degree = d
  )
}


saveRDS(ed, "saved_results/latent_space_degrees.RDS")

# ed <- readRDS("saved_results/latent_space_degrees.RDS")

eds <- ed %>%
  group_by(n) %>%
  summarise(
    mean.degree.se = sd(mean.degree),
    mean.degree = mean(mean.degree)
  ) %>%
  mutate(
    mean.degree.ratio = mean.degree / log10(n)
  )


pdf("figures/coverage_sim_latent_space_degrees.pdf", width = 3.5, height = 3)

ggplot(
  aes(x = n, y = mean.degree#,
      #ymax = mean.degree + 1.96 * mean.degree.se,
      #ymin = mean.degree - 1.96 * mean.degree.se
      ),
  data = eds
) +
  scale_x_log10(breaks = unique(ed$n)) +
  geom_line() +
  geom_point() +
  xlab("|V|") + ylab("mean degree")

dev.off()
