# using normal likelihood and simulated data

library(rstan)
library(tidyverse)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

run_sim <- function(M = 5, THETA = log(3), GAMMA = 0.5, BETA = log(2), TAU = 0.5, ...) {

  n_negative_controls <- floor(runif(M, 50, 150))
  site <- reduce(map(1:M, ~rep(.x, n_negative_controls[.x])), c)
  
  tau <- rep(.5, M)
  
  beta <- rnorm(M, BETA, TAU)
  true_b <- reduce(map(1:M, ~rnorm(n_negative_controls[.x], beta[.x], tau[.x])), c)
  se_b <- runif(sum(n_negative_controls), .2, 1)
  b <- rnorm(sum(n_negative_controls), true_b, se_b)
  
  b_i <- rnorm(M, beta, tau)
  theta_i <- rnorm(M, THETA, GAMMA)
  theta_tilde_i <- theta_i + b_i
  se_theta <- runif(M, .2, 1)
  theta_hat <- rnorm(M, theta_tilde_i, se_theta)
  
  priors <- list(
                 mean_prior_BETA = 0,
                 sd_prior_BETA = 4,
                 mean_prior_TAU = 0,
                 sd_prior_TAU = 1,
                 mean_prior_beta = c(0, 0, 0, 0, 0),
                 sd_prior_beta = c(1, 1, 1, 1, 1),
                 mean_prior_tau = c(0.5, 0.5, 0.5, 0.5, 0.5),
                 sd_prior_tau = c(1, 1, 1, 1, 1),
                 mean_prior_THETA = 0,
                 sd_prior_THETA = 4,
                 mean_prior_GAMMA = 0.5,
                 sd_prior_GAMMA = 2)
  
  dat <- list(J = sum(n_negative_controls),
              b = b,
              se_b = se_b,
              M = M,
              site = site,
              theta_hat = theta_hat,
              se_theta = se_theta)
  
  mod <- stan(here::here("stan", "NCs-multiple-sites-priors-real-effect.stan"),
              data = c(dat, priors), control = list(adapt_delta = 0.9999, max_treedepth = 15))
  
  list(
    mod_res = summary(mod, pars = c("beta", "tau", "THETA", "theta_i"), probs = c(0.025, .25, .5, .75, 0.975))$summary[,c(1, 4:8)],
    vals = c(beta, tau, THETA, theta_i),
    ests = theta_hat
  )
}

res <- replicate(100, run_sim())
res_res <- map(seq(1, 300, 3), ~pluck(res, .x))
res_true <- map(seq(2, 300, 3), ~pluck(res, .x))
is_in_int <- function(x, y, left, right) imap_lgl(y, ~between(.x, x[.y, left], x[.y, right]))
in_int_mean <- function(left, right, res_res, res_true) {
  map2(res_res, res_true, is_in_int, left = left, right = right) %>% 
    reduce(rbind) %>% 
    colMeans()
}
in_int_95 <- in_int_mean(2, 6, res_res, res_true)
in_int_50 <- in_int_mean(3, 5, res_res, res_true)

