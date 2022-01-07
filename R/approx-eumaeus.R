library(rstan)
library(tidyverse)
library(tidybayes)
library(furrr)
# options(mc.cores = parallel::detectCores()) # when running one at a time
rstan_options(auto_write = TRUE)

# not all have positive controls
# eumaeus_CohortMethod_1_21184.rds has 39 positive controls with multiple sites
# eumaeus_CohortMethod_1_21215.rds has 294 positive controls with multiple sites
eumaeus_file <- "eumaeus_HistoricalComparator_1_21215.rds"

all_dat <- read_rds(here::here("data", eumaeus_file)) %>% 
  mutate(site = as.numeric(factor(databaseId))) %>% 
  arrange(site)

# these are all the unique positive controls
positive_control_outcomes <- all_dat %>% 
  filter(!is.na(effectSize), !is.na(seLogRr)) %>% 
  select(outcomeId, negativeControlId, effectSize) %>% 
  distinct()

sim_function <- function(i, positive_control_outcomes, all_dat, 
                         min_val = log(.1), max_val = log(10), n_vals = 1000,
                         iter = 3000, chains = 1, refresh = 0,
                         pars = c("BETA", "THETA", "theta_i", "beta"), ...) {
  
  # choose a single positive control
  positive_control_i <- positive_control_outcomes[i,]
  outcome <- positive_control_i$outcomeId
  
  # but we also need to know what negative control it was based on
  NCoutcome <- positive_control_i$negativeControlId
  
  # and what the "truth" is
  effect_size <- positive_control_i$effectSize
  
  # these are the estimates (profile likelihoods) from each of the sites
  # that had that positive control
  init_ests <- all_dat %>% 
    filter(outcomeId == outcome, !is.na(seLogRr), 
           # remove any without concave likelihood within bounds
           between(logRr, min_val, max_val)) %>% 
    group_by(databaseId) %>% 
    # just use last period
    filter(periodId == max(periodId)) %>% 
    ungroup() 
  
  # skip if no estimates
  if (nrow(init_ests) == 0) return(NULL)
  # skip if only 1 site (can't really meta-analyze (stan code doesn't work because expects vector)... will do elsewhere)
  sites <- unique(init_ests$site)
  if (length(sites) < 2) return(NULL)
  
  # if the grid is not regular, interpolate first
  # note that this really just checks that there are the right number, 
  # not necessarily that they are evenly spaced, but we will assume 
  # if the right number, they are regular
  if (!all(map_lgl(init_ests$ll, ~nrow(.x) == n_vals))) {
    ests <- init_ests %>% 
      mutate(new_ll = map(ll, ~as_tibble(approx(.x$point, .x$value, xout = seq(min_val, max_val, length.out = n_vals))))) %>% 
      select(site, new_ll) %>% 
      unnest(new_ll) %>% 
      rename("point" = x, "value" = y)
  } else {
    ests <- init_ests %>% 
      select(site, ll) %>% 
      unnest(ll)
  }
  
  # this is just helping to index the log-likelihood approximations
  for_grid <- ests %>% 
    rowid_to_column() %>% 
    group_by(site) %>% 
    slice(1, n())
  
  # we want to use only the negative controls that were not used to
  # generate the positive control of interest
  # so we exclude them from the negative controls
  # as well as any sites that don't have that positive control
  NCs_all <- all_dat %>% 
    filter(is.na(effectSize),
           !is.na(seLogRr), 
           outcomeId != NCoutcome,
           site %in% sites) %>% 
    group_by(databaseId, outcomeId) %>% 
    # just use last period
    filter(periodId == max(periodId)) %>% 
    ungroup() 
  
  dat <- list(J = nrow(NCs_all), # number of negative controls overall
              b = NCs_all$logRr, # negative control point estimates
              se_b = NCs_all$seLogRr, # negative control standard errors
              M = length(sites), # number of sites
              site = as.numeric(factor(NCs_all$site)), # which site each negative control comes from (1:M)
              # this is because, e.g., site 4 might not be part of this
              # and these are used to index 1:M so can't be 1, 2, 3, 5...
              G = nrow(ests), # number of sites
              ll = ests$value, # log likelihood estimates
              vals_evaled = ests$point, # values at which the log likelihood was estimated
              ll_indices = mutate(group_by(ests, site), rowid = row_number())$rowid, # indexes the log likelihood within site
              grid_id = matrix(for_grid$rowid, ncol = 2, byrow = TRUE), # start and end indices of the log-likelihood vectors
              # recalculate in case regular but forgot to provide appropriate sequence
              n_vals = tapply(ests$point, ests$site, length), 
              min_vals = tapply(ests$point, ests$site, min),
              max_vals = tapply(ests$point, ests$site, max)) 
  
  priors <- list(mean_prior_BETA = 0,
                 sd_prior_BETA = 4,
                 mean_prior_TAU = 0,
                 sd_prior_TAU = 1,
                 mean_prior_beta = rep(0, length(sites)),
                 sd_prior_beta = rep(1, length(sites)),
                 mean_prior_tau = rep(0.5, length(sites)),
                 sd_prior_tau = rep(1, length(sites)),
                 mean_prior_THETA = 0,
                 sd_prior_THETA = 4,
                 mean_prior_GAMMA = 0.5,
                 sd_prior_GAMMA = 2)
  
  mod <- stan(here::here("stan", "NCs-multiple-sites-priors-real-effect-likelihood.stan"),
              data = c(dat, priors), iter = iter, chains = chains,
              pars = pars, refresh = refresh,
              control = list(adapt_delta = 0.9999, max_treedepth = 15))

  
  # extract the site-specific thetas (true effect) and betas (mean bias)
  little_thetas <- mod %>% 
    spread_draws(theta_i[M], beta[M]) %>% 
    ungroup() %>% 
    mutate(theta = theta_i,
           M = factor(M, labels = levels(factor(NCs_all$databaseId))))
  # and the mean of the true effect distribution and the bias distribution
  big_theta <- mod %>% 
    spread_draws(THETA, BETA) %>% 
    mutate(theta = THETA,
           beta = BETA,
           M = "Overall")
  thetas <- bind_rows(little_thetas, big_theta) %>% 
    select(-theta_i, -THETA, -BETA)
  
  # where does the true effect size fall in the posterior distribution?
  theta_percentiles <- group_by(thetas, M) %>% 
    summarise(theta_percentile = mean(exp(theta) < effect_size))
  
  to_return <- list(posterior_draws = thetas,
                    theta_percentiles = theta_percentiles,
                    effect_size = effect_size,
                    outcome_id = outcome)
  
  to_return
}

plan(multisession)
res <- future_map(1:nrow(positive_control_outcomes), sim_function, 
                    positive_control_outcomes, all_dat, iter = 2000, chains = 3)
plan(sequential)

write_rds(res, here::here("results", paste0("res_", eumaeus_file)))