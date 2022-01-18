library(rstan)
library(tidyverse)
library(tidybayes)
library(furrr)
library(progressr)
# options(mc.cores = parallel::detectCores()) # when running one at a time
rstan_options(auto_write = TRUE)

eumaeus_file <- "eumaeus_HistoricalComparator_1_21215.rds"

all_dat <- read_rds(here::here("data", eumaeus_file)) %>% 
  mutate(site = as.numeric(factor(databaseId))) %>% 
  arrange(site)

# these are all the unique positive controls
# only do those of size 2 (highly correlated with others)
positive_control_outcomes <- all_dat %>% 
  filter(effectSize == 2, !is.na(seLogRr)) %>% 
  select(outcomeId, negativeControlId, effectSize) %>% 
  distinct()

sim_function <- function(i, positive_control_outcomes, all_dat, 
                         iter = 3000, chains = 1, refresh = 0,
                         pars = c("THETA", "theta_i", "beta"), ...) {
  
  # choose a single positive control
  positive_control_i <- positive_control_outcomes[i,]
  outcome <- positive_control_i$outcomeId
  
  # but we also need to know what negative control it was based on
  NCoutcome <- positive_control_i$negativeControlId
  
  # and what the "truth" is
  effect_size <- positive_control_i$effectSize
  
  # these are the estimates (profile likelihoods) from each of the sites
  # that had that positive control
  ests <- all_dat %>% 
    filter(outcomeId == outcome, !is.na(seLogRr)) %>% 
    group_by(databaseId) %>% 
    # just use last period
    filter(periodId == max(periodId)) %>% 
    ungroup() 
  
  # skip if no estimates
  if (nrow(ests) == 0) return(NULL)
  # skip if only 1 site (can't really meta-analyze (stan code doesn't work because expects vector)... will do elsewhere)
  sites <- unique(ests$site)
  if (length(sites) < 2) return(NULL)
  
  # we want to use only the negative controls that were not used to
  # generate the positive control of interest
  # so we exclude them from the negative controls
  # as well as any sites that don't have that positive control
  NCs <- all_dat %>% 
    filter(is.na(effectSize),
           !is.na(seLogRr), 
           outcomeId != NCoutcome,
           site %in% sites) %>% 
    group_by(databaseId, outcomeId) %>% 
    # just use last period
    filter(periodId == max(periodId)) %>% 
    ungroup() 
  
  dat <- list(J = nrow(NCs), # number of negative controls overall
              M = length(sites), # number of sites
              site = as.numeric(factor(NCs$site)), # which site each negative control comes from (1:M)
              # this is because, e.g., site 4 might not be part of this
              # and these are used to index 1:M so can't be 1, 2, 3, 5...
              counterfactualOutcomesNCs = NCs$counterfactualOutcomes, 
              counterfactualDaysNCs = NCs$counterfactualDays,
              exposureOutcomesNCs = NCs$exposureOutcomes, 
              exposureDaysNCs = NCs$exposureDays,
              counterfactualOutcomesInterest = ests$counterfactualOutcomes, 
              counterfactualDaysInterest = ests$counterfactualDays,
              exposureOutcomesInterest = ests$exposureOutcomes, 
              exposureDaysInterest = ests$exposureDays
              ) 
  
  priors <- list(mean_prior_beta = rep(0, length(sites)),
                 sd_prior_beta = rep(1, length(sites)),
                 mean_prior_tau = rep(0.5, length(sites)),
                 sd_prior_tau = rep(1, length(sites)),
                 mean_prior_THETA = 0,
                 sd_prior_THETA = 4,
                 mean_prior_GAMMA = 0.5,
                 sd_prior_GAMMA = 2)
  
  mod <- stan(here::here("stan", "NCs-multiple-sites-priors-real-effect-poisson.stan"),
              data = c(dat, priors), iter = iter, chains = chains,
              pars = pars, refresh = refresh,
              control = list(adapt_delta = 0.9999, max_treedepth = 15))

  # extract the site-specific thetas (true effect) and betas (mean bias)
  little_thetas <- mod %>% 
    spread_draws(theta_i[M], beta[M]) %>% 
    ungroup() %>% 
    mutate(theta = theta_i,
           M = factor(M, labels = levels(factor(NCs$databaseId))))
  # and the mean of the true effect distribution and the bias distribution
  big_theta <- mod %>% 
    spread_draws(THETA) %>% 
    mutate(theta = THETA,
           M = "Overall")
  thetas <- bind_rows(little_thetas, big_theta) %>% 
    select(-theta_i, -THETA)
  
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
with_progress({
  p <- progressor(steps = nrow(positive_control_outcomes))
  res <- future_map(1:nrow(positive_control_outcomes), ~{
    r <- sim_function(.x, positive_control_outcomes, all_dat, iter = 2000, chains = 3)
    p()
    r
  }, .options = furrr_options(seed = TRUE))
})
plan(sequential)

write_rds(res, here::here("results", paste0("historical_", eumaeus_file)))
