library(rstan)
library(tidyverse)
library(tidybayes)
library(furrr)
# library(progressr)
# options(mc.cores = parallel::detectCores()) # when running one at a time
rstan_options(auto_write = TRUE)

eumaeus_file <- "eumaeus_HistoricalComparator_1_21215.rds"

all_dat <- read_rds(here::here("data", eumaeus_file)) %>% 
  mutate(site = as.numeric(factor(databaseId))) %>% 
  arrange(site)

# outcomes to test (both positive and negative controls)
positive_control_outcomes <- all_dat %>% 
  mutate(effectSize = ifelse(is.na(effectSize), 1, effectSize),
         negativeControlId = ifelse(is.na(negativeControlId), outcomeId, negativeControlId)) %>% 
  select(outcomeId, negativeControlId, effectSize) %>% 
  distinct()

fit_historical_meta <- function(i, t, positive_control_outcomes, all_dat, 
                                mean_prior_beta = 0, sd_prior_beta = 1,
                                mean_prior_tau = .5, sd_prior_tau = 1,
                                mean_prior_THETA = 0, sd_prior_THETA = 4,
                                mean_prior_GAMMA = 0.5, sd_prior_GAMMA = 2,
                                iter = 2000, chains = 3, refresh = 0,
                                pars = c("THETA", "theta_i", "beta"), ...) {
  
  # choose a single positive control
  positive_control_i <- positive_control_outcomes[i,]
  outcome <- positive_control_i$outcomeId
  
  # but we also need to know what negative control it was based on
  NCoutcome <- positive_control_i$negativeControlId
  
  # and what the "truth" is
  effect_size <- positive_control_i$effectSize
  
  # these are the estimates from each of the sites
  # that had that outcome
  ests <- all_dat %>% 
    filter(outcomeId == outcome,
           periodId == t)
  
  # skip if no estimates
  if (nrow(ests) == 0) return(NULL)
  sites <- unique(ests$site)
  
  # we want to use only the negative controls that were not used to
  # generate the positive control of interest
  # so we exclude them from the negative controls
  # as well as any sites that don't have that positive control
  NCs <- all_dat %>% 
    filter(is.na(effectSize),
           outcomeId != NCoutcome,
           site %in% sites) %>% 
    filter(periodId == t) %>% 
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
              counterfactualOutcomesInterest = as.array(ests$counterfactualOutcomes), 
              counterfactualDaysInterest = as.array(ests$counterfactualDays),
              exposureOutcomesInterest = as.array(ests$exposureOutcomes), 
              exposureDaysInterest = as.array(ests$exposureDays)
  )
  
  priors <- list(mean_prior_beta = as.array(rep(mean_prior_beta, length(sites))),
                 sd_prior_beta = as.array(rep(sd_prior_beta, length(sites))),
                 mean_prior_tau = as.array(rep(mean_prior_tau, length(sites))),
                 sd_prior_tau = as.array(rep(sd_prior_tau, length(sites))),
                 mean_prior_THETA = mean_prior_THETA,
                 sd_prior_THETA = sd_prior_THETA,
                 mean_prior_GAMMA = mean_prior_GAMMA,
                 sd_prior_GAMMA = sd_prior_GAMMA)
  
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
                    outcome_id = outcome,
                    period = t)
  
  to_return
}

done <- tibble(files = list.files(here::here("results", "HistoricalComparator"))) %>% 
  separate(files, into = c("sd_prior_THETA", "i", "t"), sep = "\\_", extra = "drop",
           convert = TRUE)

its <- expand_grid(i = 1:nrow(positive_control_outcomes),
                   t = unique(all_dat$periodId),
                   mean_prior_beta = 0, sd_prior_beta = 1,
                   mean_prior_tau = .5, sd_prior_tau = 1,
                   mean_prior_THETA = 0, sd_prior_THETA = 4,
                   mean_prior_GAMMA = 0.5, sd_prior_GAMMA = 2,
                   eumaeus_file = eumaeus_file) %>% 
  anti_join(done, by = c("i", "t", "sd_prior_THETA"))


plan(multisession)
future_pwalk(its,  function(i, t, eumaeus_file, mean_prior_beta, sd_prior_beta,
                            mean_prior_tau, sd_prior_tau,
                            mean_prior_THETA, sd_prior_THETA,
                            mean_prior_GAMMA, sd_prior_GAMMA, ...){
  write_rds(fit_historical_meta(i, t, positive_control_outcomes = positive_control_outcomes, all_dat = all_dat,
                                mean_prior_beta, sd_prior_beta, mean_prior_tau, sd_prior_tau,
                                mean_prior_THETA, sd_prior_THETA, mean_prior_GAMMA, sd_prior_GAMMA),
            here::here("results", "HistoricalComparator", str_glue("{sd_prior_THETA}_{i}_{t}_{eumaeus_file}")))
}, .options = furrr_options(seed = TRUE))
plan(sequential)

