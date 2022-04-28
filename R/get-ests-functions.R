metaanalysis_grid <- function(outcome, all_dat, mu_prior_sd = 10,
                              tau_prior_sd = 10, delta_prior_sd = 10,
                              gamma_prior_sd = 10, min_val = log(.1), 
                              max_val = log(10), n_vals = 1000, 
                              extreme_min_val = log(.01), extreme_max_val = log(100),
                              iter = 1000, chains = 4, adapt_delta = 0.9999, refresh = 0) {
  
  dat <- all_dat %>% 
    filter(outcomeId == outcome)
  
  vals_evaled <- seq(min_val, max_val, length.out = n_vals)
  
  # if likelihood is monotone need to extrapolate
  to_extrapolate <- map_lgl(dat$ll, ~which.max(.x$value) %in% c(1, nrow(.x)))
  # if the grid is not regular, interpolate first
  # note that this really just checks that there are the right number
  to_interpolate <- map_lgl(dat$ll, ~nrow(.x) != n_vals)
  
  if (any(to_extrapolate)) {
    warning(str_glue("Likelihood monotone; appoximating likelihood function between ",
                     "{round(extreme_min_val, 2)} and {round(extreme_max_val, 2)}. Consider changing arguments ",
                     "`extreme_min_val` and `extreme_vax_val` if that is not sufficient."))
    by <- vals_evaled[2] - vals_evaled[1]
    min_val <- extreme_min_val
    max_val <- extreme_max_val
    vals_evaled <- seq(min_val, max_val, by = by)
    n_vals <- length(vals_evaled)
    
    dat <- dat %>% 
      mutate(l_mod = map(ll, ~lm(value ~ poly(point, 4), data = .x)),
             pred = map(l_mod, predict, newdata = data.frame(point = vals_evaled)),
             new_ll_extrapolated = map(pred, ~tibble(point = vals_evaled, value = .x)),
             new_ll_interpolated = map(ll, ~as_tibble(approx(.x$point, .x$value, 
                                                             xout = vals_evaled,
                                                             rule = 2))))
    ests <- dat[to_extrapolate,] %>% 
      unnest(new_ll_extrapolated) %>% 
      select(site, point, value)
    
    if (!all(to_extrapolate)) {
      # just to avoid approximating when not necessary
      ests <- dat[!to_extrapolate,] %>% 
        unnest(new_ll_interpolated) %>% 
        rename("point" = x, "value" = y) %>% 
        select(site, point, value) %>% 
        bind_rows(ests) %>% 
        arrange(site)
    }
  } else if (any(to_interpolate)) {
    ests <- dat %>% 
      mutate(new_ll = map(ll, ~as_tibble(approx(.x$point, .x$value, 
                                                xout = vals_evaled,
                                                rule = 2)))) %>% 
      unnest(new_ll) %>% 
      rename("point" = x, "value" = y) %>% 
      select(site, point, value)
  } else {
    ests <- dat %>% 
      select(site, ll) %>% 
      unnest(ll)
  }
  
  NCs <- all_dat %>% 
    filter(outcomeId != dat$negativeControlId[1], 
           is.na(effectSize),
           databaseId %in% dat$databaseId) %>% 
    mutate(site = as.numeric(factor(site)))
  
  stan_dat <- list(M = nrow(dat),
                   ll = matrix(ests$value, ncol = nrow(dat)), # log likelihood estimates
                   vals_evaled = vals_evaled, # values at which the log likelihood was estimated
                   ll_indices = 1:n_vals, 
                   n_vals = n_vals, 
                   min_val = min_val,
                   max_val = max_val,
                   N = nrow(NCs),
                   site = NCs$site,                  
                   x = NCs$logRr,
                   s_j = NCs$seLogRr,
                   mu_prior_sd = mu_prior_sd)
  
  mod <- stan(here::here("stan", "metaanalysis-grid.stan"),
              data = stan_dat, iter = iter, chains = chains,
              control = list(adapt_delta = adapt_delta),
              refresh = refresh, init = 0)
  
  extract_mod(mod, outcome = outcome, labs = levels(factor(dat$databaseId)))
  
}


metaanalysis_normal <- function(outcome, all_dat, iter = 2500, chains = 4,
                                adapt_delta = 0.9999, refresh = 0) {
  dat <- all_dat %>% 
    filter(outcomeId == outcome)
  
  NCs <- all_dat %>% 
    filter(outcomeId != dat$negativeControlId[1], 
           is.na(effectSize),
           databaseId %in% dat$databaseId) %>% 
    mutate(site = as.numeric(factor(site)))
  
  stan_dat <- list(M = nrow(dat),
                   y = as.array(dat$logRr),
                   s_0 = as.array(dat$seLogRr),
                   N = nrow(NCs),
                   site = NCs$site,                  
                   x = NCs$logRr,
                   s_j = NCs$seLogRr)
  
  mod <- stan(here::here("stan", "metaanalysis-normal.stan"),
              data = stan_dat, iter = iter, chains = chains,
              control = list(adapt_delta = adapt_delta),
              refresh = refresh)
  
  extract_mod(mod, outcome = outcome, labs = levels(factor(dat$databaseId)))
  
}

metaanalysis_poisson <- function(outcome, all_dat, iter = 2500, chains = 4,
                                 adapt_delta = 0.9999, refresh = 0) {
  
  all_dat <- all_dat %>% 
    mutate(counterfactualExpected = exposureDays*(counterfactualOutcomes/counterfactualDays))
  
  dat <- all_dat %>% 
    filter(outcomeId == outcome)
  
  NCs <- all_dat %>% 
    filter(outcomeId != dat$negativeControlId[1], 
           is.na(effectSize),
           databaseId %in% dat$databaseId) %>% 
    mutate(site = as.numeric(factor(site)))
  
  stan_dat <- list(M = nrow(dat),
                   y = as.array(dat$exposureOutcomes),
                   y_star = as.array(dat$counterfactualExpected),
                   N = nrow(NCs),
                   site = NCs$site,                  
                   x = NCs$exposureOutcomes,
                   x_star = NCs$counterfactualExpected, 
                   zeros = rep(0, nrow(NCs)))
  
  mod <- stan(here::here("stan", "metaanalysis-poisson.stan"),
              data = stan_dat, iter = iter, chains = chains,
              control = list(adapt_delta = adapt_delta),
              refresh = refresh)
  
  extract_mod(mod, outcome = outcome, labs = levels(factor(dat$databaseId)))
  
}

extract_mod <- function(mod, outcome, labs) {
  draws <- mod %>% 
    gather_draws(theta_0[M], beta_0[M], delta[M], gamma[M], mu, tau) %>% 
    ungroup() %>% 
    mutate(outcomeId = outcome,
           databaseId = factor(M, labels = labs))
  
  summary_stats <- mod %>% 
    gather_draws(theta_0[M], beta_0[M], delta[M], gamma[M], mu, tau) %>% 
    summarise_draws() %>% 
    ungroup() %>% 
    mutate(outcomeId = outcome,
           databaseId = factor(M, labels = labs))
  
  plot <- bayesplot::mcmc_trace(mod, pars = c("mu", "tau"), regex_pars = "theta")
  
  lst(draws, summary_stats, plot)
}

