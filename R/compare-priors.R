library(rstan)
library(tidyverse)
library(tidybayes)
theme_set(theme_tidybayes())
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

source(here::here("R", "get-ests-functions.R"))

eumaeus_file <- "eumaeus_SCCS_2_211981.rds"

all_dat <- read_rds(here::here("data", eumaeus_file)) %>% 
  filter(periodId == 12) %>% 
  mutate(negativeControlId = ifelse(is.na(effectSize), outcomeId, negativeControlId),
    site = as.numeric(factor(databaseId, 
                                  levels = c("CCAE", "IBM_MDCD", "IBM_MDCR", "OptumDod", "OptumEhr")))) %>% 
  arrange(site)

outcomes <- c(23731,196347, 196625, 433716, 440367)
res_10 <- list()
res_5 <- list()
res_1 <- list()
res_01 <- list()
for (i in seq_along(outcomes)) {
  res_10[[i]] <- metaanalysis_grid(outcomes[i], all_dat = all_dat, mu_sd = 10,
                                   iter = 5000)
}
for (i in seq_along(outcomes)) {
  res_5[[i]] <- metaanalysis_grid(outcomes[i], all_dat = all_dat, mu_sd = 5,
                                  iter = 5000)
}
for (i in seq_along(outcomes)) {
  res_1[[i]] <- metaanalysis_grid(outcomes[i], all_dat = all_dat, mu_sd = 1,
                                  iter = 5000)
}
for (i in seq_along(outcomes)) {
  res_01[[i]] <- metaanalysis_grid(outcomes[i], all_dat = all_dat, mu_sd = 0.1,
                                   iter = 5000)
}

summary_stats <- map_dfr(list(res_10, res_5, res_1, res_01), 
                         ~map_dfr(.x, pluck, "summary_stats"), .id = "prior") %>% 
  mutate(prior = factor(prior, labels = c("10", "5", "1", "0.1")))

all_draws <- map_dfr(list(res_10, res_5, res_1, res_01),
                     ~map_dfr(.x, pluck, "draws"), .id = "prior") %>% 
  mutate(prior = factor(prior, labels = c("10", "5", "1", "0.1")))

all_plots <- map(list(res_10, res_5, res_1, res_01), 
                 ~map(.x, pluck, "plot"))

make_ci <- function(est, lci, uci, accuracy = 0.01) {
  glue::glue("{scales::number(est, accuracy = accuracy, big.mark = ',')}",
             " (",
             "{scales::number(lci, accuracy = accuracy, big.mark = ',')}",
             ", ",
             "{scales::number(uci, accuracy = accuracy, big.mark = ',')}",
             ")")
}

all_dat %>% 
  filter(outcomeId %in% outcomes) %>% 
  right_join({summary_stats %>% 
      filter(.variable == "theta_0")}, by = c("databaseId", "outcomeId")) %>% 
  select(method, exposureId, analysisId, periodId, outcomeId, databaseId, 
         median, q5, q95, calibratedLogRr, calibratedCi95Lb, calibratedCi95Ub,
         logRr, ci95Lb, ci95Ub, prior) %>% 
  arrange(outcomeId) %>% 
  mutate(across(c(median, q5, q95, calibratedLogRr, logRr), exp)) %>% 
  mutate(`Original estimate` = make_ci(logRr, ci95Lb, ci95Ub),
         `Calibrated estimate` = make_ci(calibratedLogRr, calibratedCi95Lb, calibratedCi95Ub),
         `My estimate` = make_ci(median, q5, q95)) %>% 
  select(method, prior, ends_with("Id"), contains("estimate")) %>% 
  gt::gt()


summary_stats %>% 
  filter(.variable == "mu") %>% 
  select(outcomeId, prior, median, q5, q95) %>% 
  mutate(across(c(median, q5, q95), exp)) %>% 
  mutate(`Overall estimate` =  make_ci(median, q5, q95)) %>% 
  arrange(outcomeId) %>% 
  select(outcomeId, prior,  `Overall estimate`) %>% 
  gt::gt()
  
