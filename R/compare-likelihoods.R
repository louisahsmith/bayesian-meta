library(rstan)
library(tidyverse)
library(tictoc)
library(tidybayes)
theme_set(theme_tidybayes())
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

source(here::here("R", "get-ests-functions.R"))

eumaeus_file <- "eumaeus_HistoricalComparator_1_21215.rds"

all_dat <- read_rds(here::here("data", eumaeus_file)) %>% 
  filter(periodId == 9) %>% 
  mutate(site = as.numeric(factor(databaseId, 
                                  levels = c("CCAE", "IBM_MDCD", "IBM_MDCR", "OptumDod", "OptumEhr")))) %>% 
  arrange(site)

outcomes <- all_dat %>% 
  filter(effectSize == 2) %>% 
  pull(outcomeId) %>% 
  unique()

res_normal <- list()
res_poisson <- list()
res_grid <- list()

tic("normal")
for (i in seq_along(outcomes)[1:5]) {
  res_normal[[i]] <- metaanalysis_normal(outcomes[i], all_dat = all_dat,
                                         iter = 2500)
}
toc() # 95.474 sec elapsed

tic("poisson")
for (i in seq_along(outcomes)[1:5]) {
  res_poisson[[i]] <- metaanalysis_poisson(outcomes[i], all_dat = all_dat,
                                          iter = 2000)
}
toc() # 138.646 sec elapsed

tic("grid")
for (i in seq_along(outcomes)[1:5]) {
  res_grid[[i]] <- metaanalysis_grid(outcomes[i], all_dat = all_dat,
                                     iter = 3500)
}
toc() # 208.663 sec elapsed # 309.317 sec elapsed (3000 to get ESS similar)


summary_stats <- map_dfr(list(res_normal, res_poisson, res_grid), 
                         ~map_dfr(.x, pluck, "summary_stats"), .id = "model") %>% 
  mutate(model = factor(model, labels = c("normal", "poisson", "grid")))

all_draws <- map_dfr(list(res_normal, res_poisson, res_grid), 
                         ~map_dfr(.x, pluck, "draws"), .id = "model") %>% 
  mutate(model = factor(model, labels = c("normal", "poisson", "grid")))

all_plots <- map(list(res_normal, res_poisson, res_grid), 
                     ~map(.x, pluck, "plot"))

summary_stats %>% 
  filter(.variable == "theta_0") %>%
  group_by(model) %>%
  summarise(across(c(mean, median, rhat, ess_bulk, ess_tail), ~mean(.x)))

# all_plots[[1]][[1]]
# all_plots[[2]][[1]]
# all_plots[[3]][[1]]
