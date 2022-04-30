library(rstan)
library(tidyverse)
library(tidybayes)
theme_set(theme_tidybayes())
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

source(here::here("R", "get-ests-functions.R"))

eumaeus_file <- "eumaeus_HistoricalComparator_1_21215.rds"

all_dat <- read_rds(here::here("data", eumaeus_file)) %>% 
  filter(periodId == 9) %>% 
  mutate(negativeControlId = ifelse(is.na(effectSize), outcomeId, negativeControlId),
    site = as.numeric(factor(databaseId, 
                                  levels = c("CCAE", "IBM_MDCD", "IBM_MDCR", "OptumDod", "OptumEhr")))) %>% 
  arrange(site)

outcomes <- c(23731, 196347, 196625, 433716, 440367)

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

write_rds(lst(res_10, res_5, res_1, res_01),
          here::here("results", "compare-priors-HistoricalComparator_1_21215.rds"))
