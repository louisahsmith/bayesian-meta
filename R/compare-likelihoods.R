library(rstan)
library(tidyverse)
library(tidybayes)
theme_set(theme_tidybayes())
library(tictoc)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

source(here::here("R", "get-ests-functions.R"))

eumaeus_file <- "eumaeus_HistoricalComparator_1_21215.rds"

all_dat <- read_rds(here::here("data", eumaeus_file)) %>% 
  filter(periodId == 9) %>% 
  mutate(negativeControlId = ifelse(is.na(negativeControlId), outcomeId, negativeControlId),
    site = as.numeric(factor(databaseId, 
                                  levels = c("CCAE", "IBM_MDCD", "IBM_MDCR", "OptumDod", "OptumEhr")))) %>% 
  arrange(site)

outcomes <- c(23731,196347, 196625, 433716, 440367)

res_normal <- list()
res_poisson <- list()
res_grid <- list()

tic("normal")
for (i in seq_along(outcomes)) {
  res_normal[[i]] <- metaanalysis_normal(outcomes[i], all_dat = all_dat,
                                         iter = 2500)
}
toc()

tic("poisson")
for (i in seq_along(outcomes)) {
  res_poisson[[i]] <- metaanalysis_poisson(outcomes[i], all_dat = all_dat,
                                          iter = 2000)
}
toc()

tic("grid")
for (i in seq_along(outcomes)) {
  res_grid[[i]] <- metaanalysis_grid(outcomes[i], all_dat = all_dat,
                                     iter = 3500)
}
toc()

write_rds(lst(res_normal, res_poisson, res_grid),
          here::here("results", "compare-likelihoods-HistoricalComparator_1_21215.rds"))

# normal: 176.279 sec elapsed
# poisson: 207.841 sec elapsed
# grid: 828.558 sec elapsed