library(rstan)
library(tidyverse)
library(tidybayes)
library(furrr)
# library(progressr)
# options(mc.cores = parallel::detectCores()) # when running one at a time
rstan_options(auto_write = TRUE)
source(here::here("R", "get-ests-functions.R"))

eumaeus_file <- "eumaeus_HistoricalComparator_1_21215.rds"

all_dat <- read_rds(here::here("data", eumaeus_file)) %>% 
  mutate(negativeControlId = ifelse(is.na(negativeControlId), outcomeId, negativeControlId),
         site = as.numeric(factor(databaseId, 
                                  levels = c("CCAE", "IBM_MDCD", "IBM_MDCR", "OptumDod", "OptumEhr")))) %>% 
  arrange(site)

# outcomes to test (both positive and negative controls)
positive_control_outcomes <- all_dat %>% 
  mutate(effectSize = ifelse(is.na(effectSize), 1, effectSize)) %>% 
  filter(effectSize %in% c(1, 1.5)) %>% 
  pull(outcomeId) %>% 
  unique()

done <- tibble(files = list.files(here::here("results", "HistoricalComparator"))) %>% 
  separate(files, into = c("i", "t"), sep = "\\_", extra = "drop",
           convert = TRUE)

its <- expand_grid(i = positive_control_outcomes,
                   t = unique(all_dat$periodId),
                   eumaeus_file = eumaeus_file) %>% 
  anti_join(done, by = c("i", "t"))

plan(multisession)
future_pwalk(its,  function(i, t, eumaeus_file, ...){
  write_rds(metaanalysis_poisson(outcome = i, 
                                 all_dat = filter(all_dat, periodId == t),
                                 iter = 1000, make_plot = FALSE),
            here::here("results", "HistoricalComparator", 
                       str_glue("{i}_{t}_{eumaeus_file}")))
}, .options = furrr_options(seed = TRUE))
plan(sequential)

