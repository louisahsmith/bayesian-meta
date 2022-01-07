library(tidyverse)
library(EmpiricalCalibration)
library(furrr)
library(progressr)

# not all have positive controls
# eumaeus_CohortMethod_1_21184.rds has 39 positive controls with multiple sites
# eumaeus_CohortMethod_1_21215.rds has 294 positive controls with multiple sites
# eumaeus_HistoricalComparator_1_21215.rds has 258 positive controls with multiple sites
eumaeus_file <- "eumaeus_CohortMethod_1_21184.rds"

all_dat <- read_rds(here::here("data", eumaeus_file)) %>% 
  mutate(site = as.numeric(factor(databaseId))) %>% 
  arrange(site)

# these are all the unique positive controls
positive_control_outcomes <- all_dat %>% 
  filter(!is.na(effectSize), !is.na(seLogRr)) %>% 
  select(outcomeId, negativeControlId, effectSize) %>% 
  distinct()

empirical_calibration_function <- function(i, positive_control_outcomes, all_dat, 
                                           min_val = log(.1), max_val = log(10), ...) {
  
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
  
  models <- NCs_all %>% 
    group_by(site) %>% 
    group_map(~convertNullToErrorModel(fitMcmcNull(.x$logRr, .x$seLogRr)))
  
  res <- imap(models, 
       ~calibrateConfidenceInterval(init_ests$logRr[.y], init_ests$seLogRr[.y], model = .x))
  
  bind_cols(bind_rows(res), bind_rows(models)) %>% 
    select(-meanSlope, -sdSlope) %>% 
    mutate(site = sites)
}
  
plan(multisession)
with_progress({
  p <- progressor(steps = nrow(positive_control_outcomes))
  res <- future_map(1:nrow(positive_control_outcomes), ~{
    r <- empirical_calibration_function(.x, positive_control_outcomes, all_dat)
    p()
    r
  }, .options = furrr_options(seed = TRUE))
})
plan(sequential)

write_rds(res, here::here("results", paste0("emp_", eumaeus_file)))



  
  
