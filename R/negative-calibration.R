library(tidyverse)
library(EmpiricalCalibration)
library(furrr)
library(progressr)

eumaeus_file <- "eumaeus_HistoricalComparator_1_21215.rds"

all_dat <- read_rds(here::here("data", eumaeus_file)) %>% 
  mutate(site = as.numeric(factor(databaseId))) %>% 
  arrange(site)

# outcomes to test (both positive and negative controls)
unique_control_outcomes <- all_dat %>% 
  mutate(effectSize = ifelse(is.na(effectSize), 1, effectSize),
         negativeControlId = ifelse(is.na(negativeControlId), outcomeId, negativeControlId)) %>% 
  select(outcomeId, negativeControlId, effectSize) %>% 
  distinct()

empirical_calibration_function <- function(i, t, unique_control_outcomes, all_dat, ...) {
  
  # choose a single control
  control_i <- unique_control_outcomes[i,]
  outcome <- control_i$outcomeId
  
  # but we also need to know what negative control it was based on
  NCoutcome <- control_i$negativeControlId
  
  # and what the "truth" is
  effect_size <- control_i$effectSize
  
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
  
  models <- NCs %>% 
    group_by(site) %>% 
    group_map(~convertNullToErrorModel(fitMcmcNull(.x$logRr, .x$seLogRr)))
  
  res <- imap(models, 
       ~calibrateConfidenceInterval(ests$logRr[.y], ests$seLogRr[.y], model = .x))
  
  bind_cols(bind_rows(res), bind_rows(models)) %>% 
    select(-meanSlope, -sdSlope) %>% 
    mutate(M = factor(sites, labels = levels(factor(NCs$databaseId))),
           outcome_id = outcome) %>% 
    rename_with(~paste0(.x, "_neg_calibrated"), c(logRr, logLb95Rr, logUb95Rr, seLogRr)) %>% 
    left_join(ests, by = c("M" = "databaseId", "outcome_id" = "outcomeId")) %>% 
    select(-(databaseName:ll))
}
  
its <- expand_grid(i = 1:nrow(unique_control_outcomes),
                   t = unique(all_dat$periodId))

plan(multisession)
with_progress({
  p <- progressor(steps = nrow(its))
  res <- future_pmap_dfr(its,  function(i, t, ...) {
    r <- empirical_calibration_function(i, t, 
                                        unique_control_outcomes = unique_control_outcomes, 
                                        all_dat = all_dat, ...)
    p()
    r
  }, .options = furrr_options(seed = TRUE))
})
plan(sequential)

write_rds(res, here::here("results", paste0("neg_", eumaeus_file)))



  
  
