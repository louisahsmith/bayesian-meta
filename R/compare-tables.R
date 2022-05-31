library(tidyverse)
library(tidybayes)
theme_set(theme_tidybayes())

make_ci <- function(est, lci, uci, accuracy = 0.01) {
  
  if (is.na(est)) return(NA_character_)
  
  glue::glue("{scales::number(est, accuracy = accuracy, big.mark = ',')}",
             " (",
             "{scales::number(lci, accuracy = accuracy, big.mark = ',')}",
             ", ",
             "{scales::number(uci, accuracy = accuracy, big.mark = ',')}",
             ")")
}

summarize_results_by_prior <- function(eumaeus_file, t, outcome, v = 1, ...) {
  
  all_dat <- read_rds(here::here("data", eumaeus_file)) %>% 
    mutate(effectSize = ifelse(is.na(effectSize), 1, effectSize)) 
  
  emp_res <- read_rds(here::here("results", paste0("emp_", eumaeus_file))) %>% 
    mutate(effectSize = round(exp(logEffectSize), 1))
  
  if (v == 1) {
    res <- read_rds(here::here("results", str_replace(eumaeus_file,
                                                      "eumaeus_",
                                                      "compare-priors-"
    )))
  } else {
    res <- read_rds(here::here("results", str_replace(eumaeus_file,
                                                      "eumaeus_",
                                                      str_glue("compare-priors{v}-")
    )))
  }
  

  
  all_draws <- map_dfr(res, ~map_dfr(.x, pluck, "draws"), .id = "prior") %>% 
    mutate(prior = factor(prior, levels = c("res_10", "res_5", "res_1", "res_01"),
                          labels = c("10", "5", 
                                     "1", "0.1")),
           databaseId = ifelse(is.na(databaseId), "Overall", as.character(databaseId)))
  
  posterior_plot <- ggplot(filter(all_draws, outcomeId == outcome), aes(y = fct_rev(databaseId), x = exp(.value))) +
    geom_vline(xintercept = 1, linetype = "dotted", color = "darkgrey") +
    stat_dist_halfeye() +
    scale_x_log10(limits = c(.1, 10)) +
    labs(x = "RR", y = NULL) +
    facet_grid(cols = vars(prior))
  
  dat_draws <- all_dat %>% 
    select(outcomeId, effectSize) %>% 
    distinct() %>% 
    right_join(all_draws, by = "outcomeId")
  
  summ <- dat_draws %>% 
    group_by(.variable, prior, databaseId, effectSize, outcomeId) %>%
    summarise(mean = mean(.value),
              median = median(.value),
              conf.lower = quantile(.value, .025),
              conf.upper = quantile(.value, .975),
              .groups = "drop")
  
  combined_res <- summ %>%
    full_join(filter(emp_res, periodId == t),
                     by = c("databaseId" = "M", "outcomeId" = "outcome_id", 
                              "effectSize" = "effectSize")) %>% 
    mutate(across(c(mean, median, conf.lower, conf.upper, logRr_redo_calibrated, 
                    logLb95Rr_redo_calibrated, logUb95Rr_redo_calibrated), 
                  exp, .names = "exp_{.col}"),
           databaseId = ifelse(is.na(databaseId), "Overall", databaseId),
           databaseId = fct_relevel(databaseId, "Overall", after = Inf)) %>% 
    filter(!is.na(.variable))
  
  all_plots <- map(res, ~map(.x, pluck, "plot"))
  example_plots <- map(all_plots, pluck, 4)
  
  tab_dat <- combined_res %>% 
    filter(.variable %in% c("theta_0", "mu")) %>% 
    select(method, exposureId, analysisId, outcomeId, databaseId, 
           median, conf.lower, conf.upper, calibratedLogRr, calibratedCi95Lb, calibratedCi95Ub,
           logRr, ci95Lb, ci95Ub, prior, exposureOutcomes) %>% 
    arrange(outcomeId) %>% 
    mutate(across(c(median, conf.lower, conf.upper, calibratedLogRr, logRr), exp)) %>% 
    rowwise() %>% 
    mutate(orig = make_ci(logRr, ci95Lb, ci95Ub),
           calibrated = make_ci(calibratedLogRr, calibratedCi95Lb, calibratedCi95Ub),
           bayesian = make_ci(median, conf.lower, conf.upper)) %>% 
    ungroup() %>% 
    select(outcomeId, databaseId, exposureOutcomes, orig, calibrated, prior, bayesian) %>% 
    mutate(across(-c(exposureOutcomes, prior), as.character)) %>% 
    arrange(outcomeId, databaseId, prior) %>% 
    mutate(across(everything(), as.character))
  
  tab <- tab_dat %>% 
    flextable::as_grouped_data(groups = "outcomeId") %>% 
    flextable::flextable() %>% 
    flextable::valign(i = 1:nrow(tab_dat), j = 1:7, valign = "top", part = "body") %>% 
    flextable::align_nottext_col(align = "center", header = FALSE) %>% 
    flextable::align_text_col(align = "left") %>% 
    flextable::align(i = 1:nrow(tab_dat), j = c(4, 5, 7), align = "center", part = "body") %>% 
    flextable::bg(i = !{count(tab_dat, outcomeId) %>% mutate(n = n + 1) %>% uncount(weights = n) %>% pull(outcomeId)} %in% unique(tab_dat$outcomeId)[seq(1, length(unique(tab_dat$outcomeId)), 2)],
                  j = 1:7, part = "body", bg = "lightgrey") %>% 
    flextable::merge_v(j = c("databaseId"),
                       target = c("databaseId", "orig", "calibrated", "exposureOutcomes")) %>% 
    flextable::set_header_labels(outcomeId = "Outcome",
                                 databaseId = "Database",
                                 exposureOutcomes = "Outcomes in exposed", 
                                 orig = "Crude estimate",
                                 calibrated = "Empirical calibration",
                                 bayesian = "Bayesian approach",
                                 prior = "Prior standard deviation"
    )
  
  if (length(unique(tab_dat$outcomeId)) == 1) {
    tab <- tab %>% 
      flextable::void(j = 1, part = "all") %>%
      flextable::width(j = 1, 0)
  }
  return(lst(posterior_plot, example_plots, tab))
}

summarize_results_by_likelihood <- function(eumaeus_file, t, outcome, ...) {
  all_dat <- read_rds(here::here("data", eumaeus_file)) %>% 
    mutate(effectSize = ifelse(is.na(effectSize), 1, effectSize)) 
  
  emp_res <- read_rds(here::here("results", paste0("emp_", eumaeus_file))) %>% 
    mutate(effectSize = round(exp(logEffectSize), 1))
  
  res <- read_rds(here::here("results", str_replace(eumaeus_file,
                                                    "eumaeus_",
                                                    "compare-likelihoods-"
  )))
  
  all_draws <- map_dfr(res, ~map_dfr(.x, pluck, "draws"), .id = "likelihood") %>% 
    mutate(likelihood = factor(likelihood, levels = c("res_normal", "res_poisson", "res_grid"),
                               labels = c("Normal", "Poisson", "Grid approximation")),
           databaseId = ifelse(is.na(databaseId), "Overall", as.character(databaseId)))
  
  posterior_plot <- ggplot(filter(all_draws, outcomeId == outcome), 
                           aes(y = fct_rev(databaseId), x = exp(.value))) +
    geom_vline(xintercept = 1, linetype = "dotted", color = "darkgrey") +
    stat_dist_halfeye() +
    scale_x_log10(limits = c(.1, 10)) +
    labs(x = "RR", y = NULL) +
    facet_grid(cols = vars(likelihood))
  
  dat_draws <- all_dat %>% 
    select(outcomeId, effectSize) %>% 
    distinct() %>% 
    right_join(all_draws, by = "outcomeId")
  
  summ <- dat_draws %>% 
    group_by(.variable, likelihood, databaseId, effectSize, outcomeId) %>%
    summarise(mean = mean(.value),
              median = median(.value),
              conf.lower = quantile(.value, .025),
              conf.upper = quantile(.value, .975),
              .groups = "drop")
  
  combined_res <- summ %>%
    full_join(filter(emp_res, periodId == t),
              by = c("databaseId" = "M", "outcomeId" = "outcome_id", 
                     "effectSize" = "effectSize")) %>% 
    mutate(across(c(mean, median, conf.lower, conf.upper, logRr_redo_calibrated, 
                    logLb95Rr_redo_calibrated, logUb95Rr_redo_calibrated), 
                  exp, .names = "exp_{.col}"),
           databaseId = ifelse(is.na(databaseId), "Overall", databaseId),
           databaseId = fct_relevel(databaseId, "Overall", after = Inf))
  
  all_plots <- map(res, ~map(.x, pluck, "plot"))
  example_plots <- map(all_plots, pluck, 4)
  
  tab <- combined_res %>% 
    filter(.variable %in% c("theta_0", "mu")) %>% 
    select(method, exposureId, analysisId, outcomeId, databaseId, 
           median, conf.lower, conf.upper, calibratedLogRr, calibratedCi95Lb, calibratedCi95Ub,
           logRr, ci95Lb, ci95Ub, likelihood, exposureOutcomes) %>% 
    arrange(outcomeId) %>% 
    mutate(across(c(median, conf.lower, conf.upper, calibratedLogRr, logRr), exp)) %>% 
    rowwise() %>% 
    mutate(orig = make_ci(logRr, ci95Lb, ci95Ub),
           calibrated = make_ci(calibratedLogRr, calibratedCi95Lb, calibratedCi95Ub),
           bayesian = make_ci(median, conf.lower, conf.upper)) %>% 
    ungroup() %>% 
    select(outcomeId, databaseId, exposureOutcomes, orig, calibrated, likelihood, bayesian) %>% 
    mutate(across(-c(exposureOutcomes, likelihood), as.character)) %>% 
    arrange(outcomeId, databaseId, likelihood) %>% 
    flextable::as_grouped_data(groups = "outcomeId") %>% 
    flextable::flextable() %>% 
    flextable::merge_v(j = c("databaseId"),
                       target = c("databaseId", "orig", "calibrated", "exposureOutcomes")) %>% 
    flextable::set_header_labels(outcomeId = "Outcome",
                                 databaseId = "Database",
                                 exposureOutcomes = "Outcomes in exposed", 
                                 orig = "Raw data",
                                 calibrated = "Empirical calibration w/ negative controls",
                                 bayesian = "Bayesian approach",
                                 likelihood = "Likelihood used"
    )
}

# sccs_priors <- summarize_results_by_prior("eumaeus_SCCS_2_211981.rds", 12, outcome = 196347)
# hc_priors <- summarize_results_by_prior("eumaeus_HistoricalComparator_1_21215.rds", 9, outcome = 196347)
# sccs_priors <- summarize_results_by_prior("eumaeus_SCCS_2_211981.rds", 12, outcome = 432513, v = 2)
# hc_priors <- summarize_results_by_prior("eumaeus_HistoricalComparator_1_21215.rds", 9, outcome = 432513, v = 2)
# hc_likelihoods <- summarize_results_by_likelihood("eumaeus_HistoricalComparator_1_21215.rds", 9)
sccs_priors <- summarize_results_by_prior("eumaeus_SCCS_1_21215.rds", 9, outcome = 432513, v = 1)

# write_rds(lst(sccs_priors, hc_priors, hc_likelihoods),
#           here::here("documents", "tables-figures-compare.rds"))
write_rds(sccs_priors,  here::here("documents", "sccs_priors-compare.rds"))
