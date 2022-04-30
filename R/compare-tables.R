library(tidyverse)

make_ci <- function(est, lci, uci, accuracy = 0.01) {
  glue::glue("{scales::number(est, accuracy = accuracy, big.mark = ',')}",
             " (",
             "{scales::number(lci, accuracy = accuracy, big.mark = ',')}",
             ", ",
             "{scales::number(uci, accuracy = accuracy, big.mark = ',')}",
             ")")
}

####################
all_dat <- read_rds(here::here("data", "eumaeus_SCCS_2_211981.rds"))
res <- read_rds(here::here("results", "compare-priors-SCCS_2_211981.rds"))
outcomes <- c(23731, 196347, 196625, 433716, 440367)

summary_stats <- map_dfr(res, ~map_dfr(.x, pluck, "summary_stats"), .id = "prior") %>% 
  mutate(prior = factor(prior, levels = c("res_10", "res_5", "res_1", "res_01"),
                        labels = c("10", "5", "1", "0.1")))

all_draws <- map_dfr(res, ~map_dfr(.x, pluck, "draws"), .id = "prior") %>% 
  mutate(prior = factor(prior, levels = c("res_10", "res_5", "res_1", "res_01"),
                        labels = c("10", "5", "1", "0.1")))

all_plots <- map(res, ~map(.x, pluck, "plot"))

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

####################
all_dat <- read_rds(here::here("data", "eumaeus_HistoricalComparator_1_21215.rds")) %>% 
  filter(periodId == 9)
res <- read_rds(here::here("results", "compare-priors-HistoricalComparator_1_21215.rds"))
outcomes <- c(23731, 196347, 196625, 433716, 440367)

summary_stats <- map_dfr(res, ~map_dfr(.x, pluck, "summary_stats"), .id = "prior") %>% 
  mutate(prior = factor(prior, levels = c("res_10", "res_5", "res_1", "res_01"),
                        labels = c("10", "5", "1", "0.1")))

all_draws <- map_dfr(res, ~map_dfr(.x, pluck, "draws"), .id = "prior") %>% 
  mutate(prior = factor(prior,  levels = c("res_10", "res_5", "res_1", "res_01"),
                        labels = c("10", "5", "1", "0.1")))

all_plots <- map(res, ~map(.x, pluck, "plot"))

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
  select(method, prior, ends_with("Id"),  contains("estimate"), -periodId) %>% 
  gt::gt()

summary_stats %>% 
  filter(.variable == "mu") %>% 
  select(outcomeId, prior, median, q5, q95) %>% 
  mutate(across(c(median, q5, q95), exp)) %>% 
  mutate(`Overall estimate` =  make_ci(median, q5, q95)) %>% 
  arrange(outcomeId) %>% 
  select(outcomeId, prior,  `Overall estimate`) %>% 
  pivot_wider(names_from = prior, values_from = `Overall estimate`) %>% 
  gt::gt()

####################
res <- read_rds(here::here("results", "compare-likelihoods-HistoricalComparator_1_21215.rds"))

summary_stats <- map_dfr(res, ~map_dfr(.x, pluck, "summary_stats"), .id = "model") %>% 
  mutate(model = factor(model, labels = c("normal", "poisson", "grid")))

all_draws <- map_dfr(res, ~map_dfr(.x, pluck, "draws"), .id = "model") %>% 
  mutate(model = factor(model, labels = c("normal", "poisson", "grid")))

all_plots <- map(res, ~map(.x, pluck, "plot"))

summary_stats %>% 
  filter(.variable == "theta_0") %>%
  select(model, M, median, q5, q95, outcomeId) %>% 
  mutate(across(c(median, q5, q95), exp),
         estimate =  make_ci(median, q5, q95)) %>% 
  select(outcomeId, model, M, estimate) %>% 
  pivot_wider(names_from = model, values_from = estimate)
