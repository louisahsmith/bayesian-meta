library(tidyverse)
library(gt)
library(tidybayes)
theme_set(theme_tidybayes())

eumaeus_file <- "eumaeus_HistoricalComparator_1_21215.rds"

dat <- read_rds(here::here("data", eumaeus_file))

emp_res <- read_rds(here::here("results", paste0("emp_", eumaeus_file)))

res <- map(list.files(here::here("results", "HistoricalComparator"), full.names = TRUE),
           read_rds) %>% 
  set_names(list.files(here::here("results", "HistoricalComparator")))

res_dat <- unnest_wider(tibble(res = res,
                               file = names(res)), 
                        res) %>% 
  separate(file, into = c("outcome", "t"), 
           sep = "\\_", extra = "drop",
           convert = TRUE)

all_stats <- unnest(res_dat, summary_stats) %>% 
  select(-draws, -plot)

long_pct <- unnest(res_dat, theta_percentiles) %>% 
  select(-posterior_draws)

my_quantile <- function(x, probs) {
  tibble(x = quantile(x, probs), probs = probs)
}
# 
# pct_res <- long_pct %>% 
#   group_by(M, effect_size, period, sd_prior_THETA) %>%
#   summarise(my_quantile(theta_percentile, seq(0, 1, .01)))
# 
#### QQ plot
# filter(pct_res, period == 9) %>% 
#   ggplot(aes(probs, x, color = M)) +
#   geom_line() +
#   facet_grid(cols = vars(effect_size),
#              rows = vars(sd_prior_THETA)) +
#   coord_equal() +
#   theme(legend.position = "top") +
#   labs(color = NULL, x = "observed", y = "theoretical", title = "Q-Q plot")

all_thetas <- unnest(res_dat, posterior_draws) %>% 
  select(-beta, -theta_percentiles)

summ <- all_thetas %>% 
  group_by(M, effect_size, period, sd_prior_THETA, outcome_id) %>%
  summarise_draws() %>% 
  filter(variable == "theta") %>% 
  ungroup()

combined_res <- summ %>%
  full_join(emp_res, by = c("M", "outcome_id", "period" = "periodId")) %>% 
  mutate(across(c(mean, median, q5, q95, logRr_redo_calibrated, 
                  logLb95Rr_redo_calibrated, logUb95Rr_redo_calibrated), exp, .names = "exp_{.col}"))

# table
# KEEP THIS
raw_res_tab <- combined_res %>%
  filter(period == 1, sd_prior_THETA == 4,
         outcome_id %in% c(10173, 10198, 10213, 378165, 10009)) %>%
  select(M, outcome_id, rr, ci95Lb, ci95Ub, exp_median, exp_q5, exp_q95, exp_logRr_redo_calibrated, exp_logLb95Rr_redo_calibrated, exp_logUb95Rr_redo_calibrated,
         calibratedRr, calibratedCi95Lb, calibratedCi95Ub) %>%
  rowwise() %>%
  mutate(min_dif = c(exp_median, exp_logRr_redo_calibrated, calibratedRr)[which.min(abs(log(c(exp_median, exp_logRr_redo_calibrated, calibratedRr)) - log(2)))],
         min_width_lower = c(exp_q5, exp_logLb95Rr_redo_calibrated, calibratedCi95Lb)[
           which.min(c(log(exp_q95) - log(exp_q5), log(exp_logUb95Rr_redo_calibrated) - log(exp_logLb95Rr_redo_calibrated), log(calibratedCi95Ub) - log(calibratedCi95Lb))
           )]) %>%
  gt(rowname_col = "M", groupname_col = "outcome_id") %>%
  cols_hide(starts_with("min_")) %>%
  fmt_number(columns = everything()) %>%
  tab_spanner(label = "Raw data",
              columns = c("rr", "ci95Lb", "ci95Ub")) %>%
  tab_spanner(label = "Bayesian approach",
              columns = c("exp_median", "exp_q5", "exp_q95")) %>%
  tab_spanner(label = "Empirical calibration w/ negative controls",
              columns = c("exp_logRr_redo_calibrated", "exp_logLb95Rr_redo_calibrated", 
                          "exp_logUb95Rr_redo_calibrated")) %>%
  tab_spanner(label = "Empirical calibration in Eumaeus",
              columns = c("calibratedRr", "calibratedCi95Lb", "calibratedCi95Ub")) %>%
  cols_label(rr = "RR",
             ci95Lb = "lower",
             ci95Ub = "upper",
             exp_median = "RR",
             exp_q5 = "lower",
             exp_q95 = "upper",
             exp_logRr_redo_calibrated = "RR",
             exp_logLb95Rr_redo_calibrated = "lower",
             exp_logUb95Rr_redo_calibrated = "upper",
             calibratedRr = "RR",
             calibratedCi95Lb = "lower",
             calibratedCi95Ub = "upper") %>%
  tab_style(style = list(cell_text(weight = "bold")),
            locations = cells_body(rows = M == "Overall")) %>%
  tab_style(style = list(cell_fill(alpha = 0.5)),
            locations = cells_body(columns = c("exp_median", "exp_q5", "exp_q95", 
                                               "calibratedRr", "calibratedCi95Lb", "calibratedCi95Ub"))) %>%
  tab_style(style = list(cell_fill(color = "#bebada", alpha = .75)),
            locations = cells_body(rows = exp_median == min_dif & M != "Overall",
                                   columns = "exp_median")) %>%
  tab_style(style = list(cell_fill(color = "#bebada", alpha = .5)),
            locations = cells_body(rows = exp_logRr_redo_calibrated == min_dif & M != "Overall",
                                   columns = "exp_logRr_redo_calibrated")) %>%
  tab_style(style = list(cell_fill(color = "#bebada", alpha = .75)),
            locations = cells_body(rows = calibratedRr == min_dif & M != "Overall",
                                   columns = "calibratedRr")) %>%
  tab_style(style = list(cell_fill(color = "#badabe", alpha = .75)),
            locations = cells_body(rows = exp_q5 == min_width_lower & M != "Overall",
                                   columns = c("exp_q5", "exp_q95"))) %>%
  tab_style(style = list(cell_fill(color = "#badabe", alpha = .5)),
            locations = cells_body(rows = exp_logLb95Rr_redo_calibrated == min_width_lower & M != "Overall",
                                   columns = c("exp_logLb95Rr_redo_calibrated", 
                                               "exp_logUb95Rr_redo_calibrated"))) %>%
  tab_style(style = list(cell_fill(color = "#badabe", alpha = .75)),
            locations = cells_body(rows = calibratedCi95Lb == min_width_lower & M != "Overall",
                                   columns = c("calibratedCi95Lb", "calibratedCi95Ub"))) %>%
  fmt_missing(columns = everything(), missing_text = "") %>%
  tab_header(title = "Point and interval estimates across methods",
             # subtitle = "With closest esimates and narrowest intervals highlighted"
             ) %>%
  cols_align(align = "center", columns = where(is.numeric)) %>%
  cols_width(where(is.numeric) ~ px(150))


outcome_samples_1 <- 378165
outcome_samples_1.5 <- 10711
outcome_samples_2 <- 10874
outcome_samples_4 <- 10044

# # when there are multiple sd's
# compare_posteriors_sd <- function(samples, all_thetas, title) {
#   filter(all_thetas, outcome_id %in% samples,
#          period == 9) %>% 
#     ggplot(aes(y = factor(sd_prior_THETA), x = exp(theta))) +
#     geom_vline(xintercept = 1, linetype = "dotted", color = "darkgrey") +
#     geom_vline(aes(xintercept = effect_size), linetype = "dashed") +
#     stat_halfeye() +
#     scale_x_log10(limits = c(.05, 20)) +
#     labs(x = "RR", y = NULL, title = title) +
#     facet_grid(cols = vars(M), rows = vars(outcome_id))
# }
# 
# # KEEP THIS (1, 4)
# compare_posteriors_plot1 <- compare_posteriors_sd(378165, all_thetas,
#                       "Posteriors for a negative control")
# 
# compare_posteriors_plot1.5 <- compare_posteriors_sd(10711, all_thetas,
#                       "Posteriors for positive control with RR = 1.5")
# 
# compare_posteriors_plot2 <- compare_posteriors_sd(10874, all_thetas,
#                       "Posteriors for positive control with RR = 2")
# 
# compare_posteriors_plot4 <- compare_posteriors_sd(10044, all_thetas,
#                       "Posteriors for positive control with RR = 4")

compare_posteriors_m <- function(samples, all_thetas, title) {
  filter(all_thetas, outcome_id %in% samples,
         period == 9) %>% 
    ggplot(aes(y = fct_rev(M), x = exp(theta))) +
    geom_vline(xintercept = 1, linetype = "dotted", color = "darkgrey") +
    geom_vline(aes(xintercept = effect_size), linetype = "dashed") +
    stat_halfeye() +
    scale_x_log10(limits = c(.05, 20)) +
    labs(x = "RR", y = NULL, title = title) +
    facet_grid(rows = vars(outcome_id))
}

# KEEP THIS (1, 4)
compare_posteriors_m_plot1 <- compare_posteriors_m(378165, all_thetas,
                                                  "Posteriors for a negative control")

compare_posteriors_m_plot1.5 <- compare_posteriors_m(10711, all_thetas,
                                                    "Posteriors for positive control with RR = 1.5")

compare_posteriors_m_plot2 <- compare_posteriors_m(10874, all_thetas,
                                                  "Posteriors for positive control with RR = 2")

compare_posteriors_m_plot4 <- compare_posteriors_m(10044, all_thetas,
                                                  "Posteriors for positive control with RR = 4")

# what kinds of questions can we ask?
# probability effect > 1

pr_gt <- filter(all_thetas, M == "Overall") %>% 
  group_by(M, effect_size, period, sd_prior_THETA, outcome_id) %>%
  summarise(pr_1 = mean(theta > 0),
            pr_1.5 = mean(theta > log(1.5)),
            pr_2 = mean(theta > log(2)),
            pr_4 = mean(theta > log(4)), .groups = "drop") %>% 
  pivot_longer(cols = starts_with("pr_"),
               names_to = c(".value", "gt"),
               names_pattern = "(.+)\\_(.+)$", 
               names_ptypes = list("pr" = numeric()))

# KEEP THIS
simple_over_time_plot <- pr_gt %>% 
  filter(outcome_id %in% outcome_samples_1,
         gt == "1") %>% 
  ggplot() +
  geom_line(aes(period, pr, col = factor(sd_prior_THETA))) +
  scale_y_continuous() +
  labs(color = "Prior standard deviation",
       x = "Month",
       y = "Probability RR > 1",
       title = "Posterior probability over time that that RR > 1")

compare_pr_gt <- function(samples, pr_gt, title) {
  pr_gt %>% 
    filter(outcome_id %in% samples) %>% 
    mutate(`prior sd` = sd_prior_THETA) %>% 
    ggplot() +
    geom_line(aes(period, pr, col = gt), size = 1) +
    scale_y_continuous(limits = c(0, 1)) +
    facet_grid(rows = vars(`prior sd`), switch = "y",
               labeller = label_both) +
    labs(color = "Pr(RR > x)",
         x = "Month",
         y = "Probability",
         title = title) +
    geom_hline(yintercept = c(.05, .5, .95), linetype = "dashed",
               color = "darkgrey")
}

# Keep this (1 and 4)
compare_pr_gt_plot1 <- compare_pr_gt(outcome_samples_1, pr_gt,
                      "Probability of effect over time\n(negative control)")

compare_pr_gt_plot1.5 <- compare_pr_gt(outcome_samples_1.5, pr_gt,
              "Probability of effect over time\n(positive control with RR = 1.5)")

compare_pr_gt_plot2 <- compare_pr_gt(outcome_samples_2, pr_gt,
              "Probability of effect over time\n(positive control with RR = 2)")

compare_pr_gt_plot4 <- compare_pr_gt(outcome_samples_4, pr_gt,
              "Probability of effect over time\n(positive control with RR = 4)")

summ_pr <- pr_gt %>% 
  filter(period == 9, sd_prior_THETA == 4) %>% 
  group_by(effect_size, gt) %>% 
  summarise(my_quantile(pr, c(.05, .5, .95)), .groups = "drop")
# 50%(probs) of outcomes with effect size effect_size said that there was at least
# an x% probability that the effect size was greater than gt

# 5%(probs) of outcomes with effect size effect_size said that there was at least
# an x% probability that the effect size was greater than gt

# summ_pr %>% 
#   pivot_wider(names_from = probs, values_from = x) %>% 
#   gt(groupname_col = "effect_size", rowname_col = "gt") %>% 
#   fmt_number(columns = everything()) 
# %>%
#   tab_spanner(label = "x% of outcomes had posterior probabilities greater than",
#               columns = c(`1`, `1.5`, `2`, `4`))


filter(all_thetas, M == "Overall") %>% 
  group_by(M, effect_size, period, sd_prior_THETA, outcome_id) %>%
  summarise(pr_1 = mean(theta > 0),
            pr_1.5 = mean(theta > log(1.5)),
            pr_2 = mean(theta > log(2)),
            pr_4 = mean(theta > log(4)), .groups = "drop") %>% 
  pivot_longer(cols = starts_with("pr_"),
               names_to = c(".value", "gt"),
               names_pattern = "(.+)\\_(.+)$", 
               names_ptypes = list("pr" = numeric()))

# 5% of posteriors
# pct_res %>% 
#   rowwise() %>% 
#   mutate(keep = any(near(probs, c(.05, .5, .95)))) %>% 
#   ungroup() %>% 
#   filter(keep, sd_prior_THETA == 4, period == 9, M == "Overall") %>% 
#   select(-keep, -M, -period, -sd_prior_THETA)

# what percentage of outcomes have > 50% probability that x > effect size
posterior_quantile <- function(theta, effect_size, probs) {
  tibble(gt = mean(exp(theta) > effect_size) < probs,
         probs = probs)
}

pct_outcomes_gt <- all_thetas %>% 
  filter(M == "Overall", sd_prior_THETA == 4, period == 9) %>% 
  group_by(outcome_id, effect_size) %>% 
  summarise(posterior_quantile(theta, effect_size, probs = seq(0, 1, .01))) %>% 
  group_by(probs, effect_size) %>% 
  summarise(pct_outcomes = mean(gt), .groups = "drop")
# 100% (pct_outcomes) of outcomes have > 1% (probs) probability that effect size > 1

# y% of outcomes have >x% posterior probabiliy that effect size > true effect size
# KEEP THIS
pct_probability_gt_plot <- pct_outcomes_gt %>% 
  mutate(`true RR` = effect_size) %>% 
  ggplot() +
  geom_step(aes(probs, pct_outcomes)) +
  facet_wrap(vars(`true RR`), labeller = label_both, scales = "free") +
  labs(x = NULL, y = NULL, 
       title = "y% of outcomes had < x% posterior probability that RR > true RR",
       subtitle = "(prior sd = 4, final time point)") +
  coord_cartesian(expand = FALSE)

# 
# combined_res %>%
#   filter(period == 1, sd_prior_THETA == 1,
#          outcome_id %in% c(10173, 10198, 10213, 378165, 10009)) %>% 
# ggplot(filter(compare_bayes, M == "Overall"),
#        aes(y = factor(outcome_id), x = exp(theta), color = method)) +
#   geom_point(position = position_dodge()) +
#   geom_errorbarh(aes(xmin = exp_q5, xmax = exp_q95), position = position_dodge()) +
#   geom_vline(aes(xintercept = effect_size), linetype = "dashed") +
#   scale_x_log10(limits = c(.1, 40)) +
#   labs(x = "RR", y = NULL)

compare_res <- combined_res %>% 
  filter(M != "Overall", sd_prior_THETA == 4, period == 9) %>% 
  select(M, outcome_id, effect_size, RR_B = exp_median, RR_E = exp_logRr_redo_calibrated, 
         RR_P = calibratedRr,
         low_B = exp_q5, hi_B = exp_q95,
         low_E = exp_logLb95Rr_redo_calibrated, 
         low_P = calibratedCi95Lb, hi_P = calibratedCi95Ub,
         hi_E = exp_logUb95Rr_redo_calibrated) %>% 
  pivot_longer(c(ends_with("_B"), ends_with("_E"), ends_with("_P")),
               names_sep = "\\_",
               names_to = c(".value", "method")) %>% 
  mutate(method = fct_recode(method, "Bayesian meta-analysis" = "B", 
                             "Empirical calibration" = "E",
                             "Eumaeus calibration" = "P")) %>% 
  rowwise() %>% 
  mutate(in_int = between(effect_size, low, hi)) %>% 
  ungroup() %>% 
  group_by(M, effect_size, method) %>% 
  mutate(o = row_number()) %>% 
  ungroup() 

# KEEP THIS
compare_calibration_tab <- compare_res %>% 
  group_by(M, method, effect_size) %>% 
  summarise(median_RR = median(RR),
            in_int = mean(in_int),
            int_width = median(log(hi) - log(low)), .groups = "drop") %>% 
  pivot_wider(names_from = method, values_from = c(median_RR, in_int, int_width)) %>% 
  gt(groupname_col = "effect_size", rowname_col = "M") %>% 
  fmt_number(where(is.numeric)) %>% 
  fmt_percent(starts_with("in_int"), decimals = 0) %>% 
  fmt_missing(columns = everything(), missing_text = "") %>% 
  tab_spanner(label = "Median RR",
              columns = starts_with("median_RR")) %>% 
  tab_spanner(label = "Interval coverage",
              columns = starts_with("in_int")) %>% 
  tab_spanner(label = "Interval width (log scale)",
              columns = starts_with("int_width")) %>% 
  cols_label(`median_RR_Bayesian meta-analysis` = "Bayesian",
             `median_RR_Empirical calibration` = "Empirical",
             `median_RR_Eumaeus calibration` = "Eumaeus",
             `in_int_Bayesian meta-analysis` = "Bayesian",
             `in_int_Empirical calibration` = "Empirical",
             `in_int_Eumaeus calibration` = "Eumaeus",
             `int_width_Bayesian meta-analysis` = "Bayesian",
             `int_width_Empirical calibration` = "Empirical",
             `int_width_Eumaeus calibration` = "Eumaeus") %>% 
  # tab_header(title = "Comparison of methods across negative and positive controls") %>% 
  cols_align(align = "center", columns = where(is.numeric)) %>% 
  cols_width(where(is.numeric) ~ px(150))

# KEEP THIS
compare_calibration_plot <- compare_res %>% 
  filter(effect_size == 1) %>% 
  filter(outcome_id %in% c(10260, 74194, 443421, 11199, 10929, 10819, 10162, 4032787, 
                           138690, 10441, 4311499, 11360, 10019, 10205, 10160, 10152, 10030, 
                           10111, 10112, 11648, 10161, 10909, 10459, 10151, 11053, 10964, 
                           11017, 10802, 11631, 10105, 10203, 192953, 11413, 10928, 
                           10131, 377877, 4239873, 10108, 10334, 10047, 196347, 765053, 
                           10056, 10245, 10839, 10695, 11163, 10264)) %>% 
  ggplot(aes(y = fct_rev(factor(outcome_id)), color = method, group = method)) +
  geom_point(aes(x = RR), position = position_dodge(width = .5)) +
  geom_errorbarh(aes(xmin = low, xmax = hi), position = position_dodge(width = .5)) +
  geom_vline(aes(xintercept = effect_size), linetype = "dashed") +
  scale_x_log10() +
  facet_grid(cols = vars(M)) +
  labs(y = NULL, color = NULL, title = "Comparison of estimates for negative controls at final timepoint across methods",
       subtitle = "(for Bayesian, prior sd = 4)") +
  theme(axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        legend.position = "top")

# compare_res %>% 
#   pivot_wider(names_from = method, values_from = RR, id_cols = c(M, o, outcome_id, effect_size)) %>% 
#   mutate(dif = log(`Bayesian meta-analysis`) - log(`Empirical calibration`)) %>% 
#   group_by(effect_size, M) %>% 
#   summarise(mean(dif))

MSE_plot <- compare_res %>% 
  group_by(method, effect_size, M) %>% 
  summarise(MSE = mean((log(RR) - log(effect_size))^2)) %>% 
  filter(!(method == "Eumaeus calibration" & M == "OptumEhr")) %>%
  ggplot() +
  aes(y = MSE, x = M, color = method) +
  geom_point() +
  facet_wrap(vars(effect_size), labeller = label_both) +
  labs(title = "Mean squared error", 
       subtitle = "Excluding Optum EHR results from Eumaeus",
       x = NULL, color = NULL) +
  geom_hline(aes(yintercept = 0))

# Let's call it positive if the interval excludes 1
finds <- compare_res %>% 
  mutate(find = !(low < 1 & hi > 1)) %>%
  group_by(M, method, effect_size) %>% 
  summarise(finds = mean(find), .groups = "drop")
  
type1_tab <- finds %>% 
  filter(effect_size == 1) %>% 
  pivot_wider(names_from = method, values_from = finds) %>% 
  gt() %>% 
  cols_hide(effect_size) %>% 
  cols_label(M = "site") %>% 
  tab_header(title = "'Type 1 error'", subtitle = "(negative control interval does not contain 1)") %>% 
  cols_align(align = "center", columns = where(is.numeric)) %>% 
  fmt_percent(columns = where(is.numeric), decimals = 1)

power_tab <- finds %>% 
  filter(effect_size != 1) %>% 
  pivot_wider(names_from = method, values_from = finds) %>% 
  gt(groupname_col = "effect_size") %>% 
  cols_label(M = "site") %>% 
  tab_header(title = "'Power' for various effect sizes", subtitle = "(interval does not contain 1)") %>% 
  cols_align(align = "center", columns = where(is.numeric)) %>% 
  fmt_percent(columns = where(is.numeric), decimals = 1)

# KEEP THIS
evaluated_tab <- filter(combined_res,
       sd_prior_THETA == 4, period == 9) %>% 
  distinct(outcome_id, effect_size) %>% 
  count(effect_size) %>% 
  mutate(`SD of prior` = "1, 4, 10",
         Months = "1-9", 
         nn = str_glue("**{n}** x 3 x 10")) %>% 
  gt() %>% 
  cols_hide(n) %>% 
  cols_label(effect_size = "RR",
             nn = "N") %>% 
  tab_header(title = "True RRs evaluated") %>% 
  fmt_markdown(columns = nn) %>% 
  cols_align(align = "center", columns = where(is.numeric)) 

sites_per_control_tab <- filter(combined_res, M != "Overall",
       sd_prior_THETA == 4, period == 9) %>% 
  count(outcome_id) %>% 
  count(n) %>% 
  gt() %>% 
  cols_label(n = "Sites",
             nn = "Controls") %>% 
  tab_header(title = "Sites per control") %>% 
  cols_align(align = "center", columns = where(is.numeric)) 

neg_control_tab <- filter(dat, is.na(effectSize)) %>% 
  group_by(databaseId, outcomeId) %>% 
  slice_head(n = 1) %>% 
  ungroup() %>% 
  count(databaseId) %>% 
  gt() %>% 
  cols_label(databaseId = "Site",
             n = "Negative controls") %>% 
  tab_header(title = "Negative controls used in calibration") %>% 
  cols_align(align = "center", columns = where(is.numeric)) 

to_keep <- lst(raw_res_tab, compare_calibration_tab, evaluated_tab,
               sites_per_control_tab, neg_control_tab, power_tab, type1_tab)

ggsave(here::here("slides", "compare_pr_gt_plot1.pdf"), 
       compare_pr_gt_plot1, height = 4*1.5, width = 6*1.5)
ggsave(here::here("slides", "compare_pr_gt_plot4.pdf"), 
       compare_pr_gt_plot4, height = 4*1.5, width = 6*1.5)
ggsave(here::here("slides", "pct_probability_gt_plot.pdf"), 
       pct_probability_gt_plot, height = 4*1.5, width = 6*1.5)
ggsave(here::here("slides", "compare_posteriors_m_plot1.pdf"), 
       compare_posteriors_m_plot1, height = 4*1.5, width = 6*1.5)
ggsave(here::here("slides", "compare_posteriors_m_plot4.pdf"), 
       compare_posteriors_m_plot4, height = 4*1.5, width = 6*1.5)
ggsave(here::here("slides", "simple_over_time_plot.pdf"), 
       simple_over_time_plot, height = 4*1.5, width = 6*1.5)
ggsave(here::here("slides", "compare_calibration_plot.pdf"), 
       compare_calibration_plot, height = 4*1.5, width = 6*1.5)
ggsave(here::here("slides", "MSE_plot.pdf"), 
       MSE_plot, height = 4*1.5, width = 6*1.5)
write_rds(
  to_keep,
  file = here::here("slides", "tabs.rds")
)
