library(tidyverse)
library(tidybayes)
theme_set(theme_tidybayes())

eumaeus_file <- "eumaeus_HistoricalComparator_1_21215.rds"

dat <- read_rds(here::here("data", eumaeus_file))

all_res <- read_rds(here::here("results", paste0("res_", eumaeus_file)))
res <- all_res[-which(sapply(all_res, is.null))]
res_t <- transpose(res)

all_res_hist <- read_rds(here::here("results", paste0("historical_", eumaeus_file)))
res_hist <- all_res_hist[-which(sapply(all_res_hist, is.null))]
res_hist_t <- transpose(res_hist)

all_res_emp <- read_rds(here::here("results", paste0("emp_", eumaeus_file)))
res_emp <- all_res_emp[-which(sapply(all_res_emp, is.null))]


# clean up results
get_theta_res <- function(res_t, ...) {
  all_draws <-  res_t$posterior_draws %>% 
    set_names(str_glue("{res_t$outcome_id}_{res_t$effect_size}")) %>% 
    bind_rows(.id = "outcome") %>% 
    separate(outcome, into = c("outcome", "effect_size"), sep = "\\_", convert = TRUE)
  
  thetas <- select(all_draws, -beta)
  betas <- select(all_draws, -theta)
  
  summ <- thetas %>% 
    group_by(M, outcome) %>%
    summarise_draws() %>% 
    filter(variable == "theta") %>% 
    ungroup()
  
  hdis <- thetas %>% 
    group_by(M, outcome) %>%
    median_hdci(theta) %>% 
    ungroup()
  
  theta_res <- full_join(summ, hdis, by = c("M", "outcome")) %>%
    left_join(bind_rows(res_emp), by = c("M", "outcome")) %>% 
    left_join(distinct(thetas, outcome, M, effect_size), by = c("M", "outcome")) %>% 
    mutate(across(c(theta, .lower, .upper, mean, median, q5, q95, logRr, logLb95Rr, logUb95Rr), exp, .names = "exp_{.col}"))
  
  theta_res
}

theta_res <- get_theta_res(res_t)
theta_res_hist <- get_theta_res(res_hist_t)

compare_bayes <- mutate(theta_res, method = "likelihood grid") %>% 
  bind_rows(mutate(theta_res_hist, method = "poisson likelihood"))

ggplot(filter(compare_bayes, M == "Overall"),
       aes(y = factor(outcome), x = exp(theta), color = method)) +
  geom_point(position = position_dodge()) +
  geom_errorbarh(aes(xmin = exp_q5, xmax = exp_q95), position = position_dodge()) +
  geom_vline(aes(xintercept = effect_size), linetype = "dashed") +
  scale_x_log10(limits = c(.1, 40)) +
  labs(x = "RR", y = NULL)

compare_bayes %>% 
  mutate(int_len = q95 - q5) %>% 
  group_by(method, M) %>% 
  summarise(mean(int_len), median(int_len)) %>% 
  pivot_wider(names_from = method, values_from = -c(method, M))

compare_res <- theta_res %>% 
  filter(M != "Overall") %>% 
  select(M, outcome, effect_size, RR_B = exp_median, RR_E = exp_logRr, 
         low_B = exp_.lower, hi_B = exp_.upper,
         low_E = exp_logLb95Rr, hi_E = exp_logUb95Rr) %>% 
  pivot_longer(c(ends_with("_B"), ends_with("_E")),
               names_sep = "\\_",
               names_to = c(".value", "method")) %>% 
  mutate(method = fct_recode(method, "Bayesian meta-analysis" = "B", 
                             "Empirical calibration" = "E")) %>% 
  rowwise() %>% 
  mutate(in_int = between(effect_size, low, hi)) %>% 
  ungroup()

compare_res %>% 
  group_by(M, method) %>% 
  summarise(median_RR = median(RR),
            in_int = mean(in_int),
            int_width = median(log(hi) - log(low)))

theta_res %>% 
  filter(M == "Overall") %>% 
  summarise(median(exp_theta))

theta_res_hist %>% 
  filter(M == "Overall") %>% 
  summarise(median(exp_theta))

compare_res %>% 
  group_by(M, effect_size) %>% 
  mutate(o = row_number()) %>% 
  ungroup() %>% 
  ggplot() +
  geom_point(aes(RR, factor(o), color = method), position = position_dodge()) +
  geom_errorbarh(aes(xmin = low, xmax = hi, y = factor(o), color = method), position = position_dodge()) +
  geom_vline(aes(xintercept = effect_size), linetype = "dashed") +
  scale_x_log10() +
  facet_grid(cols = vars(M), scales = "free", space = "free") +
  labs(y = NULL) +
  theme(axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        legend.position = "top")


ggplot(filter(thetas, M == "Overall"),
       aes(y = factor(outcome), x = exp(theta))) +
  stat_halfeye() +
  geom_vline(aes(xintercept = effect_size), linetype = "dashed") +
  scale_x_log10(limits = c(.1, 40)) +
  labs(x = "RR", y = NULL)

halfeye_plot_outcome <- function(thetas, ncol = 10, ...) {
  ggplot(thetas,
         aes(y = fct_rev(M), x = exp(theta))) +
    stat_halfeye() +
    geom_vline(aes(xintercept = effect_size), linetype = "dashed") +
    scale_x_log10(limits = c(.05, 20)) +
    labs(x = "RR", y = NULL) +
    facet_wrap(vars(outcome), ncol = ncol)
}

# # when there are multiple effect sizes, plot separately
# plots_by_outcome <- thetas %>% 
#   mutate(e = effect_size) %>% 
#   group_by(e) %>% 
#   group_map(halfeye_plot_outcome)

halfeye_plot_outcome(thetas, ncol = 12)

halfeye_plot_outcome(filter(thetas, outcome < 10055), ncol = 4)

# ggplot(all_draws, aes(x = exp(theta), y = fct_rev(M))) +
#   stat_dotsinterval(quantiles = 100) +
#   scale_x_log10() +
#   labs(x = "RR", y = NULL) +
#   facet_wrap(vars(outcome))


ggplot(filter(thetas, outcome == 10010),
       aes(y = fct_rev(M), x = exp(theta), fill = stat(x > 0))) +
  stat_halfeye() +
  geom_vline(xintercept = 1, linetype = "dashed") +
  scale_fill_manual(values = c("gray80", "skyblue"))  +
  scale_x_log10() +
  theme(legend.position = "none") +
  labs(x = "RR", y = NULL)

filter(thetas, M == "Overall") %>% 
  group_by(outcome) %>% 
  summarise(mean(theta > 0))



# res_t$theta_percentiles %>% 
#   set_names(str_glue("{res_t$outcome_id}_{res_t$effect_size}")) %>% 
#   bind_rows(.id = "outcome") %>% 
#   group_by(M) %>% 
#   summarise(median(theta_percentile))