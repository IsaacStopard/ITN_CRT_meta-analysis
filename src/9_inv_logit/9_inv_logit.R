library(cmdstanr)
library(tidyverse)
library(posterior)

orderly2::orderly_dependency("1_data_cleaning", "latest()", c("stan_data.rds"))

orderly2::orderly_shared_resource("m2.stan", "model_functions.R")

source("model_functions.R")

list2env(readRDS("stan_data.rds"), envir = .GlobalEnv)

orderly2::orderly_dependency("2_fits",
                             "latest()",
                             c("fit_2.rds")
                             )

fit_full <- readRDS(file = "fit_2.rds")

stan_model <- cmdstan_model(stan_file = "m2.stan")

data_limits <- COMBO_stan |> group_by(study, i) |> summarise(min_prev = min(BL_prev), max_prev = max(BL_prev))

target_params <- c("alpha", "alpha_i", "theta_l", "theta_li", "kappa_l", "kappa_li", "omega_l", "delta_l", "alpha_i",
                     "gamma", "kappa", "kappa_i", "omega", "delta", "e_ij", "sigma_e_r")

params <- extract_params(target_params, fit_full)

probs_in <- c(lower, 0.5, upper)

##########################################################
##### within-study difference in baseline prevalence #####
##########################################################

u_ij <- unique(COMBO_stan[,c("l", "li", "i", "study", "net", "cluster", "BL_prev")])

# log odds ratios
u_gq <- unique(COMBO_stan[,c("l", "li", "i", "study", "net")]) |> left_join(u_ij |> summarise(base_prev_mean = mean(BL_prev), .by = c("i")), by = c("i")) |> 
  mutate(pool_base_prev_mean = mean(u_ij$BL_prev))

diff_prev <- seq(-1, 1, 0.025)

u_gq_diff <- u_gq[rep(1:nrow(u_gq), length(diff_prev)), ] |> mutate(base_prev_diff = rep(diff_prev, each = nrow(u_gq))) |> 
  mutate(start_prev = base_prev_diff + base_prev_mean) |> 
  filter(start_prev >= 0 & start_prev <= 1)

u_gq_diff <- u_gq_diff |> left_join(data_limits) |> filter(start_prev >= min_prev & start_prev <= max_prev) |> select(-min_prev, -max_prev)

u_gq_diff_pool <- u_gq[rep(1:nrow(u_gq), length(diff_prev)), ] |> mutate(base_prev_diff = rep(diff_prev, each = nrow(u_gq))) |> 
  mutate(start_prev = base_prev_diff + pool_base_prev_mean) |> 
  filter(start_prev >= 0 & start_prev <= 1) |> filter(start_prev >= min(data_limits$min_prev) & start_prev <= max(data_limits$max_prev))

out_o_r_diff <- mean_time_l_or_study_fun(3, data = u_gq_diff, params, 2, TRUE) |> 
  cbind(u_gq_diff) 

out_pool_o_r_diff <- mean_time_l_or_pooled_fun(3, data = u_gq_diff_pool |> mutate(base_prev_mean = pool_base_prev_mean), params, 2, TRUE) |> 
  cbind(u_gq_diff_pool)

# relative efficacy

out_r_eff_diff <- mean_time_r_eff_study_fun(3, data = u_gq_diff, params, 2, TRUE) |> 
  cbind(u_gq_diff)

####################################
##### mean baseline prevalence #####
####################################

mean_base_prev_limits <- u_ij |> summarise(base_prev_mean = mean(BL_prev),
                                           pbo = ifelse(2 %in% unique(l), 1, 0),
                                           pp = ifelse(3 %in% unique(l), 1, 0),
                                           .by = c("i", "study")) |> 
  pivot_longer(cols = starts_with("p"), names_to = "net", values_to = "net_in") |> 
  filter(net_in == 1) |>
  summarise(
    min_bp = min(base_prev_mean),
    max_bp = max(base_prev_mean),
    .by = c("net")
  ) |> mutate(l = c(2, 3)) |> as.data.frame()

mean_prev <- c(seq(mean_base_prev_limits[1, "min_bp"], mean_base_prev_limits[1, "max_bp"], length.out = 45),
               seq(mean_base_prev_limits[2, "min_bp"], mean_base_prev_limits[2, "max_bp"], length.out = 15)) |> sort()

u_gq_mean <- unique(COMBO_stan[,c("l", "li", "i", "study", "net")]) |> mutate(base_prev_diff = 0)

u_gq_mean <- u_gq_mean[rep(1:nrow(u_gq_mean), length(mean_prev)), ] |> mutate(base_prev_mean = rep(mean_prev, each = nrow(u_gq)),
                                                                              base_prev_diff = 0) |> 
  left_join(mean_base_prev_limits |> select(-net), by = c("l")) |> 
  filter(base_prev_mean >= min_bp & base_prev_mean <= max_bp) |> 
  select(-min_bp, -max_bp)
  
out_o_r_mean <- mean_time_l_or_study_fun(3, data = u_gq_mean, params, 2, TRUE) |> cbind(u_gq_mean) 

out_pool_o_r_mean <- mean_time_l_or_pooled_fun(3, data = u_gq_mean, params, 2, TRUE) |> cbind(u_gq_mean)

# relative efficacy

out_r_eff_mean <- mean_time_r_eff_study_fun(3, data = u_gq_mean, params, 2, TRUE) |> cbind(u_gq_mean)

#############################################
##### mean prevalence with sample times #####
#############################################

u_gq_BL_prev <- unique(COMBO_stan[, c("i", "cluster", "study", "net", "l", "li", "years", "BL_prev")]) |> 
  left_join(u_ij |> summarise(base_prev_mean = mean(BL_prev), .by = c("i")))

u_gq_BL_prev <- u_gq_BL_prev[rep(1:nrow(u_gq_BL_prev), length(diff_prev)), ] |> 
  mutate(base_prev_diff = rep(diff_prev, each = nrow(u_gq_BL_prev))) |> 
  mutate(start_prev = base_prev_mean + base_prev_diff) |> 
  filter(start_prev >= 0 & start_prev <= 1)

u_gq_BL_prev_lit <- unique(u_gq_BL_prev[,c("start_prev", "base_prev_diff", "base_prev_mean", "net", "study", "i", "l", "li")]) |> mutate(lit = row_number())

u_gq_BL_prev <- u_gq_BL_prev |> left_join(
  u_gq_BL_prev_lit
  )

list2env(params$params_list, envir = environment())

set.seed(123)
e_ij_single <- matrix(rnorm(length(sigma_e_r), 0, 1), nrow = nrow(sigma_e_r)) * sigma_e_r

fit_full$expose_functions()

n_iter <- length(alpha)

logit_prob_t <- sapply(1:n_iter, function(i){
    fit_full$functions$logit_prob_fun(nrow(u_gq_BL_prev), u_gq_BL_prev$i, u_gq_BL_prev$l, u_gq_BL_prev$li, u_gq_BL_prev$i,
                                   alpha[i], alpha_i[i, ], theta_l[i, ], theta_li[i, ], u_gq_BL_prev$base_prev_mean, u_gq_BL_prev$base_prev_diff, 
                                   gamma[i], kappa[i], kappa_i[i, ], kappa_l[i, ], kappa_li[i, ], u_gq_BL_prev$years, omega[i], delta[i], omega_l[i, ], delta_l[i, ], 
                                   e_ij_single[i, ])}
  )

prob_t_mean <- rowsum(plogis(logit_prob_t), group = u_gq_BL_prev$lit) / as.numeric(table(u_gq_BL_prev$lit))
prob_t_mean <- prob_t_mean |> apply(1, quantile, probs = c(lower, 0.5, upper)) |> t() |> 
  as.data.frame() |> 
  rename("l_p" = 1, "m_p" = 2, "u_p" = 3) |> 
  bind_cols(u_gq_BL_prev_lit)

### saving all results
saveRDS(list("or_mean_3y_mean" = out_o_r_mean,
             "or_mean_3y_mean_pool" = out_pool_o_r_mean,
             "eff_mean_3y_mean" = out_r_eff_mean,
             "or_mean_3y_diff" = out_o_r_diff,
             "or_mean_3y_diff_pool" = out_pool_o_r_diff,
             "eff_mean_3y_diff" = out_r_eff_diff,
             "inv_log_BL_prev" = prob_t_mean),
        file = "m2_pred.rds")

