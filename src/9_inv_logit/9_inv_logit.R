library(rstan)
library(tidyverse)

orderly2::orderly_dependency("1_data_cleaning", "latest()", c("stan_data.rds"))

orderly2::orderly_parameters(model = NULL)

orderly2::orderly_shared_resource("m2_re.stan")

list2env(readRDS("stan_data.rds"), envir = .GlobalEnv)

orderly2::orderly_dependency("2_fits",
                             "latest(parameter:model == this:model)",
                             paste0(model, "_fit_full.rds")
                             )

fit_full <- readRDS(file = paste0(model, "_fit_full.rds"))

stan_model <- rstan::stan_model(file = paste0(model, ".stan"))

u_gq <- unique(COMBO_stan[,c("l", "li", "i", "net", "study")]) |>
  dplyr::mutate(ij = i,
                rn = dplyr::row_number())

# ##############################
# ##### model predictions ######
# ##############################

# different combinations of predictors
ucs <- unique(COMBO_stan[, c("cluster", "study", "i", "BL_prev_num", "BL_prev_denom", "BL_prev")])

ucs <- ucs |> left_join(ucs |> group_by(study) |>
               summarise(mean_BL_prev = mean(BL_prev))) |>
  mutate(diff_BL_prev = BL_prev - mean_BL_prev)

mean_prev_pooled <- mean(ucs[,"BL_prev"])

mean_BL_prev_combs <- unique(ucs$mean_BL_prev)

data_limits <- COMBO_stan |> group_by(study, i) |> summarise(min_prev = min(BL_prev), max_prev = max(BL_prev)) |>
  left_join(unique(ucs[,c("i", "mean_BL_prev")])) |> as.data.frame()

data_limits <- rbind(data_limits, data_limits |> mutate(mean_BL_prev = mean_prev_pooled,
                                                        min_prev = 0,
                                                        max_prev = 1))

diff_prev <- seq(-1, 1, 0.01)
years <- c(1, 2, 3)

bp_gq <- expand.grid("mean_prev" = mean_BL_prev_combs,
                     "diff_prev" = diff_prev,
                     "years" = years) |>
  mutate(start_prev = diff_prev + mean_prev,
         mean_prev_pooled = mean_prev_pooled,
         diff_prev_pooled = diff_prev) |>
  filter(start_prev <= 1 & start_prev >= 0)

u_gq_m <- u_gq[rep(1:nrow(u_gq), nrow(bp_gq)),] |>
  dplyr::mutate(mean_prev = rep(bp_gq[, "mean_prev"], each = nrow(u_gq)),
                start_prev = rep(bp_gq[, "start_prev"], each = nrow(u_gq)),
                diff_prev = rep(bp_gq[, "diff_prev"], each = nrow(u_gq)),
                years = rep(bp_gq[, "years"], each = nrow(u_gq)),
                mean_prev_pooled = rep(bp_gq[, "mean_prev_pooled"], each = nrow(u_gq)),
                diff_prev_pooled = rep(bp_gq[, "diff_prev_pooled"], each = nrow(u_gq))
                )

u_gq_diff_prev <- lapply(1:nrow(data_limits), FUN = function(i, u_gq_m, data_limits){

  u_gq_m |> dplyr::filter(study == !!data_limits[i, "study"],
                          start_prev >= !!data_limits[i, "min_prev"] &
                          start_prev <= !!data_limits[i, "max_prev"] &
                          round(mean_prev, digits = 3) == round(!!data_limits[i, "mean_BL_prev"], digits = 3))

  }, u_gq_m = u_gq_m, data_limits = data_limits) |> bind_rows()

mean_gq <- expand.grid("mean_prev" =
                         sort(unique(
                           c(seq(min(data_limits$mean_BL_prev), max(data_limits$mean_BL_prev), length.out = 10),
                             seq(subset(data_limits, i == 3 & min_prev != 0)$mean_BL_prev, subset(data_limits, i == 4 & min_prev != 0)$mean_BL_prev, length.out = 10))
                                  )
                              ),
                       "years" = years) |>
  mutate(diff_prev = 0,
         start_prev = mean_prev,
         mean_prev_pooled = mean_prev,
         diff_prev_pooled = diff_prev)

u_gq_mean_prev <- u_gq[rep(1:nrow(u_gq), nrow(mean_gq)),] |>
  dplyr::mutate(mean_prev = rep(mean_gq[, "mean_prev"], each = nrow(u_gq)),
                start_prev = rep(mean_gq[, "start_prev"], each = nrow(u_gq)),
                diff_prev = rep(mean_gq[, "diff_prev"], each = nrow(u_gq)),
                years = rep(mean_gq[, "years"], each = nrow(u_gq)),
                mean_prev_pooled = rep(mean_gq[, "mean_prev_pooled"], each = nrow(u_gq)),
                diff_prev_pooled = rep(mean_gq[, "diff_prev_pooled"], each = nrow(u_gq))
                )

pred_data_mean_prev <- data_in_full |>
  purrr::list_assign(gq = 1,
                     N_gq = nrow(u_gq_mean_prev),
                     N_ij_gq = length(unique(u_gq_mean_prev$ij)), # the same random effect for all the clusters in each study
                     N_ij_unq_gq = length(unique(u_gq_mean_prev$i)),
                     pmat_ij_gq = fastDummies::dummy_cols(u_gq_mean_prev$ij)[,-1] |> as.matrix(),
                     r_id_gq = seq(1, length(unique(u_gq_mean_prev$i))),
                     pmat_i_gq = fastDummies::dummy_cols(u_gq_mean_prev$i)[,-1] |> as.matrix(),
                     pmat_l_gq = fastDummies::dummy_cols(u_gq_mean_prev$l)[,-c(1, 2)] |> as.matrix(),
                     pmat_li_gq = fastDummies::dummy_cols(u_gq_mean_prev$li)[,-c(1, 2)] |> as.matrix(),

                     pmat_i_pooled_gq = matrix(0, nrow = nrow(u_gq_mean_prev), ncol = length(unique(u_gq_m$i))),
                     pmat_li_pooled_gq = matrix(0, nrow = nrow(u_gq_mean_prev), ncol = length(unique(u_gq_m$li)) - 1),
                     pmat_it_pooled_gq = matrix(0, nrow = nrow(u_gq_mean_prev), ncol = length(unique(u_gq_m$it))),

                     ij_train_gq = rep(0, length(unique(u_gq_mean_prev$ij))),
                     ij_unq_train_gq = rep(1, length(unique(u_gq_mean_prev$i))),
                     years_gq = u_gq_mean_prev$years,
                     m_prob_s_gq = u_gq_mean_prev$mean_prev,
                     d_prob_s_gq = u_gq_mean_prev$diff_prev,

                     m_prob_s_gq_pooled = u_gq_mean_prev$mean_prev_pooled,
                     d_prob_s_gq_pooled = u_gq_mean_prev$diff_prev_pooled)

pred_data_diff_prev <- data_in_full |>
  purrr::list_assign(gq = 1,
                     N_gq = nrow(u_gq_diff_prev),
                     N_ij_gq = length(unique(u_gq_diff_prev$ij)), # the same random effect for all the clusters in each study
                     N_ij_unq_gq = length(unique(u_gq_diff_prev$i)),
                     pmat_ij_gq = fastDummies::dummy_cols(u_gq_diff_prev$ij)[,-1] |> as.matrix(),
                     r_id_gq = seq(1, length(unique(u_gq_diff_prev$i))),
                     pmat_i_gq = fastDummies::dummy_cols(u_gq_diff_prev$i)[,-1] |> as.matrix(),
                     pmat_l_gq = fastDummies::dummy_cols(u_gq_diff_prev$l)[,-c(1, 2)] |> as.matrix(),
                     pmat_li_gq = fastDummies::dummy_cols(u_gq_diff_prev$li)[,-c(1, 2)] |> as.matrix(),

                     pmat_i_pooled_gq = matrix(0, nrow = nrow(u_gq_diff_prev), ncol = length(unique(u_gq_m$i))),
                     pmat_li_pooled_gq = matrix(0, nrow = nrow(u_gq_diff_prev), ncol = length(unique(u_gq_m$li)) - 1),
                     pmat_it_pooled_gq = matrix(0, nrow = nrow(u_gq_diff_prev), ncol = length(unique(u_gq_m$it))),

                     ij_train_gq = rep(0, length(unique(u_gq_diff_prev$ij))),
                     ij_unq_train_gq = rep(1, length(unique(u_gq_diff_prev$i))),
                     years_gq = u_gq_diff_prev$years,
                     m_prob_s_gq = u_gq_diff_prev$mean_prev,
                     d_prob_s_gq = u_gq_diff_prev$diff_prev,

                     m_prob_s_gq_pooled = u_gq_diff_prev$mean_prev_pooled,
                     d_prob_s_gq_pooled = u_gq_diff_prev$diff_prev_pooled)

process_gq <- function(pred_data,
                       u_gq_in){

  model_predict <- rstan::gqs(stan_model,
                              draws = as.matrix(fit_full),
                              data = pred_data,
                              seed = 123)

  inv_log_all <- rstan::extract(model_predict, "inv_logit_posterior")$inv_logit_posterior

  inv_log_all_pooled <- rstan::extract(model_predict, "inv_logit_posterior_pooled")$inv_logit_posterior_pooled

  inv_log <- inv_log_all |> apply(2, quantile, probs = c(lower, 0.5, upper)) |> t() |> cbind(u_gq_in) |> dplyr::rename(l_p = 1, m_p = 2, u_p = 3)

  inv_log_pooled <- inv_log_all_pooled |> apply(2, quantile, probs = c(lower, 0.5, upper)) |> t() |> cbind(u_gq_in) |> dplyr::rename(l_p_pooled = 1, m_p_pooled = 2, u_p_pooled = 3)

  o_r <- rstan::extract(model_predict, "o_r")$o_r |> apply(2, quantile, probs = c(lower, 0.5, upper)) |> t() |> cbind(u_gq_in) |> dplyr::rename(l_or = 1, m_or = 2, u_or = 3)

  o_r_pooled <- rstan::extract(model_predict, "o_r_pooled")$o_r_pooled |> apply(2, quantile, probs = c(lower, 0.5, upper)) |> t() |>
    cbind(u_gq_in) |>
    dplyr::rename(l_or_pooled = 1, m_or_pooled = 2, u_or_pooled = 3)

  out <- inv_log |> left_join(inv_log_pooled) |> left_join(o_r) |> left_join(o_r_pooled)

  rm(list = c("model_predict", "inv_log", "inv_log_pooled", "o_r", "o_r_pooled"))

  pred_data_p_only <- pred_data
  pred_data_p_only$pmat_l_gq[pred_data_p_only$pmat_l_gq == 1] <- 0
  pred_data_p_only$pmat_li_gq[pred_data_p_only$pmat_li_gq == 1] <- 0

  model_predict_p_only <- rstan::gqs(stan_model,
                                     draws = as.matrix(fit_full),
                                     data = pred_data_p_only,
                                     seed = 123)

  inv_log_p <- rstan::extract(model_predict_p_only, "inv_logit_posterior")$inv_logit_posterior

  inv_log_p_pooled <- rstan::extract(model_predict_p_only, "inv_logit_posterior_pooled")$inv_logit_posterior_pooled

  a_eff <- (inv_log_p - inv_log_all) |> apply(2, quantile, probs = c(lower, 0.5, upper)) |> t() |> cbind(u_gq_in) |> dplyr::rename(l_e_a = 1, m_e_a = 2, u_e_a = 3)
  a_eff_pooled <- (inv_log_p_pooled - inv_log_all_pooled) |> apply(2, quantile, probs = c(lower, 0.5, upper)) |> t() |> cbind(u_gq_in) |> dplyr::rename(l_e_a_pooled = 1, m_e_a_pooled = 2, u_e_a_pooled = 3)
  r_eff <- (1 - (inv_log_all / inv_log_p)) |> apply(2, quantile, probs = c(lower, 0.5, upper)) |> t() |> cbind(u_gq_in) |> dplyr::rename(l_e = 1, m_e = 2, u_e = 3)
  r_eff_pooled <- (1 - (inv_log_all_pooled / inv_log_p_pooled)) |> apply(2, quantile, probs = c(lower, 0.5, upper)) |> t() |> cbind(u_gq_in) |> dplyr::rename(l_e_pooled = 1, m_e_pooled = 2, u_e_pooled = 3)

  out <- left_join(out, r_eff) |> left_join(a_eff) |> left_join(a_eff_pooled) |> left_join(r_eff_pooled)

  rm(list = c("model_predict_p_only", "inv_log_p", "inv_log_p_pooled", "a_eff", "a_eff_pooled", "r_eff", "r_eff_pooled"))

  return(out)
}

pred_mean = process_gq(pred_data = pred_data_mean_prev,
                       u_gq_in = u_gq_mean_prev)

pred_diff = process_gq(pred_data = pred_data_diff_prev,
                       u_gq_in = u_gq_diff_prev)

##################################################
##### mean odds ratios between 0 and 3 years #####
##################################################

u_gq_mean_prev_or <- u_gq_mean_prev |> subset(years == 1) |>
  mutate(l_in = l - 1,
         li_in = case_when(li <= 5 ~ li - 1,
                           li == 7 ~ 5),
         mean_or_3_l = NA, mean_or_3_m = NA, mean_or_3_u = NA,
         pooled_mean_or_3_l = NA, pooled_mean_or_3_m = NA, pooled_mean_or_3_u = NA,
         mean_prev_3_l = NA, mean_prev_3_m = NA, mean_prev_3_u = NA,
         mean_eff_3_l = NA, mean_eff_3_m = NA, mean_eff_3_u = NA) |> select(-years)

u_gq_diff_prev_or <- u_gq_diff_prev |> subset(years == 1) |>
  mutate(l_in = l - 1, li_in = case_when(li <= 5 ~ li - 1,
                                         li == 7 ~ 5),
         mean_or_3_l = NA, mean_or_3_m = NA, mean_or_3_u = NA,
         pooled_mean_or_3_l = NA, pooled_mean_or_3_m = NA, pooled_mean_or_3_u = NA,
         mean_prev_3_l = NA, mean_prev_3_m = NA, mean_prev_3_u = NA,
         mean_eff_3_l = NA, mean_eff_3_m = NA, mean_eff_3_u = NA) |> select(-years)

calc_or_mean_prev <- function(u_gq_in){

  alpha <- rstan::extract(fit_full, "alpha")[[1]]
  kappa <- rstan::extract(fit_full, "kappa")[[1]]
  gamma <- rstan::extract(fit_full, "gamma")[[1]]
  delta <- rstan::extract(fit_full, "alpha")[[1]]
  omega <- rstan::extract(fit_full, "omega")[[1]]

  sigma_e_r <- rstan::extract(fit_full, "sigma_e_r_train")[[1]]

  e_ij_all <- matrix(rnorm(iter/2*4 * length(unique(u_gq_in$i)), mean = 0, sd = sigma_e_r), nrow = nrow(sigma_e_r))

  for(i in 1:nrow(u_gq_in)){

    e_ij <- e_ij_all[, u_gq_in[i,"i"]]

    alpha_i <- rstan::extract(fit_full, "alpha_i")[[1]][, u_gq_in[i, "i"]]

    if(u_gq_in[i, "l_in"] == 0){
      theta_l <- kappa_l <- omega_l <- delta_l <- rep(0, iter/2*4)
    } else{
      theta_l <- rstan::extract(fit_full, "theta_l")[[1]][, u_gq_in[i, "l_in"]]
      kappa_l <- rstan::extract(fit_full, "kappa_l")[[1]][, u_gq_in[i, "l_in"]]
      omega_l <- rstan::extract(fit_full, "omega_l")[[1]][, u_gq_in[i, "l_in"]]
      delta_l <- rstan::extract(fit_full, "delta_l")[[1]][, u_gq_in[i, "l_in"]]
    }

    theta_li <- if(u_gq_in[i, "l_in"] == 0){rep(0, iter/2*4)} else{rstan::extract(fit_full, "theta_li")[[1]][, u_gq_in[i, "li_in"]]}

    # odds ratios
    log_or_3 <- (theta_l * 3) + (theta_li * 3) + ((kappa_l * 3^2)/2) + (omega_l * u_gq_in[i, "mean_prev"] * 3) + (delta_l * u_gq_in[i, "diff_prev"] * 3)
    or_3 <- exp(log_or_3)
    mean_or_3 <- quantile(log_or_3, probs = c(lower, 0.5, upper))

    log_or_3_p <-  (theta_l * 3) + ((kappa_l * 3^2)/2) + (omega_l * u_gq_in[i, "mean_prev_pooled"] * 3) + (delta_l * u_gq_in[i, "diff_prev_pooled"] * 3)
    or_3_p <- exp(log_or_3_p)
    mean_or_3_p <- quantile(log_or_3_p, probs = c(lower, 0.5, 0.975))

    D <- kappa + kappa_l + delta * u_gq_in[i, "diff_prev"]
    N <- alpha + alpha_i + theta_l + theta_li + omega * u_gq_in[i, "mean_prev"] + gamma * u_gq_in[i, "diff_prev"] +
      omega_l * u_gq_in[i, "mean_prev"] + delta_l * u_gq_in[i, "diff_prev"] + e_ij
    mean_prev <- (log(exp(D * 3 + N) + 1) - log(exp(N) + 1)) / (D * 3)

    mean_prev_3 <- quantile(mean_prev, probs = c(lower, 0.5, upper))

    D_c <- kappa + delta * u_gq_in[i, "diff_prev"]
    N_c <- alpha + alpha_i + omega * u_gq_in[i, "mean_prev"] + gamma * u_gq_in[i, "diff_prev"] + e_ij
    mean_prev_c <- (log(exp(D_c * 3 + N_c) + 1) - log(exp(N_c) + 1)) / (D_c * 3)

    eff_3 <- (1 - (or_3 / (1 - mean_prev_c + (mean_prev_c * or_3)))) * 100
    mean_eff_3 <- quantile(eff_3, probs = c(lower, 0.5, upper))

    u_gq_in[i, "mean_or_3_l"] <- mean_or_3[1]
    u_gq_in[i, "mean_or_3_m"] <- mean_or_3[2]
    u_gq_in[i, "mean_or_3_u"] <- mean_or_3[3]

    u_gq_in[i, "pooled_mean_or_3_l"] <- mean_or_3_p[1]
    u_gq_in[i, "pooled_mean_or_3_m"] <- mean_or_3_p[2]
    u_gq_in[i, "pooled_mean_or_3_u"] <- mean_or_3_p[3]

    u_gq_in[i, "mean_prev_3_l"] <- mean_prev_3[1]
    u_gq_in[i, "mean_prev_3_m"] <- mean_prev_3[2]
    u_gq_in[i, "mean_prev_3_u"] <- mean_prev_3[3]

    u_gq_in[i, "mean_eff_3_l"] <- mean_eff_3[1]
    u_gq_in[i, "mean_eff_3_m"] <- mean_eff_3[2]
    u_gq_in[i, "mean_eff_3_u"] <- mean_eff_3[3]
  }

  return(u_gq_in)
}

u_gq_mean_prev_or_3 <- calc_or_mean_prev(u_gq_mean_prev_or)
u_gq_diff_prev_or_3 <- calc_or_mean_prev(u_gq_diff_prev_or)

#############################################
##### mean prevalence with sample times #####
#############################################

u_gq_BL_prev <- u_gq |>
  full_join(unique(COMBO_stan[,c("study", "time")]), relationship = "many-to-many") |>
  left_join(unique(ucs[,c("study", "mean_BL_prev")]))

u_gq_BL_prev <- u_gq_BL_prev[rep(1:nrow(u_gq_BL_prev), each = length(diff_prev)),] |>
  mutate(diff_prev = rep(diff_prev, nrow(u_gq_BL_prev)),
         start_prev = diff_prev + mean_BL_prev,
         diff_prev_pooled = diff_prev,
         mean_prev_pooled = mean_prev_pooled) |>
  filter(start_prev <= 1 & start_prev >= 0)

pred_data_BL_prev <- data_in_full |>
  purrr::list_assign(gq = 1,
                     N_gq = nrow(u_gq_BL_prev),
                     N_ij_gq = length(unique(u_gq_BL_prev$ij)), # the same random effect for all the clusters in each study
                     N_ij_unq_gq = length(unique(u_gq_BL_prev$i)),
                     pmat_ij_gq = fastDummies::dummy_cols(u_gq_BL_prev$ij)[,-1] |> as.matrix(),
                     r_id_gq = seq(1, length(unique(u_gq_BL_prev$i))),
                     pmat_i_gq = fastDummies::dummy_cols(u_gq_BL_prev$i)[,-1] |> as.matrix(),
                     pmat_l_gq = fastDummies::dummy_cols(u_gq_BL_prev$l)[,-c(1, 2)] |> as.matrix(),
                     pmat_li_gq = fastDummies::dummy_cols(u_gq_BL_prev$li)[,-c(1, 2)] |> as.matrix(),

                     pmat_i_pooled_gq = matrix(0, nrow = nrow(u_gq_BL_prev), ncol = length(unique(u_gq_m$i))),
                     pmat_li_pooled_gq = matrix(0, nrow = nrow(u_gq_BL_prev), ncol = length(unique(u_gq_m$li)) - 1),
                     pmat_it_pooled_gq = matrix(0, nrow = nrow(u_gq_BL_prev), ncol = length(unique(u_gq_m$it))),

                     ij_train_gq = rep(0, length(unique(u_gq_BL_prev$ij))),
                     ij_unq_train_gq = rep(1, length(unique(u_gq_BL_prev$i))),
                     years_gq = u_gq_BL_prev$time/12,
                     m_prob_s_gq = u_gq_BL_prev$mean_BL_prev,
                     d_prob_s_gq = u_gq_BL_prev$diff_prev,
                     m_prob_s_gq_pooled = u_gq_BL_prev$mean_prev_pooled,
                     d_prob_s_gq_pooled = u_gq_BL_prev$diff_prev_pooled)

model_predict_BL_prev <- rstan::gqs(stan_model,
                            draws = as.matrix(fit_full),
                            data = pred_data_BL_prev,
                            seed = 123)

inv_log_BL_prev <- rstan::extract(model_predict_BL_prev, "inv_logit_posterior")$inv_logit_posterior

u_gq_BL_prev <- u_gq_BL_prev |> left_join(unique(u_gq_BL_prev[,c("l", "i", "diff_prev_pooled", "mean_prev_pooled")]) |> mutate(group = row_number()))

inv_log_BL_prev <- apply(rowsum(t(inv_log_BL_prev), u_gq_BL_prev$group) / as.vector(table(u_gq_BL_prev$group)), 1, quantile, probs = c(lower, 0.5, upper)) |> t() |>
  cbind(u_gq_BL_prev |> select(-time) |> distinct() |> arrange(group)) |> dplyr::rename(l_p = 1, m_p = 2, u_p = 3)

### saving all results
saveRDS(list("pred_mean" = pred_mean,
             "pred_diff" = pred_diff,
             "or_mean_3y_mean" = u_gq_mean_prev_or_3,
             "or_mean_3y_diff" = u_gq_diff_prev_or_3,
             "inv_log_BL_prev" = inv_log_BL_prev),
        file = paste0(model, "_pred.rds"))

