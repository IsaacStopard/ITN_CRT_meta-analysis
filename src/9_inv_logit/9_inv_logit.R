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

diff_prev <- seq(-1, 1, 0.05)
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

saveRDS(list("pred_mean" = pred_mean,
             "pred_diff" = pred_diff),
        file = paste0(model, "_pred.rds"))

#
# filter_gq <- function(u_gq_m, plot_name, data_limits){
#
#   # subsetting so for each study the mean values are the same
#   if(plot_name != "m_plot" && plot_name != "mn_plot"){
#       u_gq_m <- lapply(1:nrow(data_limits), function(i, data_limits, u_gq_m){
#         x <- round(data_limits[i, "mean_BL_prev"], digit = 3)
#         y <- round(data_limits[i, "mean_BL_net_use"], digit = 3)
#         s <- data_limits[i, "i"]
#         return(filter(u_gq_m, round(mean_prev, digits = 3) == x &
#                         round(mean_net, digits = 3) == y &
#                         i == s)
#                )},
#         data_limits = data_limits, u_gq_m = u_gq_m) |>
#       bind_rows()
#
#       u_gq_m <- lapply(1:nrow(data_limits), function(i, data_limits, u_gq_m){
#         a <- round(data_limits[i, "min_prev"], digit = 3)
#         b <- round(data_limits[i, "min_net"], digit = 3)
#         c <- round(data_limits[i, "max_prev"], digit = 3)
#         d <- round(data_limits[i, "max_net"], digit = 3)
#         s <- data_limits[i, "i"]
#         return(filter(u_gq_m,
#                       round(start_prev, digits = 3) >= a &
#                         round(start_prev, digits = 3) <= c &
#                         round(start_net, digits = 3) >= b &
#                         round(start_net, digits = 3) <= d &
#                         i == s)
#         )},
#         data_limits = data_limits, u_gq_m = u_gq_m) |>
#         bind_rows()
#   }
# # predictions that vary the mean values
#   if(plot_name == "m_plot" | plot_name == "mn_plot"){
#     u_gq_m <- u_gq_m |> filter(round(start_prev, digits = 3) >= round(min(data_limits$mean_BL_prev), digits = 3) &
#                                  round(start_prev, digits = 3) <= round(max(data_limits$mean_BL_prev), digits = 3) &
#                                round(start_net, digits = 3) >= round(min(data_limits$mean_BL_net_use), digits = 3) &
#                                  round(start_net, digits = 3) <= round(max(data_limits$mean_BL_net_use), digits = 3)
#                                  )
#     u_gq_m <- if(plot_name == "m_plot"){
#       lapply(1:nrow(data_limits),
#              function(i, data_limits, u_gq_m){
#                y <- round(data_limits[i, "mean_BL_net_use"], digit = 3)
#                s <- data_limits[i, "i"]
#       return(
#         filter(u_gq_m,
#           round(mean_net, digits = 3) == y &
#                       i == s)
#                )
#              },
#       data_limits = data_limits, u_gq_m = u_gq_m) |>
#       bind_rows()
#
#     } else{
#       u_gq_m <- lapply(1:nrow(data_limits), function(i, data_limits, u_gq_m){
#         x <- round(data_limits[i, "mean_BL_prev"], digit = 3)
#         s <- data_limits[i, "i"]
#         return(filter(u_gq_m, round(mean_prev, digits = 3) == x &
#                         i == s)
#         )},
#         data_limits = data_limits, u_gq_m = u_gq_m) |>
#         bind_rows()
#     }
#     }
#
#   return(u_gq_m)
# }
#
# get_predictions <- function(model,
#                             u_gq,
#                             years,
#                             mean_prev,
#                             diff_prev,
#                             mean_net,
#                             diff_net,
#                             mean_prev_pooled,
#                             mean_net_pooled,
#                             data_in_full, stan_model, fit_full,
#                             plot_name, data_limits){
#
#   bp_bn_gq <- expand.grid("mean_prev" = mean_prev,
#                           "diff_prev" = diff_prev,
#                           "years" = years,
#                           "mean_net" = mean_net,
#                           "diff_net" = diff_net) |>
#           mutate(start_prev = diff_prev + mean_prev,
#                  start_net = diff_net + mean_net,
#                  mean_prev_pooled = mean_prev_pooled,
#                  diff_prev_pooled = diff_prev,
#                  mean_net_pooled = mean_net_pooled,
#                  diff_net_pooled = diff_net,
#                  start_prev_pooled = diff_prev_pooled + mean_prev_pooled,
#                  start_net_pooled = diff_net_pooled + mean_net_pooled) |>
#           filter(start_prev <= 1 & start_prev >= 0 &
#                  start_net <= 1 & start_net >= 0 &
#                  start_prev_pooled <= 1 & start_prev_pooled >= 0 &
#                  start_net_pooled <= 1 & start_net_pooled >= 0)
#
#   u_gq_m <- u_gq[rep(1:nrow(u_gq), nrow(bp_bn_gq)),] |>
#           dplyr::mutate(mean_prev = rep(bp_bn_gq[, "mean_prev"], each = nrow(u_gq)),
#                   start_prev = rep(bp_bn_gq[, "start_prev"], each = nrow(u_gq)),
#                   start_net = rep(bp_bn_gq[, "start_net"], each = nrow(u_gq)),
#                   diff_prev = rep(bp_bn_gq[, "diff_prev"], each = nrow(u_gq)),
#                   mean_net = rep(bp_bn_gq[, "mean_net"], each = nrow(u_gq)),
#                   diff_net = rep(bp_bn_gq[, "diff_net"], each = nrow(u_gq)),
#                   years = rep(bp_bn_gq[, "years"], each = nrow(u_gq)),
#                   mean_prev_pooled = rep(bp_bn_gq[, "mean_prev_pooled"], each = nrow(u_gq)),
#                   diff_prev_pooled = rep(bp_bn_gq[, "diff_prev_pooled"], each = nrow(u_gq)),
#                   mean_net_pooled = rep(bp_bn_gq[, "mean_net_pooled"], each = nrow(u_gq)),
#                   diff_net_pooled = rep(bp_bn_gq[, "diff_net_pooled"], each = nrow(u_gq)))
#
#   u_gq_m <- filter_gq(u_gq_m = u_gq_m, plot_name = plot_name, data_limits = data_limits) |>
#         mutate(
#           ij = i,
#           rn = dplyr::row_number()
#           )
#
#   pred_data <- data_in_full |>
#     purrr::list_assign(gq = 1,
#                        N_gq = nrow(u_gq_m),
#                        N_ij_gq = length(unique(u_gq_m$ij)), # the same random effect for all the clusters in each study
#                        N_ij_unq_gq = length(unique(u_gq_m$i)),
#                        pmat_ij_gq = fastDummies::dummy_cols(u_gq_m$ij)[,-1] |> as.matrix(),
#                        r_id_gq = seq(1, length(unique(u_gq_m$i))),
#                        pmat_i_gq = fastDummies::dummy_cols(u_gq_m$i)[,-1] |> as.matrix(),
#                        pmat_l_gq = fastDummies::dummy_cols(u_gq_m$l)[,-c(1, 2)] |> as.matrix(),
#                        pmat_li_gq = fastDummies::dummy_cols(u_gq_m$li)[,-c(1, 2)] |> as.matrix(),
#                        pmat_it_gq = fastDummies::dummy_cols(u_gq_m$it)[,-1] |> as.matrix(),
#
#                        pmat_i_pooled_gq = matrix(0, nrow = nrow(u_gq_m), ncol = length(unique(u_gq_m$i))),
#                        pmat_li_pooled_gq = matrix(0, nrow = nrow(u_gq_m), ncol = length(unique(u_gq_m$li)) - 1),
#                        pmat_it_pooled_gq = matrix(0, nrow = nrow(u_gq_m), ncol = length(unique(u_gq_m$it))),
#
#                        ij_train_gq = rep(0, length(unique(u_gq_m$ij))),
#                        ij_unq_train_gq = rep(1, length(unique(u_gq_m$i))),
#                        years_gq = u_gq_m$years,
#                        m_prob_s_gq = u_gq_m$mean_prev,
#                        d_prob_s_gq = u_gq_m$diff_prev,
#                        m_net_s_gq = u_gq_m$mean_net,
#                        d_net_s_gq = u_gq_m$diff_net,
#                        m_prob_s_gq_pooled = u_gq_m$mean_prev_pooled,
#                        d_prob_s_gq_pooled = u_gq_m$diff_prev_pooled,
#                        m_net_s_gq_pooled = u_gq_m$mean_net_pooled,
#                        d_net_s_gq_pooled = u_gq_m$diff_net_pooled)
#
#   model_predict <- rstan::gqs(stan_model,
#                               draws = as.matrix(fit_full),
#                               data = pred_data,
#                               seed = 123)
#
#   #inv_log_all <- rstan::extract(model_predict, "inv_logit_posterior")$inv_logit_posterior
#
#   #inv_log_all_pooled <- rstan::extract(model_predict, "inv_logit_posterior_pooled")$inv_logit_posterior_pooled
#
#   #inv_log <- inv_log_all |> apply(2, quantile, probs = c(lower, 0.5, upper)) |> t() |> cbind(u_gq_m) |> dplyr::rename(l_p = 1, m_p = 2, u_p = 3)
#
#   #inv_log_pooled <- inv_log_all_pooled |> apply(2, quantile, probs = c(lower, 0.5, upper)) |> t() |> cbind(u_gq_m) |> dplyr::rename(l_p_pooled = 1, m_p_pooled = 2, u_p_pooled = 3)
#
#   o_r <- rstan::extract(model_predict, "o_r")$o_r |> apply(2, quantile, probs = c(lower, 0.5, upper)) |> t() |> cbind(u_gq_m) |> dplyr::rename(l_or = 1, m_or = 2, u_or = 3)
#
#   o_r_pooled <- rstan::extract(model_predict, "o_r_pooled")$o_r_pooled |> apply(2, quantile, probs = c(lower, 0.5, upper)) |> t() |> cbind(u_gq_m) |> dplyr::rename(l_or_pooled = 1, m_or_pooled = 2, u_or_pooled = 3)
#
#   out <- o_r |> left_join(o_r_pooled) #inv_log |> left_join()|> left_join(inv_log_pooled)
#
#   rm(list = c("model_predict"))
#
#   #pred_data_p_only <- pred_data
#   #pred_data_p_only$pmat_l_gq[pred_data_p_only$pmat_l_gq == 1] <- 0
#   #pred_data_p_only$pmat_li_gq[pred_data_p_only$pmat_li_gq == 1] <- 0
#
#   #model_predict_p_only <- rstan::gqs(stan_model,
#   #                                   draws = as.matrix(fit_full),
#   #                                   data = pred_data_p_only,
#   #                                   seed = 123)
#
#   #inv_log_p <- rstan::extract(model_predict_p_only, "inv_logit_posterior")$inv_logit_posterior
#
#   #inv_log_p_pooled <- rstan::extract(model_predict_p_only, "inv_logit_posterior_pooled")$inv_logit_posterior_pooled
#
#   #a_eff <- (inv_log_p - inv_log_all) |> apply(2, quantile, probs = c(lower, 0.5, upper)) |> t() |> cbind(u_gq_m) |> dplyr::rename(l_e_a = 1, m_e_a = 2, u_e_a = 3)
#   #a_eff_pooled <- (inv_log_p_pooled - inv_log_all_pooled) |> apply(2, quantile, probs = c(lower, 0.5, upper)) |> t() |> cbind(u_gq_m) |> dplyr::rename(l_e_a_pooled = 1, m_e_a_pooled = 2, u_e_a_pooled = 3)
#   #r_eff <- (1 - (inv_log_all / inv_log_p)) |> apply(2, quantile, probs = c(lower, 0.5, upper)) |> t() |> cbind(u_gq_m) |> dplyr::rename(l_e = 1, m_e = 2, u_e = 3)
#   #r_eff_pooled <- (1 - (inv_log_all_pooled / inv_log_p_pooled)) |> apply(2, quantile, probs = c(lower, 0.5, upper)) |> t() |> cbind(u_gq_m) |> dplyr::rename(l_e_pooled = 1, m_e_pooled = 2, u_e_pooled = 3)
#
#   #rm(list = c("model_predict_p_only"))
#
#   #out <- left_join(out, r_eff) |> left_join(a_eff) |> left_join(a_eff_pooled) |> left_join(r_eff_pooled)
#
#   return(out |> mutate(plot_name = plot_name))
# }
#
# # different combinations of predictors
# # so each plot is done with the same mean BL prev and net use
# ucs <- unique(COMBO_stan[, c("cluster", "study", "i", "BL_prev_num", "BL_prev_denom", "net_use_num", "net_use_denom")]) |>
#   mutate(BL_prev = BL_prev_num/BL_prev_denom,
#          BL_net_use = net_use_num/net_use_denom)
#
# ucs <- ucs |> left_join(ucs |> group_by(study) |>
#               summarise(mean_BL_prev = mean(BL_prev),
#                         mean_BL_net_use = mean(BL_net_use))) |>
#   mutate(diff_BL_prev = BL_prev - mean_BL_prev,
#          diff_BL_net_use = BL_net_use - mean_BL_net_use)
#
# mean_prev_pooled <- mean(ucs[,"BL_prev"])
# mean_net_pooled <- mean(ucs[,"BL_net_use"])
#
# data_limits <- COMBO_stan |> group_by(study, i) |> summarise(min_prev = min(BL_prev), max_prev = max(BL_prev), min_net = min(net_use), max_net = max(net_use)) |>
#   left_join(unique(ucs[,c("i", "mean_BL_prev", "mean_BL_net_use")])) |> as.data.frame()
#
# # for time plots
# years_t <- c(seq(0.3, 3, 0.3), unique(COMBO_stan$time/12)) |> unique() |> sort()
# mean_prev_t <- unique(ucs$mean_BL_prev)
# diff_prev_t <- c(0)
# mean_net_t <- unique(ucs$mean_BL_net_use)
# diff_net_t <- diff_prev_t
#
# pred_t <- get_predictions(model = model,
#                           u_gq = u_gq,
#                           years = years_t,
#                           mean_prev = mean_prev_t,
#                           diff_prev = diff_prev_t,
#                           mean_net = mean_net_t,
#                           diff_net = diff_net_t,
#                           mean_prev_pooled = mean_prev_pooled,
#                           mean_net_pooled = mean_net_pooled,
#                           data_in_full = data_in_full,
#                           stan_model = stan_model, fit_full = fit_full,
#                           plot_name = "t_plot",
#                           data_limits = data_limits)
#
# # for baseline prev plots
# years_m <- c(1)
# mean_prev_m <- seq(0, 1, 0.025)
# diff_prev_m <- 0
# mean_net_m <- mean_net_t
# diff_net_m <- 0
#
# pred_m <- get_predictions(model = model, u_gq = u_gq,
#                   years = years_m,
#                   mean_prev = mean_prev_m,
#                   diff_prev = diff_prev_m,
#                   mean_net = mean_net_m,
#                   diff_net = diff_net_m,
#                   mean_prev_pooled = mean_prev_pooled,
#                   mean_net_pooled = mean_net_pooled,
#                   data_in_full = data_in_full,
#                   stan_model = stan_model, fit_full = fit_full,
#                   plot_name = "m_plot",
#                   data_limits = data_limits)
#
# years_p <- c(1)
# mean_prev_p <- mean_prev_t
# diff_prev_p <- seq(-1, 1, 0.025)
# mean_net_p <- mean_net_t
# diff_net_p <- 0
#
# pred_p <- get_predictions(model = model,
#                           u_gq = u_gq,
#                           years = years_p,
#                           mean_prev = mean_prev_p,
#                           diff_prev = diff_prev_p,
#                           mean_net = mean_net_p,
#                           diff_net = diff_net_p,
#                           mean_prev_pooled = mean_prev_pooled,
#                           mean_net_pooled = mean_net_pooled,
#                           data_in_full = data_in_full,
#                           stan_model = stan_model, fit_full = fit_full,
#                           plot_name = "p_plot",
#                           data_limits = data_limits)
#
# # for mean baseline net use plots
# years_mn <- years_p
# mean_prev_mn <- mean_prev_t
# diff_prev_mn <- 0
# mean_net_mn <- seq(0, 1.0, 0.025)
# diff_net_mn <- 0
#
# pred_mn <- get_predictions(model = model, u_gq = u_gq,
#                   years = years_mn,
#                   mean_prev = mean_prev_mn,
#                   diff_prev = diff_prev_mn,
#                   mean_net = mean_net_mn,
#                   diff_net = diff_net_mn,
#                   mean_prev_pooled = mean_prev_pooled,
#                   mean_net_pooled = mean_net_pooled,
#                   data_in_full = data_in_full,
#                   stan_model = stan_model,
#                   fit_full = fit_full,
#                   plot_name = "mn_plot",
#                   data_limits = data_limits)
#
# # for difference in baseline net use plots
# years_dn <- years_p
# mean_prev_dn <- mean_prev_t
# diff_prev_dn <- 0
# mean_net_dn <- mean_net_t
# diff_net_dn <- seq(-1, 1, 0.025)
#
# pred_dn <- get_predictions(model = model, u_gq = u_gq,
#                   years = years_dn,
#                   mean_prev = mean_prev_dn,
#                   diff_prev = diff_prev_dn,
#                   mean_net = mean_net_dn,
#                   diff_net = diff_net_dn,
#                   mean_prev_pooled = mean_prev_pooled,
#                   mean_net_pooled = mean_net_pooled,
#                   data_in_full = data_in_full,
#                   stan_model = stan_model,
#                   fit_full = fit_full,
#                   plot_name = "dn_plot",
#                   data_limits = data_limits)
#
#
# saveRDS(list("pred_t" = pred_t,
#              "pred_m" = pred_m,
#              "pred_p" = pred_p,
#              "pred_mn" = pred_mn,
#              "pred_dn" = pred_dn,
#              "data_limits" = data_limits),
#         file = paste0(model, "_pred.rds"))
#
#
