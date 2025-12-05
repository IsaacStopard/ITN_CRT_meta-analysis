library(rstan)
library(tidyverse)
library(loo)
library(patchwork)

ggplot2::theme_set(theme_bw() +
  theme(text = element_text(size = 14),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
)

##############################
##### trial GLM fits #########
##############################

orderly2::orderly_dependency("1_data_cleaning", "latest()", c("stan_data.rds"))

list2env(readRDS("stan_data.rds"), envir = .GlobalEnv)

orderly2::orderly_dependency("2_fits",
                             quote(latest(parameter:model == "m0_c")),
                             "m0_c_fit_full.rds"
)

orderly2::orderly_dependency("2_fits",
                             quote(latest(parameter:model == "m1_re")),
                             "m1_re_fit_full.rds"
)

orderly2::orderly_dependency("2_fits",
                             quote(latest(parameter:model == "m2_re")),
                             "m2_re_fit_full.rds"
)

m2_fit_full <- readRDS(file = "m2_re_fit_full.rds")

m1_fit_full <- readRDS(file = "m1_re_fit_full.rds")

m0_fit_full <- readRDS(file = "m0_c_fit_full.rds")

orderly2::orderly_dependency("9_inv_logit",
                             quote(latest(parameter:model == "m2_re")),
                             "m2_re_pred.rds")

list2env(readRDS(file = "m2_re_pred.rds"), envir = .GlobalEnv)

orderly2::orderly_shared_resource("m0_c.stan",
                                  "m1_re.stan",
                                  "m2_re.stan")

m0_stan_model <- rstan::stan_model(file = "m0_c.stan")

m1_stan_model <- rstan::stan_model(file = "m1_re.stan")

m2_stan_model <- rstan::stan_model(file = "m2_re.stan")

############################
##### model comparison #####
############################

orderly2::orderly_dependency("3_cv",
                             quote(latest(parameter:model == "m0_c")),
                             "m0_c_log_pd_kfold.rds")

orderly2::orderly_dependency("3_cv",
                             quote(latest(parameter:model == "m1_re")),
                             "m1_re_log_pd_kfold.rds")

orderly2::orderly_dependency("3_cv",
                             quote(latest(parameter:model == "m2_re")),
                             "m2_re_log_pd_kfold.rds")

m0_c_log_pd_kfold <- readRDS("m0_c_log_pd_kfold.rds")
m1_re_log_pd_kfold <- readRDS("m1_re_log_pd_kfold.rds")
m2_re_log_pd_kfold <- readRDS("m2_re_log_pd_kfold.rds")

# checking the fits
# leave one out cross validation
m0_c_elpd_kfold <- elpd(m0_c_log_pd_kfold)
m1_re_elpd_kfold <- elpd(m1_re_log_pd_kfold)
m2_re_elpd_kfold <- elpd(m2_re_log_pd_kfold)

loo_compare(m0_c_elpd_kfold, m1_re_elpd_kfold, m2_re_elpd_kfold)

# actual vs fitted plot
extract_af <- function(fit_full){
  cbind(rstan::extract(fit_full, "inv_logit_posterior")$inv_logit_posterior |>
          apply(2, quantile, probs = c(lower, 0.5, upper)) |>
          t() |> as.data.frame() |> rename(l_p = 1, m_p = 2, u_p = 3),
        COMBO_stan[,c("cluster", "time", "l", "li", "study", "net", "prev", "prev_lower", "prev_upper")])
}

m0_af <- extract_af(m0_fit_full) |> mutate(model = "model 0")
m1_af <- extract_af(m2_fit_full) |> mutate(model = "model 1")
m2_af <- extract_af(m2_fit_full) |> mutate(model = "model 2")

af_plot <- ggplot(data = rbind(m0_af, m1_af, m2_af),
                  aes(x = prev, xmin = prev_lower, xmax = prev_upper, y = m_p, ymin = l_p, ymax = u_p, colour = model)) +
  geom_errorbarh(alpha = 0.1, linewidth = 0.1) +
  geom_pointrange(alpha = 0.5, size = 0.7) +
  geom_abline(slope = 1, linetype = 2) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
  scale_x_continuous(labels = scales::percent) +
  xlab("Sampled prevalence") + ylab("Predicted prevalence") +
  scale_colour_manual(name = "", values = c("#0072B2", "#D55E00", "#009E73")) +
  guides(colour  = guide_legend(position = "inside")) +
  theme(legend.position.inside = c(0.9, 0.25))

# random effects variance
extract_sd <- function(model_fit, model_name){
  out <- rstan::extract(model_fit, "sigma_e_r_train")$sigma_e_r |> as.data.frame()
  colnames(out) <- unique(COMBO_stan[,c("i", "study")])[,2]
  return(out |> pivot_longer(cols = 1:4, names_to = "study", values_to = "sd") |> mutate(model = model_name))
}

s_all <- rbind(extract_sd(m0_fit_full, "model 0"),
               extract_sd(m1_fit_full, "model 1"),
               extract_sd(m2_fit_full, "model 2")) |>
  mutate(study_place = case_when(
    str_detect(study, "Protopopoff") ~ "Tanzania (2015)",
    str_detect(study, "Staedke") ~ "Uganda (2017)",
    str_detect(study, "Mosha") ~ "Tanzania (2019)",
    str_detect(study, "Accrombessi") ~ "Benin (2020)"
  )
  )

sd_plot <- ggplot(data = s_all, aes(x = sd, fill = model, col = model, group = interaction(study_place, model))) +
  geom_density(alpha = 0.25) +
  facet_wrap(~study_place) +
  xlab(expression(paste("Study-specific cluster random effect standard deviation (", sigma["i"], ")"))) +
  ylab("Density") +
  scale_fill_manual(name = "", values = c("#0072B2", "#D55E00", "#009E73")) +
  scale_colour_manual(name = "", values = c("#0072B2", "#D55E00", "#009E73")) +
  guides(fill  = guide_legend(position = "inside")) +
  theme(legend.position.inside = c(0.9, 0.25))

ggsave(
  file = "model_checks.pdf",
  device = "pdf",
  plot = af_plot / sd_plot  +
    plot_annotation(tag_levels = 'A'),
  width = 20, height = 35/3*2,
  units = "cm"
)

##########################
##### data wrangling #####
##########################

ucs <- unique(COMBO_stan[, c("cluster", "study", "i", "BL_prev_num", "BL_prev_denom", "BL_prev")])

COMBO_stan <- COMBO_stan |>
  left_join(ucs |> group_by(study) |>
              summarise(mean_BL_prev = mean(BL_prev))) |>
  mutate(diff_base_prev = BL_prev - mean_BL_prev,
         BL_prev_group = cut(diff_base_prev, breaks = seq(-1, 1, 0.2)),
         study_place = case_when(
           str_detect(study, "Protopopoff") ~ "Tanzania (2015)",
           str_detect(study, "Staedke") ~ "Uganda (2017)",
           str_detect(study, "Mosha") ~ "Tanzania (2019)",
           str_detect(study, "Accrombessi") ~ "Benin (2020)"
         ))

COMBO_limits_t <- COMBO_stan |> group_by(study) |>
  summarise(lt = min(time)/12, ut = max(time)/12)

COMBO_limits_prev <- unique(COMBO_stan[,c("cluster", "study", "net", "diff_base_prev")]) |> group_by(study, net) |>
  summarise(lbp = min(diff_base_prev), ubp = max(diff_base_prev))

#####################
##### GLM plots #####
#####################

plot_labels_fun <- function(df){
  out <- df |> mutate(
                      study_place = case_when(
                        str_detect(study, "Protopopoff") ~ "Tanzania (2015)",
                        str_detect(study, "Staedke") ~ "Uganda (2017)",
                        str_detect(study, "Mosha") ~ "Tanzania (2019)",
                        str_detect(study, "Accrombessi") ~ "Benin (2020)"
                        ),
                      year_plot = paste(years, "year post ITN distribution")
                      )

  out$study_place <- factor(out$study_place, levels = c("Tanzania (2015)", "Uganda (2017)", "Tanzania (2019)", "Benin (2020)"))

  return(out)
}

pred_mean <- plot_labels_fun(pred_mean)
pred_diff <- plot_labels_fun(pred_diff)

cols <- unname(palette.colors(palette = "Okabe-Ito")[c(7,2,3,8)])
cols_net <- c("blue", "aquamarine", "darkgreen")

################
##### time #####
################

round(quantile(rstan::extract(m2_fit_full, "delta")[[1]], probs = c(lower, 0.5, upper)), digits = 2)

quantile(exp(rstan::extract(m1_fit_full, "kappa")[[1]]), probs = c(lower, 0.5, upper)) |> round(digits = 2)
quantile(exp(rstan::extract(m1_fit_full, "kappa")[[1]] + rstan::extract(m1_fit_full, "kappa_l")[[1]][,1]), probs = c(lower, 0.5, upper)) |> round(digits = 2)
quantile(exp(rstan::extract(m1_fit_full, "kappa")[[1]] + rstan::extract(m1_fit_full, "kappa_l")[[1]][,2]), probs = c(lower, 0.5, upper)) |> round(digits = 2)

################################
##### baseline prevalence  #####
################################

round(quantile(rstan::extract(m2_fit_full, "gamma")[[1]] + rstan::extract(m2_fit_full, "delta")[[1]] * 1, probs = c(lower, 0.5, upper)), digits = 2)
round(quantile(rstan::extract(m2_fit_full, "gamma")[[1]] + rstan::extract(m2_fit_full, "delta")[[1]] * 1 + rstan::extract(m2_fit_full, "delta_l")[[1]][,1], probs = c(lower, 0.5, upper)), digits = 2)
round(quantile(rstan::extract(m2_fit_full, "gamma")[[1]] + rstan::extract(m2_fit_full, "delta")[[1]] * 1 + rstan::extract(m2_fit_full, "delta_l")[[1]][,2], probs = c(lower, 0.5, upper)), digits = 2)

round(quantile(rstan::extract(m2_fit_full, "delta_l")[[1]][,1], probs = c(lower, 0.5, upper)), digits = 2)
round(quantile(rstan::extract(m2_fit_full, "delta_l")[[1]][,2], probs = c(lower, 0.5, upper)), digits = 2)

round(quantile(rstan::extract(m2_fit_full, "omega")[[1]], probs = c(lower, 0.5, upper)), digits = 2)
round(quantile(rstan::extract(m2_fit_full, "omega")[[1]] +  rstan::extract(m2_fit_full, "omega_l")[[1]][,1], probs = c(lower, 0.5, upper)), digits = 2)

m2_dp_plot <- ggplot(data = pred_diff,
                     aes(x = start_prev, y = m_p, ymin = l_p, ymax = u_p, fill = net,
                         group = interaction(study_place, net, year_plot))) +
  geom_ribbon(alpha = 0.1) +
  geom_line(aes(col = net, linetype = year_plot), linewidth = 1) +
  facet_wrap(~ study_place, nrow = 1) +
  ylab("Malaria prevalence") + xlab("Baseline prevalence") +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2), labels = scales::percent) +
  scale_colour_manual(values = cols_net, name = "ITN") +
  scale_fill_manual(values = cols_net, name = "ITN") +
  scale_x_continuous(limits = c(0, 1), breaks = seq(-0, 1, 0.25), labels = scales::percent) +
  scale_linetype_manual(values = c(1, 2, 3), name = "") +
  labs(tag = "A")

m2_o_r_dp_plot <- ggplot(data = pred_diff |> filter(net != "Pyrethroid-only"),
                         aes(x = diff_prev, ymin = l_or, ymax = u_or, y = m_or, fill = study_place, group = interaction(study_place, net, year_plot))) +
  geom_ribbon(inherit.aes = FALSE,
              aes(x = diff_prev, y = m_or_pooled, ymin = l_or_pooled, ymax = u_or_pooled,
                  group = interaction(study_place, year_plot)), alpha = 0.1) +
  geom_ribbon(alpha = 0.1) +
  geom_line(aes(col = study_place, linetype = year_plot), linewidth = 1) +
  geom_hline(yintercept = 1, col = "black", linetype = 2) +
  facet_wrap(~ net, nrow = 1) +
  ylab("Odds ratio of infection in trial ITN\nclusters relative to pyrethroid-only clusters") +
  xlab("Within-study differences in baseline prevalence") +
  scale_x_continuous(limits = c(-0.55, 0.55), breaks = seq(-0.5, 0.5, 0.25), labels = scales::percent) +
  scale_colour_manual(values = cols, name = "study") +
  scale_fill_manual(values = cols, name = "study") +
  geom_line(inherit.aes = FALSE, aes(x = diff_prev, y = m_or_pooled, linetype = year_plot), linewidth = 1) +
  scale_linetype_manual(values = c(1, 2, 3), name = "") +
  labs(tag = "B")

m2_mean_BL_prev_o_r_plot <- ggplot(data = pred_mean |>
                                     left_join(COMBO_stan |> group_by(net) |> summarise(min_mean_BL_prev = min(mean_BL_prev),
                                                                                        max_mean_BL_prev = max(mean_BL_prev))) |>
                                     subset(net != "Pyrethroid-only" & mean_prev >= min_mean_BL_prev & mean_prev <= max_mean_BL_prev),
                        aes(x = mean_prev, y = m_or, ymin = l_or, ymax = u_or,
                            fill = study_place, group = interaction(study_place, net, year_plot))) +
  geom_ribbon(aes(x = mean_prev, ymin = l_or_pooled, ymax = u_or_pooled,
                  group = interaction(net, year_plot)), inherit.aes = FALSE, alpha = 0.1) +
  geom_ribbon(alpha = 0.1) +
  geom_line(aes(col = study_place, linetype = year_plot), linewidth = 1) +
  geom_line(aes(x = mean_prev, y = m_or_pooled, linetype = year_plot), inherit.aes = FALSE, linewidth = 1) +
  geom_hline(yintercept = 1, linetype = 2) +
  facet_wrap(~net, scales = "free_x") +
  ylab("Odds ratio of infection in trial ITN\nclusters relative to pyrethroid-only clusters") +
  xlab("Between study mean baseline prevalence") +
  scale_y_continuous(limits = c(0, 2), breaks = seq(0, 2.5, 0.5)) +
  scale_colour_manual(values = cols, name = "study") +
  scale_fill_manual(values = cols, name = "study") +
  scale_x_continuous(labels = scales::percent) +
  scale_linetype_manual(values = c(1, 2, 3), name = "") +
  labs(tag = "C")

ggsave(file = "regression_base_prev.pdf",
       plot = m2_dp_plot /
         (m2_o_r_dp_plot + m2_mean_BL_prev_o_r_plot +
            plot_layout(guides = "collect", axis_titles = "collect")),
       device = "pdf",
       width = 50, height = 25,
       units = "cm"
)

############################
##### parameter values #####
############################

get_quantiles <- function(model_fit, param, dim = NULL){
  out <- rstan::extract(model_fit, param)[[1]]
  if(!is.null(dim)){
    out <- out[,dim]
  }
  prob_g0 <- round(sum(out > 0) / length(out), digits = 2)
  out <- round(quantile(out, probs = c(lower, 0.5, upper)), digits = 2)
  return(paste0(out[2], " (", out[1], " - ", out[3], " 95%CI, p(>0) = ", prob_g0,")"))
}

get_params <- function(fit){
  params <- c("$\\alpha$")
  params_description <- ("pooled intercept")
  values <- c(get_quantiles(fit, "alpha", dim = NULL))
  # alpha_i
  params <- c(params, "$\\alpha_{i}$", "$\\alpha_{i}$", "$\\alpha_{i}$", "$\\alpha_i$")
  params_description <- c(params_description, "study effect (Tanzania (2015))", "study effect (Uganda (2017))", "study effect (Tanzania (2019))", "study effect (Benin (2020))")
  for(i in 1:4){
    values <- c(values, get_quantiles(fit, "alpha_i", dim = i))
  }
  # net effect
  params <- c(params, "$\\theta_{l}$")
  params_description <- c(params_description, "pyrethroid-PBO pooled treatment effect")
  values <- c(values, get_quantiles(fit, "theta_l", dim = 1))

  params_description <- c(params_description, "pyrethroid-PBO-Tanzania (2015) interaction", "pyrethroid-PBO-Uganda (2017) interaction", "pyrethroid-PBO-Tanzania (2019) interaction")
  params <- c(params, "$\\theta_{li}$", "$\\theta_{li}$", "$\\theta_{li}$")
  for(i in 1:3){
    values <- c(values, get_quantiles(fit, "theta_li", dim = i))
  }

  params <- c(params, "$\\theta_{l}$")
  params_description <- c(params_description, "pyrethroid-pyrrole pooled treatment effect")
  values <- c(values, get_quantiles(fit, "theta_l", dim = 2))

  params_description <- c(params_description, "pyrethroid-PBO-Tanzania (2019) interaction", "pyrethroid-pyrrole-Benin (2020) interaction")
  params <- c(params, "$\\theta_{li}$", "$\\theta_{li}$")
  values <- c(values, get_quantiles(fit, "theta_li", dim = 4), get_quantiles(fit, "theta_li", dim = 5))

  params_description <- c(params_description, "Tanzania (2015) cluster random effect standard deviation", "Uganda (2017) cluster random effect standard deviation",
                          "Tanzania (2019) cluster random effect standard deviation", "Benin (2020) cluster random effect standard deviation")

  params <- c(params, "$\\sigma_{i}$", "$\\sigma_{i}$", "$\\sigma_{i}$", "$\\sigma_{i}$")
  values <- c(values, get_quantiles(fit, "sigma_e_r", dim = 1), get_quantiles(fit, "sigma_e_r", dim = 2),
              get_quantiles(fit, "sigma_e_r", dim = 3), get_quantiles(fit, "sigma_e_r", dim = 4))

  params_description <- c(params_description, "pyrethroid-PBO between study random treatment effect standard deviation", "pyrethroid-pyrrole between study random treatment effect standard deviation")
  params <- c(params, "$\\tau_{l}", "$\\tau_{l}")
  values <- c(values, get_quantiles(fit, "tau_sd_li", dim = 1), get_quantiles(fit, "tau_sd_li", dim = 2))

  return(list("params" = params, "params_description" = params_description, "values" = values))
}

# model 2
m2_params <- get_params(fit = m2_fit_full)
m2_params_df <- data.frame("params" = m2_params$params, params_description = m2_params$params_description, "values" = m2_params$values)
m2_params_df <- rbind(m2_params_df,
                      data.frame(
                        "params" = c("$\\kappa$",
                                     "$\\kappa_{l}$",
                                     "$\\kappa_{l}",
                                     "$\\gamma$",
                                     "$\\gamma^{net}_{l}$",
                                     "$\\gamma^{net}_{l}$",
                                     "$\\delta$",
                                     "$\\omega$",
                                     "$\\omega^{net}_{l}$",
                                     "$\\omega^{net}_{l}$"),
                        "params_description" = c("year effect",
                                                 "pyrethroid-PBO year effect interaction",
                                                 "pyrethroid-pyrrole year effect interaction",
                                                 "within study differences in baseline prevalence effect",
                                                 "within study differences in baseline prevalence pyrethroid-PBO interaction",
                                                 "within study differences in baseline prevalence pyrethroid-pyrrole interaction",
                                                 "within study differences in baseline prevalence year interaction",
                                                 "mean study baseline prevalence effect",
                                                 "mean study baseline prevalence pyrethroid-PBO interaction",
                                                 "mean study baseline prevalence pyrethroid-pyrrole interaction"
                                                 ),
                        "values" = c(get_quantiles(model_fit = m2_fit_full, param = "kappa", dim = NULL),
                                     get_quantiles(model_fit = m2_fit_full, param = "kappa_l", dim = 1),
                                     get_quantiles(model_fit = m2_fit_full, param = "kappa_l", dim = 2),
                                     get_quantiles(model_fit = m2_fit_full, param = "gamma", dim = NULL),
                                     get_quantiles(model_fit = m2_fit_full, param = "delta_l", dim = 1),
                                     get_quantiles(model_fit = m2_fit_full, param = "delta_l", dim = 2),
                                     get_quantiles(model_fit = m2_fit_full, param = "delta", dim = NULL),
                                     get_quantiles(model_fit = m2_fit_full, param = "omega", dim = NULL),
                                     get_quantiles(model_fit = m2_fit_full, param = "omega_l", dim = 1),
                                     get_quantiles(model_fit = m2_fit_full, param = "omega_l", dim = 2)
                                     )
                        )
                      )

write.csv(m2_params_df, file = "m2_params_df.csv")

# ##########################
# ##### data wrangling #####
# ##########################
#
# COMBO_stan <- COMBO_stan |>
#   left_join(data_limits[,c("study", "i", "mean_BL_prev", "mean_BL_net_use")]) |>
#   mutate(diff_base_prev = BL_prev - mean_BL_prev,
#          diff_base_net = net_use - mean_BL_net_use,
#          BL_prev_group = cut(diff_base_prev, breaks = seq(-1, 1, 0.2)),
#          it = gsub("_", " Q", it |> str_to_title()),
#          study_place = case_when(
#            str_detect(study, "Protopopoff") ~ "Tanzania (2015)",
#            str_detect(study, "Staedke") ~ "Uganda (2017)",
#            str_detect(study, "Mosha") ~ "Tanzania (2019)",
#            str_detect(study, "Accrombessi") ~ "Benin (2020)"
#          ))
#
# COMBO_limits_t <- COMBO_stan |> group_by(study) |>
#   summarise(lt = min(time)/12, ut = max(time)/12)
#
# COMBO_limits_prev <- unique(COMBO_stan[,c("cluster", "study", "net", "diff_base_prev")]) |> group_by(study, net) |>
#   summarise(lbp = quantile(diff_base_prev, probs = 0.1), ubp = quantile(diff_base_prev, probs = 0.9))
#
# COMBO_limits_net <- unique(COMBO_stan[,c("cluster", "study", "net", "diff_base_net")]) |> group_by(study, net) |>
#   summarise(lbn = quantile(diff_base_net, probs = 0.1), ubn = quantile(diff_base_net, probs = 0.9))
#
# #####################
# ##### GLM plots #####
# #####################
#
# plot_labels_fun <- function(df){
#   out <- df |> mutate(it = gsub("_", " Q", it |> str_to_title()),
#                       study_place = case_when(
#                         str_detect(study, "Protopopoff") ~ "Tanzania (2015)",
#                         str_detect(study, "Staedke") ~ "Uganda (2017)",
#                         str_detect(study, "Mosha") ~ "Tanzania (2019)",
#                         str_detect(study, "Accrombessi") ~ "Benin (2020)"
#                       )#,
#                       #l_e_la = u_p - start_prev,
#                       #m_e_la = m_p - start_prev,
#                       #u_e_la = l_p - start_prev
#   )
#
#   out$study_place <- factor(out$study_place, levels = c("Tanzania (2015)", "Uganda (2017)", "Tanzania (2019)", "Benin (2020)"))
#
#   return(out)
# }
#
# pred_t <- plot_labels_fun(pred_t)
# pred_m <- plot_labels_fun(pred_m)
# pred_p <- plot_labels_fun(pred_p)
# pred_mn <- plot_labels_fun(pred_mn)
# pred_dn <- plot_labels_fun(pred_dn)
#
# cols <- unname(palette.colors(palette = "Okabe-Ito")[c(1,2,7,3)])
# cols_net <- c("blue", "aquamarine", "darkgreen")
#
# ######################
# ##### time plots #####
# ######################
#
# # prevalence
# m3_p_t_plot_data <- pred_t |>
#   left_join(COMBO_limits_t) |>
#   subset(it %in% c("Protopopoff Q2", "Mosha Q1", "Staedke Q1", "Accrombessi Q3") &
#            diff_net == 0 & diff_prev == 0 &
#            years >= lt & years <= ut)
#
# m3_p_t_plot <- ggplot(data = m3_p_t_plot_data) +
#   geom_ribbon(aes(x = years, ymin = l_p, ymax = u_p, y = m_p, group = interaction(study_place, net), fill = net),
#               alpha = 0.1) +
#   geom_line(linewidth = 1,
#             aes(x = years, y = m_p, group = interaction(study_place, net), col = net)) +
#   geom_boxplot(data = unique(COMBO_stan[,c("cluster", "study_place", "BL_prev")]) |>
#                  mutate(years = 0),
#                aes(x = years,
#                    y = BL_prev,
#                    group = study_place),
#                alpha = 0.1, width = 0.25,
#                position = "identity") +
#   geom_point(data = unique(COMBO_stan[,c("study_place", "mean_BL_prev")]) |> mutate(years = 0),
#              aes(x = years, y = mean_BL_prev), size = 3) +
#   facet_wrap(~ factor(study_place, levels = c("Tanzania (2018)", "Uganda (2020)", "Tanzania (2022)", "Benin (2023)")),
#              nrow = 1) +
#   scale_linetype_manual(name = "", values = c(1, 2)) +
#   ylab("Malaria prevalence") + xlab("Years post trial net distribution") +
#   scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2), labels = scales::percent_format()) +
#   scale_x_continuous(limits = c(-0.25, 3), breaks = 0:3) +
#   scale_colour_manual(values = cols_net, name = "ITN") +
#   scale_fill_manual(values = cols_net, name = "ITN")
#
# # efficacy
#
# m3_e_t_plot <- ggplot(data = m3_p_t_plot_data |> filter(net != "Pyrethroid-only"),
#                       aes(x = years, ymin = l_e, ymax = u_e, y = m_e, fill = study_place)) +
#   geom_ribbon(data = unique(subset(m3_p_t_plot_data, "it" == "Mosha Q1")[,c("l", "net", "it", "years", "l_e_pooled", "u_e_pooled", "m_e_pooled")]),
#               inherit.aes = FALSE, aes(x = years, ymin = l_e_pooled, ymax = u_e_pooled, y = m_e_pooled), alpha = 0.1) +
#   geom_line(data = unique(subset(m3_p_t_plot_data, "it" == "MoshaQ1")[,c("l", "net", "it", "years", "l_e_pooled", "u_e_pooled", "m_e_pooled")]),
#             inherit.aes = FALSE, aes(x = years, y = m_e_pooled), alpha = 0.1) +
#   geom_ribbon(alpha = 0.1) +
#   geom_line(aes(col = study_place), linewidth = 1) +
#   geom_hline(yintercept = 0, col = "black", linetype = 2) +
#   facet_wrap(~ net, nrow = 1) +
#   ylab("Relative efficacy") + xlab("Years post trial net distribution") +
#   scale_y_continuous(limits = c(-0.4, 0.75), breaks = seq(-0.25, 0.75, 0.25), labels = scales::percent_format()) +
#   scale_x_continuous(limits = c(-0.25, 3), breaks = 0:3) +
#   scale_colour_manual(values = cols, name = "study") +
#   scale_fill_manual(values = cols, name = "study")
#
# m3_e_a_t_plot <- ggplot(data = m3_p_t_plot_data |> filter(net != "Pyrethroid-only"),
#                         aes(x = years, y = m_e_a, group = interaction(study_place, net),
#                             ymin = l_e_a, ymax = u_e_a, fill = study_place)) +
#   geom_ribbon(alpha = 0.1) +
#   geom_line(aes(col = study_place), linewidth = 1) +
#   geom_hline(yintercept = 0, col = "black", linetype = 2) +
#   facet_wrap(~ net, nrow = 1) +
#   ylab("Absolute efficacy") + xlab("Years post trial net distribution") +
#   scale_y_continuous(limits = c(-0.1, 0.3), breaks = seq(-0.1, 0.3, 0.1), labels = scales::percent_format()) +
#   scale_x_continuous(limits = c(-0.25, 3), breaks = 0:3) +
#   scale_colour_manual(values = cols, name = "study") +
#   scale_fill_manual(values = cols, name = "study")
#
# # odds ratios
# m3_or_t_plot <- ggplot(data = m3_p_t_plot_data |> filter(net != "Pyrethroid-only"),
#                        aes(x = years, ymin = l_or, ymax = u_or, y = m_or, fill = study_place)) +
#   geom_ribbon(alpha = 0.1) +
#   geom_line(aes(col = study_place), linewidth = 1) +
#   geom_hline(yintercept = 1, col = "black", linetype = 2) +
#   facet_wrap(~ net, nrow = 1) +
#   ylab("Odds ratios") + xlab("Years post trial net distribution") +
#   scale_y_continuous(limits = c(0, 1.5), breaks = seq(0, 1.5, 0.25)) +
#   scale_x_continuous(limits = c(-0.25, 3), breaks = 0:3) +
#   scale_colour_manual(values = cols, name = "study") +
#   scale_fill_manual(values = cols, name = "study") +
#   geom_line(inherit.aes = FALSE, aes(x = years, y = m_or_pooled)) +
#   geom_ribbon(inherit.aes = FALSE, aes(x = years, y = m_or_pooled, ymin = l_or_pooled, ymax = u_or_pooled), alpha = 0.1)
#
# ggsave(
#   file = "regression_t_plot.pdf",
#   device = "pdf",
#   plot = m3_p_t_plot / m3_or_t_plot /
#     m3_e_t_plot / m3_e_a_t_plot +
#     plot_layout(guides = "collect",
#                 axes = "collect",
#                 axis_titles = "collect") +
#     plot_annotation(tag_levels = 'A'),
#   width = 35, height = 35,
#   units = "cm"
# )
#
# ##########################################################
# ##### baseline prevalence and baseline net use plots #####
# ##########################################################
#
# # between study differences
# # mean baseline prevalence
#
# m3_m_o_r_plot <- ggplot(data = pred_m |>
#                           subset(net != "Pyrethroid-only" &
#                                    it %in% c("Protopopoff Q2", "Mosha Q1", "Staedke Q1", "Accrombessi Q3") &
#                                    mean_prev >= min(unique(COMBO_stan$mean_BL_prev)) & mean_prev <= max(unique(COMBO_stan$mean_BL_prev))),
#                         aes(x = mean_prev, y = m_or, ymin = l_or, ymax = u_or,
#                             fill = study_place, group = interaction(study_place, net))) +
#   geom_ribbon(alpha = 0.1) +
#   geom_line(aes(col = study_place), linewidth = 1) +
#   geom_hline(yintercept = 1, linetype = 2) +
#   facet_wrap(~ net) +
#   ylab("Odds ratio") + xlab("Between study mean baseline prevalence") +
#   scale_y_continuous(limits = c(0, 2), breaks = seq(0, 2.5, 0.5)) +
#   scale_colour_manual(values = cols, name = "study") +
#   scale_fill_manual(values = cols, name = "study") +
#   scale_x_continuous(labels = scales::percent)
#
# # mean baseline net use
# m3_mn_o_r_plot <- ggplot(data = pred_mn |>
#                            subset(net != "Pyrethroid-only" &
#                                     it %in% c("Protopopoff Q2", "Mosha Q1", "Staedke Q1", "Accrombessi Q3") &
#                                     mean_net >= min(unique(COMBO_stan$mean_BL_net_use)) & mean_net <= max(unique(COMBO_stan$mean_BL_net_use))),
#                          aes(x = mean_net, y = m_or, ymin = l_or, ymax = u_or,
#                              fill = study_place, group = interaction(study_place, net))) +
#   geom_ribbon(alpha = 0.1) +
#   geom_line(aes(col = study_place), linewidth = 1) +
#   geom_hline(yintercept = 1, linetype = 2) +
#   facet_wrap(~ net) +
#   ylab("Odds ratio") + xlab("Between study mean baseline net use") +
#   scale_y_continuous(limits = c(0, 4), breaks = seq(0, 4, 1)) +
#   scale_colour_manual(values = cols, name = "study") +
#   scale_fill_manual(values = cols, name = "study") +
#   scale_x_continuous(labels = scales::percent)
#
# ggsave(file = "regression_mean_prev_net.pdf",
#        plot = m3_m_o_r_plot / m3_mn_o_r_plot +
#          plot_annotation(tag_levels = 'A') +
#          plot_layout(guides = "collect",
#                      axis_titles = "collect"),
#        device = "pdf",
#        width = 35, height = 25,
#        units = "cm"
# )
#
# ##### within-study differences in baseline prevalence
#
# m3_pred_p_data <- pred_p |>
#   subset(it %in% c("Protopopoff Q2", "Mosha Q1", "Staedke Q1", "Accrombessi Q3")) |>
#   left_join(COMBO_limits_prev) |>
#   filter(diff_prev >= lbp & diff_prev <= ubp)
#
# m3_dp_plot <- ggplot(data = m3_pred_p_data,
#                      aes(x = start_prev, y = m_p, ymin = l_p, ymax = u_p, fill = net, group = interaction(study_place, net))) +
#   geom_ribbon(alpha = 0.1) +
#   geom_line(aes(col = net), linewidth = 1) +
#   facet_wrap(~ study_place, nrow = 1) +
#   ylab("Malaria prevalence") + xlab("Baseline prevalence") +
#   scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2), labels = scales::percent) +
#   scale_colour_manual(values = cols_net, name = "ITN") +
#   scale_fill_manual(values = cols_net, name = "ITN") +
#   scale_x_continuous(limits = c(0, 1), breaks = seq(-0, 1, 0.25), labels = scales::percent)
#
# # efficacy
# m3_e_dp_plot <- ggplot(data = m3_pred_p_data |> filter(net != "Pyrethroid-only") |>
#                          left_join(COMBO_limits_prev |> filter(net == "Pyrethroid-only") |> rename(lbpp = lbp, ubpp = ubp) |>
#                                      select(study, lbpp, ubpp),
#                                    by = c("study")) |>
#                          filter(diff_prev >= lbpp & diff_prev <= ubpp),
#                        aes(x = start_prev, ymin = l_e, ymax = u_e, y = m_e, fill = study_place)) +
#   geom_ribbon(alpha = 0.1) +
#   geom_line(aes(col = study_place), linewidth = 1) +
#   geom_hline(yintercept = 0, col = "black", linetype = 2) +
#   facet_wrap(~ net, nrow = 1) +
#   ylab("Relative efficacy compared\nto pyrethroid-only ITNs") + xlab("Baseline prevalence") +
#   scale_y_continuous(limits = c(-1.1, 0.8), breaks = seq(-1, 0.75, 0.25), labels = scales::percent_format()) +
#   scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2), labels = scales::percent) +
#   scale_colour_manual(values = cols, name = "study") +
#   scale_fill_manual(values = cols, name = "study") +
#   coord_cartesian(ylim = c(-0.25, 0.8))
#
# m3_e_a_dp_plot <- ggplot(data = m3_pred_p_data |> filter(net != "Pyrethroid-only") |>
#                            left_join(COMBO_limits_prev |> filter(net == "Pyrethroid-only") |> rename(lbpp = lbp, ubpp = ubp) |>
#                                        select(study, lbpp, ubpp),
#                                      by = c("study")) |>
#                            filter(diff_prev >= lbpp & diff_prev <= ubpp),
#                          aes(x = start_prev, ymin = l_e_a, ymax = u_e_a, y = m_e_a, fill = study_place)) +
#   geom_ribbon(alpha = 0.1) +
#   geom_line(aes(col = study_place), linewidth = 1) +
#   geom_hline(yintercept = 0, col = "black", linetype = 2) +
#   facet_wrap(~ net, nrow = 1) +
#   ylab("Absolute efficacy compared\nto pyrethroid-only ITNs") + xlab("Baseline prevalence") +
#   scale_y_continuous(limits = c(-0.175, 0.4), breaks = seq(-0.1, 0.4, 0.1), labels = scales::percent_format()) +
#   scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2), labels = scales::percent) +
#   scale_colour_manual(values = cols, name = "study") +
#   scale_fill_manual(values = cols, name = "study")
#
# m3_o_r_dp_plot <- ggplot(data = m3_pred_p_data |> filter(net != "Pyrethroid-only") |>
#                            left_join(COMBO_limits_prev |> filter(net == "Pyrethroid-only") |> rename(lbpp = lbp, ubpp = ubp) |>
#                                        select(study, lbpp, ubpp),
#                                      by = c("study")) |>
#                            filter(diff_prev >= lbpp & diff_prev <= ubpp),
#                          aes(x = start_prev, ymin = l_or, ymax = u_or, y = m_or, fill = study_place)) +
#   geom_ribbon(alpha = 0.1) +
#   geom_line(aes(col = study_place), linewidth = 1) +
#   geom_hline(yintercept = 1, col = "black", linetype = 2) +
#   facet_wrap(~ net, nrow = 1) +
#   ylab("Odds ratio") + xlab("Baseline prevalence") +
#   #scale_y_continuous(limits = c(-0.175, 0.4), breaks = seq(-0.1, 0.4, 0.1), labels = scales::percent_format()) +
#   scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2), labels = scales::percent) +
#   scale_colour_manual(values = cols, name = "study") +
#   scale_fill_manual(values = cols, name = "study") +
#   geom_line(inherit.aes = FALSE, aes(x = start_prev, y = m_or_pooled)) +
#   geom_ribbon(inherit.aes = FALSE, aes(x = start_prev, y = m_or_pooled, ymin = l_or_pooled, ymax = u_or_pooled), alpha = 0.1)
#
#
# ggsave(
#   plot = m3_dp_plot /
#     m3_o_r_dp_plot /
#     m3_e_dp_plot / m3_e_a_dp_plot +
#     plot_layout(guides = "collect",
#                 axes = "collect",
#                 axis_titles = "collect") +
#     plot_annotation(tag_levels = 'A'),
#   file = "regression_diff_p_plot.pdf",
#   device = "pdf",
#   width = 35, height = 35,
#   units = "cm"
# )
#
# # odds ratios
# m3_pred_dn_data <- pred_dn |>
#   subset(it %in% c("Protopopoff Q2", "Mosha Q1", "Staedke Q1", "Accrombessi Q3")) |>
#   left_join(COMBO_limits_net) |>
#   filter(diff_net >= lbn & diff_net <= ubn) |>
#   mutate(start_net = diff_net + mean_net)
#
# m3_d_n_plot <- ggplot(data = m3_pred_dn_data,
#                       aes(x = start_net, y = m_p, ymin = l_p, ymax = u_p,
#                           fill = net, group = interaction(study_place, net))) +
#   geom_ribbon(alpha = 0.1) +
#   geom_line(aes(col = net), linewidth = 1) +
#   facet_wrap(~ study_place, nrow = 1) +
#   ylab("Malaria prevalence") + xlab("Baseline net use") +
#   scale_y_continuous(limits = c(0, 0.85), breaks = seq(0, 0.8, 0.2), labels = scales::percent) +
#   scale_colour_manual(values = cols_net, name = "ITN") +
#   scale_fill_manual(values = cols_net, name = "ITN") +
#   scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25), labels = scales::percent)
#
# m3_d_or_n_plot <- ggplot(data = m3_pred_dn_data |>
#                            filter(net != "Pyrethroid-only") |>
#                            left_join(subset(COMBO_limits_net, net == "Pyrethroid-only") |>
#                                        filter(net == "Pyrethroid-only") |> rename(lbn_p = lbn, ubn_p = ubn) |>
#                                        select(-net),
#                                      by = c("study")) |>
#                            filter(diff_net >= lbn_p & diff_net <= ubn_p),
#                          aes(x = start_net, y = m_or, ymin = l_or, ymax = u_or,
#                              fill = study_place, group = interaction(study, net))) +
#   geom_ribbon(alpha = 0.1) +
#   geom_line(aes(col = study_place), linewidth = 1) +
#   facet_wrap(~ net) +
#   geom_hline(yintercept = 1, linetype = 2) +
#   facet_wrap(~ net) +
#   ylab("Odds ratio") + xlab("Baseline net use") +
#   scale_y_continuous(limits = c(0, 1.5), breaks = seq(0, 1.5, 0.25)) +
#   scale_colour_manual(values = cols, name = "study") +
#   scale_fill_manual(values = cols, name = "study") +
#   scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25), labels = scales::percent)
#
# m3_d_e_n_plot <- ggplot(data = m3_pred_dn_data |>
#                           filter(net != "Pyrethroid-only") |>
#                           left_join(subset(COMBO_limits_net, net == "Pyrethroid-only") |>
#                                       filter(net == "Pyrethroid-only") |> rename(lbn_p = lbn, ubn_p = ubn) |>
#                                       select(-net),
#                                     by = c("study")) |>
#                           filter(diff_net >= lbn_p & diff_net <= ubn_p),
#                         aes(x = start_net, y = m_e, ymin = l_e, ymax = u_e,
#                             fill = study_place, group = interaction(study, net))) +
#   geom_ribbon(alpha = 0.1) +
#   geom_line(aes(col = study_place), linewidth = 1) +
#   facet_wrap(~ net, nrow = 1) +
#   geom_hline(yintercept = 0, linetype = 2) +
#   ylab("Relative efficacy compared\nto pyrethroid-only ITNs") + xlab("Baseline net use") +
#   scale_y_continuous(limits = c(-0.25, 1), breaks = seq(-0.25, 1, 0.25), labels = scales::percent) +
#   scale_colour_manual(values = cols, name = "study") +
#   scale_fill_manual(values = cols, name = "study") +
#   scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25), labels = scales::percent)
#
# m3_d_ea_n_plot <- ggplot(data = m3_pred_dn_data |>
#                            filter(net != "Pyrethroid-only") |>
#                            left_join(subset(COMBO_limits_net, net == "Pyrethroid-only") |>
#                                        filter(net == "Pyrethroid-only") |> rename(lbn_p = lbn, ubn_p = ubn) |>
#                                        select(-net),
#                                      by = c("study")) |>
#                            filter(diff_net >= lbn_p & diff_net <= ubn_p),
#                          aes(x = start_net, y = m_e_a, ymin = l_e_a, ymax = u_e_a,
#                              fill = study_place, group = interaction(study, net))) +
#   geom_ribbon(alpha = 0.1) +
#   geom_line(aes(col = study_place), linewidth = 1) +
#   facet_wrap(~ net, nrow = 1) +
#   geom_hline(yintercept = 0, linetype = 2) +
#   ylab("Absolute efficacy compared\nto pyrethroid-only ITNs") + xlab("Baseline net use") +
#   scale_y_continuous(limits = c(-0.05, 0.3), breaks = seq(0, 0.3, 0.1), labels = scales::percent) +
#   scale_colour_manual(values = cols, name = "study") +
#   scale_fill_manual(values = cols, name = "study") +
#   scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25), labels = scales::percent)
#
# ggsave(
#   plot = m3_d_n_plot / m3_d_or_n_plot /
#     m3_d_e_n_plot / m3_d_ea_n_plot +
#     plot_layout(guides = "collect",
#                 axes = "collect",
#                 axis_titles = "collect") +
#     plot_annotation(tag_levels = 'A'),
#   file = "regression_diff_n_plot.pdf",
#   device = "pdf",
#   width = 35, height = 35,
#   units = "cm"
# )
#
# ############################
# ##### model comparison #####
# ############################
#
# orderly2::orderly_dependency("3_cv",
#                              quote(latest(parameter:model == "m0_c")),
#                              "m0_c_log_pd_kfold.rds")
#
# orderly2::orderly_dependency("3_cv",
#                              quote(latest(parameter:model == "m3_re")),
#                              "m3_re_log_pd_kfold.rds")
#
# m0_c_log_pd_kfold <- readRDS("m0_c_log_pd_kfold.rds")
# m3_re_log_pd_kfold <- readRDS("m3_re_log_pd_kfold.rds")
#
# # checking the fits
# # leave one out cross validation
# m0_c_elpd_kfold <- elpd(m0_c_log_pd_kfold)
# m3_re_elpd_kfold <- elpd(m3_re_log_pd_kfold)
#
# loo_compare(m0_c_elpd_kfold, m3_re_elpd_kfold)
#
# # actual vs fitted plot
# extract_af <- function(fit_full){
#   cbind(rstan::extract(fit_full, "inv_logit_posterior")$inv_logit_posterior |>
#           apply(2, quantile, probs = c(lower, 0.5, upper)) |>
#           t() |> as.data.frame() |> rename(l_p = 1, m_p = 2, u_p = 3),
#         COMBO_stan[,c("cluster", "time", "seasonality", "l", "li", "study", "net", "prev", "prev_lower", "prev_upper")])
# }
#
# m0_af <- extract_af(m0_fit_full) |> mutate(model = "standard")
# m3_af <- extract_af(fit_full) |> mutate(model = "covariate-adjusted")
#
# af_plot <- ggplot(data = rbind(m0_af, m3_af),
#                   aes(x = prev, xmin = prev_lower, xmax = prev_upper, y = m_p, ymin = l_p, ymax = u_p, fill = model)) +
#   geom_errorbarh(alpha = 0.1, linewidth = 0.1) +
#   geom_pointrange(shape = 21, alpha = 0.5, size = 0.7) +
#   geom_abline(slope = 1, linetype = 2) +
#   scale_y_continuous(labels = scales::percent, limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
#   scale_x_continuous(labels = scales::percent) +
#   xlab("Sampled prevalence") + ylab("Predicted prevalence") +
#   scale_fill_manual(name = "", values = c("#0072B2", "#D55E00")) +
#   guides(fill  = guide_legend(position = "inside")) +
#   theme(legend.position.inside = c(0.9, 0.25))
#
# extract_yrep <- function(fit_full, model_name){
#   rstan::extract(fit_full, "y_rep")$y_rep |> t() |>
#     as.data.frame() |>
#     cbind(COMBO_stan[,c("cluster", "time", "seasonality", "l", "i", "li", "study", "net", "it")]) |>
#     pivot_longer(cols = 1:16000, names_to = "iter", values_to = "yrep") |>
#     mutate(model = model_name,
#            iter = as.numeric(gsub("V", "", iter)))
# }
#
# m0_yrep <- extract_yrep(m0_fit_full, model_name = "standard")
# m3_yrep <- extract_yrep(fit_full, model_name = "covariate-adjusted")
#
# i_s <- sample(1:16000, 50, replace = FALSE)
# yrep_plot <- rbind(m0_yrep |> filter(iter %in% i_s),
#                    m3_yrep |> filter(iter %in% i_s))
#
# ppd_plot <- ggplot() +
#   geom_histogram(data = COMBO_stan, aes(x = prev_num), position = "identity", binwidth = 1, col = "black") +
#   geom_histogram(data = yrep_plot,
#                  aes(x = yrep, group = interaction(iter, model), col = model, fill = model), linewidth = 0.05,
#                  alpha = 0.025,
#                  position = "identity",
#                  binwidth = 1) +
#   scale_colour_manual(name = "", values = c("#0072B2", "#D55E00")) +
#   scale_fill_manual(name = "", values = c("#0072B2", "#D55E00")) +
#   xlab("Number positive for malaria") +
#   ylab("Frequency") +
#   facet_wrap(~model) +
#   theme(legend.position = "none")
#
# # random effects variance
# extract_sd <- function(model_fit, model_name){
#   out <- rstan::extract(model_fit, "sigma_e_r_train")$sigma_e_r |> as.data.frame()
#   colnames(out) <- unique(COMBO_stan[,c("i", "study")])[,2]
#   return(out |> pivot_longer(cols = 1:4, names_to = "study", values_to = "sd") |> mutate(model = model_name))
# }
#
# s_all <- rbind(extract_sd(m0_fit_full, "standard"),
#                extract_sd(fit_full, "covariate-adjusted")) |>
#   mutate(study_place = case_when(
#     str_detect(study, "Protopopoff") ~ "Tanzania (2018)",
#     str_detect(study, "Staedke") ~ "Uganda (2020)",
#     str_detect(study, "Mosha") ~ "Tanzania (2022)",
#     str_detect(study, "Accrombessi") ~ "Benin (2023)"
#   )
#   )
#
# sd_plot <- ggplot(data = s_all, aes(x = sd, fill = model, group = interaction(study_place, model))) +
#   geom_density(alpha = 0.5) +
#   facet_wrap(~study_place) +
#   xlab(expression(paste("Study-specific cluster random effect standard deviation (", sigma["i"], ")"))) +
#   ylab("Density") +
#   scale_fill_manual(name = "", values = c("#0072B2", "#D55E00")) +
#   guides(fill  = guide_legend(position = "inside")) +
#   theme(legend.position.inside = c(0.9, 0.25))
#
# ggsave(
#   file = "model_checks.pdf",
#   device = "pdf",
#   plot = af_plot / ppd_plot / sd_plot  +
#     plot_annotation(tag_levels = 'A'),
#   width = 20, height = 35,
#   units = "cm"
# )
#
#
#
# params <- get_params(fit = fit_full)$params
# values <- get_params(fit = fit_full)$values
#
# # seasonality
# params <- c(params, "Benin (2023) Q3", "Tanzania (2022) Q1", "Tanzania (2022) Q3", "Tanzania (2018) Q2", "Tanzania (2018) Q4",
#             "Uganda (2020) Q1", "Uganda (2020) Q2", "Uganda (2020) Q3", "Uganda (2020) Q4")
#
# for(i in 1:9){
#   values <- c(values, get_quantiles(fit_full, "beta_it", dim = i))
# }
#
# params <- c(params, "time (years)")
# values <- c(values, get_quantiles(fit_full, "kappa", dim = NULL))
#
# #
# params <- c(params, "pyrethroid-PBO-time interaction", "pyrethroid-pyrrole-time interaction")
# values <- c(values, get_quantiles(fit_full, "kappa_l", dim = 1), get_quantiles(fit_full, "kappa_l", dim = 2))
#
# params <- c(params, "pooled within-study differences in baseline prevalence", "pooled mean study baseline prevalence",
#             "within-study differences in baseline prevalence-time interaction")
# values <- c(values, get_quantiles(fit_full, "gamma", dim = NULL), get_quantiles(fit_full, "omega", dim = NULL), get_quantiles(fit_full, "delta", dim = NULL))
#
# params <- c(params, "pooled within-study differences in baseline ITN use", "pooled mean study baseline ITN use")
# values <- c(values, get_quantiles(fit_full, "tau_d", dim = NULL), get_quantiles(fit_full, "tau_m", dim = NULL))
#
# params <- c(params, "pyrethroid-PBO-within study differences in baseline prevalence interaction", "pyrethroid-pyrrole-within study differences in baseline prevalence interaction")
# values <- c(values, get_quantiles(fit_full, "delta_l", dim = 1), get_quantiles(fit_full, "delta_l", dim = 2))
#
# params <- c(params, "pyrethroid-PBO-mean study baseline prevalence interaction", "pyrethroid-pyrrole-mean study baseline prevalence interaction")
# values <- c(values, get_quantiles(fit_full, "omega_l", dim = 1), get_quantiles(fit_full, "omega_l", dim = 2))
#
# params <- c(params, "pyrethroid-PBO-within study differences in baseline ITN use interaction", "pyrethroid-pyrrole-within study differences in baseline ITN use interaction")
# values <- c(values, get_quantiles(fit_full, "tau_ld", dim = 1), get_quantiles(fit_full, "tau_ld", dim = 2))
#
# params <- c(params, "pyrethroid-PBO-mean study baseline ITN use interaction", "pyrethroid-pyrrole-mean study baseline ITN use interaction")
# values <- c(values, get_quantiles(fit_full, "tau_lm", dim = 1), get_quantiles(fit_full, "tau_lm", dim = 2))
#
# write.csv(data.frame("params" = params, "values" = values), file = "m3_params.csv")
