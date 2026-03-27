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

loo_comp <- round(loo_compare(m0_c_elpd_kfold, m1_re_elpd_kfold, m2_re_elpd_kfold), digits = 1)

elpd_table <- data.frame("model" = c("model 0", "model 1", "model 2"),
                         "elpd" = c(loo_comp["model1", "elpd"], loo_comp["model2", "elpd"], loo_comp["model3", "elpd"]),
                         "se_elpd" = c(loo_comp["model1", "se_elpd"], loo_comp["model2", "se_elpd"], loo_comp["model3", "se_elpd"]),
                         "elpd_diff" = c(loo_comp["model1", "elpd_diff"], loo_comp["model2", "elpd_diff"], loo_comp["model3", "elpd_diff"]))

write.csv(elpd_table, "elpd_table.csv")

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

or_mean_3y_diff <- or_mean_3y_diff |> mutate(study_place = case_when(str_detect(study, "Protopopoff") ~ "Tanzania (2015)",
                                                                     str_detect(study, "Staedke") ~ "Uganda (2017)",
                                                                     str_detect(study, "Mosha") ~ "Tanzania (2019)",
                                                                     str_detect(study, "Accrombessi") ~ "Benin (2020)"))
or_mean_3y_diff$study_place <- factor(or_mean_3y_diff$study_place, levels = c("Tanzania (2015)", "Uganda (2017)", "Tanzania (2019)", "Benin (2020)"))

or_mean_3y_mean <- or_mean_3y_mean |> mutate(study_place = case_when(str_detect(study, "Protopopoff") ~ "Tanzania (2015)",
                                                                     str_detect(study, "Staedke") ~ "Uganda (2017)",
                                                                     str_detect(study, "Mosha") ~ "Tanzania (2019)",
                                                                     str_detect(study, "Accrombessi") ~ "Benin (2020)"))

or_mean_3y_mean$study_place <- factor(or_mean_3y_mean$study_place, levels = c("Tanzania (2015)", "Uganda (2017)", "Tanzania (2019)", "Benin (2020)"))

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

m2_dp_plot <- ggplot(data = or_mean_3y_diff, #pred_diff,
                     aes(x = start_prev, y = mean_prev_3_m, ymin = mean_prev_3_l, ymax = mean_prev_3_u, # y = m_p, ymin = l_p, ymax = u_p,
                         fill = net, group = interaction(study_place, net))) + #, year_plot
  geom_ribbon(alpha = 0.1) +
  geom_line(aes(col = net), linewidth = 1) +#, linetype = year_plot
  facet_wrap(~ study_place, nrow = 1) +
  ylab("Malaria prevalence") + xlab("Baseline prevalence") +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2), labels = scales::percent) +
  scale_colour_manual(values = cols_net, name = "ITN") +
  scale_fill_manual(values = cols_net, name = "ITN") +
  scale_x_continuous(limits = c(0, 1), breaks = seq(-0, 1, 0.25), labels = scales::percent) +
  #scale_linetype_manual(values = c(1, 2, 3), name = "") +
  labs(tag = "A")

m2_o_r_dp_plot <- ggplot(data = or_mean_3y_diff |> filter(net != "Pyrethroid-only"), #pred_diff
                         aes(x = diff_prev, y = mean_or_3_m, ymin = mean_or_3_l, ymax = mean_or_3_u, #ymin = l_or, ymax = u_or, y = m_or,
                             fill = study_place, group = interaction(study_place, net))) +#, year_plot
  geom_ribbon(inherit.aes = FALSE,
              aes(x = diff_prev, y = pooled_mean_or_3_m, ymin = pooled_mean_or_3_l, ymax = pooled_mean_or_3_u, #y = m_or_pooled, ymin = l_or_pooled, ymax = u_or_pooled,
                  group = interaction(study_place)), alpha = 0.1) + #, year_plot
  geom_ribbon(alpha = 0.1) +
  geom_line(aes(col = study_place), linewidth = 1) +#, linetype = year_plot
  geom_hline(yintercept = 0, col = "black", linetype = 2) +
  facet_wrap(~ net, nrow = 1) +
  ylab("Log odds ratio of infection in trial ITN\nclusters relative to pyrethroid-only clusters") +
  xlab("Within-study differences in baseline prevalence") +
  scale_x_continuous(limits = c(-0.55, 0.55), breaks = seq(-0.5, 0.5, 0.25), labels = scales::percent) +
  scale_colour_manual(values = cols, name = "study") +
  scale_fill_manual(values = cols, name = "study") +
  geom_line(inherit.aes = FALSE, aes(x = diff_prev, y = pooled_mean_or_3_m), linewidth = 1) +#, linetype = year_plot # m_or_pooled
  #scale_linetype_manual(values = c(1, 2, 3), name = "") +
  scale_y_continuous(breaks = seq(-4, 2, 1)) +
  labs(tag = "B")

m2_eff_dp_plot <- ggplot(data = or_mean_3y_diff |> filter(net != "Pyrethroid-only"),
                         aes(x = diff_prev, y = mean_eff_3_m / 100, ymin = mean_eff_3_l / 100, ymax = mean_eff_3_u / 100,
                             fill = study_place, group = interaction(study_place, net))) +
  geom_ribbon(alpha = 0.1) +
  geom_line(aes(col = study_place), linewidth = 1) +
  geom_hline(yintercept = 0, col = "black", linetype = 2) +
  facet_wrap(~ net, nrow = 1) +
  ylab("Relative efficacy of trial ITN clusters\nrelative to pyrethroid-only clusters") +
  xlab("Within-study differences in baseline prevalence") +
  scale_x_continuous(limits = c(-0.55, 0.55), breaks = seq(-0.5, 0.5, 0.25), labels = scales::percent) +
  scale_y_continuous(labels = scales::percent, breaks = seq(-2, 1, 0.5)) +
  scale_colour_manual(values = cols, name = "study") +
  scale_fill_manual(values = cols, name = "study") +
  coord_cartesian(ylim = c(-2, 1)) +
  labs(tag = "C")

m2_mean_BL_prev_o_r_plot <- ggplot(data = or_mean_3y_mean |> #pred_mean
                                     left_join(COMBO_stan |> group_by(net) |> summarise(min_mean_BL_prev = min(mean_BL_prev),
                                                                                        max_mean_BL_prev = max(mean_BL_prev))) |>
                                     subset(net != "Pyrethroid-only" & mean_prev >= min_mean_BL_prev & mean_prev <= max_mean_BL_prev),
                        aes(x = mean_prev, y = mean_or_3_m, ymin = mean_or_3_l, ymax = mean_or_3_u,
                            fill = study_place, group = interaction(study_place, net))) + #, year_plot
  geom_ribbon(aes(x = mean_prev, ymin = pooled_mean_or_3_l, ymax = pooled_mean_or_3_u,
                  group = interaction(net)), inherit.aes = FALSE, alpha = 0.1) +
  geom_ribbon(alpha = 0.1) +
  geom_line(aes(col = study_place), linewidth = 1) +
  geom_line(aes(x = mean_prev, y = pooled_mean_or_3_m), inherit.aes = FALSE, linewidth = 1) +
  geom_hline(yintercept = 0, linetype = 2) +
  facet_wrap(~net, scales = "free_x") +
  ylab("Log odds ratio of infection in trial ITN\nclusters relative to pyrethroid-only clusters") +
  xlab("Between study mean baseline prevalence") +
  scale_y_continuous(limits = c(-5, 2), breaks = seq(-5, 2, 1)) +
  scale_colour_manual(values = cols, name = "study") +
  scale_fill_manual(values = cols, name = "study") +
  scale_x_continuous(labels = scales::percent) +
  labs(tag = "A")

m2_mean_BL_prev_eff_plot <- ggplot(data = or_mean_3y_mean |>
                                     left_join(COMBO_stan |> group_by(net) |> summarise(min_mean_BL_prev = min(mean_BL_prev),
                                                            max_mean_BL_prev = max(mean_BL_prev))) |>
         subset(net != "Pyrethroid-only" & mean_prev >= min_mean_BL_prev & mean_prev <= max_mean_BL_prev),
       aes(x = mean_prev, y = mean_eff_3_m/100, ymin = mean_eff_3_l/100, ymax = mean_eff_3_u/100,
           fill = study_place, group = interaction(study_place, net))) + #, year_plot
  geom_ribbon(alpha = 0.1) +
  geom_line(aes(col = study_place), linewidth = 1) +
  geom_hline(yintercept = 0, linetype = 2) +
  facet_wrap(~net, scales = "free_x") +
  ylab("Relative efficacy of trial ITN clusters\nrelative to pyrethroid-only clusters") +
  xlab("Between study mean baseline prevalence") +
  scale_y_continuous(breaks = seq(-1.25, 1, 0.25), labels = scales::percent) +
  scale_colour_manual(values = cols, name = "study") +
  scale_fill_manual(values = cols, name = "study") +
  scale_x_continuous(labels = scales::percent) +
  labs(tag = "B")

ggsave(file = "regression_diff_base_prev.pdf",
       plot = m2_dp_plot /
         (m2_o_r_dp_plot + m2_eff_dp_plot +
            plot_layout(guides = "collect")),
       device = "pdf",
       width = 45, height = 20,
       units = "cm"
)

ggsave(file = "regression_mean_base_prev.pdf",
       plot = m2_mean_BL_prev_o_r_plot /
         m2_mean_BL_prev_eff_plot  +
         plot_layout(guides = "collect"),
       device = "pdf",
       width = 27.5, height = 20,
       units = "cm")

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
