library(cmdstanr)
library(posterior)
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

orderly2::orderly_dependency("2_fits", "latest()", c("fit_0.rds", "fit_1.rds", "fit_2.rds"))

m2_fit_full <- readRDS(file = "fit_2.rds")
m1_fit_full <- readRDS(file = "fit_1.rds")
m0_fit_full <- readRDS(file = "fit_0.rds")

orderly2::orderly_dependency("9_inv_logit",
                             "latest()",
                             "m2_pred.rds")

list2env(readRDS(file = "m2_pred.rds"), envir = .GlobalEnv)

orderly2::orderly_shared_resource("m0.stan",
                                  "m1.stan",
                                  "m2.stan",
                                  "model_functions.R")

source("model_functions.R")

m0_stan_model <- cmdstan_model(stan_file = "m0.stan")
m1_stan_model <- cmdstan_model(stan_file = "m1.stan")
m2_stan_model <- cmdstan_model(stan_file = "m2.stan")

############################
##### model comparison #####
############################

orderly2::orderly_dependency("3_cv",
                             quote(latest(parameter:model == "m0")),
                             "m0_log_pd_kfold.rds")

orderly2::orderly_dependency("3_cv",
                             quote(latest(parameter:model == "m1")),
                             "m1_log_pd_kfold.rds")

orderly2::orderly_dependency("3_cv",
                             quote(latest(parameter:model == "m2")),
                             "m2_log_pd_kfold.rds")

m0_log_pd_kfold <- readRDS("m0_log_pd_kfold.rds")
m1_log_pd_kfold <- readRDS("m1_log_pd_kfold.rds")
m2_log_pd_kfold <- readRDS("m2_log_pd_kfold.rds")

# checking the fits
# leave one out cross validation
m0_elpd_kfold <- elpd(m0_log_pd_kfold)
m1_elpd_kfold <- elpd(m1_log_pd_kfold)
m2_elpd_kfold <- elpd(m2_log_pd_kfold)

loo_comp <- round(loo_compare(m0_elpd_kfold, m1_elpd_kfold, m2_elpd_kfold), digits = 1)

elpd_table <- data.frame("model" = c("model 0", "model 1", "model 2"),
                         "elpd" = c(loo_comp["model1", "elpd"], loo_comp["model2", "elpd"], loo_comp["model3", "elpd"]),
                         "se_elpd" = c(loo_comp["model1", "se_elpd"], loo_comp["model2", "se_elpd"], loo_comp["model3", "se_elpd"]),
                         "elpd_diff" = c(loo_comp["model1", "elpd_diff"], loo_comp["model2", "elpd_diff"], loo_comp["model3", "elpd_diff"]))

write.csv(elpd_table, "elpd_table.csv")

# actual vs fitted plot
extract_af <- function(fit_full){
  cbind(fit_full$draws("inv_logit_lp") |> as_draws_matrix() |>
          apply(2, quantile, probs = c(lower, 0.5, upper)) |>
          t() |> as.data.frame() |> rename(l_p = 1, m_p = 2, u_p = 3),
        COMBO_stan[,c("cluster", "time", "l", "li", "study", "net", "prev", "prev_lower", "prev_upper")])
}

m0_af <- extract_af(m0_fit_full) |> mutate(model = "model 0")
m1_af <- extract_af(m2_fit_full) |> mutate(model = "model 1")
m2_af <- extract_af(m2_fit_full) |> mutate(model = "model 2")

af_plot <- ggplot(data = rbind(m0_af, m1_af, m2_af),
                  aes(x = prev, xmin = prev_lower, xmax = prev_upper, y = m_p, ymin = l_p, ymax = u_p, colour = model)) +
  geom_errorbar(alpha = 0.1, linewidth = 0.15, orientation = "y") +
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
  out <- model_fit$draws("sigma_e_r") |> as_draws_matrix() |> as.data.frame()
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

s_all$study_place = factor(s_all$study_place, levels = c("Tanzania (2015)", "Uganda (2017)", "Tanzania (2019)", "Benin (2020)"))

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
              summarise(base_prev_mean = mean(BL_prev))) |>
  mutate(base_prev_diff = BL_prev - base_prev_mean,
         study_place = case_when(
           str_detect(study, "Protopopoff") ~ "Tanzania (2015)",
           str_detect(study, "Staedke") ~ "Uganda (2017)",
           str_detect(study, "Mosha") ~ "Tanzania (2019)",
           str_detect(study, "Accrombessi") ~ "Benin (2020)"
         ))

COMBO_limits_t <- COMBO_stan |> group_by(study) |>
  summarise(lt = min(time)/12, ut = max(time)/12)

COMBO_limits_prev <- unique(COMBO_stan[,c("cluster", "study", "net", "base_prev_diff")]) |> group_by(study, net) |>
  summarise(lbp = min(base_prev_diff), ubp = max(base_prev_diff))

COMBO_limits_BL_prev <- unique(COMBO_stan[,c("cluster", "study", "net", "BL_prev")]) |> group_by(study, net) |>
  summarise(lbp = min(BL_prev), ubp = max(BL_prev))

COMBO_mean_prev <- COMBO_stan |> group_by(cluster, net, BL_prev, study_place) |> summarise(mean_prev = mean(prev))

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
                        )
                      )

  out$study_place <- factor(out$study_place, levels = c("Tanzania (2015)", "Uganda (2017)", "Tanzania (2019)", "Benin (2020)"))

  return(out)
}

or_mean_3y_diff <- plot_labels_fun(or_mean_3y_diff) |> left_join(COMBO_limits_prev) |>
  filter(base_prev_diff >= lbp & base_prev_diff <= ubp)
or_mean_3y_diff_pool <- plot_labels_fun(or_mean_3y_diff_pool) |> left_join(COMBO_limits_prev) |>
  filter(base_prev_diff >= lbp & base_prev_diff <= ubp)
or_mean_3y_mean <- plot_labels_fun(or_mean_3y_mean)
or_mean_3y_mean_pool <- plot_labels_fun(or_mean_3y_mean_pool)
inv_log_BL_prev <- plot_labels_fun(inv_log_BL_prev) |> left_join(COMBO_limits_BL_prev) |>
  filter(start_prev >= lbp & start_prev <= ubp)

eff_mean_3y_diff <- plot_labels_fun(eff_mean_3y_diff) |> left_join(COMBO_limits_prev) |>
  filter(base_prev_diff >= lbp & base_prev_diff <= ubp)

eff_mean_3y_mean <- plot_labels_fun(eff_mean_3y_mean)

cols <- unname(palette.colors(palette = "Okabe-Ito")[c(7,2,3,8)])
cols_net <- c("blue", "aquamarine", "darkgreen")

COMBO_mean_prev$study_place <- factor(COMBO_mean_prev$study_place, levels = c("Tanzania (2015)", "Uganda (2017)", "Tanzania (2019)", "Benin (2020)"))

m2_dp_plot <- ggplot(data = inv_log_BL_prev, #pred_diff,
                     aes(x = start_prev, y = m_p, ymin = l_p, ymax = u_p,
                         fill = net, group = interaction(study_place, net))) + #, year_plot
  geom_abline(intercept = 0, slope = 1, linetype = 2, linewidth = 0.75, col = "grey") +
  geom_ribbon(alpha = 0.1) +
  geom_line(aes(col = net), linewidth = 1) +#, linetype = year_plot
  facet_wrap(~ study_place, nrow = 1) +
  geom_point(data = COMBO_mean_prev, aes(x = BL_prev, y = mean_prev, col = net), inherit.aes = FALSE) +
  ylab("Malaria prevalence") + xlab("Baseline prevalence") +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2), labels = scales::percent) +
  scale_colour_manual(values = cols_net, name = "ITN") +
  scale_fill_manual(values = cols_net, name = "ITN") +
  scale_x_continuous(limits = c(0, 1), breaks = seq(-0, 1, 0.25), labels = scales::percent) +
  #scale_linetype_manual(values = c(1, 2, 3), name = "") +
  labs(tag = "A")

m2_o_r_dp_plot <- ggplot(data = or_mean_3y_diff |> filter(net != "Pyrethroid-only"), #pred_diff
                         aes(x = base_prev_diff, y = med, ymin = low, ymax = up, #ymin = l_or, ymax = u_or, y = m_or,
                             fill = study_place, group = interaction(study_place, net))) +#, year_plot
  geom_ribbon(inherit.aes = FALSE,
              data = subset(or_mean_3y_diff_pool, l != 1 ),
              aes(x = base_prev_diff, y = med, ymin = low, ymax = up, #y = m_or_pooled, ymin = l_or_pooled, ymax = u_or_pooled,
                  group = interaction(study_place)), alpha = 0.1) + #, year_plot
  geom_ribbon(alpha = 0.1) +
  geom_line(aes(col = study_place), linewidth = 1) +#, linetype = year_plot
  geom_hline(yintercept = 0, col = "grey", linetype = 2, linewidth = 0.75) +
  facet_wrap(~ net, nrow = 1) +
  ylab("Log odds ratio of infection in trial ITN\nclusters relative to pyrethroid-only clusters") +
  xlab("Within-study differences in baseline prevalence") +
  scale_x_continuous(limits = c(-0.55, 0.55), breaks = seq(-0.5, 0.5, 0.25), labels = scales::percent) +
  scale_colour_manual(values = cols, name = "study") +
  scale_fill_manual(values = cols, name = "study") +
  geom_line(inherit.aes = FALSE, data = subset(or_mean_3y_diff_pool, l != 1),
            aes(x = base_prev_diff, y = med), linewidth = 1) +#, linetype = year_plot # m_or_pooled
  #scale_linetype_manual(values = c(1, 2, 3), name = "") +
  scale_y_continuous(breaks = seq(-4, 2, 1)) +
  labs(tag = "B")

m2_eff_dp_plot <- ggplot(data = eff_mean_3y_diff |> filter(net != "Pyrethroid-only"),
                         aes(x = base_prev_diff, y = med, ymin = low, ymax = up,
                             fill = study_place, group = interaction(study_place, net))) +
  geom_ribbon(alpha = 0.1) +
  geom_line(aes(col = study_place), linewidth = 1) +
  geom_hline(yintercept = 0, col = "grey", linetype = 2, linewidth = 0.75) +
  facet_wrap(~ net, nrow = 1) +
  ylab("Relative efficacy of trial ITN clusters\nrelative to pyrethroid-only clusters") +
  xlab("Within-study differences in baseline prevalence") +
  scale_x_continuous(limits = c(-0.55, 0.55), breaks = seq(-0.5, 0.5, 0.25), labels = scales::percent) +
  scale_y_continuous(labels = scales::percent, breaks = seq(-2, 1, 0.25)) +
  scale_colour_manual(values = cols, name = "study") +
  scale_fill_manual(values = cols, name = "study") +
  coord_cartesian(ylim = c(-1, 1)) +
  labs(tag = "C")

m2_mean_BL_prev_o_r_plot <- ggplot(data = or_mean_3y_mean |>
                                     subset(net != "Pyrethroid-only"),
                        aes(x = base_prev_mean, y = med, ymin = low, ymax = up,
                            fill = study_place, group = interaction(study_place, net))) + #, year_plot
  geom_ribbon(data = subset(or_mean_3y_mean_pool, net != "Pyrethroid-only"), 
              aes(x = base_prev_mean, ymin = low, ymax = up,
                  group = interaction(net)), inherit.aes = FALSE, alpha = 0.1) +
  geom_ribbon(alpha = 0.1) +
  geom_line(aes(col = study_place), linewidth = 1) +
  geom_line(data = subset(or_mean_3y_mean_pool, net != "Pyrethroid-only"), 
            aes(x = base_prev_mean, y = med), inherit.aes = FALSE, linewidth = 1) +
  geom_hline(yintercept = 0, linetype = 2, linewidth = 0.75, col = "grey") +
  facet_wrap(~net, scales = "free_x") +
  ylab("Log odds ratio of infection in trial ITN\nclusters relative to pyrethroid-only clusters") +
  xlab("Between study mean baseline prevalence") +
  scale_y_continuous(limits = c(-2, 1), breaks = seq(-2, 1, 0.5)) +
  scale_colour_manual(values = cols, name = "study") +
  scale_fill_manual(values = cols, name = "study") +
  scale_x_continuous(labels = scales::percent) +
  labs(tag = "A")

m2_mean_BL_prev_eff_plot <- ggplot(data = eff_mean_3y_mean |> subset(net != "Pyrethroid-only"),
       aes(x = base_prev_mean, y = med, ymin = low, ymax = up,
           fill = study_place, group = interaction(study_place, net))) + #, year_plot
  geom_ribbon(alpha = 0.1) +
  geom_line(aes(col = study_place), linewidth = 1) +
  geom_hline(yintercept = 0, linetype = 2, linewidth = 0.75, col = "grey") +
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

##########################################
##### univariable baseline prevalence #####
##########################################

COMBO_stan_BL_prev <- COMBO_stan |> group_by(study_place, net, BL_prev, cluster) |>
  summarise(prev_num = sum(prev_num), prev_denom = sum(prev_denom)) |>
  mutate(prev = prev_num / prev_denom)

COMBO_stan_BL_prev$study_place <- factor(COMBO_stan_BL_prev$study_place, levels = c("Tanzania (2015)", "Uganda (2017)", "Tanzania (2019)", "Benin (2020)"))

glm(cbind(prev_num, prev_denom - prev_num) ~ BL_prev, family = "binomial", data = COMBO_stan_BL_prev)

uni_BL_prev_plot <- ggplot(data = COMBO_stan_BL_prev, aes(x = BL_prev, y = prev, group = interaction(study_place, net), col = net, fill = net)) +
  geom_abline(intercept = 0, slope = 1, linetype = 2, linewidth = 0.75, col = "grey") +
  geom_point() +
  geom_smooth(se = FALSE, formula = y ~ x, method = "glm", method.args = list(family = "binomial"), aes(weight = prev_denom)) +
  facet_wrap(~study_place, nrow = 2) +
  ylab("Mean malaria prevalence post ITN distribution") + xlab("Baseline prevalence") +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2), labels = scales::percent) +
  scale_colour_manual(values = cols_net, name = "ITN") +
  scale_fill_manual(values = cols_net, name = "ITN") +
  scale_x_continuous(limits = c(0, 1), breaks = seq(-0, 1, 0.25), labels = scales::percent)

ggsave(file = "univariable_BL_prev_plots.pdf",
       plot = uni_BL_prev_plot,
       device = "pdf",
       width = 27.5, height = 17.5,
       units = "cm")

############################
##### parameter values #####
############################

# model 2
m2_params <- get_params(fit = m2_fit_full)
m2_params_df <- data.frame("params" = m2_params$params, params_description = m2_params$params_description, "values" = m2_params$values)

m2_params_df <- rbind(m2_params_df,
                      data.frame(
                        "params" = c("$\\kappa^{\\text{pooled}}$", "$\\kappa_{1}$", "$\\kappa_{2}", "$\\kappa_{3}$", "$\\kappa_{4}",
                                     "$\\kappa^{\\text{net}}_{1}$", "$\\kappa^{\\text{net}}_{2}$", "$\\kappa^{\\text{study}}_{2}$",
                                     "$\\kappa^{\\text{study}}_{3}$", "$\\kappa^{\\text{study}}_{4}$", "$\\kappa^{\\text{study}}_{5}$",
                                     "$\\kappa^{\\text{study}}_{6}$",
                                     "$\\gamma$",
                                     "$\\gamma^{\\text{net}}_{1}$",
                                     "$\\gamma^{\\text{net}}_{2}$",
                                     "$\\delta$",
                                     "$\\omega$",
                                     "$\\omega^{\\text{net}}_{1}$",
                                     "$\\omega^{\\text{net}}_{2}$"),
                        "params_description" = c("pooled year effect",
                                                 "Tanzania (2015) year effect",
                                                 "Uganda (2017) year effect",
                                                 "Tanzania (2019) year effect",
                                                 "Benin (2020) year effect",
                                                 "pooled pyrethroid-PBO year effect interaction", 
                                                 "pooled pyrethroid-pyrrole year effect interaction",
                                                 "Tanzania (2015) pyrethroid-PBO year effect interaction",
                                                 "Uganda (2017) pyrethroid-PBO year effect interaction",
                                                 "Tanzania (2019) pyrethroid-PBO year effect interaction",
                                                 "Tanzania (2019) pyrethroid-pyrrole year effect interaction",
                                                 "Benin (2020) pyrethroid-pyrrole year effect interaction",
                                                 "within study differences in baseline prevalence effect",
                                                 "within study differences in baseline prevalence pyrethroid-PBO interaction",
                                                 "within study differences in baseline prevalence pyrethroid-pyrrole interaction",
                                                 "within study differences in baseline prevalence year interaction",
                                                 "mean study baseline prevalence effect",
                                                 "mean study baseline prevalence pyrethroid-PBO interaction",
                                                 "mean study baseline prevalence pyrethroid-pyrrole interaction"
                                                 ),
                        "values" = c(get_quantiles(model_fit = m1_fit_full, param = "kappa", dim = NULL),
                                     get_quantiles(model_fit = m1_fit_full, param = "kappa_i", dim = 1),
                                     get_quantiles(model_fit = m1_fit_full, param = "kappa_i", dim = 2),
                                     get_quantiles(model_fit = m1_fit_full, param = "kappa_i", dim = 3),
                                     get_quantiles(model_fit = m1_fit_full, param = "kappa_i", dim = 4),
                                     get_quantiles(model_fit = m1_fit_full, param = "kappa_l", dim = 1),
                                     get_quantiles(model_fit = m1_fit_full, param = "kappa_l", dim = 2),
                                     get_quantiles(model_fit = m1_fit_full, param = "kappa_li", dim = 2),
                                     get_quantiles(model_fit = m1_fit_full, param = "kappa_li", dim = 3),
                                     get_quantiles(model_fit = m1_fit_full, param = "kappa_li", dim = 4),
                                     get_quantiles(model_fit = m1_fit_full, param = "kappa_li", dim = 5),
                                     get_quantiles(model_fit = m1_fit_full, param = "kappa_li", dim = 6),
                                     get_quantiles(model_fit = m2_fit_full, param = "gamma", dim = NULL),
                                     get_quantiles(model_fit = m2_fit_full, param = "delta_l_raw", dim = 1),
                                     get_quantiles(model_fit = m2_fit_full, param = "delta_l_raw", dim = 2),
                                     get_quantiles(model_fit = m2_fit_full, param = "delta", dim = NULL),
                                     get_quantiles(model_fit = m2_fit_full, param = "omega", dim = NULL),
                                     get_quantiles(model_fit = m2_fit_full, param = "omega_l_raw", dim = 1),
                                     get_quantiles(model_fit = m2_fit_full, param = "omega_l_raw", dim = 2)
                                     )
                        )
                      )

write.csv(m2_params_df, file = "m2_params_df.csv")
