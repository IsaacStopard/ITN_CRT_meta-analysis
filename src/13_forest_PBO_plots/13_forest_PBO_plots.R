library(cmdstanr)
library(posterior)
library(tidyverse)
library(loo)
library(patchwork)

ggplot2::theme_set(theme_bw() +
                     theme(text = element_text(size = 17),
                           legend.text = element_text(size = 11),
                           legend.title = element_text(size = 11),
                           panel.background = element_blank(),
                           panel.grid.major = element_blank(),
                           panel.grid.minor = element_blank(),
                           axis.text = element_text(colour = "black", size = 17),
                           strip.text = element_text(size = 17))
)

orderly2::orderly_shared_resource("m1.stan", "model_functions.R")

source("model_functions.R")

orderly2::orderly_dependency("12_forest_PBO_fit",
                             "latest()",
                             c("m1_re_PBO_subset_fit_full.rds",
                               "COMBO_stan_PBO_subset.rds",
                               "data_in_PBO_subset.rds")
                             )

m1_fit_full <- readRDS(file = "m1_re_PBO_subset_fit_full.rds")

m1_stan_model <- cmdstan_model(stan_file = "m1.stan")

COMBO_stan <- readRDS(file = "COMBO_stan_PBO_subset.rds") |> mutate(study_place = case_when(
  str_detect(study, "Protopopoff") ~ "Tanzania (2015)",
  str_detect(study, "Staedke") ~ "Uganda (2017)",
  str_detect(study, "Mosha") ~ "Tanzania (2019)",
  str_detect(study, "Accrombessi") ~ "Benin (2020)"))

COMBO_stan$study_place <- factor(COMBO_stan$study_place, levels = c("Tanzania (2015)", "Uganda (2017)", "Tanzania (2019)", "Benin (2020)"))

data_in_full <- readRDS("data_in_PBO_subset.rds")

iter <- 3000
warmup <- floor(iter/2)

lower <- 0.025
upper <- 0.975

###########################
##### m1 forest plots #####
###########################

# parameters
target_params <- c("theta_l", "theta_li", "kappa_l", "kappa_li", "alpha_i", "alpha",
                   "kappa", "kappa_i", "e_ij", "sigma_e_r")

params_model_1 <- extract_params(target_params, m1_fit_full)

years <- seq(0.1, 3, 0.1)

u_li_pred <- unique(COMBO_stan[,c("l", "li", "i", "net", "study_place")]) |> arrange(li)

u_li_all <- u_li_pred[rep(1:nrow(u_li_pred), length(years)),] |> mutate(years = rep(years, each = nrow(u_li_pred)), ij = i)

m1_fit_full$expose_functions()

list2env(params_model_1$params_list, envir = environment())

# extrapolating to new cluster for each study
e_ij_single <- matrix(rnorm(length(sigma_e_r), 0, 1), nrow = nrow(sigma_e_r)) * sigma_e_r

n_iter <- nrow(sigma_e_r)

logit_prob_t <- sapply(1:n_iter, function(i){
  m1_fit_full$functions$logit_prob_fun(nrow(u_li_all), u_li_all$i, u_li_all$l, u_li_all$li, u_li_all$i,
                                       alpha[i], alpha_i[i, ], theta_l[i, ], theta_li[i, ], kappa[i], kappa_i[i, ], 
                                       kappa_l[i, ], kappa_li[i, ], u_li_all$years, e_ij_single[i, ])}
) |> plogis()

probs_in <- c(lower, 0.5, upper)

u_li_all <- cbind(
  cbind(
    cbind(u_li_all,
          l_or_pooled_fun(u_li_all, params_model_1, model = 1, quantile_out = TRUE) |> rename("o_r_pooled_l" = 1, "o_r_pooled_m" = 2, "o_r_pooled_u" = 3)
    ),
    l_or_study_fun(u_li_all, params_model_1, model = 1, quantile_out = TRUE) |> rename("o_r_l" = 1, "o_r_m" = 2, "o_r_u" = 3)
  ),
  apply(logit_prob_t, 1, quantile, probs = c(lower, 0.5, upper)) |>
    t() |> as.data.frame() |> rename("l_p" = 1, "m_p" = 2, "u_p" = 3)
)


# colours
cols <- unname(palette.colors(palette = "Okabe-Ito")[c(7,3,8)])
cols_net <- c("blue", "aquamarine", "darkgreen")

m1_or_plot <- ggplot(data = u_li_all |> subset(l != 1)) +
  geom_ribbon(aes(x = years, ymin = o_r_l, ymax = o_r_u, fill = study_place), alpha = 0.25) +
  geom_line(aes(x = years, y = o_r_m, col = study_place), linewidth = 1.5) +
  geom_ribbon(aes(x = years, ymin = o_r_pooled_l, ymax = o_r_pooled_u), alpha = 0.25) +
  geom_line(aes(x = years, y = o_r_pooled_m), linewidth = 1.5) +
  facet_wrap(~net) +
  ylab("Log odds ratio of infection in trial ITN\nclusters relative to pyrethroid-only clusters") +
  scale_fill_manual(values = cols, name = "") +
  scale_colour_manual(values = cols, name = "") +
  geom_hline(yintercept = 0, linetype = 2) +
  theme(legend.position = c(0.925, 0.925)) +
  scale_y_continuous(limits = c(-1.5, 0.5), breaks = seq(-1.5, 0.5, 0.5))

######################################################
##### calculating the times the upper ORs pass 1 #####
######################################################
u_li_all |> group_by(study_place, i) |> filter(o_r_u > 0) |> group_by(study_place) |> filter(years == min(years))

u_li_all |> group_by(study_place, i) |> filter(o_r_pooled_u > 0) |> group_by(study_place) |> filter(years == min(years))

############################
##### time odds ratios #####
############################

u_li_pred |> cbind(apply(
  exp(
    replicate(nrow(u_li_pred), as.vector(as_draws_matrix(m1_fit_full$draws("kappa")))) + 
    as_draws_matrix(m1_fit_full$draws("kappa_i"))[,u_li_pred$i] +
    as_draws_matrix(m1_fit_full$draws("kappa_l"))[,u_li_pred$l] +
    as_draws_matrix(m1_fit_full$draws("kappa_li"))[,u_li_pred$li]
    ), 2, quantile,
  probs = c(lower, 0.5, upper)) |> round(digits = 2) |> t())

########################################
##### calculating the mean log ORs #####
########################################

u_li_comp <- u_li_pred

u_li_comp_m1_3 <- u_li_comp |> mutate(model = "model 1", timespan = "3-years", model_timespan = "model 1, 3-years")
u_li_comp_m1_2 <- u_li_comp |> mutate(model = "model 1", timespan = "2-years", model_timespan = "model 1, 2-years")

# calculating the means over 2 and 3 years

u_li_comp_m1_3 <- cbind(u_li_comp_m1_3, 
                        mean_time_l_or_study_fun(3, u_li_comp_m1_3, params_model_1, model = 1, quantile_out = TRUE) |> rename("l_o_r" = 1, "m_o_r" = 2, "u_o_r" = 3)
) |> 
  cbind(mean_time_l_or_pooled_fun(3, u_li_comp_m1_3, params_model_1, model = 1, quantile_out = TRUE) |> rename("l_o_r_pooled" = 1, "m_o_r_pooled" = 2, "u_o_r_pooled" = 3))

u_li_comp_m1_2 <- cbind(u_li_comp_m1_2, 
                        mean_time_l_or_study_fun(2, u_li_comp_m1_2, params_model_1, model = 1, quantile_out = TRUE) |> rename("l_o_r" = 1, "m_o_r" = 2, "u_o_r" = 3)
) |> 
  cbind(mean_time_l_or_pooled_fun(2, u_li_comp_m1_2, params_model_1, model = 1, quantile_out = TRUE) |> rename("l_o_r_pooled" = 1, "m_o_r_pooled" = 2, "u_o_r_pooled" = 3))

u_li_pooled_m1_3 <- u_li_comp_m1_3[c(5, 7), c("net", "l_o_r_pooled", "m_o_r_pooled", "u_o_r_pooled", "model_timespan")] |> mutate(study_place = "Pooled")
u_li_pooled_m1_2 <- u_li_comp_m1_2[c(5, 7), c("net", "l_o_r_pooled", "m_o_r_pooled", "u_o_r_pooled", "model_timespan")] |> mutate(study_place = "Pooled")

u_li_comp <- rbind(u_li_comp_m1_2, u_li_comp_m1_3)

u_li_comp_pbo <- subset(u_li_comp, net == "Pyrethroid-PBO")
u_li_comp_pbo$study_place <- factor(u_li_comp_pbo$study_place, levels = c("Pooled", "Tanzania (2015)", "Uganda (2017)", "Tanzania (2019)"))
u_li_comp_pp <- subset(u_li_comp, net == "Pyrethroid-pyrrole")
u_li_comp_pp$study_place <- factor(u_li_comp_pp$study_place, levels = c("Pooled", "Tanzania (2019)", "Benin (2020)"))

# checking the log OR for Tanzania (2019)
subset(u_li_comp_pbo, study_place == "Tanzania (2019)")

log_o_r_plot_m1_subset_y2_y3 <- ggplot(data = u_li_comp_pbo)  +
  scale_y_discrete(drop = FALSE) +
  geom_point(aes(x = m_o_r, y = study_place, col = model_timespan), position = position_dodge(width = 0.25), size = 5, alpha = 0.75) +
  geom_errorbar(aes(xmin = l_o_r, xmax = u_o_r, y = study_place, col = model_timespan), position = position_dodge(width = 0.25), height = 0, alpha = 0.75, orientation = "y") +
  geom_polygon(data = data.frame(x = unname(unlist(c(u_li_pooled_m1_3[1, c(2, 3, 4, 3)]))), y = c(1 + 0.1, 1.1 + 0.1, 1 + 0.1, 0.9 + 0.1), model_timespan = u_li_pooled_m1_3[1, "model_timespan"]),
               aes(x = x, y = y, fill = model_timespan),
               inherit.aes = FALSE, alpha = 0.5) +
  geom_polygon(data = data.frame(x = unname(unlist(c(u_li_pooled_m1_2[1, c(2, 3, 4, 3)]))), y = c(1 + 0.05, 1.1 + 0.05, 1 + 0.05, 0.9 + 0.05), model_timespan = u_li_pooled_m1_2[1, "model_timespan"]),
               aes(x = x, y = y, fill = model_timespan),
               inherit.aes = FALSE, alpha = 0.5) +
  scale_x_continuous(limits = c(-1.5, 0.25), breaks = seq(-1.5, 0.5, 0.5)) +
  xlab("Log odds ratio") + ylab("") +
  geom_vline(xintercept = 0, linetype = 2) +
  theme(axis.ticks.y = element_blank(),
        panel.border = element_blank(),
        axis.line.x = element_line(colour = "black"),
        legend.position = "none") +
  scale_colour_manual(values = c("red", "skyblue"), name = "") +
  scale_fill_manual(values = c("red", "skyblue"), name = "") +
  facet_wrap(~net) + labs(tag = "B") +
  ggplot(data = u_li_comp_pp)  +
  scale_y_discrete(drop = FALSE) +
  geom_point(aes(x = m_o_r, y = study_place, col = model_timespan), position = position_dodge(width = 0.25), size = 5, alpha = 0.75) +
  geom_errorbar(aes(xmin = l_o_r, xmax = u_o_r, y = study_place, col = model_timespan), position = position_dodge(width = 0.25), height = 0, alpha = 0.75, orientation = "y") +
  geom_polygon(data = data.frame(x = unname(unlist(c(u_li_pooled_m1_3[2, c(2, 3, 4, 3)]))), y = c(1 + 0.1, 1.1 + 0.1, 1 + 0.1, 0.9 + 0.1), model_timespan = u_li_pooled_m1_3[2, "model_timespan"]),
               aes(x = x, y = y, fill = model_timespan),
               inherit.aes = FALSE, alpha = 0.5) +
  geom_polygon(data = data.frame(x = unname(unlist(c(u_li_pooled_m1_2[2, c(2, 3, 4, 3)]))), y = c(1 + 0.05, 1.1 + 0.05, 1 + 0.05, 0.9 + 0.05), model_timespan = u_li_pooled_m1_2[2, "model_timespan"]),
               aes(x = x, y = y, fill = model_timespan),
               inherit.aes = FALSE, alpha = 0.5) +
  scale_x_continuous(limits = c(-1.5, 0.25), breaks = seq(-1.5, 0.5, 0.5)) +
  xlab("Log odds ratio") + ylab("") +
  geom_vline(xintercept = 0, linetype = 2) +
  theme(axis.ticks.y = element_blank(),
        panel.border = element_blank(),
        axis.line.x = element_line(colour = "black"),
        legend.position = c(0.175, 0.925)) +
  scale_colour_manual(values = c("red", "skyblue"), name = "") +
  scale_fill_manual(values = c("red", "skyblue"), name = "") +
  facet_wrap(~net)

ggsave(file = "m1_forest_plot_comb_subset.pdf",
       device = "pdf",
       plot = (m1_or_plot + labs(tag = "A")) + (log_o_r_plot_m1_subset_y2_y3),
       width = 70, height = 20,
       units = "cm"
)

##################
##### params #####
##################

get_params_less <- function(fit){
  params <- c("$\\alpha^{\\text{pooled}}$")
  params_description <- ("pooled intercept")
  values <- c(get_quantiles(fit, "alpha", dim = NULL))
  # alpha_i
  params <- c(params, "$\\alpha_{i}$", "$\\alpha_{i}$", "$\\alpha_{i}$")
  params_description <- c(params_description, "study effect (Uganda (2017))", "study effect (Tanzania (2019))", "study effect (Benin (2020))")
  for(i in 1:3){
    values <- c(values, get_quantiles(fit, "alpha_i", dim = i))
  }
  # net effect
  params <- c(params, "$\\theta^{\\text{net}}_{1}$")
  params_description <- c(params_description, "pyrethroid-PBO pooled treatment effect")
  values <- c(values, get_quantiles(fit, "theta_l", dim = 2))
  
  params_description <- c(params_description, "pyrethroid-PBO-Uganda (2017) interaction", "pyrethroid-PBO-Tanzania (2019) interaction")
  params <- c(params, "$\\theta^{\\text{study}}_{2}$", "$\\theta^{\\text{study}}_{3}$")
  for(i in 2:3){
    values <- c(values, get_quantiles(fit, "theta_li", dim = i))
  }
  
  params <- c(params, "$\\theta^{\\text{net}}_{2}$")
  params_description <- c(params_description, "pyrethroid-pyrrole pooled treatment effect")
  values <- c(values, get_quantiles(fit, "theta_l", dim = 2))
  
  params_description <- c(params_description, "pyrethroid-pyrrole-Tanzania (2019) interaction", "pyrethroid-pyrrole-Benin (2020) interaction")
  params <- c(params, "$\\theta^{\\text{study}}_{4}$", "$\\theta^{\\text{study}}_{5}$")
  values <- c(values, get_quantiles(fit, "theta_li", dim = 4), get_quantiles(fit, "theta_li", dim = 5))
  
  params_description <- c(params_description, "Uganda (2017) cluster random effect standard deviation",
                          "Tanzania (2019) cluster random effect standard deviation", "Benin (2020) cluster random effect standard deviation")
  
  params <- c(params, "$\\sigma_{1}$", "$\\sigma_{2}$", "$\\sigma_{3}$")
  values <- c(values, get_quantiles(fit, "sigma_e_r", dim = 1), get_quantiles(fit, "sigma_e_r", dim = 2),
              get_quantiles(fit, "sigma_e_r", dim = 3))
  
  params_description <- c(params_description, "pyrethroid-PBO between study random treatment effect standard deviation", "pyrethroid-pyrrole between study random treatment effect standard deviation")
  params <- c(params, "$\\tau_{1}", "$\\tau_{2}")
  values <- c(values, get_quantiles(fit, "tau_sd_net_li", dim = 1), get_quantiles(fit, "tau_sd_net_li", dim = 2))
  
  return(list("params" = params, "params_description" = params_description, "values" = values))
}

m1_params <- get_params_less(fit = m1_fit_full)
m1_params_df <- data.frame("params" = m1_params$params, params_description = m1_params$params_description, "values" = m1_params$values)
m1_params_df <- rbind(m1_params_df,
                      data.frame(
                        "params" = c("$\\kappa^{\\text{pooled}}$", "$\\kappa_{1}$", "$\\kappa_{2}", "$\\kappa_{3}$",
                                     "$\\kappa^{\\text{net}}_{1}$", "$\\kappa^{\\text{net}}_{2}$", "$\\kappa^{\\text{study}}_{2}$",
                                     "$\\kappa^{\\text{study}}_{3}$", "$\\kappa^{\\text{study}}_{4}$", "$\\kappa^{\\text{study}}_{5}$"),
                        "params_description" = c("pooled year effect",
                                                 "Uganda (2017) year effect",
                                                 "Tanzania (2019) year effect",
                                                 "Benin (2020) year effect",
                                                 "pooled pyrethroid-PBO year effect interaction", 
                                                 "pooled pyrethroid-pyrrole year effect interaction",
                                                 "Uganda (2017) pyrethroid-PBO year effect interaction",
                                                 "Tanzania (2019) pyrethroid-PBO year effect interaction",
                                                 "Tanzania (2019) pyrethroid-pyrrole year effect interaction",
                                                 "Benin (2020) pyrethroid-pyrrole year effect interaction"),
                        "values" = c(get_quantiles(model_fit = m1_fit_full, param = "kappa", dim = NULL),
                                     get_quantiles(model_fit = m1_fit_full, param = "kappa_i", dim = 1),
                                     get_quantiles(model_fit = m1_fit_full, param = "kappa_i", dim = 2),
                                     get_quantiles(model_fit = m1_fit_full, param = "kappa_i", dim = 3),
                                     get_quantiles(model_fit = m1_fit_full, param = "kappa_l", dim = 1),
                                     get_quantiles(model_fit = m1_fit_full, param = "kappa_l", dim = 2),
                                     get_quantiles(model_fit = m1_fit_full, param = "kappa_li", dim = 2),
                                     get_quantiles(model_fit = m1_fit_full, param = "kappa_li", dim = 3),
                                     get_quantiles(model_fit = m1_fit_full, param = "kappa_li", dim = 4),
                                     get_quantiles(model_fit = m1_fit_full, param = "kappa_li", dim = 5))
                        ))

write.csv(m1_params_df, file = "m1_params_df_subset.csv")

