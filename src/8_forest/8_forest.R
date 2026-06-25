library(cmdstanr)
library(tidyverse)
library(loo)
library(patchwork)
library(posterior)
library(bayesplot)

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

orderly2::orderly_dependency("1_data_cleaning", "latest()", c("stan_data.rds"))

orderly2::orderly_dependency("4_time_m0_c_fits", quote(latest(parameter:model == "m0")), 
                             c("m0_time_fits_data.rds", "m0_c_fit_1.rds", "m0_c_fit_2.rds",
                               "m0_c_fit_3.rds", "m0_c_fit_4.rds"))

orderly2::orderly_shared_resource("m0.stan",
                                  "m1.stan",
                                  "model_functions.R")

source("model_functions.R")

orderly2::orderly_dependency("2_fits", "latest()", c("fit_0.rds", "fit_1.rds"))

list2env(readRDS("stan_data.rds"), envir = .GlobalEnv)

m0_fit_full <- readRDS(file = "fit_0.rds")
m1_fit_full <- readRDS(file = "fit_1.rds")

m0_fit_1 <- readRDS(file = "m0_c_fit_1.rds")
m0_fit_2 <- readRDS(file = "m0_c_fit_2.rds")
m0_fit_3 <- readRDS(file = "m0_c_fit_3.rds")
m0_fit_4 <- readRDS(file = "m0_c_fit_4.rds")

m0_stan_model <- cmdstan_model(stan_file = "m0.stan", cpp_options = list(stan_threads = TRUE))

m1_stan_model <- cmdstan_model(stan_file = "m1.stan", cpp_options = list(stan_threads = TRUE))

list2env(readRDS(file = "m0_time_fits_data.rds"), envir = .GlobalEnv)

ucs <- unique(COMBO_stan[, c("cluster", "study", "i", "BL_prev_num", "BL_prev_denom", "BL_prev")])

ucs <- ucs |> left_join(ucs |> group_by(study) |>
                          summarise(mean_BL_prev = mean(BL_prev))) |>
  mutate(diff_BL_prev = BL_prev - mean_BL_prev)

COMBO_stan <- COMBO_stan |>
  mutate(study_place = case_when(
           str_detect(study, "Protopopoff") ~ "Tanzania (2015)",
           str_detect(study, "Staedke") ~ "Uganda (2017)",
           str_detect(study, "Mosha") ~ "Tanzania (2019)",
           str_detect(study, "Accrombessi") ~ "Benin (2020)"
         )) |>
  left_join(ucs)

COMBO_stan$study_place = factor(COMBO_stan$study_place, levels = c("Tanzania (2015)", "Uganda (2017)", "Tanzania (2019)", "Benin (2020)"))

########################
##### model checks #####
########################

# pars_in <- c("^theta_li", "^theta_l", "tau_sd_net_li", "alpha", "^alpha_i")
# 
# m0_posterior_draws <- as_draws_array(m0_fit_full)
# m0_1_posterior_draws <- as_draws_array(m0_fit_1)
# m0_2_posterior_draws <- as_draws_array(m0_fit_2)
# m0_3_posterior_draws <- as_draws_array(m0_fit_3)
# m0_4_posterior_draws <- as_draws_array(m0_fit_4)
# 
# mcmc_pairs(m0_posterior_draws, regex_pars = pars_in, off_diag_fun = "hex")
# mcmc_trace(m0_posterior_draws, regex_pars = pars_in)
# 
# mcmc_pairs(m0_1_posterior_draws, regex_pars = pars_in, off_diag_fun = "hex")
# mcmc_trace(m0_1_posterior_draws, regex_pars = pars_in)
# 
# mcmc_pairs(m0_2_posterior_draws, regex_pars = pars_in, off_diag_fun = "hex")
# mcmc_trace(m0_2_posterior_draws, regex_pars = pars_in)
# 
# mcmc_pairs(m0_3_posterior_draws, regex_pars = pars_in, off_diag_fun = "hex")
# mcmc_trace(m0_3_posterior_draws, regex_pars = pars_in)
# 
# mcmc_pairs(m0_4_posterior_draws, regex_pars = pars_in, off_diag_fun = "hex")
# mcmc_trace(m0_4_posterior_draws, regex_pars = pars_in)

###############################
##### model 0 forest plot #####
###############################

u_li <- unique(COMBO_stan[,c("l", "li", "i", "net", "study")]) |> arrange(li) |> filter(l!=1) |>
  mutate(study_place = case_when(
    str_detect(study, "Protopopoff") ~ "Tanzania (2015)",
    str_detect(study, "Staedke") ~ "Uganda (2017)",
    str_detect(study, "Mosha") ~ "Tanzania (2019)",
    str_detect(study, "Accrombessi") ~ "Benin (2020)"))

gen_forest <- function(fit_full, u_li, title, COMBO_stan){


  COMBO_text <- COMBO_stan |> mutate(study_place = case_when(
                         str_detect(study, "Protopopoff") ~ "Tanzania (2015)",
                         str_detect(study, "Staedke") ~ "Uganda (2017)",
                         str_detect(study, "Mosha") ~ "Tanzania (2019)",
                         str_detect(study, "Accrombessi") ~ "Benin (2020)")) |>
  group_by(study_place, net) |>
    summarise(positive = as.character(sum(prev_num)),
              negative = as.character(sum(prev_denom) - sum(prev_num)))

  for(i in 1:nrow(u_li)){
  q <- quantile(as_draws_matrix(fit_full$draws("theta_l"))[, u_li[i, "l"]] + as_draws_matrix(fit_full$draws("theta_li"))[, u_li[i, "li"]], probs = c(lower, 0.5, upper))
  u_li[i, "low"] <- q[1]
  u_li[i, "med"] <- q[2]
  u_li[i, "up"] <- q[3]
  }

  u_li[,"text_label"] <- paste(round(u_li[,"med"], digits = 2), " [", round(u_li[,"low"], digits = 2), ", ", round(u_li[,"up"], digits = 2),"]")

  u_li_pbo <- subset(u_li, net == "Pyrethroid-PBO")
  u_li_pbo$study_place <- factor(u_li_pbo$study_place, levels = c("Pooled", "Tanzania (2015)", "Uganda (2017)", "Tanzania (2019)"))

  u_li_pp <- subset(u_li, net == "Pyrethroid-pyrrole")
  u_li_pp$study_place <- factor(u_li_pp$study_place, levels = c("Pooled", "Tanzania (2019)", "Benin (2020)"))

  u_li_pooled <- data.frame("net" = c("Pyrethroid-PBO", "Pyrethroid-pyrrole"),
                          "low" = c(quantile(as_draws_matrix(fit_full$draws("theta_l_raw"))[, 1], probs = c(lower)),
                                    quantile(as_draws_matrix(fit_full$draws("theta_l_raw"))[, 2], probs = c(lower))),
                          "med" = c(quantile(as_draws_matrix(fit_full$draws("theta_l_raw"))[, 1], probs = c(0.5)),
                                    quantile(as_draws_matrix(fit_full$draws("theta_l_raw"))[, 2], probs = c(0.5))),
                          "up" = c(quantile(as_draws_matrix(fit_full$draws("theta_l_raw"))[, 1], probs = c(upper)),
                                   quantile(as_draws_matrix(fit_full$draws("theta_l_raw"))[, 2], probs = c(upper))),
                          "study_place" = rep("Pooled", 2)
                          )

  u_li_pooled[,"text_label"] <- paste(round(u_li_pooled[,"med"], digits = 2), " [", round(u_li_pooled[,"low"], digits = 2), ", ", round(u_li_pooled[,"up"], digits = 2),"]")

  m0_tau <- round(apply(as_draws_matrix(fit_full$draws("tau_sd_net_li")), 2, quantile, probs = c(lower, 0.5, upper)), digits = 2)

m0_forest <- ggplot(data = u_li_pbo,
             aes(x = med, y = study_place, xmin = low, xmax = up)) +
  scale_y_discrete(drop = FALSE) +
  geom_polygon(data = data.frame(x = unname(unlist(c(u_li_pooled[1, c(2, 3, 4, 3)]))), y = c(1, 1.1, 1, 0.9)),
               aes(x = x, y = y), fill = "white", col = "black",
               inherit.aes = FALSE) +
  geom_point(size = 3) +
  geom_errorbar(height = 0, orientation = "y") +
  scale_x_continuous(limits = c(-3.25, 1.25), breaks = seq(-1.5, 0.5, 0.5)) +
  xlab("Log odds ratio") + ylab("") +
  geom_vline(xintercept = 0, linetype = 2) +
  theme(axis.ticks.y = element_blank(),
        panel.border = element_blank(),
        axis.line.x = element_line(colour = "black")) +
  facet_wrap(~net) +
  ggtitle(title) +
  geom_text(aes(x = 0.5, y = study_place, hjust = 0, label = text_label), size = 3) +
  geom_text(aes(x = 0.5, y = 1, hjust = 0, label = u_li_pooled[1, "text_label"]), size = 3) +
  geom_text(aes(x = 0.5, y = nlevels(u_li_pbo[,"study_place"]) + 0.25, label = "Log OR (95% CI)", hjust = 0), size = 3) +
  geom_text(data = subset(COMBO_text, net == "Pyrethroid-only" & study_place != "Benin (2020)") |> select(-net), aes(x = -2.25, y = study_place, label = positive),
            hjust = 0, inherit.aes = FALSE) +
  geom_text(data = subset(COMBO_text, net == "Pyrethroid-only" & study_place != "Benin (2020)") |> select(-net), aes(x = -1.75, y = study_place, label = negative),
            hjust = 0, inherit.aes = FALSE) +
  geom_text(aes(x = -2, y = nlevels(u_li_pbo[,"study_place"]) + 0.35, label = "Control"), hjust = 0) +
  geom_text(aes(x = -2.25, y = nlevels(u_li_pbo[,"study_place"]) + 0.2, label = "+"), hjust = 0) +
  geom_text(aes(x = -1.75, y = nlevels(u_li_pbo[,"study_place"]) + 0.2, label = "-"), hjust = 0) +
  geom_text(data = subset(COMBO_text, net == "Pyrethroid-PBO" & study_place != "Benin (2020)") |> select(-net), aes(x = -3.25, y = study_place, label = positive),
            hjust = 0, inherit.aes = FALSE) +
  geom_text(data = subset(COMBO_text, net == "Pyrethroid-PBO" & study_place != "Benin (2020)") |> select(-net), aes(x = -2.75, y = study_place, label = negative),
            hjust = 0, inherit.aes = FALSE) +
  geom_text(aes(x = -3, y = nlevels(u_li_pbo[,"study_place"]) + 0.35, label = "Treatment"), hjust = 0) +
  geom_text(aes(x = -3.25, y = nlevels(u_li_pbo[,"study_place"]) + 0.2, label = "+"), hjust = 0) +
  geom_text(aes(x = -2.75, y = nlevels(u_li_pbo[,"study_place"]) + 0.2, label = "-"), hjust = 0) +
  annotate("text", x = -2.75, y = 0.5, label = paste("tau^2 == ", m0_tau[2,1], " (", m0_tau[1,1], " - ", m0_tau[3,1], ")"), parse = TRUE, size = 5) +

ggplot(data = u_li_pp,
         aes(x = med, y = study_place, xmin = low, xmax = up)) +
  scale_y_discrete(drop = FALSE) +
  geom_polygon(data = data.frame(x = unname(unlist(c(u_li_pooled[2, c(2, 3, 4, 3)]))), y = c(1, 1.1, 1, 0.9)),
               aes(x = x, y = y), fill = "white", col = "black",
               inherit.aes = FALSE) +
  geom_point(size = 3) +
  geom_errorbar(height = 0, orientation = "y") + scale_x_continuous(limits = c(-3.25, 1.25), breaks = seq(-1.5, 0.5, 0.5)) +
  xlab("Log odds ratio") + ylab("") +
  geom_vline(xintercept = 0, linetype = 2) +
  theme(axis.ticks.y = element_blank(),
        panel.border = element_blank(),
        axis.line.x = element_line(colour = "black")) +
  facet_wrap(~net) +
  geom_text(aes(x = 0.5, y = study_place, hjust = 0, label = text_label), size = 3) +
  geom_text(aes(x = 0.5, y = 1, hjust = 0, label = u_li_pooled[2, "text_label"]), size = 3) +
  geom_text(aes(x = 0.5, y = nlevels(u_li_pp[,"study_place"]) + 0.25, label = "Log OR (95% CI)", hjust = 0), size = 3) +
  geom_text(data = subset(COMBO_text, net == "Pyrethroid-only" & study_place %in% c("Benin (2020)", "Tanzania (2019)")) |> select(-net), aes(x = -2.25, y = study_place, label = positive),
            hjust = 0, inherit.aes = FALSE) +
  geom_text(data = subset(COMBO_text, net == "Pyrethroid-only" & study_place %in% c("Benin (2020)", "Tanzania (2019)")) |> select(-net), aes(x = -1.75, y = study_place, label = negative),
            hjust = 0, inherit.aes = FALSE) +
  geom_text(aes(x = -2, y = nlevels(u_li_pp[,"study_place"]) + 0.35, label = "Control"), hjust = 0) +
  geom_text(aes(x = -2.25, y = nlevels(u_li_pp[,"study_place"]) + 0.2, label = "+"), hjust = 0) +
  geom_text(aes(x = -1.75, y = nlevels(u_li_pp[,"study_place"]) + 0.2, label = "-"), hjust = 0) +
  geom_text(data = subset(COMBO_text, net == "Pyrethroid-pyrrole" & study_place  %in% c("Benin (2020)", "Tanzania (2019)")) |> select(-net), aes(x = -3.25, y = study_place, label = positive),
            hjust = 0, inherit.aes = FALSE) +
  geom_text(data = subset(COMBO_text, net == "Pyrethroid-pyrrole" & study_place  %in% c("Benin (2020)", "Tanzania (2019)")) |> select(-net), aes(x = -2.75, y = study_place, label = negative),
            hjust = 0, inherit.aes = FALSE) +
  geom_text(aes(x = -3, y = nlevels(u_li_pp[,"study_place"]) + 0.35, label = "Treatment"), hjust = 0) +
  geom_text(aes(x = -3.25, y = nlevels(u_li_pp[,"study_place"]) + 0.2, label = "+"), hjust = 0) +
  geom_text(aes(x = -2.75, y = nlevels(u_li_pp[,"study_place"]) + 0.2, label = "-"), hjust = 0) +
  annotate("text", x = -2.75, y = 0.5, label = paste("tau^2 == ", m0_tau[2,2], " (", m0_tau[1,2], " - ", m0_tau[3,2], ")"), parse = TRUE, size = 5)


return(list("u_li" = u_li, "u_li_pooled" = u_li_pooled, "tau_sd_net_li" = m0_tau, "plot" = m0_forest))
}

m0_forest_all <- gen_forest(fit_full = m0_fit_full, u_li = u_li, COMBO_stan = COMBO_stan, title = "B: [0-3] years")

m1_forest <- gen_forest(fit_full = m0_fit_1, u_li = u_li, COMBO_stan = COMBO_stan_1, title = "A: [0-1] years")
m2_forest <- gen_forest(fit_full = m0_fit_2, u_li = u_li, COMBO_stan = COMBO_stan_2, title = "B: (1-2] years")
m3_forest <- gen_forest(fit_full = m0_fit_3, u_li = u_li, COMBO_stan = COMBO_stan_3, title = "C: (2-3] years")
m4_forest <- gen_forest(fit_full = m0_fit_4, u_li = u_li, COMBO_stan = COMBO_stan_4, title = "A: [0-2] years")

ggsave(file = "time_stratified_forest_plots.pdf",
       device = "pdf",
       plot =  m1_forest$plot /
         m2_forest$plot /
         m3_forest$plot +
  plot_layout(ncol = 1),
  width = 50, height = 35,
  units = "cm"
)

ggsave(
  file = "m0_forest_plot.pdf",
  device = "pdf",
  plot = m4_forest$plot / m0_forest_all$plot +
    plot_layout(ncol = 1),
  width = 50, height = (35/3)*2,
  units = "cm"
)

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
cols <- unname(palette.colors(palette = "Okabe-Ito")[c(7,2,3,8)])
cols_net <- c("blue", "aquamarine", "darkgreen")

year_prev_plot <- ggplot(data = u_li_all) +
  geom_pointrange(data = COMBO_stan, aes(x = time/12, y = prev, ymin = prev_lower, ymax = prev_upper, fill = net, col = net),
                  shape = 1, alpha = 0.5, size = 0.25, inherit.aes = FALSE, position = position_jitter(width = 0.075)) +
  geom_ribbon(aes(x = years, ymin = l_p, ymax = u_p, y = m_p,
                  group = interaction(study_place, net), fill = net),
              alpha = 0.1) +
  geom_line(linewidth = 1,
            aes(x = years, y = m_p, group = interaction(study_place, net), col = net)) +
  geom_boxplot(data = unique(COMBO_stan[,c("cluster", "study_place", "BL_prev")]) |>
                 mutate(years = 0),
               aes(x = years,
                   y = BL_prev,
                   group = study_place),
               alpha = 0.1, width = 0.25,
               position = "identity") +
  geom_point(data = unique(COMBO_stan[,c("study_place", "mean_BL_prev")]) |> mutate(years = 0),
             aes(x = years, y = mean_BL_prev), size = 3) +
  facet_wrap(~ factor(study_place, levels = c("Tanzania (2015)", "Uganda (2017)", "Tanzania (2019)", "Benin (2020)")),
             nrow = 1) +
  scale_linetype_manual(name = "", values = c(1, 2)) +
  ylab("Malaria prevalence") + xlab("Years post trial net distribution") +
  scale_y_continuous(limits = c(-0.01, 1.01), breaks = seq(0, 1, 0.2), labels = scales::percent_format()) +
  scale_x_continuous(limits = c(-0.25, 3.1), breaks = 0:3) +
  scale_colour_manual(values = cols_net, name = "ITN") +
  scale_fill_manual(values = cols_net, name = "ITN")

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
  theme(legend.position = c(0.625, 0.875)) +
  scale_y_continuous(limits = c(-1.5, 0.5), breaks = seq(-1.5, 0.5, 0.5))

######################################################
##### calculating the times the upper ORs pass 1 #####
######################################################
u_li_all |> group_by(study_place, i) |> filter(o_r_u > 0) |> group_by(study_place) |> filter(years == min(years))

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

u_li_comp <- u_li

u_li_comp_m0_3 <- u_li_comp |> mutate(model = "model 0", timespan = "3-years", model_timespan = "model 0, 3-years")
u_li_comp_m1_3 <- u_li_comp |> mutate(model = "model 1", timespan = "3-years", model_timespan = "model 1, 3-years")
u_li_comp_m0_2 <- u_li_comp |> mutate(model = "model 0", timespan = "2-years", model_timespan = "model 0, 2-years")
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

u_li_comp_m0_2 <- cbind(u_li_comp_m0_2,
                        apply(
                          as_draws_matrix(m0_fit_4$draws("theta_l"))[,u_li_comp_m0_2$l] + as_draws_matrix(m0_fit_4$draws("theta_li"))[,u_li_comp_m0_2$li], 
                              2, quantile, probs = probs_in) |> t() |> as.data.frame() |> rename("l_o_r" = 1, "m_o_r" = 2, "u_o_r" = 3)
                        ) |> 
  cbind(apply(
    as_draws_matrix(m0_fit_4$draws("theta_l"))[,u_li_comp_m0_2$l], 
    2, quantile, probs = probs_in) |> t() |> as.data.frame() |> rename("l_o_r_pooled" = 1, "m_o_r_pooled" = 2, "u_o_r_pooled" = 3)
  )

u_li_comp_m0_3 <- cbind(u_li_comp_m0_3,
                        apply(
                          as_draws_matrix(m0_fit_full$draws("theta_l"))[,u_li_comp_m0_3$l] + as_draws_matrix(m0_fit_full$draws("theta_li"))[,u_li_comp_m0_3$li], 
                          2, quantile, probs = probs_in) |> t() |> as.data.frame() |> rename("l_o_r" = 1, "m_o_r" = 2, "u_o_r" = 3)
) |> 
  cbind(apply(
    as_draws_matrix(m0_fit_full$draws("theta_l"))[,u_li_comp_m0_3$l], 
    2, quantile, probs = probs_in) |> t() |> as.data.frame() |> rename("l_o_r_pooled" = 1, "m_o_r_pooled" = 2, "u_o_r_pooled" = 3))

u_li_pooled_m0_3 <- u_li_comp_m0_3[c(1, 4), c("net", "l_o_r_pooled", "m_o_r_pooled", "u_o_r_pooled", "model_timespan")] |> mutate(study_place = "Pooled")
u_li_pooled_m1_3 <- u_li_comp_m1_3[c(1, 4), c("net", "l_o_r_pooled", "m_o_r_pooled", "u_o_r_pooled", "model_timespan")] |> mutate(study_place = "Pooled")
u_li_pooled_m0_2 <- u_li_comp_m0_2[c(1, 4), c("net", "l_o_r_pooled", "m_o_r_pooled", "u_o_r_pooled", "model_timespan")] |> mutate(study_place = "Pooled")
u_li_pooled_m1_2 <- u_li_comp_m1_2[c(1, 4), c("net", "l_o_r_pooled", "m_o_r_pooled", "u_o_r_pooled", "model_timespan")] |> mutate(study_place = "Pooled")

u_li_comp <- rbind(u_li_comp_m0_2, u_li_comp_m0_3, u_li_comp_m1_2, u_li_comp_m1_3)

u_li_comp_pbo <- subset(u_li_comp, net == "Pyrethroid-PBO")
u_li_comp_pbo$study_place <- factor(u_li_comp_pbo$study_place, levels = c("Pooled", "Tanzania (2015)", "Uganda (2017)", "Tanzania (2019)"))
u_li_comp_pp <- subset(u_li_comp, net == "Pyrethroid-pyrrole")
u_li_comp_pp$study_place <- factor(u_li_comp_pp$study_place, levels = c("Pooled", "Tanzania (2019)", "Benin (2020)"))

# checking the log OR for Tanzania (2019)
subset(u_li_comp_pbo, study_place == "Tanzania (2019)")

log_o_r_plot_m1_m0_y2_y3 <- ggplot(data = u_li_comp_pbo)  +
  scale_y_discrete(drop = FALSE) +
  geom_point(aes(x = m_o_r, y = study_place, col = model_timespan), position = position_dodge(width = 0.25), size = 5, alpha = 0.75) +
  geom_errorbar(aes(xmin = l_o_r, xmax = u_o_r, y = study_place, col = model_timespan), position = position_dodge(width = 0.25), height = 0, alpha = 0.75, orientation = "y") +
  geom_polygon(data = data.frame(x = unname(unlist(c(u_li_pooled_m1_3[1, c(2, 3, 4, 3)]))), y = c(1 + 0.1, 1.1 + 0.1, 1 + 0.1, 0.9 + 0.1), model_timespan = u_li_pooled_m1_3[1, "model_timespan"]),
               aes(x = x, y = y, fill = model_timespan),
               inherit.aes = FALSE, alpha = 0.5) +
  geom_polygon(data = data.frame(x = unname(unlist(c(u_li_pooled_m0_3[1, c(2, 3, 4, 3)]))), y = c(1 - 0.05, 1.1 - 0.05, 1 - 0.05, 0.9 - 0.05), model_timespan = u_li_pooled_m0_3[1, "model_timespan"]),
               aes(x = x, y = y, fill = model_timespan), alpha = 0.5,
               inherit.aes = FALSE)  +
  geom_polygon(data = data.frame(x = unname(unlist(c(u_li_pooled_m1_2[1, c(2, 3, 4, 3)]))), y = c(1 + 0.05, 1.1 + 0.05, 1 + 0.05, 0.9 + 0.05), model_timespan = u_li_pooled_m1_2[1, "model_timespan"]),
               aes(x = x, y = y, fill = model_timespan),
               inherit.aes = FALSE, alpha = 0.5) +
  geom_polygon(data = data.frame(x = unname(unlist(c(u_li_pooled_m0_2[1, c(2, 3, 4, 3)]))), y = c(1 - 0.1, 1.1 - 0.1, 1 - 0.1, 0.9 - 0.1), model_timespan = u_li_pooled_m0_2[1, "model_timespan"]),
               aes(x = x, y = y, fill = model_timespan), alpha = 0.5,
               inherit.aes = FALSE)  +
  scale_x_continuous(limits = c(-1.5, 0.25), breaks = seq(-1.5, 0.5, 0.5)) +
  xlab("Log odds ratio") + ylab("") +
  geom_vline(xintercept = 0, linetype = 2) +
  theme(axis.ticks.y = element_blank(),
        panel.border = element_blank(),
        axis.line.x = element_line(colour = "black"),
        legend.position = "none") +
  scale_colour_manual(values = c("darkred", "darkblue", "red", "skyblue"), name = "") +
  scale_fill_manual(values = c("darkred", "darkblue", "red", "skyblue"), name = "") +
  facet_wrap(~net) + labs(tag = "C") +
  ggplot(data = u_li_comp_pp)  +
  scale_y_discrete(drop = FALSE) +
  geom_point(aes(x = m_o_r, y = study_place, col = model_timespan), position = position_dodge(width = 0.25), size = 5, alpha = 0.75) +
  geom_errorbar(aes(xmin = l_o_r, xmax = u_o_r, y = study_place, col = model_timespan), position = position_dodge(width = 0.25), height = 0, alpha = 0.75, orientation = "y") +
  geom_polygon(data = data.frame(x = unname(unlist(c(u_li_pooled_m1_3[2, c(2, 3, 4, 3)]))), y = c(1 + 0.1, 1.1 + 0.1, 1 + 0.1, 0.9 + 0.1), model_timespan = u_li_pooled_m1_3[2, "model_timespan"]),
               aes(x = x, y = y, fill = model_timespan),
               inherit.aes = FALSE, alpha = 0.5) +
  geom_polygon(data = data.frame(x = unname(unlist(c(u_li_pooled_m0_3[2, c(2, 3, 4, 3)]))), y = c(1 - 0.05, 1.1 - 0.05, 1 - 0.05, 0.9 - 0.05), model_timespan = u_li_pooled_m0_3[2, "model_timespan"]),
               aes(x = x, y = y, fill = model_timespan), alpha = 0.5,
               inherit.aes = FALSE)  +
  geom_polygon(data = data.frame(x = unname(unlist(c(u_li_pooled_m1_2[2, c(2, 3, 4, 3)]))), y = c(1 + 0.05, 1.1 + 0.05, 1 + 0.05, 0.9 + 0.05), model_timespan = u_li_pooled_m1_2[2, "model_timespan"]),
               aes(x = x, y = y, fill = model_timespan),
               inherit.aes = FALSE, alpha = 0.5) +
  geom_polygon(data = data.frame(x = unname(unlist(c(u_li_pooled_m0_2[2, c(2, 3, 4, 3)]))), y = c(1 - 0.1, 1.1 - 0.1, 1 - 0.1, 0.9 - 0.1), model_timespan = u_li_pooled_m0_2[2, "model_timespan"]),
               aes(x = x, y = y, fill = model_timespan), alpha = 0.5,
               inherit.aes = FALSE) +
  scale_x_continuous(limits = c(-1.5, 0.25), breaks = seq(-1.5, 0.5, 0.5)) +
  xlab("Log odds ratio") + ylab("") +
  geom_vline(xintercept = 0, linetype = 2) +
  theme(axis.ticks.y = element_blank(),
        panel.border = element_blank(),
        axis.line.x = element_line(colour = "black"),
        legend.position = c(0.175, 0.925)) +
  scale_colour_manual(values = c("darkred", "darkblue", "red", "skyblue"), name = "") +
  scale_fill_manual(values = c("darkred", "darkblue", "red", "skyblue"), name = "") +
  facet_wrap(~net)

ggsave(file = "m0_m1_forest_plot_comb.pdf",
       device = "pdf",
  plot = (year_prev_plot + labs(tag = "A")) / ((m1_or_plot + labs(tag = "B")) + log_o_r_plot_m1_m0_y2_y3),
  width = 70, height = 40,
  units = "cm"
  )

###########################################
##### extracting the parameter values #####
###########################################

m0_0_3_params <- get_params(fit = m0_fit_full)
m0_0_2_params <- get_params(fit = m0_fit_4)
m0_0_1_params <- get_params(fit = m0_fit_1)
m0_1_2_params <- get_params(fit = m0_fit_2)
m0_2_3_params <- get_params(fit = m0_fit_3)

m0_params_0_3_df <- data.frame("params" = m0_0_3_params$params, "params_description" = m0_0_3_params$params_description, "values" = m0_0_3_params$values)
m0_params_0_2_df <- data.frame("params" = m0_0_2_params$params,"params_description" = m0_0_2_params$params_description, "values" = m0_0_2_params$values)
m0_params_0_1_df <- data.frame("params" = m0_0_1_params$params,"params_description" = m0_0_1_params$params_description, "values" = m0_0_1_params$values)
m0_params_1_2_df <- data.frame("params" = m0_1_2_params$params,"params_description" = m0_1_2_params$params_description, "values" = m0_1_2_params$values)
m0_params_2_3_df <- data.frame("params" = m0_2_3_params$params,"params_description" = m0_2_3_params$params_description, "values" = m0_2_3_params$values)

write.csv(m0_params_0_3_df, file = "m0_params_0_3_df.csv")
write.csv(m0_params_0_2_df, file = "m0_params_0_2_df.csv")
write.csv(m0_params_0_1_df, file = "m0_params_0_1_df.csv")
write.csv(m0_params_1_2_df, file = "m0_params_1_2_df.csv")
write.csv(m0_params_2_3_df, file = "m0_params_2_3_df.csv")

# model 1
m1_params <- get_params(fit = m1_fit_full)
m1_params_df <- data.frame("params" = m1_params$params, params_description = m1_params$params_description, "values" = m1_params$values)
m1_params_df <- rbind(m1_params_df,
                   data.frame(
                     "params" = c("$\\kappa^{\\text{pooled}}$", "$\\kappa_{1}$", "$\\kappa_{2}", "$\\kappa_{3}$", "$\\kappa_{4}",
                                  "$\\kappa^{\\text{net}}_{1}$", "$\\kappa^{\\text{net}}_{2}$", "$\\kappa^{\\text{study}}_{2}$",
                                  "$\\kappa^{\\text{study}}_{3}$", "$\\kappa^{\\text{study}}_{4}$", "$\\kappa^{\\text{study}}_{5}$",
                                  "$\\kappa^{\\text{study}}_{6}$"),
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
                                              "Benin (2020) pyrethroid-pyrrole year effect interaction"),
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
                                  get_quantiles(model_fit = m1_fit_full, param = "kappa_li", dim = 6)
                                  )))

write.csv(m1_params_df, file = "m1_params_df.csv")
