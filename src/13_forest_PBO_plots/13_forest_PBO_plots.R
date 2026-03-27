library(rstan)
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

orderly2::orderly_shared_resource("m1_re.stan")

orderly2::orderly_dependency("12_forest_PBO_fit",
                             "latest()",
                             c("m1_re_PBO_subset_fit_full.rds",
                               "COMBO_stan_PBO_subset.rds",
                               "data_in_PBO_subset.rds")
                             )

m1_fit_full <- readRDS(file = "m1_re_PBO_subset_fit_full.rds")

m1_stan_model <- rstan::stan_model(file = "m1_re.stan")

COMBO_stan <- readRDS(file = "COMBO_stan_PBO_subset.rds")

data_in_full <- readRDS("data_in_PBO_subset.rds")

iter <- 8000
warmup <- floor(iter/2)

lower <- 0.025
upper <- 0.975

###########################
##### m1 forest plots #####
###########################

years <- seq(0, 3, 0.1)

u_li_pred <- unique(COMBO_stan[,c("l", "li", "i", "net", "study")]) |> arrange(li) |>
  mutate(study_place = case_when(
    str_detect(study, "Protopopoff") ~ "Tanzania (2015)",
    str_detect(study, "Staedke") ~ "Uganda (2017)",
    str_detect(study, "Mosha") ~ "Tanzania (2019)",
    str_detect(study, "Accrombessi") ~ "Benin (2020)"))

u_li_all <- u_li_pred[rep(1:nrow(u_li_pred), length(years)),] |> mutate(years = rep(years, each = nrow(u_li_pred)), ij = i)

pred_data <- data_in_full |>
  purrr::list_assign(gq = 1,
                     N_gq = nrow(u_li_all),
                     N_ij_gq = length(unique(u_li_all$ij)), # the same random effect for all the clusters in each study
                     N_ij_unq_gq = length(unique(u_li_all$i)),
                     pmat_ij_gq = fastDummies::dummy_cols(u_li_all$ij)[,-1] |> as.matrix(),
                     r_id_gq = seq(1, length(unique(u_li_all$i))),
                     pmat_i_gq = fastDummies::dummy_cols(u_li_all$i)[,-1] |> as.matrix(),
                     pmat_l_gq = fastDummies::dummy_cols(u_li_all$l)[,-c(1, 2)] |> as.matrix(),
                     pmat_li_gq = fastDummies::dummy_cols(u_li_all$li)[,-c(1, 2)] |> as.matrix(),

                     pmat_i_pooled_gq = matrix(0, nrow = nrow(u_li_all), ncol = length(unique(u_li_all$i))),
                     pmat_li_pooled_gq = matrix(0, nrow = nrow(u_li_all), ncol = length(unique(u_li_all$li)) - 1),

                     ij_train_gq = rep(0, length(unique(u_li_all$ij))),
                     ij_unq_train_gq = rep(1, length(unique(u_li_all$i))),
                     years_gq = u_li_all$years)

model_predict <- rstan::gqs(m1_stan_model,
                            draws = as.matrix(m1_fit_full),
                            data = pred_data,
                            seed = 123)

u_li_all <- cbind(
  cbind(
    cbind(u_li_all,
          apply(rstan::extract(model_predict, "o_r_pooled")$o_r_pooled, 2, quantile, probs = c(lower, 0.5, upper)) |>
            t() |> as.data.frame() |> rename("o_r_pooled_l" = 1, "o_r_pooled_m" = 2, "o_r_pooled_u" = 3)
    ),
    apply(rstan::extract(model_predict, "o_r")$o_r, 2, quantile, probs = c(lower, 0.5, upper)) |>
      t() |> as.data.frame() |> rename("o_r_l" = 1, "o_r_m" = 2, "o_r_u" = 3)
  ),
  apply(rstan::extract(model_predict, "inv_logit_posterior")$inv_logit_posterior, 2, quantile, probs = c(lower, 0.5, upper)) |>
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
  ylab("Odds ratio of infection in trial ITN\nclusters relative to pyrethroid-only clusters") +
  scale_fill_manual(values = cols, name = "") +
  scale_colour_manual(values = cols, name = "") +
  geom_hline(yintercept = 1, linetype = 2) +
  theme(legend.position = c(0.925, 0.925))

######################################################
##### calculating the times the upper ORs pass 1 #####
######################################################
u_li_all |> group_by(study_place, i) |> filter(o_r_u > 1) |> group_by(study_place) |> filter(years == min(years))

u_li_all |> group_by(study_place, i) |> filter(o_r_pooled_u > 1) |> group_by(study_place) |> filter(years == min(years))

############################
##### time odds ratios #####
############################

quantile(exp(rstan::extract(m1_fit_full, "kappa")[[1]]), probs = c(lower, 0.5, upper)) |> round(digits = 2)
quantile(exp(rstan::extract(m1_fit_full, "kappa")[[1]] + rstan::extract(m1_fit_full, "kappa_l")[[1]][,1]), probs = c(lower, 0.5, upper)) |> round(digits = 2)
quantile(exp(rstan::extract(m1_fit_full, "kappa")[[1]] + rstan::extract(m1_fit_full, "kappa_l")[[1]][,2]), probs = c(lower, 0.5, upper)) |> round(digits = 2)

########################################
##### calculating the mean log ORs #####
########################################

u_li_comp <- unique(COMBO_stan[,c("l", "li", "i", "net", "study")]) |> arrange(li) |> filter(l!=1) |>
  mutate(l_in = l - 1,
         li_in = row_number(),
         study_place = case_when(
           str_detect(study, "Protopopoff") ~ "Tanzania (2015)",
           str_detect(study, "Staedke") ~ "Uganda (2017)",
           str_detect(study, "Mosha") ~ "Tanzania (2019)",
           str_detect(study, "Accrombessi") ~ "Benin (2020)"))

u_li_comp_m1_3 <- u_li_comp |> mutate(model = "model 1", timespan = "3-years", model_timespan = "model 1, 3-years")
u_li_comp_m1_2 <- u_li_comp |> mutate(model = "model 1", timespan = "2-years", model_timespan = "model 1, 2-years")

# calculating the means over 2 and 3 years
for(i in 1:nrow(u_li_comp)){
  a <- rstan::extract(m1_fit_full, "theta_l")[[1]][, u_li_comp[i, "l_in"]]
  b <- rstan::extract(m1_fit_full, "theta_li")[[1]][, u_li_comp[i, "li_in"]]
  c <- rstan::extract(m1_fit_full, "kappa_l")[[1]][, u_li_comp[i, "l_in"]]

  # 3 years
  mean_o_r_3 <- quantile((a * 3 + b * 3 + (c * 3^2)/2)/3, probs = c(lower, 0.5, upper))
  mean_o_r_pooled_3 <- quantile((a * 3 + (c * 3^2)/2)/3, probs = c(lower, 0.5, upper))

  mean_o_r_2 <- quantile((a * 2 + b * 2 + (c * 2^2)/2)/2, probs = c(lower, 0.5, upper))
  mean_o_r_pooled_2 <- quantile((a * 2 + (c * 2^2)/2)/2, probs = c(lower, 0.5, upper))

  u_li_comp_m1_3[i, "l_o_r"] <- mean_o_r_3[1]
  u_li_comp_m1_3[i, "m_o_r"] <- mean_o_r_3[2]
  u_li_comp_m1_3[i, "u_o_r"] <- mean_o_r_3[3]

  u_li_comp_m1_3[i, "l_o_r_pooled"] <- mean_o_r_pooled_3[1]
  u_li_comp_m1_3[i, "m_o_r_pooled"] <- mean_o_r_pooled_3[2]
  u_li_comp_m1_3[i, "u_o_r_pooled"] <- mean_o_r_pooled_3[3]

  u_li_comp_m1_2[i, "l_o_r"] <- mean_o_r_2[1]
  u_li_comp_m1_2[i, "m_o_r"] <- mean_o_r_2[2]
  u_li_comp_m1_2[i, "u_o_r"] <- mean_o_r_2[3]

  u_li_comp_m1_2[i, "l_o_r_pooled"] <- mean_o_r_pooled_2[1]
  u_li_comp_m1_2[i, "m_o_r_pooled"] <- mean_o_r_pooled_2[2]
  u_li_comp_m1_2[i, "u_o_r_pooled"] <- mean_o_r_pooled_2[3]

}

u_li_pooled_m1_3 <- u_li_comp_m1_3[c(1, 4), c("net", "l_o_r_pooled", "m_o_r_pooled", "u_o_r_pooled", "model_timespan")] |> mutate(study_place = "Pooled")
u_li_pooled_m1_2 <- u_li_comp_m1_2[c(1, 4), c("net", "l_o_r_pooled", "m_o_r_pooled", "u_o_r_pooled", "model_timespan")] |> mutate(study_place = "Pooled")

u_li_comp <- rbind(u_li_comp_m1_2, u_li_comp_m1_3)

u_li_comp_pbo <- subset(u_li_comp, net == "Pyrethroid-PBO")
u_li_comp_pbo$study_place <- factor(u_li_comp_pbo$study_place, levels = c("Pooled", "Uganda (2017)", "Tanzania (2019)"))
u_li_comp_pp <- subset(u_li_comp, net == "Pyrethroid-pyrrole")
u_li_comp_pp$study_place <- factor(u_li_comp_pp$study_place, levels = c("Pooled", "Tanzania (2019)", "Benin (2020)"))

u_li_comp_pbo

# checking the log OR for Tanzania (2019)

log_o_r_plot_m1_subset_y2_y3 <- ggplot(data = u_li_comp_pbo)  +
  scale_y_discrete(drop = FALSE) +
  geom_point(aes(x = m_o_r, y = study_place, col = model_timespan), position = position_dodge(width = 0.25), size = 5, alpha = 0.75) +
  geom_errorbarh(aes(xmin = l_o_r, xmax = u_o_r, y = study_place, col = model_timespan), position = position_dodge(width = 0.25), height = 0, alpha = 0.75) +
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
  geom_errorbarh(aes(xmin = l_o_r, xmax = u_o_r, y = study_place, col = model_timespan), position = position_dodge(width = 0.25), height = 0, alpha = 0.75) +
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
  params <- c(params, "$\\alpha_{i}$", "$\\alpha_{i}$", "$\\alpha_{i}$")
  params_description <- c(params_description, "study effect (Uganda (2017))", "study effect (Tanzania (2019))", "study effect (Benin (2020))")
  for(i in 1:3){
    values <- c(values, get_quantiles(fit, "alpha_i", dim = i))
  }
  # net effect
  params <- c(params, "$\\theta_{l}$")
  params_description <- c(params_description, "pyrethroid-PBO pooled treatment effect")
  values <- c(values, get_quantiles(fit, "theta_l", dim = 1))

  params_description <- c(params_description, "pyrethroid-PBO-Uganda (2017) interaction", "pyrethroid-PBO-Tanzania (2019) interaction")
  params <- c(params, "$\\theta_{li}$", "$\\theta_{li}$")
  for(i in 1:2){
    values <- c(values, get_quantiles(fit, "theta_li", dim = i))
  }

  params <- c(params, "$\\theta_{l}$")
  params_description <- c(params_description, "pyrethroid-pyrrole pooled treatment effect")
  values <- c(values, get_quantiles(fit, "theta_l", dim = 2))

  params_description <- c(params_description, "pyrethroid-PBO-Tanzania (2019) interaction", "pyrethroid-pyrrole-Benin (2020) interaction")
  params <- c(params, "$\\theta_{li}$", "$\\theta_{li}$")
  values <- c(values, get_quantiles(fit, "theta_li", dim = 3), get_quantiles(fit, "theta_li", dim = 4))

  params_description <- c(params_description, "Uganda (2017) cluster random effect standard deviation",
                          "Tanzania (2019) cluster random effect standard deviation", "Benin (2020) cluster random effect standard deviation")

  params <- c(params, "$\\sigma_{i}$", "$\\sigma_{i}$", "$\\sigma_{i}$")
  values <- c(values, get_quantiles(fit, "sigma_e_r", dim = 1), get_quantiles(fit, "sigma_e_r", dim = 2),
              get_quantiles(fit, "sigma_e_r", dim = 3))

  params_description <- c(params_description, "pyrethroid-PBO between study random treatment effect standard deviation", "pyrethroid-pyrrole between study random treatment effect standard deviation")
  params <- c(params, "$\\tau_{l}", "$\\tau_{l}")
  values <- c(values, get_quantiles(fit, "tau_sd_li", dim = 1), get_quantiles(fit, "tau_sd_li", dim = 2))

  return(list("params" = params, "params_description" = params_description, "values" = values))
}

m1_params <- get_params(fit = m1_fit_full)
m1_params_df <- data.frame("params" = m1_params$params, params_description = m1_params$params_description, "values" = m1_params$values)
m1_params_df <- rbind(m1_params_df,
                      data.frame(
                        "params" = c("$\\kappa$", "$\\kappa_{l}$", "$\\kappa_{l}"),
                        "params_description" = c("year effect", "pyrethroid-PBO year effect interaction", "pyrethroid-pyrrole year effect interaction"),
                        "values" = c(get_quantiles(model_fit = m1_fit_full, param = "kappa", dim = NULL),
                                     get_quantiles(model_fit = m1_fit_full, param = "kappa_l", dim = 1),
                                     get_quantiles(model_fit = m1_fit_full, param = "kappa_l", dim = 2))))

write.csv(m1_params_df, file = "m1_params_df_subset.csv")

