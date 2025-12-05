##############################
##### malariasim results #####
##############################

library(tidyverse)
library(malariasimulation)
library(patchwork)

ggplot2::theme_set(theme_bw() +
                     theme(text = element_text(size = 14),
                           legend.text = element_text(size = 10),
                           legend.title = element_text(size = 10),
                           panel.background = element_blank(),
                           panel.grid.major = element_blank(),
                           panel.grid.minor = element_blank())
                   )


orderly2::orderly_dependency("6_malariasim", "latest()", c("sim_results.rds"))

msim_results <- readRDS(file = "sim_results.rds")

model_runs <- msim_results$model_runs

year <- 365
human_population <- 17500

retention <- 5 * year
rate <- 1 - exp(-1/retention)

vals <- msim_results$vals |> dplyr::mutate(EIR = msim_results$EIR |> unlist() |> as.vector(),
                                           base_net_m = base_net * (exp(-0 * 365 * rate) - exp(-3 * 365 * rate)) / (rate * (3 * 365 - 0 * 365)),
                                           trial_net_m = trial_net * (exp(-0 * 365 * rate) - exp(-3 * 365 * rate)) / (rate * (3 * 365 - 0 * 365)))

# calibrated EIR plot
cali_EIR_plot <- ggplot(data = vals |> filter(net_type == "pp" & trial_net_m == 0.75),
                        aes(x = base_prev, y = EIR, fill = factor(base_net_m))) +
  geom_point(size = 4, shape = 21) +
  scale_fill_brewer(name = "Mean baseline net\nuse over three years") +

  ylab("Calibrated entomological inoculation rate (EIR)") +
  xlab("Baseline prevalence") +
  scale_x_continuous(limits = c(0.05, 0.55), breaks = seq(0.1, 0.5, 0.1), labels = scales::percent) +
  scale_y_continuous(breaks = seq(0, 1000, 200), limits = c(0, 1000))

ggsave(file = "cali_EIR_plot.pdf",
       plot = cali_EIR_plot,
       dpi = "retina",
       width = 20, height = 10,
       units = "cm")

model_runs <- lapply(1:length(model_runs),
                     function(i, model_runs){
                       model_runs[[i]] <- model_runs[[i]] |>
                         filter(timestep >= floor(8.75 * year) & timestep <= 12 * year) |>
                         mutate(prev = n_detect_lm_0_36500 / n_age_0_36500,
                                n_inc_clinical_pp = n_inc_clinical_0_36500 / n_age_0_36500,
                                base_prev = round(vals$base_prev[i], digits = 2),
                                base_net = round(vals$base_net[i], digits = 2),
                                trial_net = round(vals$trial_net[i], digits = 2),
                                net_type = vals$net_type[i],
                                timestep = timestep - floor(9 * year),
                                years = timestep / year,
                                year_floor = floor(years),
                                EIR_pp_py = EIR_gamb / human_population * year
                         )},
                     model_runs = model_runs) |>
  bind_rows()


model_runs <- left_join(model_runs,
                        subset(model_runs, net_type == "pyr") |>
                          select(prev, base_prev, base_net, trial_net, timestep, EIR_pp_py, n_inc_clinical_pp) |>
                          rename(p_prev = prev,
                                 p_EIR_pp_py = EIR_pp_py,
                                 p_n_inc_clinical_pp = n_inc_clinical_pp),
                        by = c("base_prev", "base_net", "trial_net", "timestep")) |>
  mutate(efficacy = (p_prev - prev) / p_prev,
         abs_efficacy = p_prev - prev,
         efficacy_EIR = (p_EIR_pp_py - EIR_pp_py) / p_EIR_pp_py,
         abs_efficacy_EIR = p_EIR_pp_py - EIR_pp_py,
         base_net_m = round(base_net * (exp(-0 * 365 * rate) - exp(-3 * 365 * rate)) / (rate * (3 * 365 - 0 * 365)), digits = 2),
         trial_net_m = round(trial_net * (exp(-0 * 365 * rate) - exp(-3 * 365 * rate)) / (rate * (3 * 365 - 0 * 365)), digits = 2),
         base_prev_plot = paste0("Baseline prevalence: ", base_prev * 100, "%"),
         trial_net_m_plot = paste0(trial_net_m*100, "%"),
         base_net_m_plot = paste0(base_net_m*100, "%"),
         net_type = ifelse(net_type == "pp", "Pyrethroid-pyrrole", "Pyrethroid-only"))

### baseline prevalence plots
# prevalence
prev_plots_all <- ggplot(data = model_runs |>
         filter(trial_net_m == 0.5 & round(base_prev, digits = 2) %in% round(seq(0.1, 0.55, 0.15), digits = 2)),
       aes(x = years, y = prev,
           col = base_net_m_plot, group = interaction(base_net_m_plot, net_type),
           linetype = net_type)) +
  geom_line(linewidth = 0.75) +
  ylab("Simulated prevalence") +
  xlab("Years post trial ITN distribution") +
  scale_linetype_manual(values = c(3, 1), name = "Simulated\nprevalence ITN") +
  scale_colour_brewer(name = "Baseline\nnet use", palette = "BrBG") +
  scale_y_continuous(labels = scales::percent, limits = c(0, 0.7), breaks = seq(0, 0.6, 0.2)) +
  geom_vline(xintercept = 0, linetype = 2) +
  facet_wrap(~ base_prev_plot, nrow = 1) +

  ggplot(data = model_runs |>
           filter(trial_net_m == 0.5 & net_type == "Pyrethroid-pyrrole" & round(base_prev, digits = 2) %in% round(seq(0.1, 0.55, 0.15), digits = 2)),
         aes(x = years, y = efficacy,
             col = base_net_m_plot, group = base_net_m_plot)) +
  geom_line(linewidth = 0.5) +
  scale_colour_brewer(name = "Baseline\nnet use", palette = "BrBG") +
  ylab("Relative efficacy") + xlab("Years post trial ITN distribution") +
  scale_y_continuous(limits = c(-0.1, 1), labels = scales::percent) +
  geom_vline(xintercept = 0, linetype = 2) +
  facet_wrap(~ base_prev_plot, nrow = 1) +

ggplot(data = model_runs |>
           filter(trial_net_m == 0.5 & net_type == "Pyrethroid-pyrrole" & round(base_prev, digits = 2) %in% round(seq(0.1, 0.55, 0.15), digits = 2)),
         aes(x = years, y = abs_efficacy,
             col = base_net_m_plot, group = base_net_m_plot)) +
  geom_line(linewidth = 0.5) +
  scale_colour_brewer(name = "Baseline\nnet use", palette = "BrBG") +
  ylab("Absolute efficacy") + xlab("Years post trial ITN distribution") +
  scale_y_continuous(limits = c(-0.1, 0.125), labels = scales::percent) +
  geom_vline(xintercept = 0, linetype = 2) +
  facet_wrap(~ base_prev_plot, nrow = 1) +
  plot_layout(nrow = 3, axes = "collect") +
  plot_annotation(tag_levels = c("A"))

# baseline prevalence plots subset
prev_time_plot_p <- ggplot(data = model_runs |>
                             filter(trial_net_m == 0.5 & base_net == 0 & base_prev %in% c(0.1, 0.25, 0.4, 0.55)),
                           aes(x = years, y = prev,
                               col = base_prev, group = interaction(base_prev, net_type),
                               linetype = net_type)) +
  geom_line(linewidth = 0.4) +
  ylab("Simulated prevalence") +
  xlab("Years post trial ITN distribution") +
  scale_linetype_manual(values = c(3, 1), name = "Simulated\nprevalence ITN") +
  scale_colour_distiller(name = "Baseline\nprevalence", palette = "BrBG", labels = scales::percent) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 0.7), breaks = seq(0, 0.6, 0.2)) +
  geom_vline(xintercept = 0, linetype = 2)

abs_eff_plot_p <- ggplot(data = model_runs |>
                           filter(trial_net_m == 0.5 & base_net == 0 & net_type == "Pyrethroid-pyrrole"),
                         aes(x = years, y = abs_efficacy,
                             col = base_prev, group = base_prev_plot)) +
  geom_line(linewidth = 0.4) +
  scale_colour_distiller(name = "Baseline\nprevalence", palette = "BrBG", labels = scales::percent) +
  ylab("Absolute efficacy") + xlab("Years post trial ITN distribution") +
  scale_y_continuous(limits = c(-0.1, 0.25), labels = scales::percent) +
  geom_vline(xintercept = 0, linetype = 2)

eff_plot_p <- ggplot(data = model_runs |>
                       filter(trial_net_m == 0.5 & base_net == 0 & net_type == "Pyrethroid-pyrrole"),
                     aes(x = years, y = efficacy,
                         col = base_prev, group = base_prev_plot)) +
  geom_line(linewidth = 0.4) +
  scale_colour_distiller(name = "Baseline\nprevalence", palette = "BrBG", labels = scales::percent) +
  ylab("Relative efficacy") + xlab("Years post trial ITN distribution") +
  scale_y_continuous(limits = c(-0.1, 1), labels = scales::percent) +
  geom_vline(xintercept = 0, linetype = 2)

# baseline net use plots
time_plot_n <- ggplot(data = model_runs |>
                        filter(trial_net_m == 0.5 & base_prev %in% c(0.1, 0.5)),
                      aes(x = years, y = prev,
                          col = base_net_m, group = interaction(base_net_m, net_type),
                          linetype = net_type)) +
  geom_line(linewidth = 0.4) +
  scale_linetype_manual(values = c(3, 1), name = "Simulated\nprevalence ITN") +
  ylab("Simulated prevalence") +
  xlab("Years post trial ITN distribution") +
  scale_y_continuous(labels = scales::percent, limits = c(0, 0.7), breaks = seq(0, 0.6, 0.2)) +
  geom_vline(xintercept = 0, linetype = 2) +
  scale_colour_distiller(name = "Mean baseline\nnet use", palette = "RdBu", labels = scales::percent) +
  facet_wrap(~base_prev_plot, nrow = 2)

eff_plot_n <- ggplot(data = model_runs |> filter(trial_net_m == 0.5 &
                                                   net_type == "Pyrethroid-pyrrole" &
                                                   base_prev %in% c(0.1, 0.5)),
                     aes(x = years, y = efficacy,
                         col = base_net_m, group = factor(base_net_m))) +
  geom_line(linewidth = 0.4) +
  scale_colour_distiller(name = "Mean baseline\nnet use", palette = "RdBu", labels = scales::percent) +
  ylab("Relative efficacy") + xlab("Years post trial ITN distribution") +
  scale_y_continuous(limits = c(-0.1, 1), labels = scales::percent) +
  geom_vline(xintercept = 0, linetype = 2) +
  facet_wrap(~base_prev_plot, nrow = 2)

abs_eff_plot_n <- ggplot(data = model_runs |> filter(trial_net_m == 0.5 &
                                                       net_type == "Pyrethroid-pyrrole" &
                                                       base_prev %in% c(0.1, 0.5)),
                         aes(x = years, y = abs_efficacy,
                             col = base_net_m, group = factor(base_net_m))) +
  geom_line(linewidth = 0.4) +
  scale_colour_distiller(name = "Mean baseline\nnet use", palette = "RdBu", labels = scales::percent) +
  ylab("Absolute efficacy") + xlab("Years post trial ITN distribution") +
  scale_y_continuous(limits = c(-0.1, 0.25), labels = scales::percent) +
  geom_vline(xintercept = 0, linetype = 2) +
  facet_wrap(~base_prev_plot, nrow = 2)


ggsave(
  file = "malariasim_plots_prev.pdf",
  plot = (prev_time_plot_p +
    eff_plot_p  + abs_eff_plot_p) / (time_plot_n +
    eff_plot_n + abs_eff_plot_n) +
    plot_layout(guides = "collect", heights = c(1, 2)) +
    plot_annotation(tag_levels = 'A'),
  dpi = "retina",
  width = 35, height = 22.5,
  units = "cm"
)

# EIR
EIR_time_plot_p <- ggplot(data = model_runs |>
                            filter(trial_net_m == 0.5 & base_net == 0 & base_prev %in% c(0.1, 0.25, 0.4, 0.55)),
                          aes(x = years, y = EIR_pp_py,
                              col = base_prev, group = interaction(base_prev, net_type),
                              linetype = net_type)) +
  geom_line(linewidth = 1) +
  ylab("Simulated EIR\nper person per year") +
  xlab("Years post trial ITN distribution") +
  scale_linetype_manual(values = c(3, 1), name = "Simulated\nprevalence ITN") +
  scale_colour_distiller(name = "Baseline\nprevalence", palette = "BrBG", labels = scales::percent) +
  scale_y_sqrt(limits = c(0, 550), breaks = seq(0, 500, 100)) +
  geom_vline(xintercept = 0, linetype = 2)

EIR_abs_eff_plot_p <- ggplot(data = model_runs |>
                               filter(trial_net_m == 0.5 & base_net == 0 & net_type == "Pyrethroid-pyrrole"),
                             aes(x = years, y = abs_efficacy_EIR,
                                 col = base_prev, group = base_prev_plot)) +
  geom_line(linewidth = 1) +
  scale_colour_distiller(name = "Baseline\nprevalence", palette = "BrBG", labels = scales::percent) +
  ylab("EIR absolute efficacy") + xlab("Years post trial ITN distribution") +
  scale_y_continuous(limits = c(-10, 105), breaks = seq(0, 100, 25)) +
  geom_vline(xintercept = 0, linetype = 2)

EIR_eff_plot_p <- ggplot(data = model_runs |>
                           filter(trial_net_m == 0.5 & base_net == 0 & net_type == "Pyrethroid-pyrrole"),
                         aes(x = years, y = efficacy_EIR,
                             col = base_prev, group = base_prev_plot)) +
  geom_line(linewidth = 1) +
  scale_colour_distiller(name = "Baseline\nprevalence", palette = "BrBG", labels = scales::percent) +
  ylab("EIR relative efficacy") + xlab("Years post trial ITN distribution") +
  scale_y_continuous(limits = c(-0.1, 1), labels = scales::percent) +
  geom_vline(xintercept = 0, linetype = 2)

EIR_time_plot_n <- ggplot(data = model_runs |>
                        filter(trial_net_m == 0.5 & base_prev %in% c(0.1, 0.5)),
                      aes(x = years, y = EIR_pp_py,
                          col = base_net_m, group = interaction(base_net_m, net_type),
                          linetype = net_type)) +
  geom_line(linewidth = 1) +
  scale_linetype_manual(values = c(3, 1), name = "Simulated\nprevalence ITN") +
  ylab("Simulated EIR\nper person per year") +
  xlab("Years post trial ITN distribution") +
  scale_y_sqrt(limits = c(0, 140), breaks = seq(0, 125, 25)) +
  geom_vline(xintercept = 0, linetype = 2) +
  scale_colour_distiller(name = "Mean baseline\nnet use", palette = "RdBu", labels = scales::percent) +
  facet_wrap(~base_prev_plot, nrow = 2)

EIR_eff_plot_n <- ggplot(data = model_runs |> filter(trial_net_m == 0.5 &
                                                   net_type == "Pyrethroid-pyrrole" &
                                                   base_prev %in% c(0.1, 0.5)),
                     aes(x = years, y = efficacy_EIR,
                         col = base_net_m, group = factor(base_net_m))) +
  geom_line(linewidth = 1) +
  scale_colour_distiller(name = "Mean baseline\nnet use", palette = "RdBu", labels = scales::percent) +
  ylab("EIR relative efficacy") + xlab("Years post trial ITN distribution") +
  scale_y_continuous(limits = c(-0.1, 1), labels = scales::percent) +
  geom_vline(xintercept = 0, linetype = 2) +
  facet_wrap(~base_prev_plot, nrow = 2)

EIR_abs_eff_plot_n <- ggplot(data = model_runs |> filter(trial_net_m == 0.5 &
                                                       net_type == "Pyrethroid-pyrrole" &
                                                       base_prev %in% c(0.1, 0.5)),
                         aes(x = years, y = abs_efficacy_EIR,
                             col = base_net_m, group = factor(base_net_m))) +
  geom_line(linewidth = 1) +
  scale_colour_distiller(name = "Mean baseline\nnet use", palette = "RdBu", labels = scales::percent) +
  ylab("EIR absolute efficacy") + xlab("Years post trial ITN distribution") +
  #scale_y_continuous(limits = c(-5, 50)) +
  geom_vline(xintercept = 0, linetype = 2) +
  facet_wrap(~base_prev_plot, nrow = 2)

ggsave(
  file = "malariasim_plots_EIR.pdf",
  plot = (EIR_time_plot_p +
    EIR_eff_plot_p + EIR_abs_eff_plot_p) / (EIR_time_plot_n +
    EIR_eff_plot_n + EIR_abs_eff_plot_n) +
    plot_layout(guides = "collect", heights = c(1, 2)) +
    plot_annotation(tag_levels = 'A'),
  dpi = "retina",
  width = 35, height = 22.5,
  units = "cm"
)

