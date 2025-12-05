library(tidyverse); library(rio); library(skimr); library(GGally); library(cowplot); library(ggsignif); library(patchwork)

#setwd("/Volumes/malaria/Isaac/trials")

orderly2::orderly_dependency("1_data_cleaning", "latest()",
                             c("stan_data.rds"))

list2env(readRDS("stan_data.rds"), envir = .GlobalEnv)

###### NEW PAIRS PLOT, FUNCTION TO REPAEAT ACROSS TIMEPOINTS

unc_bin <- function(u, x, n){
  qbeta(u, 1 + x, 1 + n - x)
}

##### data #####
pairs_data <- COMBO_stan |>
  mutate(time = time,
         time_6m = cut(time,
                       breaks = c(0,6,12,18,24,30,36),
                       labels = c("≤6 m", "7-12 m", "13-18 m", "19-24 m", "25-30 m", "31-36 m")),
         i = factor(i,
                    levels = c(1,2,3,4),
                    labels = c("Protopopoff", "Staedke", "Mosha", "Accrombessi")),
         prev = prev_num / prev_denom,
         prev_num = prev_num,
         prev_denom = prev_denom,
         BL_prev = BL_prev_num/BL_prev_denom,
         BL_prev_l = unc_bin(lower, BL_prev_num, BL_prev_denom),
         BL_prev_u = unc_bin(upper, BL_prev_num, BL_prev_denom),
         BL_prev_num = BL_prev_num,
         BL_prev_denom = BL_prev_denom,
         net_use = net_use_num/net_use_denom,
         net_use_num = net_use_num,
         net_use_l = unc_bin(lower, net_use_num, net_use_denom),
         net_use_u = unc_bin(upper, net_use_num, net_use_denom),
         net_use_denom = net_use_denom,
         tr_net_use = tr_net_use_num/tr_net_use_denom,
         tr_net_use_l = unc_bin(lower, tr_net_use_num, tr_net_use_denom),
         tr_net_use_u = unc_bin(upper, tr_net_use_num, tr_net_use_denom),
         tr_net_use_num = tr_net_use_num,
         tr_net_use_denom = tr_net_use_denom,
         study_place = factor(case_when(
           str_detect(study, "Protopopoff") ~ "Tanzania (2015)",
           str_detect(study, "Staedke") ~ "Uganda (2017)",
           str_detect(study, "Mosha") ~ "Tanzania (2019)",
           str_detect(study, "Accrombessi") ~ "Benin (2020)"
         ), levels = c("Tanzania (2015)", "Uganda (2017)", "Tanzania (2019)", "Benin (2020)")),
         .keep = "none")

# colours
cols <- unname(palette.colors(palette = "Okabe-Ito")[c(7,2,3,8)])
cols_net <- c("blue", "aquamarine", "darkgreen")

##### plotting functions #####
pairs_theme <-
  theme(panel.border = element_rect(fill = NA, linewidth = 1),
        axis.line = element_blank(),
        plot.title = element_text(face = "bold"),
        strip.placement = "outside",
        strip.clip = "off",
        strip.text.x = element_text(vjust = +3.5),
        strip.text.y = element_text(vjust = -2),
        axis.text.x = element_text(angle = 90, vjust = 0.4),
        strip.background = element_blank(),
        plot.margin = unit(c(0.25, 0.25, 0.25, 0.25), "cm"),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())


pairs_plot_month <- function(months){

  df <- subset(pairs_data, time_6m == months)

lower_plots_custom <- function(data,
                               mapping,
                               mapping_eb,
                               data_smooth,
                               mapping_smooth,
                               formula,
                               ...) {


  ggplot(...) +
    geom_pointrange(data = data, mapping = mapping, alpha = 0.5, size = 0.15, linewidth = 0.025) +
    geom_errorbarh(data = data, mapping = mapping_eb, linewidth = 0.025) +
    scale_y_continuous(labels=scales::percent, limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
    scale_x_continuous(labels=scales::percent, limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
    pairs_theme +
    scale_colour_manual(values = c("Tanzania (2015)" = cols[1], "Uganda (2017)" = cols[2], "Tanzania (2019)" = cols[3], "Benin (2020)" = cols[4]), name = "")  +
    scale_fill_manual(values = c("Tanzania (2015)" = cols[1], "Uganda (2017)" = cols[2], "Tanzania (2019)" = cols[3], "Benin (2020)" = cols[4]), name = "") +
    scale_shape_manual(values = c("Tanzania (2015)" = 15, "Uganda (2017)" = 16, "Tanzania (2019)" = 17, "Benin (2020)" = 18), name = "") +
    geom_smooth(data = data_smooth,
                method = "glm",
                mapping = mapping_smooth,
                method.args = list(family = "binomial"),
                formula = formula,
                se = FALSE,
                linewidth = 0.75) +
    theme(legend.position = "none", text = element_text(size = 14))
}



bp_nu_plot <- lower_plots_custom(data = df,
                   mapping = aes(x = BL_prev, xmin = BL_prev_l, xmax = BL_prev_u,
                                 y = net_use, ymin = net_use_l, ymax = net_use_u,
                                 col = study_place, fill = study_place, shape = study_place, group = study_place
                                 ),
                   mapping_eb = aes(xmin = BL_prev_l, xmax = BL_prev_u, y = net_use, col = study_place),
                   data_smooth = df[,c("BL_prev", "net_use_num", "net_use_denom", "study_place")] %>%
                     rowwise() %>%
                     reframe(y = rep(c(1, 0), times = c(net_use_num, net_use_denom - net_use_num)), BL_prev = BL_prev, study_place = study_place),
                   mapping_smooth = aes(x = BL_prev, y = y, col = study_place, fill = study_place, group = study_place),
                   formula = y ~ x
                   ) +
  xlab("Baseline prevalence") + ylab("Baseline net use") +
  theme(text = element_text(size = 14))

bp_tnu_plot <- lower_plots_custom(data = df,
                                 mapping = aes(x = BL_prev, xmin = BL_prev_l, xmax = BL_prev_u,
                                               y = tr_net_use, ymin = tr_net_use_l, ymax = tr_net_use_u,
                                               col = study_place, fill = study_place, shape = study_place, group = study_place
                                 ),
                                 mapping_eb = aes(xmin = BL_prev_l, xmax = BL_prev_u, y = tr_net_use, col = study_place),
                                 data_smooth = df[,c("BL_prev", "tr_net_use_num", "tr_net_use_denom", "study_place")] %>%
                                   rowwise() %>%
                                   reframe(y = rep(c(1, 0), times = c(tr_net_use_num, tr_net_use_denom - tr_net_use_num)), BL_prev = BL_prev, study_place = study_place),
                                 mapping_smooth = aes(x = BL_prev, y = y, col = study_place, fill = study_place, group = study_place),
                                 formula = y ~ x
) +
  xlab("Baseline prevalence") + ylab("Trial net use") +
  theme(text = element_text(size = 14))

nu_tnu_plot <- lower_plots_custom(data = df,
                                  mapping = aes(x = net_use, xmin = net_use_l, xmax = net_use_u,
                                                y = tr_net_use, ymin = tr_net_use_l, ymax = tr_net_use_u,
                                                col = study_place, fill = study_place, shape = study_place, group = study_place
                                  ),
                                  mapping_eb = aes(xmin = net_use_l, xmax = net_use_u, y = tr_net_use, col = study_place),
                                  data_smooth = df[,c("net_use", "tr_net_use_num", "tr_net_use_denom", "study_place")] %>%
                                    rowwise() %>%
                                    reframe(y = rep(c(1, 0), times = c(tr_net_use_num, tr_net_use_denom - tr_net_use_num)), net_use = net_use,
                                            study_place = study_place),
                                  mapping_smooth = aes(x = net_use, y = y, col = study_place, fill = study_place, group = study_place),
                                  formula = y ~ x
) +
  xlab("Baseline net use") + ylab("Trial net use") +
  theme(text = element_text(size = 14))

diag_plots <- function(data,
                       mapping,
                       ...) {

  ggplot(data = data,
         mapping = mapping,
         ...) +
    geom_histogram(binwidth = 0.05, boundary = 0, closed = "left", ..., alpha = 0.85) +
    scale_x_continuous(labels=scales::percent, limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
    pairs_theme +
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          text = element_text(size = 14)) +
    scale_colour_manual(values = c("Tanzania (2015)" = cols[1], "Uganda (2017)" = cols[2], "Tanzania (2019)" = cols[3], "Benin (2020)" = cols[4]), name = "")  +
    scale_fill_manual(values = c("Tanzania (2015)" = cols[1], "Uganda (2017)" = cols[2], "Tanzania (2019)" = cols[3], "Benin (2020)" = cols[4]), name = "")
}

bnu_hist_plot <- diag_plots(data = df,
           mapping = aes(x = net_use, col = study_place, fill = study_place, shape = study_place, group = study_place
           )) +
  xlab("Baseline net use") +
  theme(axis.title.y = element_blank())

tnu_hist_plot <- diag_plots(data = df,
                            mapping = aes(x = tr_net_use, col = study_place, fill = study_place, shape = study_place, group = study_place
                            )) +
  xlab("Trial net use") +
  theme(axis.title.y = element_blank())

bp_hist_plot <- diag_plots(data = df,
                            mapping = aes(x = BL_prev, col = study_place, fill = study_place, shape = study_place, group = study_place
                            )) +
  xlab("Baseline prevalence") +
  ylab("Baseline prevalence")

out <- if(months != "13-18 m"){
  (bp_hist_plot + theme(legend.position = "none")) + plot_spacer() + plot_spacer() + (bp_nu_plot + theme(legend.position = "none")) +
    (bnu_hist_plot + theme(legend.position = "none")) + plot_spacer() +
    (bp_tnu_plot + theme(legend.position = "none")) + (nu_tnu_plot + theme(legend.position = "none")) + (tnu_hist_plot + theme(legend.position = "none")) +
    plot_layout(guides = "collect", axes = "collect", axis_titles = "collect")

} else{
  bp_hist_plot + plot_spacer() + plot_spacer() + bp_nu_plot + bnu_hist_plot + plot_spacer() +
  bp_tnu_plot + nu_tnu_plot + tnu_hist_plot +
    plot_layout(guides = "collect", axes = "collect", axis_titles = "collect")
}

return(out)
}

p1 <- cowplot::plot_grid(
  cowplot::plot_grid(pairs_plot_month("≤6 m"),
                         pairs_plot_month("7-12 m"),
                         pairs_plot_month("19-24 m"),
                         pairs_plot_month("25-30 m"),
                         pairs_plot_month("31-36 m"),
                         nrow = 3,
                         labels = c("A (0-6 months)","B (7-12 months)","C (19-24 months)","D (25-30 months)","E (31-36 months)"),
                     hjust = -3.1
                         ),
                         cowplot::get_legend(pairs_plot_month(months = "13-18 m")),
                         rel_widths = c(1, 0.1),
                         nrow = 1)

ggsave("pairs_plot_supplement.pdf",
       plot = p1,
       device = "pdf",
       width = 40, height = 50,
       units = "cm",
       bg="white")

## GANTT CHART OF TRIAL SURVEY TIMEPOINTS
df_gantt <- pairs_data %>%
  select(author = i,
         months = time,
         time_6m) %>%
  distinct() %>%
  mutate(author = factor(author, levels = c("Accrombessi","Mosha","Staedke","Protopopoff"),
                         labels = c("Benin (2020)", "Tanzania (2019)", "Uganda (2017)", "Tanzania (2015)")))

gantt_plot <- ggplot(df_gantt) +
  annotate("rect", xmin = 0, xmax = 6.5, ymin = 0, ymax = Inf, alpha = 0.35, fill = "lightgrey", col = "grey") +
  annotate("rect", xmin = 6.5, xmax = 12.5, ymin = 0, ymax = Inf, alpha = 0.35, fill = "lightgrey", col = "grey") +
  annotate("rect", xmin = 12.5, xmax = 18.5, ymin = 0, ymax = Inf, alpha = 0.7, fill = "darkgrey", col = "darkgrey") +
  annotate("rect", xmin = 18.5, xmax = 24.5, ymin = 0, ymax = Inf, alpha = 0.35, fill = "lightgrey", col = "darkgrey") +
  annotate("rect", xmin = 24.5, xmax = 30.5, ymin = 0, ymax = Inf, alpha = 0.35, fill = "lightgrey", col = "grey") +
  annotate("rect", xmin = 30.5, xmax = 36.5, ymin = 0, ymax = Inf, alpha = 0.35, fill = "lightgrey", col = "grey") +
  geom_vline(xintercept = 0, linetype = "dashed", col = "black") +
  geom_tile(aes(x = months, y = author, fill = author), width = 1, height = 0.7) +
  theme_classic() +
  theme(panel.border = element_rect(fill = NA, linewidth = 1),
        axis.line = element_blank(),
        legend.position = "none",
        #plot.margin = margin(unit(c(3,3,3,10),"cm")),
        axis.title.y = element_blank(),
        text = element_text(size = 14)) +
  scale_x_continuous(breaks = seq(0,36,6), limits = c(0,37)) +
  labs(x = "Months since net distribution") +
  scale_fill_manual(values = c("Tanzania (2015)" = cols[1], "Uganda (2017)" = cols[2], "Tanzania (2019)" = cols[3], "Benin (2020)" = cols[4]))

### differences in baseline

arm_comp_plot <- function(a = "BL_prev_num", b = "BL_prev_denom", label = "Baseline prevalence"){

  s <- function(...){as.symbol(...)}

  return(COMBO_stan |>
    mutate(author = factor(author, levels = c("Accrombessi","Mosha","Staedke","Protopopoff"),
                           labels = c("Benin (2020)", "Tanzania (2019)", "Uganda (2017)", "Tanzania (2015)")),
           l = factor(l, levels = c(1, 2, 3, 4), labels = c("Pyrethroid-only", "Pyrethroid-PBO", "Pyrethroid-pyrrole", "Pyrethroid-pyriproxyfen"))) |>
    select(cluster, author, l, !!s(a), !!s(b)) |>
    distinct() |>
  mutate(y = !!s(a) / !!s(b),
         y_lower = unc_bin(0.025, !!s(a), !!s(b)),
         y_upper = unc_bin(0.975, !!s(a), !!s(b))) |>
    mutate(label = label) |> select(-c(!!s(a), !!s(b))))
}

comp_plot_fun <- function(df, xlab_s){
  ggplot(data = df,
       aes(y = author,
           x = y,
           fill = factor(l),
           group = interaction(author, l))) +
  geom_boxplot(alpha = 0.9) +
  pairs_theme +
  scale_x_continuous(labels = scales::percent, limits = c(-0.025, 1.025), breaks = seq(0, 1, 0.25)) +
  theme(legend.title = element_blank(),
        legend.text = element_text(size = 14),
        strip.text = element_text(face = "bold", size = 10)) + ylab("Study") +
  xlab(xlab_s) +
  scale_fill_manual(values = cols_net) +
    guides(fill = guide_legend(override.aes = list(size = 0.25)))
}

comp_plot <-  comp_plot_fun(arm_comp_plot(a = "BL_prev_num", b = "BL_prev_denom", label = "Baseline prevalence"),
                            "Baseline prevalence") +

  comp_plot_fun(arm_comp_plot(a = "net_use_num", b = "net_use_denom", label = "Baseline net use"), "Baseline net use") +
  comp_plot_fun(COMBO_stan |>
           mutate(author = factor(author, levels = c("Accrombessi","Mosha","Staedke","Protopopoff"),
                                  labels = c("Benin (2020)", "Tanzania (2019)", "Uganda (2017)", "Tanzania (2015)")),
                  l = factor(l, levels = c(1, 2, 3, 4), labels = c("Pyrethroid-only", "Pyrethroid-PBO", "Pyrethroid-pyrrole", "Pyrethroid-pyriproxyfen")),
                  y = tr_net_use_num/tr_net_use_denom, y_lower = unc_bin(0.025, tr_net_use_num, tr_net_use_denom),
                  y_upper = unc_bin(0.975, tr_net_use_num, tr_net_use_denom), label = "Trial net use"
           ) |>
           select(cluster, author, l, y, y_lower, y_upper, label),
        "Trial net use") +
  plot_layout(guides = "collect",
              axes = "collect",
              axis_titles = "collect") +
  guide_area()


# main plot
p2 <- cowplot::plot_grid(gantt_plot,
                         cowplot::plot_grid(comp_plot, pairs_plot_month(months = "13-18 m"),
                         nrow = 1,
                         labels = c("B","C")),
                         labels = c("A",""),
                         nrow = 2,
                         rel_heights = c(0.3,1))

ggsave("pairs_plot_plus_gantt_main.pdf",
       plot = p2,
       device = "pdf",
       width = 45,
       height = 25,
       units = "cm",
       bg="white")

# lower_plots <- function(data,
#                         mapping,
#                         ...) {
#
#   ggplot(data = data,
#          mapping = mapping,
#          ...) +
#     geom_pointrange(alpha = 0.5, size = 1.25, ...) +
#     scale_y_continuous(labels=scales::percent, limits = c(0, 1)) +
#     scale_x_continuous(labels=scales::percent, limits = c(0, 1))
# }

# diag_plots <- function(data,
#                        mapping,
#                        ...) {
#
#   ggplot(data = data,
#          mapping = mapping,
#          ...) +
#     geom_histogram(binwidth = 0.05, boundary = 0, closed = "left", ..., alpha = 0.85) +
#     scale_x_continuous(labels=scales::percent, limits = c(0, 1)) +
#     theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
# }

# pairs_plot_month <- function(months) {
#   ggpairs(filter(pairs_data, time_6m == months),
#           columns = c(which(colnames(pairs_data) %in% c("BL_prev", "net_use", "trial_net_use"))),
#           upper = NULL,
#           diag = list(continuous = diag_plots),
#           lower = list(continuous = lower_plots),
#           switch = "both",
#           aes(color = factor(i), fill = factor(i), shape = factor(i)),
#           columnLabels = c("Baseline\nprevalence",
#                             "Baseline\nITN use",
#                             "Trial\nITN use")
#           ) +
#     pairs_theme +
#     scale_color_manual(values = c("Protopopoff" = cols[1], "Staedke" = cols[2], "Mosha" = cols[3], "Accrombessi" = cols[4])) +
#     scale_fill_manual(values = c("Protopopoff" = cols[1], "Staedke" = cols[2], "Mosha" = cols[3], "Accrombessi" = cols[4])) +
#     scale_shape_manual(values = c("Protopopoff" = 15, "Staedke" = 16, "Mosha" = 17, "Accrombessi" = 18))
# }

# p1 <- cowplot::plot_grid(ggmatrix_gtable(pairs_plot_month("≤6 m") + labs(title = "0-6 months")),
#                    ggmatrix_gtable(pairs_plot_month("7-12 m") + labs(title = "7-12 months")),
#                    ggmatrix_gtable(pairs_plot_month("13-18 m") + labs(title = "13-18 months")),
#                    ggmatrix_gtable(pairs_plot_month("19-24 m") + labs(title = "19-24 months")),
#                    ggmatrix_gtable(pairs_plot_month("25-30 m") + labs(title = "25-30 months")),
#                    ggmatrix_gtable(pairs_plot_month("31-36 m") + labs(title = "31-36 months")),
#                    nrow = 3, labels = c("A","B","C","D","E","F")
#                    )
