library(tidyverse)
library(metafor)
library(lme4)

############################################################################
###### frequentist individual participant data two-stage meta-analsyis #####
############################################################################

##### data

orderly2::orderly_dependency("1_data_cleaning", "latest()",
                             c("stan_data.rds"))

list2env(readRDS("stan_data.rds"), envir = .GlobalEnv)

u_ij <- unique(COMBO_stan[,c("i", "cluster")]) |> dplyr::mutate(ij = dplyr::row_number())

COMBO_stan <- COMBO_stan |> dplyr::left_join(u_ij) |>
  mutate(study = case_when(study == "Mosha et al" ~ "Tanzania (2019)",
                           study == "Accrombessi et al" ~ "Benin (2020)",
                           study == "Staedke et al" ~ "Uganda (2017)",
                           study == "Protopopoff et al" ~ "Tanzania (2015)"))

COMBO_stan$study <- factor(COMBO_stan$study, levels = c("Tanzania (2015)", "Uganda (2017)", "Tanzania (2019)", "Benin (2020)"))

data <- COMBO_stan |> group_by(ij, study, net, l, i, cluster) |> summarise(pos = sum(prev_num), test = sum(prev_denom)) |>
  mutate(cluster = as.character(cluster))

P_data <- subset(data, study == "Tanzania (2015)")
S_data <- subset(data, study == "Uganda (2017)")
M_data <- subset(data, study == "Tanzania (2019)")
A_data <- subset(data, study == "Benin (2020)")

# tests
nrow(P_data) + nrow(S_data) + nrow(M_data) + nrow(A_data) == nrow(data)
sum(data$test) == sum(COMBO_stan$prev_denom)

### model fits
# Protopopoff et al
P_fit <- glmer(cbind(pos, test - pos) ~ factor(l) + (1|cluster), data = P_data, family = "binomial")
S_fit <- glmer(cbind(pos, test - pos) ~ factor(l) + (1|cluster), data = S_data, family = "binomial")
M_fit <- glmer(cbind(pos, test - pos) ~ factor(l) + (1|cluster), data = M_data, family = "binomial")
A_fit <- glmer(cbind(pos, test - pos) ~ factor(l) + (1|cluster), data = A_data, family = "binomial")

effect_data <- data.frame(study = c("Tanzania (2015)", "Uganda (2017)", "Tanzania (2019)", "Tanzania (2019)", "Benin (2020)"),
                          net = c("Pyrethroid-PBO", "Pyrethroid-PBO", "Pyrethroid-PBO", "Pyrethroid-pyrrole", "Pyrethroid-pyrrole"),
                          yi = c(summary(P_fit)$coef[2,"Estimate"], summary(S_fit)$coef[2,"Estimate"],
                                 summary(M_fit)$coef[2,"Estimate"], summary(M_fit)$coef[3,"Estimate"], summary(A_fit)$coef[2,"Estimate"]),
                          se = c(summary(P_fit)$coef[2,"Std. Error"], summary(S_fit)$coef[2,"Std. Error"],
                                 summary(M_fit)$coef[2,"Std. Error"], summary(M_fit)$coef[3,"Std. Error"],  summary(A_fit)$coef[2,"Std. Error"]))

data_plot <- data |> group_by(study, net) |> summarise(pos = sum(pos), neg = sum(test) - sum(pos))

data_plot_pyr <- subset(data_plot, net == "Pyrethroid-only")
data_plot_pp <- subset(data_plot, net == "Pyrethroid-pyrrole")
data_plot_pbo <- subset(data_plot, net == "Pyrethroid-PBO")

### meta-analysis random effects model
rma_fit <- rma.uni(yi = subset(effect_data, net == "Pyrethroid-PBO")$yi,
               sei = subset(effect_data, net == "Pyrethroid-PBO")$se,
               method = "REML",
               slab = subset(effect_data, net == "Pyrethroid-PBO")$study)

rma_fit_pp <- rma(yi = subset(effect_data, net == "Pyrethroid-pyrrole")$yi,
                  sei = subset(effect_data, net == "Pyrethroid-pyrrole")$se,
                  method = "REML",
                  slab = subset(effect_data, net == "Pyrethroid-pyrrole")$study)


pdf(file = "two_stage_forest_plot_pbo.pdf",
  width = 15, height = 10)
forest.rma(rma_fit,
           ilab = cbind(data_plot_pbo$pos, data_plot_pbo$neg, subset(data_plot_pyr, study != "Benin (2020)")$pos, subset(data_plot_pyr, study != "Benin (2020)")$neg),
           ilab.lab = c("+", "-", "+", "-"))
text(c(-2.75, -1.95), rma_fit$k+2.5, c("Pyrethroid-PBO", "Pyrethroid-only"), cex=0.75, font=2)
text(-0.75, 5.5, pos=4, cex=0.75, bquote(
  paste("RE Model for All Studies (Q = ",
        .(formatC(rma_fit$QE, digits=2, format="f")),
        ", df = ", .(rma_fit$k - rma_fit$p),", p = ",
        .(formatC(rma_fit$QEp, digits=2, format="f")),
        "; ", I^2, " = ",
        .(formatC(rma_fit$I2, digits=1, format="f")),
        "%)", "; ", tau^2 ==
          .(formatC(rma_fit$tau2, digits=2, format="f")))))
dev.off()

pdf(file = "two_stage_forest_plot_pp.pdf",
    width = 15, height = 10)
forest.rma(rma_fit_pp,
           ilab = cbind(data_plot_pp$pos, data_plot_pp$neg, subset(data_plot_pyr, study %in% c("Tanzania (2019)", "Benin (2020)"))$pos,
                        subset(data_plot_pyr, study %in% c("Tanzania (2019)", "Benin (2020)"))$neg),
           ilab.lab = c("+", "-", "+", "-"))
text(c(-1.75, -1.3), rma_fit_pp$k+2.5, c("Pyrethroid-pyrrole", "Pyrethroid-only"), cex=0.75, font=2)
text(-0.75, 4.5, pos=4, cex=0.75, bquote(
  paste("RE Model for All Studies (Q = ",
        .(formatC(rma_fit_pp$QE, digits=2, format="f")),
        ", df = ", .(rma_fit_pp$k - rma_fit_pp$p),", p = ",
        .(formatC(rma_fit_pp$QEp, digits=2, format="f")),
        "; ", I^2, " = ",
        .(formatC(rma_fit_pp$I2, digits=1, format="f")),
        "%)", "; ", tau^2 ==
          .(formatC(rma_fit_pp$tau2, digits=2, format="f")))))
dev.off()
