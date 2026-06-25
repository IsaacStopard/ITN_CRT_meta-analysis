library(cmdstanr)
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

orderly2::orderly_dependency("1_data_cleaning", "latest()", c("stan_data.rds"))

orderly2::orderly_shared_resource("m1.stan")

orderly2::orderly_dependency("1_data_cleaning", "latest()",
                             c("stan_data.rds"))

list2env(readRDS("stan_data.rds"), envir = .GlobalEnv)

COMBO_stan_subset <- subset(COMBO_stan, study != "Protopopoff et al") |>
  mutate(
  # study
  i = case_when(
    str_detect(author, "Staedke") ~ 1,
    str_detect(author, "Mosha") ~ 2,
    str_detect(author, "Accrombessi") ~ 3
  ),
  # net study
  li = case_when(
    l == 1 ~ 1, # if pyrethroid only nets then it is coded as 1, these should not be included in pmat, so just the intercept is used for these clusters
    i == 1 & l == 2 ~ 2,
    i == 2 & l == 2 ~ 3,
    i == 2 & l == 3 ~ 4,
    i == 3 & l == 3 ~ 5,
    .default = NA
  ))

data_in_subset <- extract_data_in(COMBO_stan_subset, train_folds = 1:K_folds)

saveRDS(COMBO_stan_subset, file = "COMBO_stan_PBO_subset.rds")
saveRDS(data_in_subset, file = "data_in_PBO_subset.rds")

adapt <- 0.9999
m_td <- 12
step_in <- 0.25

stan_model_in <- cmdstan_model(stan_file = "m1.stan", cpp_options = list(stan_threads = TRUE))

fit_full <- stan_model_in$sample(data = data_in_subset,
                                 iter_sampling = iter - warmup,
                                 iter_warmup = warmup,
                                 max_treedepth = m_td, 
                                 adapt_delta = adapt, 
                                 step_size = step_in,
                                 parallel_chains = 4,
                                 threads_per_chain = 1,
                                 seed = 12345)

fit_full$save_object(file = "m1_re_PBO_subset_fit_full.rds")

