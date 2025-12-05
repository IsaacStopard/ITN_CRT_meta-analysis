library(rstan)
library(tidyverse)
library(fastDummies)
library(loo)

orderly2::orderly_dependency("1_data_cleaning", "latest()",
                             c("stan_data.rds"))

orderly2::orderly_parameters(model = NULL)

orderly2::orderly_shared_resource("m0_c_FE.stan", "m0_c.stan")

list2env(readRDS("stan_data.rds"), envir = .GlobalEnv)

adapt <- 0.999
m_td <- 11.25
step_in <- 0.75

stan_model_in <- rstan::stan_model(file = paste0(model, ".stan"))

################
##### data #####
################

COMBO_stan_1 <- subset(COMBO_stan, time <= 12)
COMBO_stan_2 <- subset(COMBO_stan, time <= 24 & time > 12)
COMBO_stan_3 <- subset(COMBO_stan, time <= 36 & time > 24)
COMBO_stan_4 <- subset(COMBO_stan, time >= 0 & time <= 24)

nrow(COMBO_stan) == (nrow(COMBO_stan_1) + nrow(COMBO_stan_2) + nrow(COMBO_stan_3))

data_in_1 <- extract_data_in(COMBO_stan = COMBO_stan_1, train_inds = 1:nrow(COMBO_stan_1))
data_in_2 <- extract_data_in(COMBO_stan = COMBO_stan_2, train_inds = 1:nrow(COMBO_stan_2))
data_in_3 <- extract_data_in(COMBO_stan = COMBO_stan_3, train_inds = 1:nrow(COMBO_stan_3))
data_in_4 <- extract_data_in(COMBO_stan = COMBO_stan_4, train_inds = 1:nrow(COMBO_stan_4))

##############################
##### first time segment #####
##############################

m0_c_fit_1 <- tryCatch(
  {rstan::sampling(stan_model_in,
                            data = data_in_1,
                            iter = iter,
                            warmup = warmup,
                            control = list(max_treedepth = m_td, adapt_delta = adapt, stepsize = step_in),
                            cores = 4,
                            chains = 4,
                            seed = 123)},
      warning = function(w){
        return(NA)
      })

###############################
##### second time segment #####
###############################

m0_c_fit_2 <- tryCatch(
  {rstan::sampling(stan_model_in,
                              data = data_in_2,
                              iter = iter,
                              warmup = warmup,
                              control = list(max_treedepth = m_td, adapt_delta = adapt, stepsize = step_in),
                              cores = 4,
                              chains = 4,
                              seed = 123)},
  warning = function(w){
    return(NA)
  })

##############################
##### third time segment #####
##############################

m0_c_fit_3 <- tryCatch(
  {rstan::sampling(stan_model_in,
                              data = data_in_3,
                              iter = iter,
                              warmup = warmup,
                              control = list(max_treedepth = m_td, adapt_delta = adapt, stepsize = step_in),
                              cores = 4,
                              chains = 4,
                              seed = 123)},
  warning = function(w){
    return(NA)
  })

##############################
##### fourth time segment #####
##############################

m0_c_fit_4 <- tryCatch(
  {rstan::sampling(stan_model_in,
                   data = data_in_4,
                   iter = iter,
                   warmup = warmup,
                   control = list(max_treedepth = m_td, adapt_delta = adapt, stepsize = step_in),
                   cores = 4,
                   chains = 4,
                   seed = 123)},
  warning = function(w){
    return(NA)
  })


saveRDS(list("COMBO_stan_1" = COMBO_stan_1,
             "COMBO_stan_2" = COMBO_stan_2,
             "COMBO_stan_3" = COMBO_stan_3,
             "COMBO_stan_4" = COMBO_stan_3,
             "data_in_1" = data_in_1,
             "data_in_2" = data_in_2,
             "data_in_3" = data_in_3,
             "data_in_4" = data_in_4,
             "m0_c_fit_1" = m0_c_fit_1,
             "m0_c_fit_2" = m0_c_fit_2,
             "m0_c_fit_3" = m0_c_fit_3,
             "m0_c_fit_4" = m0_c_fit_4
             ),
        file = paste0(model, "_time_fits.rds"))
