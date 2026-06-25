library(cmdstanr)
library(tidyverse)
library(loo)

orderly2::orderly_dependency("1_data_cleaning", "latest()",
                             c("stan_data.rds"))

pars <- orderly2::orderly_parameters(model = NULL)
model <- pars$model

orderly2::orderly_shared_resource("m0.stan")

list2env(readRDS("stan_data.rds"), envir = .GlobalEnv)

stan_model_in <- cmdstan_model(stan_file = paste0(model, ".stan"), cpp_options = list(stan_threads = TRUE))

################
##### data #####
################

COMBO_stan_1 <- subset(COMBO_stan, time <= 12)
COMBO_stan_2 <- subset(COMBO_stan, time <= 24 & time > 12)
COMBO_stan_3 <- subset(COMBO_stan, time <= 36 & time > 24)
COMBO_stan_4 <- subset(COMBO_stan, time >= 0 & time <= 24)

nrow(COMBO_stan) == (nrow(COMBO_stan_1) + nrow(COMBO_stan_2) + nrow(COMBO_stan_3))

data_in_1 <- extract_data_in(COMBO_stan_1, train_folds = 1:K_folds)
data_in_2 <- extract_data_in(COMBO_stan_2, train_folds = 1:K_folds)
data_in_3 <- extract_data_in(COMBO_stan_3, train_folds = 1:K_folds)
data_in_4 <- extract_data_in(COMBO_stan_4, train_folds = 1:K_folds)

fit_model_m0 <- function(data_in, stan_model_in){
  
  adapt <- 0.999
  m_td <- 12
  step_in <- 0.75
  
  fit_model <- function(seed){
    stan_model_in$sample(data = data_in,
                              iter_sampling = iter - warmup,
                              iter_warmup = warmup,
                              max_treedepth = m_td, 
                              adapt_delta = adapt, 
                              step_size = step_in,
                              parallel_chains = 4,
                              threads_per_chain = 1,
                              seed = seed)
  }
  
  fit <- fit_model(seed = 123)
  
  diagnostics <- fit$diagnostic_summary()
  
  if(any(diagnostics$num_divergent > 0) || any(diagnostics$ebfmi < 0.3)){
    
    fit <- fit_model(seed = 1234)
    
    diagnostics <- fit$diagnostic_summary()
    
    if(any(diagnostics$num_divergent > 0) || any(diagnostics$ebfmi < 0.3)){
      fit <- fit_model(seed = 12345)
      
      diagnostics <- fit$diagnostic_summary()
      
      if(any(diagnostics$num_divergent > 0) || any(diagnostics$ebfmi < 0.3)){
        
        fit <- NA
        
      }
    }
  }
  
  return(fit)
}

##############################
##### first time segment #####
##############################

m0_c_fit_1 <- fit_model_m0(data_in_1, stan_model_in)

if(!is.environment(m0_c_fit_1)){
  stop("model not fit without warnings")
}

m0_c_fit_1$save_object(file = "m0_c_fit_1.rds")

###############################
##### second time segment #####
###############################

m0_c_fit_2 <- fit_model_m0(data_in_2, stan_model_in)

if(!is.environment(m0_c_fit_2)){
  stop("model not fit without warnings")
}

m0_c_fit_2$save_object(file = "m0_c_fit_2.rds")

##############################
##### third time segment #####
##############################

m0_c_fit_3 <- fit_model_m0(data_in_3, stan_model_in)

if(!is.environment(m0_c_fit_3)){
  stop("model not fit without warnings")
}

m0_c_fit_3$save_object(file = "m0_c_fit_3.rds")

##############################
##### fourth time segment #####
##############################

m0_c_fit_4 <- fit_model_m0(data_in_4, stan_model_in)

if(!is.environment(m0_c_fit_4)){
  stop("model not fit without warnings")
}

m0_c_fit_4$save_object(file = "m0_c_fit_4.rds")

saveRDS(list("COMBO_stan_1" = COMBO_stan_1,
             "COMBO_stan_2" = COMBO_stan_2,
             "COMBO_stan_3" = COMBO_stan_3,
             "COMBO_stan_4" = COMBO_stan_3,
             "data_in_1" = data_in_1,
             "data_in_2" = data_in_2,
             "data_in_3" = data_in_3,
             "data_in_4" = data_in_4
             ),
        file = paste0(model, "_time_fits_data.rds"))
