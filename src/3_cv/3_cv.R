library(cmdstanr)
library(tidyverse)
library(loo)
library(posterior)

orderly2::orderly_dependency("1_data_cleaning", "latest()",
                             c("stan_data.rds"))

orderly2::orderly_shared_resource("m0.stan",
                                  "m1.stan",
                                  "m2.stan")

list2env(readRDS("stan_data.rds"), envir = .GlobalEnv)

pars <- orderly2::orderly_parameters(model = NULL)

model <- pars$model

###################################
##### k-fold cross validation #####
###################################

# https://discourse.mc-stan.org/t/k-fold-cv-with-observation-level-random-effects/35244

run_cv <- function(COMBO_stan, model, stan_model_in){
  log_pd_kfold <- matrix(nrow = (iter - warmup) * 4, ncol = nrow(COMBO_stan))

  K_folds <- length(unique(COMBO_stan[,"fold"]))

  for(i in 1:K_folds){

    train_inds <- which(COMBO_stan$fold != i)
    test_inds <- which(COMBO_stan$fold == i)
    
    train_folds <- unique(COMBO_stan$fold)[!(unique(COMBO_stan$fold) %in% i)]

    data_in_train <- extract_data_in(COMBO_stan, train_folds = train_folds)

    adapt <- if(model %in% c("m0", "m1")){0.9999} else{0.99}
    m_td <- if(model %in% c("m0", "m1")){15} else {12}
    step_in <- if(model %in% c("m0", "m1")){0.25} else{1}

    seed_in <- if(model == "m0" & i %in% c(3, 4) | model == "m1" & i %in% c(1)){
      1234
    } else if(model == "m1" & i %in% c(2)){
        123456
    } else{
        123
      }

    fit <- stan_model_in$sample(data = data_in_train,
                                seed = seed_in,
                                iter_sampling = (iter - warmup),
                                iter_warmup = warmup,
                                parallel_chains = 4,
                                threads_per_chain = 1,
                                init = 2,
                                max_treedepth = m_td, 
                                adapt_delta = adapt, 
                                step_size = step_in
                                )
    
    diagnostics <- fit$diagnostic_summary()
    
    if(any(diagnostics$num_divergent > 0) || any(diagnostics$ebfmi < 0.3)){
      fit <- NA
    }

    log_pd_kfold[, test_inds] <- fit$draws(variables = "log_lik") |> as_draws_matrix()

    rm(list = c("data_in_train", "fit"))
  }

  return(log_pd_kfold)
}

stan_model_in <- cmdstan_model(stan_file = paste0(model,".stan"), cpp_options = list(stan_threads = TRUE))

log_pd_kfold <- run_cv(COMBO_stan = COMBO_stan, model = model, stan_model_in = stan_model_in)

saveRDS(log_pd_kfold, paste0(model, "_log_pd_kfold.rds"))

