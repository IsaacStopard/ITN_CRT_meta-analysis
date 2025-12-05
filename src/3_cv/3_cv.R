library(rstan)
library(tidyverse)
library(fastDummies)
library(loo)

orderly2::orderly_dependency("1_data_cleaning", "latest()",
                             c("stan_data.rds"))

orderly2::orderly_parameters(model = NULL)

orderly2::orderly_shared_resource("m0_c.stan",
                                  "m0_c_FE.stan",
                                  "m1_re.stan",
                                  "m2_re.stan",
                                  "m3_re.stan",
                                  "m3_FE.stan")

list2env(readRDS("stan_data.rds"), envir = .GlobalEnv)

###################################
##### k-fold cross validation #####
###################################

# https://discourse.mc-stan.org/t/k-fold-cv-with-observation-level-random-effects/35244

run_cv <- function(COMBO_stan, model_i, stan_model_in){
  log_pd_kfold <- matrix(nrow = (iter - warmup) * 4, ncol = nrow(COMBO_stan))

  K_folds <- length(unique(COMBO_stan[,"fold"]))

  for(i in 1:K_folds){

    train_inds <- which(COMBO_stan$fold != i)
    test_inds <- which(COMBO_stan$fold == i)

    data_in_train <- extract_data_in(COMBO_stan, train_inds = train_inds)

    #adapt <- if(model == "m4" | model == "m3"){0.9}else{0.8}
    #m_td <- if(model == "m0"){14}else{12}
    adapt <- if(model == "m3_FE"){0.8} else if(model == "m0_c"){0.9999} else if(model == "m1_re"){0.99999} else{0.99}
    m_td <- if(model == "m0_c" | model == "m0_c_FE" | model == "m3_FE" | model == "m1_re"){11.25}else{14}
    step_in <- if(model == "m0_c"){0.5}else if( model == "m1_re"){0.25}else{1}

    set.seed(123)

    fit <- tryCatch(
      {rstan::sampling(stan_model_in,
                           data = data_in_train,
                           iter = iter,
                           chains = 4,
                           warmup = warmup,
                           control = list(max_treedepth = m_td, adapt_delta = adapt, stepsize = step_in),
                           cores = 4)
        },
      warning = function(w){
        return(NA)
      }
                    )

    log_pd_kfold[, test_inds] <- rstan::extract(fit, "log_lik")$log_lik[,test_inds]

    rm(list = c("data_in_train", "fit"))
  }

  return(log_pd_kfold)
}

stan_model_in <- rstan::stan_model(file = paste0(model,".stan"))

log_pd_kfold <- run_cv(COMBO_stan = COMBO_stan, model_i = model, stan_model_in = stan_model_in)

saveRDS(log_pd_kfold, paste0(model, "_log_pd_kfold.rds"))
