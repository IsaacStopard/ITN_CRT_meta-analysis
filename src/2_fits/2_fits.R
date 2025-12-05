library(rstan)
library(tidyverse)
library(fastDummies)

orderly2::orderly_dependency("1_data_cleaning", "latest()",
                             c("stan_data.rds"))

orderly2::orderly_shared_resource("m0_c.stan",
                                  "m0_c_FE.stan",
                                  "m1_re.stan",
                                  "m2_re.stan",
                                  "m3_re.stan",
                                  "m3_FE.stan")

orderly2::orderly_parameters(model = NULL)

list2env(readRDS("stan_data.rds"), envir = .GlobalEnv)

adapt <- if(model == "m3_FE"){0.8}else if(model == "m0_c" | model == "m1_re"){0.999}else{0.99}
m_td <- if(model == "m0_c" | model == "m0_c_FE" | model == "m3_FE"){11.5}else{14}
step_in <- if(model == "m0_c" | model == "m1_re"){0.25}else{1}

stan_model_in <- rstan::stan_model(file = paste0(model, ".stan"))

fit_full <- rstan::sampling(stan_model_in,
                            data = data_in_full,
                            iter = iter,
                            warmup = warmup,
                            control = list(max_treedepth = m_td, adapt_delta = adapt, stepsize = step_in),
                            cores = 4,
                            chains = 4,
                            seed = 123)

saveRDS(fit_full, file = paste0(model, "_fit_full.rds"))
