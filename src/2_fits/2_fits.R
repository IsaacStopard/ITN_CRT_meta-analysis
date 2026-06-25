library(cmdstanr)
library(tidyverse)

orderly2::orderly_dependency("1_data_cleaning", "latest()",
                             c("stan_data.rds"))

orderly2::orderly_shared_resource("m0.stan",
                                  "m1.stan",
                                  "m2.stan")

list2env(readRDS("stan_data.rds"), envir = .GlobalEnv)

model_0 <- cmdstan_model(stan_file = "m0.stan", cpp_options = list(stan_threads = TRUE))
model_1 <- cmdstan_model(stan_file = "m1.stan", cpp_options = list(stan_threads = TRUE))
model_2 <- cmdstan_model(stan_file = "m2.stan", cpp_options = list(stan_threads = TRUE))

fit_0 <- model_0$sample(data = data_in_full,
                        seed = 123,
                        iter_sampling = iter - warmup,
                        iter_warmup = warmup,
                        chains = 4,
                        init = 2,
                        parallel_chains = 4,
                        threads_per_chain = 1,
                        adapt_delta = 0.9999,
                        max_treedepth = 15,
                        step_size = 0.25)

fit_0$save_object(file = "fit_0.rds")

fit_1 <- model_1$sample(data = data_in_full,
                        seed = 123,
                        iter_sampling = iter - warmup,
                        iter_warmup = warmup,
                        chains = 4,
                        init = 2,
                        parallel_chains = 4,
                        threads_per_chain = 1,
                        adapt_delta = 0.9999,
                        max_treedepth = 15,
                        step_size = 0.25)

fit_1$save_object(file = "fit_1.rds")

fit_2 <- model_2$sample(data = data_in_full,
                        seed = 123,
                        iter_sampling = iter - warmup,
                        iter_warmup = warmup,
                        chains = 4,
                        init = 2,
                        parallel_chains = 4,
                        threads_per_chain = 1,
                        adapt_delta = 0.99,
                        max_treedepth = 12)

fit_2$save_object(file = "fit_2.rds")

