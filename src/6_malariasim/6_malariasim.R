library(malariasimulation); library(tidyverse); library(cali); library(parallel);

##### default values #####
year <- 365
sim_length <- 15 * year
res <- 0.5

simparams <- malariasimulation::get_parameters(
  list(human_population = 17500,
       prevalence_rendering_min_ages = c(0) * 365,
       prevalence_rendering_max_ages = c(100) * 365,
       incidence_rendering_min_ages = c(0) * 365,
       incidence_rendering_max_ages = c(100) * 365,
       clinical_incidence_rendering_min_ages = c(0) * 365,
       clinical_incidence_rendering_max_ages = c(100) * 365,
       severe_incidence_rendering_min_ages = c(0) * 365,
       severe_incidence_rendering_max_ages = c(100) * 365,
       model_seasonality = FALSE,
       individual_mosquitoes = FALSE,
       bednets = TRUE,
       smc = FALSE)
)

dist_step <- 3
trial_time <- 9

bednet_timesteps <- c(seq(0, 6, dist_step), trial_time) * year

retention <- 5 * year
rate <- 1 - exp(-1/retention)
base_net <- seq(0, 0.75, 0.25) * rate * (3 * 365 - 0 * 365) / (exp(-0 * 365 * rate) - exp(-3 * 365 * rate))

# bed net parameters
gen_median <- function(.df_w){
  return(.df_w |> dplyr::group_by(resistance, bioassay_mortality) |>
           dplyr::summarise(dn0_med = median(dn0),
                            rn0_med = median(rn0),
                            gamman_med = median(gamman)
           ) |>
           ungroup() |>
           as.data.frame()
  )
}

# pyrethroid-only
dat_res_pyr <- readRDS("pyrethroid_binomial_uncertainty 1") |>
  dplyr::rename(rn0 = rn_pyr, gamman = mean_duration) |>
  dplyr::mutate(bioassay_mortality = 1 - resistance) |>
  select(dn0, rn0, gamman, resistance, bioassay_mortality) |>
  gen_median()

# pyrethroid-pyrrole
dat_res_pp <- readRDS("pyrrole_binomial_uncertainty_using_pyrethroid_dn0_for_mn_durability 1") |>
  dplyr::rename(rn0 = rn_pbo, gamman = mn_dur, resistance = resistance.x) |>
  dplyr::mutate(bioassay_mortality = 1 - resistance) |>
  select(dn0, rn0, gamman, resistance, bioassay_mortality) |> gen_median()

dn0_pyr <- dat_res_pyr[dat_res_pyr$resistance == res, "dn0_med"][[1]]
rn0_pyr <- dat_res_pyr[dat_res_pyr$resistance == res, "rn0_med"][[1]]
gamman_pyr <- dat_res_pyr[dat_res_pyr$resistance == res, "gamman_med"][[1]]

dn0_pp <- dat_res_pp[dat_res_pp$resistance == res, "dn0_med"][[1]]
rn0_pp <- dat_res_pp[dat_res_pp$resistance == res, "rn0_med"][[1]]
gamman_pp <- dat_res_pp[dat_res_pp$resistance == res, "gamman_med"][[1]]

vals_pyr <- expand.grid(base_prev = seq(0.1, 0.55, 0.05),
                    base_net = base_net,
                    trial_net = base_net[-which(base_net == 0)]) |>
  mutate(net_type = "pyr")

vals_pp <- vals_pyr |> mutate(net_type = "pp")

#####################
##### functions #####
#####################

get_params <- function(i,
                       vals,
                       simparams,
                       bednet_timesteps,
                       retention,
                       dn0_b = dn0_pyr,
                       dn0_t = dn0_pp,
                       rn0_b = rn0_pyr,
                       rn0_t = rn0_pp,
                       gamman_b = gamman_pyr,
                       gamman_t = gamman_pp
                       ){

  base_net <- vals[i, "base_net"]
  trial_net <- vals[i, "trial_net"]

  n_b_t <- length(bednet_timesteps) - 1

  rn0_t_in <- if(vals[i, "net_type"] == "pyr"){rn0_b}else{rn0_t}
  dn0_t_in <- if(vals[i, "net_type"] == "pyr"){dn0_b}else{dn0_t}
  gamman_t_in <- if(vals[i, "net_type"] == "pyr"){gamman_b}else{gamman_t}

  bednetparams <- malariasimulation::set_bednets(
    simparams,
    timesteps = bednet_timesteps,
    coverages = c(rep(base_net, n_b_t), trial_net),  # Each round is distributed to 50% of the population.
    retention = retention, # Nets are kept on average 5 years
    dn0 = matrix(c(rep(dn0_b, n_b_t), dn0_t_in), nrow = n_b_t + 1, ncol = 1), # Matrix of death probabilities for each mosquito species over time
    rn = matrix(c(rep(rn0_b, n_b_t), rn0_t_in), nrow = n_b_t + 1, ncol = 1), # Matrix of repelling probabilities for each mosquito species over time
    rnm = matrix(c(rep(.24, n_b_t), .24), nrow = n_b_t + 1, ncol = 1), # Matrix of minimum repelling probabilities for each mosquito species over time
    gamman = (c(rep(gamman_b/log(2), n_b_t), gamman_t_in/log(2)) * 365) # Vector of bed net half-lives for each distribution timestep
    )

  return(bednetparams)
}

summary_pfpr <- function(.output){
  return(
    mean(.output[(9 * 365 - 1),'n_detect_lm_0_36500'] / .output[(9 * 365 - 1),'n_age_0_36500']) # baseline prevalence calibrated the day before the trial starts
         )
}

calibrate_model <- function(i, vals, params, trial_time, year, summary_pfpr){
  params[[i]]$timesteps <- trial_time * year

  set.seed(12345)
  out <- cali::calibrate(parameters = params[[i]],
                   target = vals[i, "base_prev"],
                   eq_prevalence = vals[i, "base_prev"],
                   summary_function = summary_pfpr,
                   eir_limits = c(0, 1500),
                   human_population = 10000,
                   max_attempts = 15,
                   correlations_int = 'bednets',
                   correlations_round_rho =  1
  )
  return(out)
}

run_model_pre_trial <- function(i,
                                EIR_values,
                                year,
                                trial_time,
                                params){


  EIR <- EIR_values[[i]]

  params[[i]] <- malariasimulation::set_equilibrium(parameters = params[[i]], init_EIR = EIR)
  correlations_in <- get_correlation_parameters(params[[i]])
  correlations_in$inter_round_rho('bednets', 1) # if one nets are given to people who already have nets

  set.seed(12345)
  out <- malariasimulation::run_resumable_simulation(timesteps = trial_time * year - 1,
                                                     parameters = params[[i]],
                                                     correlations = correlations_in,
                                                     initial_state = NULL,
                                                     restore_random_state = FALSE)
  return(out)
}

run_model_post_trial <- function(i,
                                 EIR_values,
                                 sim_length,
                                 params,
                                 init_states){

  EIR <- EIR_values[[i]]

  params[[i]] <- malariasimulation::set_equilibrium(parameters = params[[i]], init_EIR = EIR)
  correlations_in <- get_correlation_parameters(params[[i]])
  correlations_in$inter_round_rho('bednets', 1) # if one nets are given to people who already have nets

  set.seed(12345)

  out <- malariasimulation::run_resumable_simulation(timesteps = sim_length,
                                                     parameters = params[[i]],
                                                     correlations = correlations_in,
                                                     initial_state = init_states[[i]]$state,
                                                     restore_random_state = TRUE)

  return(rbind(init_states[[i]]$data, out$data))
}

###############
##### run #####
###############

params_pyr <- lapply(seq(1, nrow(vals_pyr)),
                 get_params,
                 vals = vals_pyr, simparams = simparams,
                 bednet_timesteps = bednet_timesteps, retention = retention)

params_pp <- lapply(seq(1, nrow(vals_pp)),
                     get_params,
                     vals = vals_pp, simparams = simparams,
                     bednet_timesteps = bednet_timesteps, retention = retention)

cores <- 5
cl <- parallel::makeCluster(cores)
parallel::clusterEvalQ(cl,{
  library(tidyr)
  library(dplyr)
  library(malariasimulation)
  library(cali)
}
)

parallel::clusterExport(cl,
                        varlist = c("vals_pyr",
                                    "params_pyr",
                                    "trial_time",
                                    "year",
                                    "calibrate_model",
                                    "summary_pfpr")
                        )

EIR_values <- parallel::parLapply(cl,
                                  seq(1, nrow(vals_pyr)),
                                  calibrate_model,
                                  vals = vals_pyr,
                                  params = params_pyr,
                                  trial_time = trial_time, year = year,
                                  summary_pfpr = summary_pfpr
                                  )

parallel::stopCluster(cl)


#
cl <- parallel::makeCluster(cores)
parallel::clusterEvalQ(cl,{
  library(tidyr)
  library(dplyr)
  library(malariasimulation)
}
)

parallel::clusterExport(cl,
                        varlist = c("vals_pyr",
                                    "params_pyr",
                                    "trial_time",
                                    "year",
                                    "run_model_pre_trial",
                                    "EIR_values")
)

init_states <- parLapply(cl,
                        seq(1, nrow(vals_pyr)),
                        run_model_pre_trial,
                        EIR_values = EIR_values,
                        trial_time = trial_time,
                        year = year,
                        params = params_pyr)

parallel::stopCluster(cl)

cl <- parallel::makeCluster(cores)
parallel::clusterEvalQ(cl,{
  library(tidyr)
  library(dplyr)
  library(malariasimulation)
}
)

parallel::clusterExport(cl,
                        varlist = c("vals_pyr",
                                    "params_pyr",
                                    "vals_pp",
                                    "params_pp",
                                    "sim_length",
                                    "run_model_post_trial",
                                    "EIR_values",
                                    "init_states")
)

model_runs_pyr <- parLapply(cl,
                         seq(1, nrow(vals_pyr)),
                         run_model_post_trial,
                         EIR_values = EIR_values,
                         sim_length = sim_length,
                         params = params_pyr,
                         init_states = init_states)

model_runs_pp <- parLapply(cl,
                            seq(1, nrow(vals_pp)),
                            run_model_post_trial,
                            EIR_values = EIR_values,
                            sim_length = sim_length,
                            params = params_pp,
                            init_states = init_states)

parallel::stopCluster(cl)

saveRDS(list("model_runs" = append(model_runs_pyr, model_runs_pp),
             "params" = append(params_pyr, params_pp),
             "EIR_values" = append(EIR_values, EIR_values),
             "vals" = rbind(vals_pyr, vals_pp)
             ),
        file = "sim_results.rds")




