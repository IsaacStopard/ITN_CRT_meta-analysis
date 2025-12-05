# calibration functions with correlation parameters
# modified from cali package https://github.com/mrc-ide/cali/tree/main

print_update <- function(eir, objective){
  output_print <- data.frame(
    EIR = signif(eir, 2),
    Objective = signif(objective, 2)
  )
  print(knitr::kable(output_print, format = "simple"))
  cat("\n")
}

check_fit <- function(proposed_eir, parameters, target, summary_function, correlations){
  parameters <- parameters |>
    malariasimulation::set_equilibrium(
      init_EIR = proposed_eir
    )
  simulation <- malariasimulation::run_simulation(
    timesteps = parameters$timesteps,
    parameters = parameters,
    correlations = correlations # added this
  )
  output <- summary_function(simulation)

  output_print <- data.frame(
    Current = signif(output, 2),
    Target = signif(target, 2)
  )
  print(knitr::kable(output_print, format = "simple"))
  cat("\n")
  elimination <- cali::check_elimination(output, target)
  objective <- ifelse(elimination, NA, sum(output - target))
  return(objective)
}

calibrate <- function(parameters, target, summary_function, eq_prevalence,
                      eq_ft = 0, human_population = c(1000, 10000, 100000),
                      eir_limits = c(0, 1500), max_attempts = 10, correlations){

  eir <- rep(0, 2)
  objective <- rep(NA, 2)
  parameters$human_population <- human_population[1]

  # Estimate starting guess EIR with equilibrium solution
  message("Initialising EIR\n")
  eir[1] <- cali::get_eq_eir(
    target_pfpr = eq_prevalence,
    ft = eq_ft
  )
  if(eir[1] < eir_limits[1]){
    eir[1] <- eir_limits[1] + 1
  }
  if(eir[1] > eir_limits[2]){
    eir[1] <- eir_limits[2] - 1
  }

  # Ensure the starting EIR is valid and obtain objective evaluation
  min_eir <- 0
  attempts <- 0
  message("Slice sampling EIR, side 1 \n")
  while(is.na(objective[1])){
    message("Attempt ", attempts + 1, " of ", max_attempts)

    objective[1] <- check_fit(
      proposed_eir = eir[1],
      parameters = parameters,
      target = target,
      summary_function = summary_function,
      correlations = correlations
    )

    cali::print_update(eir, objective)
    attempts <- attempts + 1

    # In the case that our starting EIR leads to unwanted elimination:
    ## First approach is to increase the human population if available.
    ## Second approach is to increase EIR and update our lower bound on EIR.
    if(is.na(objective[1])){
      if(parameters$human_population < max(human_population)){
        parameters$human_population <- human_population[which(human_population == parameters$human_population) + 1]
        message("Increasing human population due to elimination")
        message("Running with new population size: ", parameters$human_population, "\n")
      } else {
        parameters$human_population <- human_population[1]
        min_eir <- eir[1]
        eir[1] <- cali::proposal(
          current_eir = eir[1],
          limits = eir_limits,
          direction = "increase",
          step = 2
        )
      }
    }
    if(attempts > max_attempts){
      stop("Failure due to max attempts reached before successful run")
    }
  }

  # Estimate the second EIR
  message("Slice sampling EIR, side 2 \n")
  direction <- ifelse(objective[1] > 0, "decrease", "increase")
  if(direction == "decrease"){
    eir_limits[1] <- min_eir
  }
  estimated_eir <- NA
  eir[2] <- eir[1]
  repeat{
    message("Attempt ", attempts + 1, " of ", max_attempts)
    # Propose new EIR
    eir[2] <- cali::proposal(
      current_eir = eir[2],
      limits = eir_limits,
      direction = direction
    )
    # Evaluate objective
    objective[2] <- check_fit(
      proposed_eir = eir[2],
      parameters = parameters,
      target = target,
      summary_function = summary_function,
      correlations = correlations
    )
    cali::print_update(eir, objective)
    # Stop or update
    if(!is.na(objective[2])){
      if(cali::good_move(objective, direction)){
        estimated_eir <- cali::linear_interpolate(eir, objective)
        message("Success")
        break
      } else{
        # If we haven't moved far enough, we can update eir1 and fit1
        eir[1] <- eir[2]
        objective[1] <- objective[2]
      }
    } else {
      if(direction == "decrease"){
        eir_limits[1] <- eir[2]
        direction <- "increase"
      }
    }
    attempts <- attempts + 1
    if(attempts > max_attempts){
      estimated_eir <- mean(eir, na.rm = TRUE)
      message("Terminating as max attempts reached")
      break
    }
  }
  return(estimated_eir)
}
