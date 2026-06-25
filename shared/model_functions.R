####################################################
##### function to extract the parameter values #####
####################################################

extract_params <- function(target_params, fit){
  
  draws_mat <- fit$draws(variables = target_params) |> as_draws_matrix()
  
  params_list <- lapply(target_params, function(p){
    subset_draws(draws_mat, variable = p)
  }
  )
  
  names(params_list) <- target_params
  
  out <- list("draws_mat" = draws_mat, "params_list" = params_list)
  
  return(out)
}

#######################
##### odds ratios #####
#######################

l_or_study_fun <- function(data, params, model = 2, quantile_out = TRUE){
  
  list2env(params$params_list, envir = environment())
  
  o_r <- t(theta_l[, data$l]) + 
    t(theta_li[, data$li]) + 
    t(kappa_l[, data$l] + kappa_li[, data$li]) * data$years
  
  if(model == 2){
    o_r <- o_r + t(omega_l[, data$l]) * data$base_prev_mean + t(delta_l[, data$l]) * data$base_prev_diff
  }
  
  if(quantile_out){
    o_r <- t(apply(o_r, 1, quantile, probs = probs_in)) |> as.data.frame() |> rename("low" = 1, "med" = 2, "up" = 3)
    rownames(o_r) <- NULL
  }
  
  return(o_r)
}

l_or_pooled_fun <- function(data, params, model = 2, quantile_out = TRUE){
  
  list2env(params$params_list, envir = environment())
  
  o_r <- t(theta_l[, data$l]) + 
    t(kappa_l[, data$l]) * data$years 
  
  if(model == 2){
    o_r <- o_r + t(omega_l[, data$l]) * data$base_prev_mean + t(delta_l[, data$l]) * data$base_prev_diff
  }
  
  if(quantile_out){
    o_r <- t(apply(o_r, 1, quantile, probs = probs_in)) |> as.data.frame() |> rename("low" = 1, "med" = 2, "up" = 3)
    rownames(o_r) <- NULL
  }
  
  return(o_r)
  
}

mean_time_l_or_study_fun <- function(year_max, data, params, model = 2, quantile_out = TRUE){
  
  list2env(params$params_list, envir = environment())
  
  o_r <-(
    t(theta_l[, data$l] + theta_li[, data$li]) * year_max + 
      t(kappa_l[, data$l] + kappa_li[, data$li]) * year_max^2 / 2)
  
  if(model == 2){
    o_r <- o_r + t(omega_l[, data$l]) * data$base_prev_mean * year_max + t(delta_l[, data$l]) * data$base_prev_diff * year_max
  }
  
  o_r <- (1 / year_max) * o_r
  
  if(quantile_out){
    o_r <- t(apply(o_r, 1, quantile, probs = probs_in)) |> as.data.frame() |> rename("low" = 1, "med" = 2, "up" = 3)
    rownames(o_r) <- NULL
  }
  
  return(o_r) 
}

mean_time_l_or_pooled_fun <- function(year_max, data, params, model = 2, quantile_out = TRUE){
  
  list2env(params$params_list, envir = environment())
  
  o_r <- (
    t(theta_l[, data$l]) * year_max + 
      t(kappa_l[, data$l]) * year_max^2 / 2
  )
  
  if(model == 2){
    o_r <- o_r + 
      t(omega_l[, data$l]) * data$base_prev_mean * year_max + 
      t(delta_l[, data$l]) * data$base_prev_diff * year_max
  }
  
  o_r <- (1 / year_max) * o_r
  
  if(quantile_out){
    o_r <- t(apply(o_r, 1, quantile, probs = probs_in)) |> as.data.frame() |> rename("low" = 1, "med" = 2, "up" = 3)
    rownames(o_r) <- NULL
  }
  
  return(o_r) 
  
}

#######################
##### prevalence ######
#######################

mean_time_prev_study_fun <- function(year_max, data, params, model = 2, quantile_out){
  list2env(params$params_list, envir = environment())
  
  set.seed(123)
  e_ij_single <- matrix(rnorm(length(sigma_e_r), 0, 1), nrow = nrow(sigma_e_r)) * sigma_e_r
  
  one_indices <- rep(1, nrow(data))
  
  u <- alpha[, one_indices] + alpha_i[, data[, "i"]] + e_ij_single[, data[,"i"]] + theta_l[, data[,"l"]] + theta_li[, data[,"li"]] + 
    e_ij_single[, data[, "i"]]
  
  v <- kappa[, one_indices] + kappa_i[, data[, "i"]] + kappa_l[, data[, "l"]] + kappa_li[, data[, "li"]]
  
  if(model == 2){
    u <- u + outer(as.vector(omega), data[, "base_prev_mean"], `*`) + outer(as.vector(gamma), data[, "base_prev_diff"], `*`) + 
      omega_l[, data[, "l"]] * data[, "base_prev_mean"] +
      delta_l[, data[, "l"]] * data[, "base_prev_diff"]
    
    v <- v + outer(as.vector(delta), data[, "base_prev_diff"], `*`)
  }
  
  mean_prev <- (log(exp(v * year_max + u) + 1) - log(exp(u) + 1)) / (v * year_max)
  
  if(quantile_out == TRUE){
    mean_prev <- apply(mean_prev, 2, quantile, probs = c(probs_in)) |> t() |> as.data.frame() |> rename("low" = 1, "med" = 2, "up" = 3)
    rownames(mean_prev) <- NULL
  }
  
  return(mean_prev)
  
}

mean_time_r_eff_study_fun <- function(year_max, data, params, model, quantile_out = TRUE){
  
  # odds ratios
  log_o_r <- mean_time_l_or_study_fun(year_max, data, params = params, model, FALSE)
  o_r <- exp(log_o_r)
  
  # prevalence
  data_c <- data |> mutate(l = 1, li = 1)
  
  mean_prev_c <- mean_time_prev_study_fun(year_max, data_c, params, model, FALSE) |> t()
  
  r_eff <- (1 - (o_r / (1 - mean_prev_c + (mean_prev_c * o_r))))
  
  if(quantile_out == TRUE){
    r_eff <- apply(r_eff, 1, quantile, probs = probs_in) |> t() |> as.data.frame() |> rename("low" = 1, "med" = 2, "up" = 3)
    rownames(r_eff) <- NULL
  }
  
  return(r_eff)
}

############################
##### parameter values #####
############################

get_quantiles <- function(model_fit, param, dim = NULL){
  out <- as_draws_matrix(model_fit$draws(param))
  if(!is.null(dim)){
    out <- out[,dim]
  }
  prob_g0 <- round(sum(out > 0) / length(out), digits = 2)
  out <- round(quantile(out, probs = c(lower, 0.5, upper)), digits = 2)
  return(paste0(out[2], " (", out[1], " - ", out[3], " 95%CI, p(>0) = ", prob_g0,")"))
}

# model 0
# alpha
get_params <- function(fit){
  params <- c("$\\alpha^{\\text{pooled}}$")
  params_description <- ("pooled intercept")
  values <- c(get_quantiles(fit, "alpha", dim = NULL))
  # alpha_i
  params <- c(params, "$\\alpha_{i}$", "$\\alpha_{i}$", "$\\alpha_{i}$", "$\\alpha_i$")
  params_description <- c(params_description, "study effect (Tanzania (2015))", "study effect (Uganda (2017))", "study effect (Tanzania (2019))", "study effect (Benin (2020))")
  for(i in 1:4){
    values <- c(values, get_quantiles(fit, "alpha_i", dim = i))
  }
  # net effect
  params <- c(params, "$\\theta^{\\text{net}}_{1}$")
  params_description <- c(params_description, "pyrethroid-PBO pooled treatment effect")
  values <- c(values, get_quantiles(fit, "theta_l", dim = 2))
  
  params_description <- c(params_description, "pyrethroid-PBO-Tanzania (2015) interaction", "pyrethroid-PBO-Uganda (2017) interaction", "pyrethroid-PBO-Tanzania (2019) interaction")
  params <- c(params, "$\\theta^{\\text{study}}_{2}$", "$\\theta^{\\text{study}}_{3}$", "$\\theta^{\\text{study}}_{4}$")
  for(i in 2:4){
    values <- c(values, get_quantiles(fit, "theta_li", dim = i))
  }
  
  params <- c(params, "$\\theta^{\\text{net}}_{2}$")
  params_description <- c(params_description, "pyrethroid-pyrrole pooled treatment effect")
  values <- c(values, get_quantiles(fit, "theta_l", dim = 3))
  
  params_description <- c(params_description, "pyrethroid-pyrrole-Tanzania (2019) interaction", "pyrethroid-pyrrole-Benin (2020) interaction")
  params <- c(params, "$\\theta^{\\text{study}}_{5}$", "$\\theta^{\\text{study}}_{6}$")
  values <- c(values, get_quantiles(fit, "theta_li", dim = 5), get_quantiles(fit, "theta_li", dim = 6))
  
  params_description <- c(params_description, "Tanzania (2015) cluster random effect standard deviation", "Uganda (2017) cluster random effect standard deviation",
                          "Tanzania (2019) cluster random effect standard deviation", "Benin (2020) cluster random effect standard deviation")
  
  params <- c(params, "$\\sigma_{1}$", "$\\sigma_{2}$", "$\\sigma_{3}$", "$\\sigma_{4}$")
  values <- c(values, get_quantiles(fit, "sigma_e_r", dim = 1), get_quantiles(fit, "sigma_e_r", dim = 2),
              get_quantiles(fit, "sigma_e_r", dim = 3), get_quantiles(fit, "sigma_e_r", dim = 4))
  
  params_description <- c(params_description, "pyrethroid-PBO between study random treatment effect standard deviation", "pyrethroid-pyrrole between study random treatment effect standard deviation")
  params <- c(params, "$\\tau_{1}", "$\\tau_{2}")
  values <- c(values, get_quantiles(fit, "tau_sd_net_li", dim = 1), get_quantiles(fit, "tau_sd_net_li", dim = 2))
  
  return(list("params" = params, "params_description" = params_description, "values" = values))
}
