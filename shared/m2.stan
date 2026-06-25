
functions{
  
  vector logit_prob_fun(int N, array[] int index_i, array[] int index_l, array[] int index_li, array[] int index_ij,
                        real alpha, vector alpha_i, vector theta_l, vector theta_li, vector base_prev_mean, vector base_prev_diff, 
                        real gamma, real kappa, vector kappa_i, vector kappa_l, vector kappa_li, vector years, real omega, real delta, vector omega_l, vector delta_l, vector e_ij){
                          
                          vector[N] logit_prob = alpha + // intercept
                          alpha_i[index_i] + // study-specific intercept
                          theta_l[index_l] + // pooled treatment effect
                          theta_li[index_li] + // study-specific treatment effect
                          (kappa + kappa_i[index_i]) .* years + // reference net time effect
                          (kappa_l[index_l] + kappa_li[index_li]) .* years + // dual AI net-time interaction effect
                          gamma * base_prev_diff + // within-study baseline prevalence effect for reference net
                          delta * base_prev_diff .* years + // within-study baseline prevalence-time interaction effect
                          delta_l[index_l] .* base_prev_diff + // within-study baseline prevalence-dual AI
                          omega * base_prev_mean + // mean study level baseline prevalence effect
                          omega_l[index_l] .* base_prev_mean +  // mean study level baseline prevalence-dual AI net interaction effect
                          e_ij[index_ij]; // cluster random intercept
                          
                        return logit_prob;
  }
  
  real partial_sum(array[] int pos_slice, int start, int end, array[] int test, vector global_logit_prob){
                       
                       return binomial_logit_lpmf(pos_slice | test[start:end], global_logit_prob[start:end]);
  }
}

data{
  // numbers
  int<lower=0> N; // number of data points
  int<lower=0> N_i; // number of unique CRTs - the intercept is study specific so all the studies must be included.
  int<lower=0> N_i_pbo;
  int<lower=0> N_i_pp;
  int<lower=0> N_l; // number of unique net types (excluding pyrethroid-only)
  int<lower=0> N_li; // number of unique combinations of net types (excluding pyrethroid-only) and study
  int<lower=0> N_ij; // number of unique combinations of study and cluster in the actual dataset
  
  array[N] int<lower=1, upper=N_i> index_i; // matrix of study predictors
  array[N] int<lower=1, upper=N_l> index_l; // matrix of net type predictors - pyrethroid only is always coded as 1
  array[N] int<lower=1, upper=N_li> index_li; // matrix of the nested study and net type predictors - if the net type is pyrethroid only then it is always coded as 1
  array[N] int<lower=1, upper=N_ij> index_ij; // matrix of the nested study and cluster predictors
  
  vector<lower=0>[N] years;
  
  // testing and training data
  int<lower=1, upper = N> N_train;
  array[N_train] int<lower=0> train_inds; // indices of training data
  int<lower=0, upper=(N-1)> N_test;
  array[N_test] int<lower=0> test_inds; // indices of training data
  
  // indexing on the cluster-level effects
  array[N_ij] int<lower=0, upper=N_ij> ij_index_i; // study index for each random effect predictor, previously r_id
  array[N_ij] int<lower=-1> ij_index_ij_train; // index to link the cluster ij to it's position in the training dataset, previously r_id_ij
  
  array[N_ij] int<lower=0, upper=1> ij_train_bin;
  
  int N_ij_train;
  
  // samples
  array[N] int<lower=0> pos; // number positive per cluster per survey
  array[N] int<lower=0> test; // number tested per cluster per survey
  
  // baseline prevalence samples
  array[N_ij] int<lower=0> pos_base_prev; // number malaria positive at baseline
  array[N_ij] int<lower=0> test_base_prev; // number tested at baseline
  
  real<lower=0> prior_sd;
  vector<lower=0, upper=1>[N_i] study_mean_base_prev;
  array[N_ij_train] int<lower=0, upper=N_i> ij_train_index_i;
  
}

transformed data{
  
  // used for calculating the uncertainty error model variables
  vector<lower=0>[N_i] count_base_prev_mean = rep_vector(0, N_i);
  
  array[N_train] int pos_train_data = pos[train_inds];
  array[N_train] int test_train_data = test[train_inds];
  array[N_train] int index_i_train_data = index_i[train_inds]; 
  array[N_train] int index_l_train_data = index_l[train_inds]; 
  array[N_train] int index_li_train_data = index_li[train_inds]; 
  array[N_train] int index_ij_train_data = index_ij[train_inds];
  vector[N_train] years_train_data = years[train_inds];
  
  for(i in 1:N_ij){
    if(ij_train_bin[i] == 1){
      count_base_prev_mean[ij_index_i[i]] += 1;
      }
      }
      
  real s_i = sqrt(N_i / (N_i - 1.0));
  real s_i_pp = sqrt(N_i_pp / (N_i_pp - 1.0));
  real s_i_pbo = sqrt(N_i_pbo / (N_i_pbo - 1.0));
}

// The parameters accepted by the model. Our model
// accepts two parameters 'mu' and 'sigma'.
parameters{
  // intercept parameters
  real alpha; // study specific intercept parameter
  sum_to_zero_vector[N_i] alpha_i; // study specific intercept parameter

  // net treatment effect
  vector[(N_l-1)] theta_l_raw;
  
  // time effect
  real kappa;
  
  vector[(N_l-1)] kappa_l_raw;
  
  // study-specific net random effect
  sum_to_zero_vector[N_i_pbo] e_raw_net_li_pbo;
  sum_to_zero_vector[N_i_pp] e_raw_net_li_pp;
  vector<lower=0>[(N_l-1)] tau_sd_net_li;
  
  // study-specific net-time random effect
  sum_to_zero_vector[N_i_pbo] e_raw_time_li_pbo;
  sum_to_zero_vector[N_i_pp] e_raw_time_li_pp;
  vector<lower=0>[(N_l-1)] tau_sd_time_li;
  
  sum_to_zero_vector[N_i] e_raw_time_i;
  real<lower=0> tau_sd_time_i;
  
  // random effects
  vector[N_ij_train] e_raw_train; // cluster specific intercept
  vector<lower=0>[N_i] sigma_e_r; // study specific standard deviation for the cluster specific intercepts

  // baseline prevalence parameters
  real gamma;
  real omega;
  vector[(N_l-1)] omega_l_raw;
  real delta;
  vector[(N_l-1)] delta_l_raw;

  vector<lower=0, upper=1>[N_ij_train] base_prev_train;

  
  vector<lower=0, upper=1>[N_i] BL_prev_mean;
  vector<lower=0>[N_i] BL_prev_conc;
}

transformed parameters{
  
  // fixed effects
  // adding zero for pyrethroid-only
  vector[N_l] theta_l = append_row(0, theta_l_raw); 
  vector[N_l] kappa_l = append_row(0, kappa_l_raw);
  vector[N_l] omega_l = append_row(0, omega_l_raw);
  vector[N_l] delta_l = append_row(0, delta_l_raw);
  
  // random effects
  vector[N_ij] e_raw = rep_vector(0, N_ij); // no random effect for testing data
  vector[N_ij] e_ij;
  vector[N_i] kappa_i = e_raw_time_i * tau_sd_time_i * s_i;
  
  // uncertainty error model
  vector<lower=0, upper=1>[N_ij] base_prev = rep_vector(0, N_ij);
  vector<lower=-1, upper=1>[N_ij] base_prev_diff_ij; // difference
  vector<lower=0, upper=1>[N_i] base_prev_mean_i = rep_vector(0, N_i); // mean base_prev
  
  vector<lower=0, upper=1>[N] base_prev_mean; // mean base_prev
  vector<lower=-1, upper=1>[N] base_prev_diff; // difference
  
  vector<lower=0>[N_i] total_base_prev_mean = rep_vector(0, N_i);
  vector[N_li] theta_li = append_row(0, append_row(e_raw_net_li_pbo * tau_sd_net_li[1] * s_i_pbo, e_raw_net_li_pp * tau_sd_net_li[2] * s_i_pp));
  vector[N_li] kappa_li = append_row(0, append_row(e_raw_time_li_pbo * tau_sd_time_li[1] * s_i_pbo, e_raw_time_li_pp * tau_sd_time_li[2] * s_i_pp));
  
  vector<lower=0>[N_i] a = BL_prev_mean .* BL_prev_conc;
  vector<lower=0>[N_i] b = BL_prev_conc .* (1 - BL_prev_mean);
  
  // filling in the missing probabilities
  // calculating the base prevalence for each cluster
  // if the cluster is in the training dataset then it is set at 0
  // the same is done for the cluster-level random intercept
  
  for(i in 1:N_ij){
    if(ij_train_bin[i] == 1){
      e_raw[i] = e_raw_train[ij_index_ij_train[i]];
      base_prev[i] = base_prev_train[ij_index_ij_train[i]];
      total_base_prev_mean[ij_index_i[i]] += base_prev[i];
    } else{
      base_prev[i] = pos_base_prev[i] * 1.0 / test_base_prev[i];
    }
  }
  
  for(i in 1:N_i){
    base_prev_mean_i[i] = total_base_prev_mean[i] / count_base_prev_mean[i];
  }
  
  base_prev_diff_ij = base_prev - base_prev_mean_i[ij_index_i];
  base_prev_diff = base_prev_diff_ij[index_ij];
  base_prev_mean = base_prev_mean_i[index_i];
  
  e_ij = sigma_e_r[ij_index_i] .* e_raw; // equivalent to e_[i] ~ normal(0, sigma_e_r[[ij_index_i][i]])
    
}

model{
  
  vector[N_train] base_prev_diff_train_data = base_prev_diff[train_inds];
  vector[N_train] base_prev_mean_train_data = base_prev_mean[train_inds];
  
  int grainsize = 1;
  
  // intercept parameters
  alpha ~ normal(0, prior_sd);
  alpha_i ~ normal(0, prior_sd * s_i);

  sigma_e_r ~ exponential(1);
  e_raw_train ~ std_normal();

  tau_sd_net_li ~ exponential(1);
  e_raw_net_li_pbo ~ std_normal();
  e_raw_net_li_pp ~ std_normal();
  
  tau_sd_time_li ~ exponential(1);
  e_raw_time_li_pbo ~ std_normal();
  e_raw_time_li_pp ~ std_normal();
  
  tau_sd_time_i ~ exponential(1);
  e_raw_time_i ~ std_normal();

  theta_l_raw ~ normal(0, prior_sd);

  kappa ~ normal(0, prior_sd);
  kappa_l_raw ~ normal(0, prior_sd);

  gamma ~ normal(0, prior_sd);
  omega ~ normal(0, prior_sd);
  omega_l_raw ~ normal(0, prior_sd);
  delta ~ normal(0, prior_sd);
  delta_l_raw ~ normal(0, prior_sd);

  base_prev_train ~ beta(a[ij_train_index_i], b[ij_train_index_i]);

  BL_prev_mean ~ beta(study_mean_base_prev * 25, (1 - study_mean_base_prev) * 25);
  BL_prev_conc ~ exponential(0.1);
  
  for(i in 1:N_ij){
    if(ij_train_bin[i] == 1){
      pos_base_prev[i] ~ binomial(test_base_prev[i], base_prev_train[ij_index_ij_train[i]]);
    }
  }
  
  vector[N_train] global_logit_prob_train = logit_prob_fun(N_train, index_i_train_data, index_l_train_data, 
                                                           index_li_train_data, index_ij_train_data, alpha, 
                                                           alpha_i, theta_l, theta_li, base_prev_mean_train_data, 
                                                           base_prev_diff_train_data, gamma, kappa, kappa_i, kappa_l, kappa_li, years_train_data, 
                                                           omega, delta, omega_l, delta_l, e_ij);

  // only incrementing the log-likelihood for the training data
  target += reduce_sum(partial_sum, pos_train_data, grainsize, test_train_data, global_logit_prob_train);

}

generated quantities{
  // posterior predictive checks
  vector[N] logit_lp;
  vector[N] inv_logit_lp;
  array[N] int y_rep;
  array[N_test] real log_lik;
  vector[N_ij] e_ij_full;
  
  for(i in 1:N_ij){
    e_ij_full[i] = ij_train_bin[i] == 0 ? normal_rng(0.0, 1.0) * sigma_e_r[ij_index_i[i]] : e_ij[i];
  }
  
  logit_lp = logit_prob_fun(N, index_i, index_l, index_li, index_ij,
                            alpha, alpha_i, theta_l, theta_li, base_prev_mean, base_prev_diff, 
                            gamma, kappa, kappa_i, kappa_l, kappa_li, years, omega, delta, omega_l, delta_l, e_ij_full);
  
  inv_logit_lp = inv_logit(logit_lp);
  
  if(N_test > 0){
    
    for(i in 1:N_test){
      log_lik[i] = binomial_logit_lpmf(pos[test_inds[i]] | test[test_inds[i]], logit_lp[test_inds[i]]);
    }
  }
  
  y_rep = binomial_rng(test, inv_logit_lp);

}

