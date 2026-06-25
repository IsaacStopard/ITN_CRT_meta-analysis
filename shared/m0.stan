
functions{
  
  vector logit_prob_fun(int N, array[] int index_i, array[] int index_l, array[] int index_li, array[] int index_ij,
                        real alpha, vector alpha_i, vector theta_l, vector theta_li, vector e_ij){
                          
                          vector[N] logit_prob = alpha + // intercept
                          alpha_i[index_i] + // study-specific intercept
                          theta_l[index_l] + // pooled treatment effect
                          theta_li[index_li] + // study-specific treatment effect
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
  
  real<lower=0> prior_sd;

}

transformed data{
  
  // used for calculating the uncertainty error model variables
  
  array[N_train] int pos_train_data = pos[train_inds];
  array[N_train] int test_train_data = test[train_inds];
  array[N_train] int index_i_train_data = index_i[train_inds]; 
  array[N_train] int index_l_train_data = index_l[train_inds]; 
  array[N_train] int index_li_train_data = index_li[train_inds]; 
  array[N_train] int index_ij_train_data = index_ij[train_inds];
  
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
  
  // study-specific net random effect
  sum_to_zero_vector[N_i_pbo] e_raw_net_li_pbo;
  sum_to_zero_vector[N_i_pp] e_raw_net_li_pp;
  vector<lower=0>[(N_l-1)] tau_sd_net_li;
  
  // random effects
  vector[N_ij_train] e_raw_train; // cluster specific intercept
  vector<lower=0>[N_i] sigma_e_r; // study specific standard deviation for the cluster specific intercepts
  
}

transformed parameters{
  
  // fixed effects
  // adding zero for pyrethroid-only
  vector[N_l] theta_l = append_row(0, theta_l_raw); 
  
  // random effects
  vector[N_ij] e_raw = rep_vector(0, N_ij); // no random effect for testing data
  vector[N_ij] e_ij;
  
  // uncertainty error model
  vector[N_li] theta_li = append_row(0, append_row(e_raw_net_li_pbo * tau_sd_net_li[1] * s_i_pbo, e_raw_net_li_pp * tau_sd_net_li[2] * s_i_pp));
  
  for(i in 1:N_ij){
    if(ij_train_bin[i] == 1){
      e_raw[i] = e_raw_train[ij_index_ij_train[i]];
    }
  }
  
  e_ij = sigma_e_r[ij_index_i] .* e_raw; // equivalent to e_[i] ~ normal(0, sigma_e_r[[ij_index_i][i]])
    
}

model{
  
  int grainsize = 1;
  
  // intercept parameters
  alpha ~ normal(0, prior_sd);
  alpha_i ~ normal(0, prior_sd * s_i);

  sigma_e_r ~ exponential(1);
  e_raw_train ~ std_normal();

  tau_sd_net_li ~ exponential(1);
  e_raw_net_li_pbo ~ std_normal();
  e_raw_net_li_pp ~ std_normal();

  theta_l_raw ~ normal(0, prior_sd);
  
  vector[N_train] global_logit_prob_train = logit_prob_fun(N_train, index_i_train_data, index_l_train_data, 
                                                           index_li_train_data, index_ij_train_data, alpha, 
                                                           alpha_i, theta_l, theta_li, 
                                                           e_ij);

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
                            alpha, alpha_i, theta_l, theta_li, e_ij_full);
  
  inv_logit_lp = inv_logit(logit_lp);
  
  if(N_test > 0){
    
    for(i in 1:N_test){
      log_lik[i] = binomial_logit_lpmf(pos[test_inds[i]] | test[test_inds[i]], logit_lp[test_inds[i]]);
    }
  }
  
  y_rep = binomial_rng(test, inv_logit_lp);

}

