
functions{
  vector logit_prob_fun_no_cl(int N, matrix pmat_i, matrix pmat_l, matrix pmat_li,
                              real alpha, vector alpha_i, vector theta_l, vector theta_li){

                          vector[N] logit_prob;
                          vector[N] theta_ij;

                          // bed net efficacy parameters

                          theta_ij = (pmat_l * theta_l) + (pmat_li * theta_li);

                          logit_prob = alpha + (pmat_i * alpha_i) + theta_ij;

                          return logit_prob;
  }

  vector random_effects_fun(int N, int N_ij, array[] int r_id, matrix pmat_ij, vector sigma_e_r, vector e_raw){
    vector[N_ij] e_ij;
    vector[N] e_out;

    for(i in 1:N_ij){
      e_ij[i] = sigma_e_r[r_id[i]] * e_raw[i]; // equivalent to e_[i] ~ normal(0, sigma_e_r[[r_id][i]])
      }

      e_out = pmat_ij * e_ij;
      return e_out;
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
  // design matrices for study, net and cluster
  matrix[N, N_i] pmat_i; // matrix of study predictors
  matrix[N, N_l] pmat_l; // matrix of net type predictors - pyrethroid only is always coded as 0
  matrix[N, N_li] pmat_li; // matrix of the nested study and net type predictors - if the net type is pyrethroid only then it is always coded as 0
  matrix[N, N_ij] pmat_ij; // matrix of the nested study and cluster predictors

  // outcome
  array[N] int<lower=0> pos; // number positive per cluster per survey
  array[N] int<lower=0> test; // number tested per cluster per survey

  real<lower=0> prior_sd;

  // indexing on the random effects: group index for each random effect predictor

  array[N_ij] int<lower=0, upper=N_ij> r_id;

  array[N_ij] int<lower=0> ij_train; // boolean indicating whether this value is fitted - should be 1 if so
  array[N_ij] int<lower=-1> r_id_ij;

  int<lower=0> N_ij_unq; // unique number of cluster and study combinations (including pyrethroid-only nets) - for the random effects
  array[N_ij_unq] int<lower=0> ij_unq_train; // boolean indicating whether this value is fitted - should be 1 if so
  array[N_ij_unq] int<lower=-1> r_id_ij_unq; // index of ij_unq_train in ij_unq - set as zero is study is missing (when ij_unq_train == 0)

  int<lower=0> N_ij_train;
  // number of sigma parameters that are fitted
  int<lower=0> N_ij_unq_train; // can be equal to N_ij_unq

  int<lower=0, upper = N> N_train;
  array[N_train] int<lower=0> train_inds; // indices of training data

  // generated quantities
  int gq; // indicator as to whether just for generated quantities
  int N_gq;
  int N_ij_gq;
  int N_ij_unq_gq;
  matrix[N_gq, N_ij_gq] pmat_ij_gq;
  array[N_ij_gq] int<lower=0, upper=N_ij_gq> r_id_gq;
  // design matrices for study, net and cluster
  matrix[N_gq, N_i] pmat_i_gq; // matrix of study predictors
  matrix[N_gq, N_l] pmat_l_gq; // matrix of net type predictors - pyrethroid only is always coded as 0
  matrix[N_gq, N_li] pmat_li_gq; // matrix of the nested study and net type predictors - if the net type is pyrethroid only then it is always coded as 0

  matrix[N_gq, N_i] pmat_i_pooled_gq;
  matrix[N_gq, N_li] pmat_li_pooled_gq;

  array[N_ij_gq] int<lower=0> ij_train_gq; // boolean indicating whether this value is fitted - should be 1 if so
  array[N_ij_unq_gq] int<lower=0> ij_unq_train_gq; // boolean indicating whether this value is fitted - should be 1 if so

}

// The parameters accepted by the model. Our model
// accepts two parameters 'mu' and 'sigma'.
parameters{
  // intercept parameters
  real alpha;
  sum_to_zero_vector[N_i] alpha_i; // study specific intercept parameter

  // net treatment effect
  vector[N_l] theta_l;

  //sum_to_zero_vector[N_i_pbo] theta_li_pbo;
  //sum_to_zero_vector[N_i_pp] theta_li_pp;
  sum_to_zero_vector[N_i_pbo] e_raw_li_pbo;
  sum_to_zero_vector[N_i_pp] e_raw_li_pp;
  vector<lower=0>[N_l] tau_sd_li;

  // random effects
  vector[N_ij_train] e_raw_train; // cluster specific intercept
  vector<lower=0>[N_ij_unq_train] sigma_e_r_train; // study specific standard deviation for the cluster specific intercepts
}

transformed parameters{
  vector[N] mu_ij;
  vector[N] z_ij;
  vector[N_ij_unq] sigma_e_r;
  vector[N_ij] e_raw;
  vector[N_li] theta_li = append_row(e_raw_li_pbo * tau_sd_li[1], e_raw_li_pp * tau_sd_li[2]);

  // mean predictions
  mu_ij = logit_prob_fun_no_cl(N, pmat_i, pmat_l, pmat_li, alpha, alpha_i, theta_l, theta_li);

  // cluster level random effects
  // filling in the missing values
  for(i in 1:N_ij_unq){
    if(ij_unq_train[i] == 0){
      sigma_e_r[i] = 0; // no random effect for testing data
    } else{
      sigma_e_r[i] = sigma_e_r_train[r_id_ij_unq[i]];
    }
  }

  for(i in 1:N_ij){
    if(ij_train[i] == 0){
      e_raw[i] = 0; // no random effect for testing data
    } else{
      e_raw[i] = e_raw_train[r_id_ij[i]];
    }
  }

  z_ij = random_effects_fun(N, N_ij, r_id, pmat_ij, sigma_e_r, e_raw);

}

model{

  // intercept parameters
  alpha ~ normal(0, prior_sd);
  alpha_i ~ normal(0, prior_sd);

  sigma_e_r_train ~ exponential(1);
  e_raw_train ~ std_normal();

  tau_sd_li ~ exponential(1);
  e_raw_li_pbo ~ std_normal();
  e_raw_li_pp ~ std_normal();

  theta_l ~ normal(0, prior_sd);

  // only incrementing the log-likelihood for the training data
  pos[train_inds] ~ binomial_logit(test[train_inds], (mu_ij[train_inds] + z_ij[train_inds]));
}

generated quantities{
  // posterior predictive checks
  vector[N] log_lik;
  array[N] int y_rep;
  // odds ratios
  vector<lower=0, upper=1>[N_gq] inv_logit_mu;
  vector<lower=0, upper=1>[N_gq] inv_logit_posterior;
  vector<lower=0, upper=1>[N_gq] inv_logit_posterior_pooled;
  real z_ij_pooled = normal_rng(0, 1) * exponential_rng(1);

  vector[N_gq] mu = gq == 1? logit_prob_fun_no_cl(N_gq, pmat_i_gq, pmat_l_gq, pmat_li_gq, alpha, alpha_i, theta_l, theta_li) : mu_ij;
  vector[N_gq] mu_pooled = logit_prob_fun_no_cl(N_gq, pmat_i_pooled_gq, pmat_l_gq, pmat_li_pooled_gq, alpha, alpha_i, theta_l, theta_li);

   // odds ratios
  vector[N_gq] o_r = exp((pmat_l_gq * theta_l) + (pmat_li_gq * theta_li));

  vector[N_gq] o_r_pooled = exp((pmat_l_gq * theta_l));

  // random effects
  vector[N_gq] z_ij_full;
  vector[N_ij_gq] e_raw_full;
  vector[N_ij_unq_gq] sigma_e_r_full;

  for(i in 1:N_ij_unq_gq){
    sigma_e_r_full[i] = ij_unq_train_gq[i] == 0 ? exponential_rng(1) : sigma_e_r[i]; // no random effect for testing data
  }

  // for the predictions ij train needs to be set to zero for all values
  for(i in 1:N_ij_gq){
    e_raw_full[i] = ij_train_gq[i] == 0 ? normal_rng(0, 1) : e_raw[i]; // no random effect for testing data
  }

  z_ij_full = random_effects_fun(N_gq, N_ij_gq, r_id_gq, pmat_ij_gq, sigma_e_r_full, e_raw_full);

  inv_logit_mu = inv_logit(mu);
  inv_logit_posterior = inv_logit(mu + z_ij_full);
  inv_logit_posterior_pooled = inv_logit(mu_pooled + z_ij_pooled);

  if(gq == 0){

    //posterior predictive checks
   for(i in 1:N){
       log_lik[i] = binomial_logit_lpmf(pos[i] | test[i], mu[i] + z_ij_full[i]);
    }

    y_rep = binomial_rng(test, inv_logit_posterior);
  } else{
    log_lik = rep_vector(positive_infinity(), N);
    for(i in 1:N){
       y_rep[i] = -1;
    }
  }
}

