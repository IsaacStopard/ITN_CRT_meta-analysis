
functions{
  vector logit_prob_fun_no_cl(int N, matrix pmat_i, matrix pmat_l, matrix pmat_li, matrix pmat_it,
                        real alpha, vector alpha_i, vector theta_l, vector theta_li, vector beta_it, vector m_prob_s, vector d_prob_s,
                        real gamma, real kappa, vector kappa_l, vector years, real omega, real delta, vector m_net_s, vector d_net_s, real tau_m, real tau_d,
                        vector omega_l, vector delta_l, vector tau_lm, vector tau_ld){

                          vector[N] logit_prob;
                          vector[N] theta_ij;

                          // bed net efficacy parameters

                          theta_ij = (pmat_l * theta_l) + (pmat_li * theta_li);

                          logit_prob = alpha + // intercept
                          (pmat_i * alpha_i) + // study-specific intercept
                          theta_ij + // study-specific treatment effect
                          (pmat_it * beta_it) + kappa .* years + (pmat_l * kappa_l) .* years +
                          gamma .* d_prob_s + omega .* m_prob_s + delta .* d_prob_s .* years + tau_m .* m_net_s + tau_d .* d_net_s +
                          (pmat_l * omega_l) .* m_prob_s + (pmat_l * delta_l) .* d_prob_s + (pmat_l * tau_lm) .* m_net_s + (pmat_l * tau_ld) .* d_net_s;

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
  int<lower=0> N_it; // number of unique combinations of study and time of year
  int<lower=0> N_ij; // number of unique combinations of study and cluster in the actual dataset
  // design matrices for study, net and cluster
  matrix[N, N_i] pmat_i; // matrix of study predictors
  matrix[N, N_l] pmat_l; // matrix of net type predictors - pyrethroid only is always coded as 0
  matrix[N, N_li] pmat_li; // matrix of the nested study and net type predictors - if the net type is pyrethroid only then it is always coded as 0
  matrix[N, N_ij] pmat_ij; // matrix of the nested study and cluster predictors
  matrix[N, N_it] pmat_it; // matrix of the nested study and cluster predictors

  vector<lower=0>[N] years;

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

  array[N_ij] int<lower=0> pos_s; // number malaria positive at baseline
  array[N_ij] int<lower=0> test_s; // number tested at baseline

  array[N_ij] int<lower=0> base_net_pos; // number using net at baseline
  array[N_ij] int<lower=0> base_net_test; // number tested at baseline

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
  matrix[N_gq, N_it] pmat_it_gq; // matrix of the nested study and cluster predictors
  array[N_ij_gq] int<lower=0> ij_train_gq; // boolean indicating whether this value is fitted - should be 1 if so
  array[N_ij_unq_gq] int<lower=0> ij_unq_train_gq; // boolean indicating whether this value is fitted - should be 1 if so

  vector<lower=0, upper=1>[N_gq] m_prob_s_gq;
  vector<lower=-1, upper=1>[N_gq] d_prob_s_gq;

  vector<lower=0, upper=1>[N_gq] m_prob_s_gq_pooled;
  vector<lower=-1, upper=1>[N_gq] d_prob_s_gq_pooled;

  matrix[N_gq, N_i] pmat_i_pooled_gq;
  matrix[N_gq, N_li] pmat_li_pooled_gq;
  matrix[N_gq, N_it] pmat_it_pooled_gq;

  vector<lower=0>[N_gq] years_gq;

  vector<lower=0, upper = 1>[N_gq] m_net_s_gq;
  vector<lower=-1, upper = 1>[N_gq] d_net_s_gq;

  vector<lower=0, upper=1>[N_gq] m_net_s_gq_pooled;
  vector<lower=-1, upper=1>[N_gq] d_net_s_gq_pooled;

}

// The parameters accepted by the model. Our model
// accepts two parameters 'mu' and 'sigma'.
parameters{
  // intercept parameters
  real alpha; // study specific intercept parameter
  sum_to_zero_vector[N_i] alpha_i; // study specific intercept parameter

  // net treatment effect
  vector[N_l] theta_l;

  vector[N_it] beta_it;

  real kappa;
  vector[N_l] kappa_l;

  sum_to_zero_vector[N_i_pbo] e_raw_li_pbo;
  sum_to_zero_vector[N_i_pp] e_raw_li_pp;
  vector<lower=0>[N_l] tau_sd_li;

  real tau_m;
  real tau_d;
  vector[N_l] tau_lm;
  vector[N_l] tau_ld;

  // random effects
  vector[N_ij_train] e_raw_train; // cluster specific intercept
  vector<lower=0>[N_ij_unq_train] sigma_e_r_train; // study specific standard deviation for the cluster specific intercepts

  // baseline prevalence parameters
  real gamma;
  real omega;
  vector[N_l] omega_l;
  real delta;
  vector[N_l] delta_l;

  vector<lower=0, upper=1>[N_ij_train] prob_s_train;
  vector<lower=0, upper=1>[N_ij_train] base_net_train;

  real<lower=0> a;
  real<lower=0> b;
  real<lower=0> c;
  real<lower=0> d;
}

transformed parameters{
  vector[N] mu_ij;
  vector[N] z_ij;
  vector[N_ij_unq] sigma_e_r;
  vector[N_ij] e_raw;

  vector<lower=0, upper=1>[N_ij] prob_s;

  vector<lower=0, upper=1>[N_i] m_prob_s_i; // mean prob_s
  vector<lower=-1, upper=1>[N_ij] d_prob_s_ij; // difference

  vector<lower=0, upper=1>[N] m_prob_s; // mean prob_s
  vector<lower=-1, upper=1>[N] d_prob_s; // difference

  vector<lower=0>[N_i] count_m_prob = rep_vector(0, N_i);
  vector<lower=0>[N_i] total_m_prob = rep_vector(0, N_i);

  vector<lower=0, upper=1>[N_ij] base_net;

  vector<lower=0, upper=1>[N_i] m_net_s_i; // mean prob_s
  vector<lower=-1, upper=1>[N_ij] d_net_s_ij; // difference

  vector<lower=0, upper=1>[N] m_net_s; // mean prob_s
  vector<lower=-1, upper=1>[N] d_net_s; // difference

  vector<lower=0>[N_i] count_m_net = rep_vector(0, N_i);
  vector<lower=0>[N_i] total_m_net = rep_vector(0, N_i);

  vector[N_li] theta_li = append_row(e_raw_li_pbo * tau_sd_li[1], e_raw_li_pp * tau_sd_li[2]);

  // filling in the missing probabilities
  for(i in 1:N_ij){
    if(ij_train[i] == 1){
      e_raw[i] = e_raw_train[r_id_ij[i]];

      prob_s[i] = prob_s_train[r_id_ij[i]];
      count_m_prob[r_id[i]] += 1;
      total_m_prob[r_id[i]] += prob_s[i];

      base_net[i] = base_net_train[r_id_ij[i]];
      count_m_net[r_id[i]] += 1;
      total_m_net[r_id[i]] += base_net[i];

    } else{
      prob_s[i] = 0;
      base_net[i] = 0; // filling missing clusters with zero
      e_raw[i] = 0; // no random effect for testing data
    }
  }

  for(i in 1:N_i){
    if(ij_unq_train[i] == 1){
      m_prob_s_i[i] = total_m_prob[i] / count_m_prob[i];
      m_net_s_i[i] = total_m_net[i] / count_m_net[i];
    } else{
      m_prob_s_i[i] = 1;
      m_net_s_i[i] = 1;
    }
  }

  for(i in 1:N_ij){
     d_prob_s_ij[i] = prob_s[i] - m_prob_s_i[r_id[i]];
     d_net_s_ij[i] = base_net[i] - m_net_s_i[r_id[i]];
  }

  d_prob_s = pmat_ij * d_prob_s_ij;
  m_prob_s = pmat_i * m_prob_s_i;

  d_net_s = pmat_ij * d_net_s_ij;
  m_net_s = pmat_i * m_net_s_i;

  // cluster level random effects
  // filling in the missing values
  for(i in 1:N_ij_unq){
    if(ij_unq_train[i] == 0){
      sigma_e_r[i] = 0; // no random effect for testing data
    } else{
      sigma_e_r[i] = sigma_e_r_train[r_id_ij_unq[i]];
    }
  }

  // mean predictions
  mu_ij = logit_prob_fun_no_cl(N, pmat_i, pmat_l, pmat_li, pmat_it, alpha, alpha_i, theta_l, theta_li, beta_it, m_prob_s, d_prob_s, gamma, kappa, kappa_l, years, omega, delta, m_net_s, d_net_s, tau_m, tau_d, omega_l, delta_l, tau_lm, tau_ld);
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

  kappa ~ normal(0, prior_sd);
  kappa_l ~ normal(0, prior_sd);

  beta_it ~ normal(0, prior_sd);
  gamma ~ normal(0, prior_sd);
  omega ~ normal(0, prior_sd);
  omega_l ~ normal(0, prior_sd);
  delta ~ normal(0, prior_sd);
  delta_l ~ normal(0, prior_sd);
  tau_m ~ normal(0, prior_sd);
  tau_d ~ normal(0, prior_sd);
  tau_lm ~ normal(0, prior_sd);
  tau_ld ~ normal(0, prior_sd);

  prob_s_train ~ beta(a, b);

  base_net_train ~ beta(c, d);

  a ~ gamma(1.5, 1.5);
  b ~ gamma(1.5, 1.5);
  c ~ gamma(1.5, 1.5);
  d ~ gamma(1.5, 1.5);

  for(i in 1:N_ij){
    if(ij_train[i] == 1){
      pos_s[i] ~ binomial(test_s[i], prob_s_train[r_id_ij[i]]);
      base_net_pos[i] ~ binomial(base_net_test[i], base_net_train[r_id_ij[i]]);
    }
  }

  // only incrementing the log-likelihood for the training data
  pos[train_inds] ~ binomial_logit(test[train_inds], (mu_ij[train_inds] + z_ij[train_inds]));

}

generated quantities{
  // posterior predictive checks
  vector[N] log_lik;
  array[N] int y_rep;

  // vector<lower=0, upper=1>[N_gq] inv_logit_mu;
  vector<lower=0, upper=1>[N_gq] inv_logit_posterior;
  vector<lower=0, upper=1>[N_gq] inv_logit_posterior_pooled;

  vector[N_gq] mu = gq == 1? logit_prob_fun_no_cl(N_gq, pmat_i_gq, pmat_l_gq, pmat_li_gq, pmat_it_gq, alpha, alpha_i, theta_l, theta_li, beta_it, m_prob_s_gq, d_prob_s_gq, gamma, kappa, kappa_l, years_gq, omega, delta, m_net_s_gq, d_net_s_gq, tau_m, tau_d, omega_l, delta_l, tau_lm, tau_ld) : mu_ij;
  vector[N_gq] mu_pooled = logit_prob_fun_no_cl(N_gq, pmat_i_pooled_gq, pmat_l_gq, pmat_li_pooled_gq, pmat_it_pooled_gq, alpha, alpha_i, theta_l, theta_li, beta_it, m_prob_s_gq_pooled, d_prob_s_gq_pooled, gamma, kappa, kappa_l, years_gq, omega, delta, m_net_s_gq_pooled, d_net_s_gq_pooled, tau_m, tau_d, omega_l, delta_l, tau_lm, tau_ld);

   // odds ratios
  vector[N_gq] o_r = exp((pmat_l_gq * theta_l) + (pmat_li_gq * theta_li) + (pmat_l_gq * kappa_l) .* years_gq + (pmat_l_gq * omega_l) .* m_prob_s_gq + (pmat_l_gq * delta_l) .* d_prob_s_gq + (pmat_l_gq * tau_lm) .* m_net_s_gq + (pmat_l_gq * tau_ld) .* d_net_s_gq);

  vector[N_gq] o_r_pooled = exp((pmat_l_gq * theta_l) + (pmat_l_gq * kappa_l) .* years_gq + (pmat_l_gq * omega_l) .* m_prob_s_gq_pooled + (pmat_l_gq * delta_l) .* d_prob_s_gq_pooled + (pmat_l_gq * tau_lm) .* m_net_s_gq_pooled + (pmat_l_gq * tau_ld) .* d_net_s_gq_pooled);
  // random effects
  vector[N_gq] z_ij_full;
  vector[N_ij_gq] e_raw_full;
  vector[N_ij_unq_gq] sigma_e_r_full;
  real z_ij_pooled = normal_rng(0, 1) * exponential_rng(1);

  for(i in 1:N_ij_unq_gq){
    sigma_e_r_full[i] = ij_unq_train_gq[i] == 0 ? exponential_rng(1) : sigma_e_r[i]; // no random effect for testing data
  }

  // for the predictions ij train needs to be set to zero for all values
  for(i in 1:N_ij_gq){
    e_raw_full[i] = ij_train_gq[i] == 0 ? normal_rng(0, 1) : e_raw[i]; // no random effect for testing data
  }

  z_ij_full = random_effects_fun(N_gq, N_ij_gq, r_id_gq, pmat_ij_gq, sigma_e_r_full, e_raw_full);

  // inv_logit_mu = inv_logit(mu);
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

