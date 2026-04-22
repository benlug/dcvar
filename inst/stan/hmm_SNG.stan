// HMM Copula VAR Model with Skew-Normal-Gaussian Margins

functions {
#include functions/gaussian_copula_uv.stan
#include functions/var_residuals.stan
#include functions/hmm_algorithms.stan
}

data {
  int<lower=2> n_time;
  int<lower=2> D;
  matrix[n_time, D] Y;
  int<lower=2> K;

  real<lower=0> sigma_mu_prior;
  real<lower=0> sigma_phi_prior;
  real<lower=0> kappa;
  real<lower=0> alpha_off;
  real<lower=0> z_rho_prior_sd;
}

transformed data {
  int n_time_eff = n_time - 1;
  real SQRT_2_OVER_PI = sqrt(2.0 / pi());
  array[K] vector[K] dirichlet_prior;
  for (k in 1:K) {
    for (j in 1:K) {
      dirichlet_prior[k][j] = (k == j) ? kappa : alpha_off;
    }
  }
}

parameters {
  vector[D] mu;
  matrix[D, D] Phi;
  vector<lower=0>[D] omega;
  vector<lower=-1, upper=1>[D] delta;
  ordered[K] z_rho;
  simplex[K] pi0;
  array[K] simplex[K] A;
}

transformed parameters {
  vector[K] rho_state;
  matrix[n_time_eff, D] eps;
  vector[D] sn_alpha;  // renamed to avoid conflict with alpha_off
  vector[D] xi;
  matrix[n_time_eff, K] obs_ll;
  matrix[K, K] log_A;
  vector[K] log_pi0;
  matrix[n_time_eff, K] log_alpha;

  sn_alpha = delta ./ sqrt(1 - square(delta));
  xi = -omega .* (delta * SQRT_2_OVER_PI);

  for (k in 1:K) rho_state[k] = tanh(z_rho[k]);

  for (j in 1:K) {
    log_pi0[j] = log(pi0[j]);
    for (k in 1:K) log_A[j, k] = log(A[j][k]);
  }

  eps = compute_var_residuals(Y, mu, Phi, n_time_eff, D);

  for (t in 1:n_time_eff) {
    real marginal_ll = 0;
    vector[2] u_vec;
    for (i in 1:D) {
      marginal_ll += skew_normal_lpdf(eps[t, i] | xi[i], omega[i], sn_alpha[i]);
      u_vec[i] = skew_normal_cdf(eps[t, i] | xi[i], omega[i], sn_alpha[i]);
    }
    for (k in 1:K) {
      obs_ll[t, k] = marginal_ll + gaussian_copula_uv_lpdf(u_vec | rho_state[k]);
    }
  }

  log_alpha = hmm_forward(obs_ll, log_A, log_pi0, n_time_eff, K);
}

model {
  mu ~ normal(0, sigma_mu_prior);
  to_vector(Phi) ~ normal(0, sigma_phi_prior);
  z_rho ~ normal(0, z_rho_prior_sd);
  pi0 ~ dirichlet(rep_vector(1.0, K));
  for (k in 1:K) A[k] ~ dirichlet(dirichlet_prior[k]);
  omega ~ normal(0, 1);
  delta ~ normal(0, 0.5);

  target += log_sum_exp(to_vector(log_alpha[n_time_eff, ]));
}

generated quantities {
  matrix[n_time_eff, K] gamma;
  array[n_time_eff] int viterbi_state;
  vector[n_time_eff] rho_hmm;
  vector[n_time_eff] log_lik;
  matrix[n_time_eff, D] eps_rep;

  gamma = hmm_state_posteriors(log_alpha, obs_ll, log_A, n_time_eff, K);

  viterbi_state = hmm_viterbi(obs_ll, log_A, log_pi0, n_time_eff, K);

  rho_hmm = hmm_rho_average(gamma, rho_state, n_time_eff, K);

  log_lik = hmm_log_lik(log_alpha, n_time_eff, K);

  // NOTE: eps_rep contains copula-level z-scores, not skew-normal residuals,
  // because Stan lacks a skew-normal inverse CDF. plot_ppc() rejects
  // skew-normal fits until replicated residuals are available on the
  // margin scale.
  for (t in 1:n_time_eff) {
    real z1_rep = std_normal_rng();
    real z2_rep = rho_hmm[t] * z1_rep + sqrt(1 - square(rho_hmm[t])) * std_normal_rng();
    eps_rep[t, 1] = z1_rep;
    eps_rep[t, 2] = z2_rep;
  }
}
