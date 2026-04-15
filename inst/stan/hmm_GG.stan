// HMM Copula VAR Model with Gamma-Gaussian Margins

functions {
#include functions/gaussian_copula_uv.stan
#include functions/var_residuals.stan
#include functions/hmm_algorithms.stan
}

data {
  int<lower=2> T;
  int<lower=2> D;
  matrix[T, D] Y;
  int<lower=2> K;
  vector[D] skew_direction;

  real<lower=0> sigma_mu_prior;
  real<lower=0> sigma_phi_prior;
  real<lower=0> kappa;
  real<lower=0> alpha_off;
  real<lower=0> z_rho_prior_sd;
}

transformed data {
  int T_eff = T - 1;
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
  vector[D] eta;
  real<lower=0> shape_gam;
  ordered[K] z_rho;
  simplex[K] pi0;
  array[K] simplex[K] A;
}

transformed parameters {
  vector[K] rho_state;
  matrix[T_eff, D] eps;
  matrix[T_eff, K] obs_ll;
  matrix[K, K] log_A;
  vector[K] log_pi0;
  matrix[T_eff, K] log_alpha;

  for (k in 1:K) rho_state[k] = tanh(z_rho[k]);

  for (j in 1:K) {
    log_pi0[j] = log(pi0[j]);
    for (k in 1:K) log_A[j, k] = log(A[j][k]);
  }

  eps = compute_var_residuals(Y, mu, Phi, T_eff, D);

  {
    real sqrt_shape = sqrt(shape_gam);
    real sigma_eps = 1e-9;
    vector[D] mean_lb;
    vector[D] mean_gam;
    vector[D] sigma_gam;
    vector[D] rate_gam;

    for (i in 1:D) {
      real m = -skew_direction[i] * eps[1, i];
      for (t in 2:T_eff) m = fmax(m, -skew_direction[i] * eps[t, i]);
      mean_lb[i] = fmax(m, 0);
      mean_gam[i] = mean_lb[i] + exp(eta[i]) + sigma_eps;
      sigma_gam[i] = mean_gam[i] / sqrt_shape;
    }
    rate_gam = shape_gam ./ mean_gam;

    for (t in 1:T_eff) {
      real marginal_ll = 0;
      vector[2] u_vec;
      for (i in 1:D) {
        real x_shifted = mean_gam[i] + skew_direction[i] * eps[t, i];
        marginal_ll += gamma_lpdf(x_shifted | shape_gam, rate_gam[i]);
        u_vec[i] = gamma_cdf(x_shifted | shape_gam, rate_gam[i]);
        if (skew_direction[i] < 0) u_vec[i] = 1.0 - u_vec[i];
      }
      for (k in 1:K) {
        obs_ll[t, k] = marginal_ll + gaussian_copula_uv_lpdf(u_vec | rho_state[k]);
      }
    }
  }

  log_alpha = hmm_forward(obs_ll, log_A, log_pi0, T_eff, K);
}

model {
  mu ~ normal(0, sigma_mu_prior);
  to_vector(Phi) ~ normal(0, sigma_phi_prior);
  z_rho ~ normal(0, z_rho_prior_sd);
  pi0 ~ dirichlet(rep_vector(1.0, K));
  for (k in 1:K) A[k] ~ dirichlet(dirichlet_prior[k]);
  shape_gam ~ lognormal(log(1), 0.5);

  // Induced prior on sigma_gam
  {
    real sqrt_shape = sqrt(shape_gam);
    real sigma_eps = 1e-9;
    vector[D] mean_lb;
    vector[D] mean_gam;
    for (i in 1:D) {
      real m = -skew_direction[i] * eps[1, i];
      for (t in 2:T_eff) m = fmax(m, -skew_direction[i] * eps[t, i]);
      mean_lb[i] = fmax(m, 0);
      mean_gam[i] = mean_lb[i] + exp(eta[i]) + sigma_eps;
      target += lognormal_lpdf(mean_gam[i] | 0, 0.5) + eta[i];
    }
  }

  target += log_sum_exp(to_vector(log_alpha[T_eff, ]));
}

generated quantities {
  matrix[T_eff, K] gamma;
  array[T_eff] int viterbi_state;
  vector[T_eff] rho_hmm;
  vector[T_eff] log_lik;
  matrix[T_eff, D] eps_rep;
  vector[D] sigma_gam;
  vector[D] b_gq;
  vector[D] rate_gam;

  gamma = hmm_state_posteriors(log_alpha, obs_ll, log_A, T_eff, K);

  viterbi_state = hmm_viterbi(obs_ll, log_A, log_pi0, T_eff, K);

  rho_hmm = hmm_rho_average(gamma, rho_state, T_eff, K);

  log_lik = hmm_log_lik(log_alpha, T_eff, K);

  {
    real sqrt_shape = sqrt(shape_gam);
    real sigma_eps = 1e-9;
    vector[D] mean_gam;
    for (i in 1:D) {
      real m = -skew_direction[i] * eps[1, i];
      for (t in 2:T_eff) m = fmax(m, -skew_direction[i] * eps[t, i]);
      b_gq[i] = fmax(m, 0);
      mean_gam[i] = b_gq[i] + exp(eta[i]) + sigma_eps;
      sigma_gam[i] = mean_gam[i] / sqrt_shape;
      rate_gam[i] = shape_gam / mean_gam[i];
    }
  }

  // NOTE: eps_rep contains copula-level z-scores, not gamma residuals,
  // because Stan lacks a gamma inverse CDF. plot_ppc() rejects gamma
  // fits until replicated residuals are available on the margin scale.
  for (t in 1:T_eff) {
    real z1_rep = std_normal_rng();
    real z2_rep = rho_hmm[t] * z1_rep + sqrt(1 - square(rho_hmm[t])) * std_normal_rng();
    eps_rep[t, 1] = z1_rep;
    eps_rep[t, 2] = z2_rep;
  }
}
