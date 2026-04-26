// Naive SEM score VAR(1) with Exponential margins and Gaussian copula
// Observed y is the T x 2 matrix of row-mean factor scores.

functions {
#include functions/gaussian_copula_uv.stan
}

data {
  int<lower=1> n_time;
  matrix[n_time, 2] y;
  vector[2] skew_direction;

  // Prior hyperparameters
  real<lower=0> prior_mu_sd;
  real<lower=0> prior_phi_sd;
  real<lower=0> prior_sigma_sd;
  real<lower=0> prior_rho_sd;
}

parameters {
  vector[2] mu;
  real<lower=-0.99, upper=0.99> phi11;
  real<lower=-0.99, upper=0.99> phi12;
  real<lower=-0.99, upper=0.99> phi21;
  real<lower=-0.99, upper=0.99> phi22;
  vector[2] eta;
  real rho_raw;
}

transformed parameters {
  matrix[2, 2] Phi_T;
  matrix[n_time, 2] eps;
  vector[2] sigma_lb;
  vector[2] sigma_exp;
  vector[2] rate_exp;
  real rho;

  rho = 0.97 * tanh(rho_raw);
  Phi_T[1, 1] = phi11; Phi_T[2, 1] = phi12;
  Phi_T[1, 2] = phi21; Phi_T[2, 2] = phi22;

  {
    matrix[n_time, 2] ylag = rep_matrix(0.0, n_time, 2);
    if (n_time > 1) {
      ylag[2:n_time, ] = y[1:(n_time - 1), ];
    }
    eps = y - (rep_matrix(mu', n_time) + ylag * Phi_T);
  }

  for (i in 1:2) {
    real m = -skew_direction[i] * eps[1, i];
    for (t in 2:n_time) {
      m = fmax(m, -skew_direction[i] * eps[t, i]);
    }
    sigma_lb[i] = fmax(m, 0);
    sigma_exp[i] = sigma_lb[i] + exp(eta[i]) + 1e-9;
    rate_exp[i] = 1.0 / sigma_exp[i];
  }
}

model {
  mu ~ normal(0, prior_mu_sd);
  phi11 ~ normal(0, prior_phi_sd);
  phi12 ~ normal(0, prior_phi_sd);
  phi21 ~ normal(0, prior_phi_sd);
  phi22 ~ normal(0, prior_phi_sd);
  rho_raw ~ normal(0, prior_rho_sd);

  for (i in 1:2) {
    target += lognormal_lpdf(sigma_exp[i] | 0, prior_sigma_sd) + eta[i];
  }

  for (t in 1:n_time) {
    vector[2] u;
    for (i in 1:2) {
      real x_shifted = sigma_exp[i] + skew_direction[i] * eps[t, i];
      target += exponential_lpdf(x_shifted | rate_exp[i]);
      u[i] = exponential_cdf(x_shifted | rate_exp[i]);
      if (skew_direction[i] < 0) {
        u[i] = 1.0 - u[i];
      }
    }
    target += gaussian_copula_uv_lpdf(u | rho);
  }
}

generated quantities {
  matrix[2, 2] Phi;
  vector[n_time] log_lik;
  vector[2] b_gq;
  vector[2] slack;

  Phi[1, 1] = phi11; Phi[1, 2] = phi12;
  Phi[2, 1] = phi21; Phi[2, 2] = phi22;

  for (i in 1:2) {
    b_gq[i] = sigma_lb[i];
    slack[i] = exp(eta[i]);
  }

  for (t in 1:n_time) {
    vector[2] u;
    log_lik[t] = 0;
    for (i in 1:2) {
      real x_shifted = sigma_exp[i] + skew_direction[i] * eps[t, i];
      log_lik[t] += exponential_lpdf(x_shifted | rate_exp[i]);
      u[i] = exponential_cdf(x_shifted | rate_exp[i]);
      if (skew_direction[i] < 0) {
        u[i] = 1.0 - u[i];
      }
    }
    log_lik[t] += gaussian_copula_uv_lpdf(u | rho);
  }
}
