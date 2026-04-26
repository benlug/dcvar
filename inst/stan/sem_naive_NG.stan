// Naive SEM score VAR(1) with Normal margins and Gaussian copula
// Observed y is the T x 2 matrix of row-mean factor scores.

functions {
#include functions/gaussian_copula.stan
}

data {
  int<lower=1> n_time;
  matrix[n_time, 2] y;

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
  vector<lower=0>[2] sigma;
  real rho_raw;
}

transformed parameters {
  matrix[2, 2] Phi_T;
  matrix[n_time, 2] eps;
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
}

model {
  mu ~ normal(0, prior_mu_sd);
  phi11 ~ normal(0, prior_phi_sd);
  phi12 ~ normal(0, prior_phi_sd);
  phi21 ~ normal(0, prior_phi_sd);
  phi22 ~ normal(0, prior_phi_sd);
  rho_raw ~ normal(0, prior_rho_sd);
  sigma ~ lognormal(0, prior_sigma_sd);

  {
    real log_sigma_sum = log(sigma[1]) + log(sigma[2]);
    for (t in 1:n_time) {
      vector[2] z;
      z[1] = eps[t, 1] / sigma[1];
      z[2] = eps[t, 2] / sigma[2];
      target += std_normal_lpdf(z[1]) + std_normal_lpdf(z[2]) - log_sigma_sum;
      target += gaussian_copula_z_lpdf(z | rho);
    }
  }
}

generated quantities {
  matrix[2, 2] Phi;
  vector[n_time] log_lik;

  Phi[1, 1] = phi11; Phi[1, 2] = phi12;
  Phi[2, 1] = phi21; Phi[2, 2] = phi22;

  {
    real log_sigma_sum = log(sigma[1]) + log(sigma[2]);
    for (t in 1:n_time) {
      vector[2] z;
      z[1] = eps[t, 1] / sigma[1];
      z[2] = eps[t, 2] / sigma[2];
      log_lik[t] = std_normal_lpdf(z[1]) + std_normal_lpdf(z[2]) - log_sigma_sum
                   + gaussian_copula_z_lpdf(z | rho);
    }
  }
}
