// Covariate Dynamic Copula VAR Model without residual drift.
// Gaussian innovations with Fisher-z correlation explained by observed
// covariates only: zeta_i = beta_0 + x_{i+1}' beta.

functions {
#include functions/gaussian_copula.stan
#include functions/var_residuals.stan
}

data {
  int<lower=2> n_time;
  int<lower=2> D;
  matrix[n_time, D] Y;

  real<lower=0> sigma_mu_prior;
  real<lower=0> sigma_phi_prior;
  real<lower=0> sigma_eps_prior;
  real<lower=0> rho_init_prior_sd;

  int<lower=1> P;
  matrix[n_time, P] X;
  real<lower=0> sigma_beta_prior;
}

transformed data {
  int n_time_eff = n_time - 1;
}

parameters {
  vector[D] mu;
  matrix[D, D] Phi;
  vector<lower=0.01>[D] sigma_eps;

  real beta_0;
  vector[P] beta;
}

transformed parameters {
  matrix[n_time_eff, D] eps = compute_var_residuals(Y, mu, Phi, n_time_eff, D);
  matrix[n_time_eff, D] eps_std;
  vector[n_time_eff] z_rho;
  vector[n_time_eff] rho;

  for (t in 1:n_time_eff) {
    z_rho[t] = beta_0 + dot_product(to_vector(X[t + 1, ]), beta);
    rho[t] = tanh(z_rho[t]);
    eps_std[t, ] = eps[t, ] ./ sigma_eps';
  }
}

model {
  mu ~ normal(0, sigma_mu_prior);
  to_vector(Phi) ~ normal(0, sigma_phi_prior);
  sigma_eps ~ exponential(1.0 / sigma_eps_prior);

  beta_0 ~ normal(0, rho_init_prior_sd);
  beta ~ normal(0, sigma_beta_prior);

  for (t in 1:n_time_eff) {
    for (d in 1:D) {
      target += normal_lpdf(eps[t, d] | 0, sigma_eps[d]);
    }
    target += gaussian_copula_z_lpdf(to_vector(eps_std[t, ]) | rho[t]);
  }
}

generated quantities {
  vector[n_time_eff] log_lik;
  matrix[n_time_eff, D] eps_rep;

  for (t in 1:n_time_eff) {
    log_lik[t] = 0;
    for (d in 1:D) {
      log_lik[t] += normal_lpdf(eps[t, d] | 0, sigma_eps[d]);
    }
    log_lik[t] += gaussian_copula_z_lpdf(to_vector(eps_std[t, ]) | rho[t]);

    real z1_rep = std_normal_rng();
    real z2_rep = rho[t] * z1_rep + sqrt(1 - square(rho[t])) * std_normal_rng();
    eps_rep[t, 1] = z1_rep * sigma_eps[1];
    eps_rep[t, 2] = z2_rep * sigma_eps[2];
  }
}
