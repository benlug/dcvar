// Constant Copula VAR Model with Normal margins and Clayton copula

functions {
#include functions/clayton_copula.stan
#include functions/var_residuals.stan
}

data {
  int<lower=2> n_time;
  int<lower=2> D;
  matrix[n_time, D] Y;

  // Prior hyperparameters
  real<lower=0> sigma_mu_prior;
  real<lower=0> sigma_phi_prior;
  real<lower=0> sigma_eps_prior;
}

transformed data {
  int n_time_eff = n_time - 1;
}

parameters {
  vector[D] mu;
  matrix[D, D] Phi;
  vector<lower=0.01>[D] sigma_eps;
  real<lower=1e-6> theta;
}

transformed parameters {
  matrix[n_time_eff, D] eps = compute_var_residuals(Y, mu, Phi, n_time_eff, D);
  matrix[n_time_eff, D] eps_std;

  for (t in 1:n_time_eff) {
    eps_std[t, ] = eps[t, ] ./ sigma_eps';
  }
}

model {
  mu ~ normal(0, sigma_mu_prior);
  to_vector(Phi) ~ normal(0, sigma_phi_prior);
  sigma_eps ~ normal(0, sigma_eps_prior);
  theta ~ lognormal(0, 1);

  for (t in 1:n_time_eff) {
    for (d in 1:D) {
      target += normal_lpdf(eps[t, d] | 0, sigma_eps[d]);
    }
    target += clayton_copula_z_lpdf(to_vector(eps_std[t, ]) | theta);
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
    log_lik[t] += clayton_copula_z_lpdf(to_vector(eps_std[t, ]) | theta);

    {
      real u = uniform_rng(0, 1);
      real q = uniform_rng(0, 1);
      real v = pow(1 + pow(u, -theta) * (pow(q, -theta / (1 + theta)) - 1), -1 / theta);
      eps_rep[t, 1] = inv_Phi(u) * sigma_eps[1];
      eps_rep[t, 2] = inv_Phi(v) * sigma_eps[2];
    }
  }
}
