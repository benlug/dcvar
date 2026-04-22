// Constant Copula VAR Model
// VAR(1) with time-invariant Gaussian copula dependence (single fixed rho)
// Baseline model for comparison against the time-varying DC-VAR
// Author: Working implementation for research project
// Date: 2026

functions {
#include functions/gaussian_copula.stan
#include functions/var_residuals.stan
}

data {
  int<lower=2> n_time;                    // Number of time points
  int<lower=2> D;                    // Number of variables (typically 2)
  matrix[n_time, D] Y;                    // Observed data (n_time x D)

  // Prior hyperparameters
  real<lower=0> sigma_mu_prior;      // Prior SD for intercepts
  real<lower=0> sigma_phi_prior;     // Prior SD for VAR coefficients
  real<lower=0> sigma_eps_prior;     // Prior SD for innovation SDs
  real<lower=0> z_rho_prior_sd;      // Prior SD for rho (on z/tanh scale)
}

transformed data {
  int n_time_eff = n_time - 1;  // Effective time points for VAR (exclude first)
}

parameters {
  // VAR parameters (constant over time)
  vector[D] mu;                      // Intercepts
  matrix[D, D] Phi;                  // VAR(1) coefficient matrix (unconstrained: stationarity
                                     // is not enforced, allowing unit roots / near-unit-root
                                     // dynamics. The normal prior provides soft regularization.)
  vector<lower=0.01>[D] sigma_eps;   // Marginal innovation SDs (lower bound prevents numerical issues)

  // Constant copula parameter
  real z_rho;                        // Copula correlation on Fisher-z scale
}

transformed parameters {
  // Residuals from VAR
  matrix[n_time_eff, D] eps = compute_var_residuals(Y, mu, Phi, n_time_eff, D);

  // Standardized residuals (z-scores for copula)
  matrix[n_time_eff, D] eps_std;

  // Constant rho (on original scale) via tanh: maps (-inf, inf) to (-1, 1)
  // Mathematically equivalent to inv_fisher_z but numerically stable for large |z|
  real rho = tanh(z_rho);

  // Standardize residuals (z-scores used directly in copula)
  for (t in 1:n_time_eff) {
    eps_std[t, ] = eps[t, ] ./ sigma_eps';
  }
}

model {
  // Priors for VAR parameters
  mu ~ normal(0, sigma_mu_prior);
  to_vector(Phi) ~ normal(0, sigma_phi_prior);

  // Exponential priors for scale parameters (avoid half-normal mass near zero)
  sigma_eps ~ exponential(1.0 / sigma_eps_prior);

  // Prior for constant rho (on Fisher-z scale)
  z_rho ~ normal(0, z_rho_prior_sd);

  // Likelihood: Gaussian copula for residuals at each time point
  for (t in 1:n_time_eff) {
    // Marginal likelihoods (normal)
    for (d in 1:D) {
      target += normal_lpdf(eps[t, d] | 0, sigma_eps[d]);
    }

    // Copula contribution (using z-scores directly, constant rho)
    target += gaussian_copula_z_lpdf(to_vector(eps_std[t, ]) | rho);
  }
}

generated quantities {
  // Log-likelihood for model comparison (LOO-CV)
  vector[n_time_eff] log_lik;

  // Posterior predictive checks
  matrix[n_time_eff, D] eps_rep;

  for (t in 1:n_time_eff) {
    // Log-likelihood
    log_lik[t] = 0;
    for (d in 1:D) {
      log_lik[t] += normal_lpdf(eps[t, d] | 0, sigma_eps[d]);
    }
    log_lik[t] += gaussian_copula_z_lpdf(to_vector(eps_std[t, ]) | rho);

    // Generate replicated data (for posterior predictive checks)
    // Sample from bivariate normal with constant correlation rho
    {
      real z1_rep = std_normal_rng();
      real z2_rep = rho * z1_rep + sqrt(1 - square(rho)) * std_normal_rng();
      eps_rep[t, 1] = z1_rep * sigma_eps[1];
      eps_rep[t, 2] = z2_rep * sigma_eps[2];
    }
  }
}
