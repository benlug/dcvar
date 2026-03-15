// Dynamic Copula VAR Model (DC-VAR) - Non-Centered Parameterization
// Time-varying Gaussian copula dependence with constant VAR coefficients
// Uses non-centered parameterization to avoid funnel geometry and divergences
// Author: Working implementation for research project
// Date: 2026

functions {
#include functions/gaussian_copula.stan
}

data {
  int<lower=2> T;                    // Number of time points
  int<lower=2> D;                    // Number of variables (typically 2)
  matrix[T, D] Y;                    // Observed data (T x D)

  // Prior hyperparameters
  real<lower=0> sigma_mu_prior;      // Prior SD for intercepts
  real<lower=0> sigma_phi_prior;     // Prior SD for VAR coefficients
  real<lower=0> sigma_eps_prior;     // Prior SD for innovation SDs
  real<lower=0> sigma_omega_prior;   // Prior SD for rho process innovation
  real<lower=0> rho_init_prior_sd;   // Prior SD for initial rho (on z-scale)
}

transformed data {
  int T_eff = T - 1;  // Effective time points for VAR (exclude first)
}

parameters {
  // VAR parameters (constant over time)
  vector[D] mu;                      // Intercepts
  matrix[D, D] Phi;                  // VAR(1) coefficient matrix (unconstrained: stationarity
                                     // is not enforced, allowing unit roots / near-unit-root
                                     // dynamics. The normal prior provides soft regularization.)
  vector<lower=0.01>[D] sigma_eps;   // Marginal innovation SDs (lower bound prevents numerical issues)

  // Time-varying copula parameter (NON-CENTERED)
  real z_rho_init;                   // Initial value of z-transformed rho
  real<lower=0.001> sigma_omega;     // Innovation SD for rho process (lower bound prevents numerical issues)
  vector[T_eff] omega_raw;           // RAW innovations (std_normal)
}

transformed parameters {
  // Residuals from VAR
  matrix[T_eff, D] eps;

  // Standardized residuals (z-scores for copula)
  matrix[T_eff, D] eps_std;

  // Time-varying rho (on z-scale and original scale)
  vector[T_eff] z_rho;
  vector[T_eff] rho;

  // Compute VAR residuals
  for (t in 1:T_eff) {
    // VAR(1): y_t = mu + Phi * (y_{t-1} - mu) + eps_t
    vector[D] y_prev = to_vector(Y[t, ]);
    vector[D] y_curr = to_vector(Y[t + 1, ]);
    vector[D] y_hat = mu + Phi * (y_prev - mu);

    eps[t, ] = to_row_vector(y_curr - y_hat);
  }

  // NON-CENTERED PARAMETERIZATION for random walk
  // Instead of: z_rho[t] = z_rho[t-1] + sigma_omega * omega[t]
  // We use: z_rho[t] = z_rho_init + sigma_omega * cumsum(omega_raw)
  // This avoids the funnel geometry when sigma_omega is small
  {
    real cumsum_omega = 0;
    for (t in 1:T_eff) {
      cumsum_omega += omega_raw[t];
      z_rho[t] = z_rho_init + sigma_omega * cumsum_omega;
    }
  }

  // Transform to correlation scale via tanh: maps (-inf, inf) to (-1, 1)
  // Mathematically equivalent to inv_fisher_z but numerically stable for large |z|
  for (t in 1:T_eff) {
    rho[t] = tanh(z_rho[t]);
  }

  // Standardize residuals (z-scores used directly in copula, avoiding
  // the Phi -> clamp -> inv_Phi roundtrip that truncates tail information)
  for (t in 1:T_eff) {
    eps_std[t, ] = eps[t, ] ./ sigma_eps';
  }
}

model {
  // Priors for VAR parameters
  mu ~ normal(0, sigma_mu_prior);
  to_vector(Phi) ~ normal(0, sigma_phi_prior);

  // Exponential priors for scale parameters (avoid half-normal mass near zero)
  sigma_eps ~ exponential(1.0 / sigma_eps_prior);

  // Priors for rho process
  z_rho_init ~ normal(0, rho_init_prior_sd);
  // Use exponential prior for sigma_omega (proper prior for scale parameter)
  // This is more appropriate than normal(0, .) which allows negative values
  sigma_omega ~ exponential(1.0 / sigma_omega_prior);

  // NON-CENTERED: sample raw innovations from standard normal
  omega_raw ~ std_normal();

  // Likelihood: Gaussian copula for residuals at each time point
  for (t in 1:T_eff) {
    // Marginal likelihoods (normal)
    for (d in 1:D) {
      target += normal_lpdf(eps[t, d] | 0, sigma_eps[d]);
    }

    // Copula contribution (using z-scores directly)
    target += gaussian_copula_z_lpdf(to_vector(eps_std[t, ]) | rho[t]);
  }
}

generated quantities {
  // Log-likelihood for model comparison
  vector[T_eff] log_lik;

  // Posterior predictive checks
  matrix[T_eff, D] eps_rep;

  for (t in 1:T_eff) {
    // Log-likelihood
    log_lik[t] = 0;
    for (d in 1:D) {
      log_lik[t] += normal_lpdf(eps[t, d] | 0, sigma_eps[d]);
    }
    log_lik[t] += gaussian_copula_z_lpdf(to_vector(eps_std[t, ]) | rho[t]);

    // Generate replicated data (for posterior predictive checks)
    // Sample from bivariate normal with correlation rho[t]
    real z1_rep = std_normal_rng();
    real z2_rep = rho[t] * z1_rep + sqrt(1 - square(rho[t])) * std_normal_rng();
    eps_rep[t, 1] = z1_rep * sigma_eps[1];
    eps_rep[t, 2] = z2_rep * sigma_eps[2];
  }
}
