// SEM Copula VAR(1) with fixed measurement model
// Normal-Gaussian margins, latent innovations as parameters
// rho = 0.97 * tanh(rho_raw) to avoid boundary singularity

functions {
#include functions/gaussian_copula.stan
}

data {
  int<lower=1> T;
  int<lower=1> J;
  matrix[T, 2 * J] y;              // indicators: y11..y1J y21..y2J
  vector[J] lambda;                 // fixed factor loadings
  real<lower=0> sigma_e;            // fixed measurement error SD

  // Prior hyperparameters
  real<lower=0> prior_mu_sd;
  real<lower=0> prior_phi_sd;
  real<lower=0> prior_sigma_sd;
  real<lower=0> prior_rho_sd;
}

parameters {
  vector[2] mu;
  // VAR coefficient bounds: each phi is restricted to (-0.99, 0.99) to prevent
  // the sampler from approaching the non-stationary boundary of the VAR process.
  // This is more restrictive than other dcVar models (dcVar, dcVar_constant,
  // dcVar_hmm), where Phi elements are unconstrained. As a consequence, models
  // with very strong autoregressive or cross-lag effects near +-1 cannot be
  // estimated with the SEM variant.
  real<lower=-0.99, upper=0.99> phi11;
  real<lower=-0.99, upper=0.99> phi12;
  real<lower=-0.99, upper=0.99> phi21;
  real<lower=-0.99, upper=0.99> phi22;
  vector<lower=0>[2] sigma;         // innovation SDs
  real rho_raw;                     // unconstrained; rho = 0.97 * tanh(rho_raw)
  matrix[T, 2] zeta;               // latent innovations
}

transformed parameters {
  matrix[2, 2] B;
  matrix[T, 2] state;
  real rho;

  // Map rho_raw in (-Inf, Inf) to rho in (-0.97, 0.97) via scaled tanh.
  // The 0.97 scaling prevents boundary singularity in the Gaussian copula
  // density at rho = +-1, where the log-density diverges. This is slightly
  // more restrictive than other dcVar models, which use tanh without scaling
  // and allow rho to approach +-1 more closely.
  rho = 0.97 * tanh(rho_raw);
  B[1, 1] = phi11; B[1, 2] = phi12;
  B[2, 1] = phi21; B[2, 2] = phi22;

  // Latent state recursion with x_0 = 0. The SEM simulator mirrors this
  // convention directly rather than returning post-burn-in states.
  {
    vector[2] s;
    for (t in 1:T) {
      vector[2] zt = to_vector(zeta[t]);
      s = (t == 1) ? (mu + zt) : (mu + B * s + zt);
      state[t, 1] = s[1];
      state[t, 2] = s[2];
    }
  }
}

model {
  // Priors
  mu ~ normal(0, prior_mu_sd);
  phi11 ~ normal(0, prior_phi_sd);
  phi12 ~ normal(0, prior_phi_sd);
  phi21 ~ normal(0, prior_phi_sd);
  phi22 ~ normal(0, prior_phi_sd);
  rho_raw ~ normal(0, prior_rho_sd);
  sigma ~ lognormal(0, prior_sigma_sd);

  // Innovation density: Normal margins + Gaussian copula (z-score form)
  {
    real log_sigma_sum = log(sigma[1]) + log(sigma[2]);
    for (t in 1:T) {
      vector[2] z;
      z[1] = zeta[t, 1] / sigma[1];
      z[2] = zeta[t, 2] / sigma[2];
      target += std_normal_lpdf(z[1]) + std_normal_lpdf(z[2]) - log_sigma_sum;
      target += gaussian_copula_z_lpdf(z | rho);
    }
  }

  // Measurement model (independent Normal errors)
  for (t in 1:T) {
    for (j in 1:J) {
      y[t, j]     ~ normal(lambda[j] * state[t, 1], sigma_e);
      y[t, J + j] ~ normal(lambda[j] * state[t, 2], sigma_e);
    }
  }
}

generated quantities {
  matrix[2, 2] Phi;
  vector[T] log_lik;

  Phi[1, 1] = phi11; Phi[1, 2] = phi12;
  Phi[2, 1] = phi21; Phi[2, 2] = phi22;

  {
    real log_sigma_sum = log(sigma[1]) + log(sigma[2]);
    for (t in 1:T) {
      vector[2] z;
      z[1] = zeta[t, 1] / sigma[1];
      z[2] = zeta[t, 2] / sigma[2];

      log_lik[t] = std_normal_lpdf(z[1]) + std_normal_lpdf(z[2]) - log_sigma_sum
                  + gaussian_copula_z_lpdf(z | rho);

      for (j in 1:J) {
        log_lik[t] += normal_lpdf(y[t, j] | lambda[j] * state[t, 1], sigma_e);
        log_lik[t] += normal_lpdf(y[t, J + j] | lambda[j] * state[t, 2], sigma_e);
      }
    }
  }
}
