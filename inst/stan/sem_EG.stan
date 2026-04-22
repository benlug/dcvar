// SEM Copula VAR(1) with fixed measurement model
// Exponential-Gaussian latent innovation margins, latent innovations as parameters
// rho = 0.97 * tanh(rho_raw) to avoid boundary singularity

functions {
#include functions/gaussian_copula_uv.stan
}

data {
  int<lower=1> n_time;
  int<lower=1> J;
  matrix[n_time, 2 * J] y;              // indicators: y11..y1J y21..y2J
  vector[J] lambda;                // fixed factor loadings
  real<lower=0> sigma_e;           // fixed measurement error SD
  vector[2] skew_direction;        // +1 / -1 exponential skew directions

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
  vector[2] eta;                   // unconstrained: sigma_exp = lb + exp(eta)
  real rho_raw;                    // unconstrained; rho = 0.97 * tanh(rho_raw)
  matrix[n_time, 2] zeta;               // latent innovations
}

transformed parameters {
  matrix[2, 2] B;
  matrix[n_time, 2] state;
  vector[2] sigma_lb;
  vector[2] sigma_exp;
  vector[2] rate_exp;
  real rho;

  rho = 0.97 * tanh(rho_raw);
  B[1, 1] = phi11; B[1, 2] = phi12;
  B[2, 1] = phi21; B[2, 2] = phi22;

  for (i in 1:2) {
    real m = -skew_direction[i] * zeta[1, i];
    for (t in 2:n_time) {
      m = fmax(m, -skew_direction[i] * zeta[t, i]);
    }
    sigma_lb[i] = fmax(m, 0);
    sigma_exp[i] = sigma_lb[i] + exp(eta[i]) + 1e-9;
    rate_exp[i] = 1.0 / sigma_exp[i];
  }

  {
    vector[2] s;
    for (t in 1:n_time) {
      vector[2] zt = to_vector(zeta[t]);
      s = (t == 1) ? (mu + zt) : (mu + B * s + zt);
      state[t, 1] = s[1];
      state[t, 2] = s[2];
    }
  }
}

model {
  mu ~ normal(0, prior_mu_sd);
  phi11 ~ normal(0, prior_phi_sd);
  phi12 ~ normal(0, prior_phi_sd);
  phi21 ~ normal(0, prior_phi_sd);
  phi22 ~ normal(0, prior_phi_sd);
  rho_raw ~ normal(0, prior_rho_sd);

  // Induced prior on sigma_exp via change-of-variables.
  for (i in 1:2) {
    target += lognormal_lpdf(sigma_exp[i] | 0, prior_sigma_sd) + eta[i];
  }

  // Latent innovation density: shifted exponential margins + Gaussian copula.
  for (t in 1:n_time) {
    vector[2] u_vec;
    for (i in 1:2) {
      real x_shifted = sigma_exp[i] + skew_direction[i] * zeta[t, i];
      target += exponential_lpdf(x_shifted | rate_exp[i]);
      u_vec[i] = exponential_cdf(x_shifted | rate_exp[i]);
      if (skew_direction[i] < 0) {
        u_vec[i] = 1.0 - u_vec[i];
      }
    }
    target += gaussian_copula_uv_lpdf(u_vec | rho);
  }

  // Measurement model (independent Normal errors)
  for (t in 1:n_time) {
    for (j in 1:J) {
      y[t, j]     ~ normal(lambda[j] * state[t, 1], sigma_e);
      y[t, J + j] ~ normal(lambda[j] * state[t, 2], sigma_e);
    }
  }
}

generated quantities {
  matrix[2, 2] Phi;
  vector[n_time] log_lik;

  Phi[1, 1] = phi11; Phi[1, 2] = phi12;
  Phi[2, 1] = phi21; Phi[2, 2] = phi22;

  for (t in 1:n_time) {
    vector[2] u_vec;
    log_lik[t] = 0;

    for (i in 1:2) {
      real x_shifted = sigma_exp[i] + skew_direction[i] * zeta[t, i];
      log_lik[t] += exponential_lpdf(x_shifted | rate_exp[i]);
      u_vec[i] = exponential_cdf(x_shifted | rate_exp[i]);
      if (skew_direction[i] < 0) {
        u_vec[i] = 1.0 - u_vec[i];
      }
    }
    log_lik[t] += gaussian_copula_uv_lpdf(u_vec | rho);

    for (j in 1:J) {
      log_lik[t] += normal_lpdf(y[t, j] | lambda[j] * state[t, 1], sigma_e);
      log_lik[t] += normal_lpdf(y[t, J + j] | lambda[j] * state[t, 2], sigma_e);
    }
  }
}
