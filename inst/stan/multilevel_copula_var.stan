// Multilevel Copula VAR(1) with random Phi_i and global rho
// Normal-Gaussian margins, person-mean centered data (no intercept)
// Stationarity is monitored in generated quantities but not imposed.

functions {
#include functions/gaussian_copula.stan
}

data {
  int<lower=1> N;                    // Number of units
  int<lower=2> T;                    // Time points per unit
  array[N] matrix[T, 2] y;          // Centered data: y[i] is T x 2

  // Prior hyperparameters
  real<lower=0> prior_phi_bar_sd;
  real<lower=0> prior_tau_phi_scale;
  real<lower=0> prior_sigma_sd;
  real<lower=0> prior_rho_sd;
}

parameters {
  // Population-level
  vector[4] phi_bar;                 // Population mean: (phi11, phi12, phi21, phi22)
  vector<lower=0>[4] tau_phi;        // Population SD per element

  // Unit-level (non-centered)
  matrix[N, 4] z_phi;               // Standardized deviates

  // Global parameters
  vector<lower=0>[2] sigma;         // Innovation SDs
  real<lower=-1, upper=1> rho;      // Global copula correlation
}

transformed parameters {
  matrix[N, 4] phi_unit;
  array[N] matrix[2, 2] Phi_T;      // Transposed for row-vector multiplication

  for (i in 1:N) {
    for (k in 1:4) {
      phi_unit[i, k] = phi_bar[k] + tau_phi[k] * z_phi[i, k];
    }
    Phi_T[i][1, 1] = phi_unit[i, 1];  // phi11
    Phi_T[i][2, 1] = phi_unit[i, 2];  // phi12
    Phi_T[i][1, 2] = phi_unit[i, 3];  // phi21
    Phi_T[i][2, 2] = phi_unit[i, 4];  // phi22
  }
}

model {
  // Hyperpriors
  phi_bar ~ normal(0, prior_phi_bar_sd);
  tau_phi ~ student_t(3, 0, prior_tau_phi_scale);

  // Unit-level latents
  to_vector(z_phi) ~ std_normal();

  // Global priors
  sigma ~ normal(0, prior_sigma_sd);
  rho ~ normal(0, prior_rho_sd);

  // Likelihood
  {
    real log_sigma_sum = log(sigma[1]) + log(sigma[2]);

    for (i in 1:N) {
      for (t in 2:T) {
        row_vector[2] pred = y[i][t - 1] * Phi_T[i];
        row_vector[2] res = y[i][t] - pred;
        row_vector[2] z = res ./ sigma';

        target += std_normal_lpdf(z[1]) + std_normal_lpdf(z[2]) - log_sigma_sum;
        target += gaussian_copula_z_lpdf(to_vector(z) | rho);
      }
    }
  }
}

generated quantities {
  // Unit-specific Phi matrices (standard orientation)
  array[N] matrix[2, 2] Phi;

  // Stationarity monitoring
  vector[N] spectral_radius;
  int<lower=0> n_nonstationary = 0;

  // Per-unit log-likelihood (for LOO)
  vector[N] log_lik;

  for (i in 1:N) {
    Phi[i] = Phi_T[i]';

    // Eigenvalue computation for 2x2
    real tr = Phi[i][1, 1] + Phi[i][2, 2];
    real det_phi = Phi[i][1, 1] * Phi[i][2, 2] - Phi[i][1, 2] * Phi[i][2, 1];
    real disc = square(tr) - 4 * det_phi;

    if (disc >= 0) {
      real sqrt_disc = sqrt(disc);
      real lam1 = 0.5 * (tr + sqrt_disc);
      real lam2 = 0.5 * (tr - sqrt_disc);
      spectral_radius[i] = fmax(abs(lam1), abs(lam2));
    } else {
      spectral_radius[i] = sqrt(abs(det_phi));
    }

    if (spectral_radius[i] >= 1.0) n_nonstationary += 1;

    // Log-likelihood for this unit
    log_lik[i] = 0;
    {
      real log_sigma_sum = log(sigma[1]) + log(sigma[2]);
      for (t in 2:T) {
        row_vector[2] pred = y[i][t - 1] * Phi_T[i];
        row_vector[2] res = y[i][t] - pred;
        row_vector[2] z = res ./ sigma';
        log_lik[i] += std_normal_lpdf(z[1]) + std_normal_lpdf(z[2]) - log_sigma_sum;
        log_lik[i] += gaussian_copula_z_lpdf(to_vector(z) | rho);
      }
    }
  }
}
