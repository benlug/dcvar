// Multilevel Copula VAR(1) with per-variable (mixed) margins
// Random per-unit Phi_i (partial pooling) and a single global Gaussian copula
// correlation, where each variable can use its own marginal family. Person-mean
// centered data (no intercept). The per-dimension marginal dispatch mirrors the
// single-level mixed models; the exponential / gamma feasibility bound is taken
// globally across all units (as in multilevel_EG.stan).
//
// Family codes (per dimension): 1 = normal, 2 = exponential, 3 = skew_normal,
// 4 = gamma. All per-family parameters are declared unconditionally and a
// dimension that does not use a parameter samples it from its (proper) prior.

functions {
#include functions/gaussian_copula_uv.stan
}

data {
  int<lower=1> N;                              // Number of units
  int<lower=2> n_time;                         // Time points per unit
  array[N] matrix[n_time, 2] y;                // Centered data: y[i] is n_time x 2
  array[2] int<lower=1, upper=4> family;       // Per-dimension margin family code
  vector[2] skew_direction;                    // Consulted only for exp/gamma dims

  // Prior hyperparameters
  real<lower=0> prior_phi_bar_sd;
  real<lower=0> prior_tau_phi_scale;
  real<lower=0> prior_sigma_sd;
  real<lower=0> prior_rho_sd;
}

transformed data {
  int n_time_eff = n_time - 1;
  real SQRT_2_OVER_PI = sqrt(2.0 / pi());
}

parameters {
  vector[4] phi_bar;
  vector<lower=0>[4] tau_phi;
  matrix[N, 4] z_phi;
  real<lower=-1, upper=1> rho;

  // Union of per-dimension marginal parameters
  vector<lower=0.01>[2] sigma_eps;             // normal innovation SDs
  vector[2] eta;                               // exponential / gamma scale (unconstrained)
  vector<lower=0>[2] omega;                    // skew_normal scale
  vector<lower=-1, upper=1>[2] delta;          // skew_normal CP skewness
  vector<lower=0>[2] shape_gam;                // gamma shape (per dimension)
}

transformed parameters {
  matrix[N, 4] phi_unit;
  array[N] matrix[2, 2] Phi_T;

  // Derived per-dimension scale parameters (filled on their own family's dims).
  vector[2] sigma_exp = rep_vector(0, 2);
  vector[2] rate_exp = rep_vector(0, 2);
  vector[2] mean_gam = rep_vector(0, 2);
  vector[2] rate_gam = rep_vector(0, 2);
  vector[2] sigma_gam = rep_vector(0, 2);

  // Skew-normal CP -> DP transform (used where family == 3)
  vector[2] alpha = delta ./ sqrt(1 - square(delta));
  vector[2] xi = -omega .* (delta * SQRT_2_OVER_PI);

  for (i in 1:N) {
    for (k in 1:4) phi_unit[i, k] = phi_bar[k] + tau_phi[k] * z_phi[i, k];
    Phi_T[i][1, 1] = phi_unit[i, 1];
    Phi_T[i][2, 1] = phi_unit[i, 2];
    Phi_T[i][1, 2] = phi_unit[i, 3];
    Phi_T[i][2, 2] = phi_unit[i, 4];
  }

  // Global feasibility bound (over all units/times) for exponential/gamma dims
  for (j in 1:2) {
    if (family[j] == 2 || family[j] == 4) {
      real m = negative_infinity();
      for (i in 1:N) {
        for (t in 2:n_time) {
          row_vector[2] res = y[i][t] - y[i][t - 1] * Phi_T[i];
          m = fmax(m, -skew_direction[j] * res[j]);
        }
      }
      real lb = fmax(m, 0);
      if (family[j] == 2) {
        sigma_exp[j] = lb + exp(eta[j]) + 1e-9;
        rate_exp[j] = 1.0 / sigma_exp[j];
      } else {
        mean_gam[j] = lb + exp(eta[j]) + 1e-9;
        rate_gam[j] = shape_gam[j] / mean_gam[j];
        sigma_gam[j] = mean_gam[j] / sqrt(shape_gam[j]);
      }
    }
  }
}

model {
  // Hyperpriors and global priors
  phi_bar ~ normal(0, prior_phi_bar_sd);
  tau_phi ~ student_t(3, 0, prior_tau_phi_scale);
  to_vector(z_phi) ~ std_normal();
  rho ~ normal(0, prior_rho_sd);

  // Union-parameter priors kept proper unconditionally
  sigma_eps ~ normal(0, prior_sigma_sd);
  omega ~ normal(0, 1);
  delta ~ normal(0, 0.5);
  shape_gam ~ lognormal(log(1), 0.5);

  // Per-dimension induced scale priors / eta priors
  for (j in 1:2) {
    if (family[j] == 2) {
      target += lognormal_lpdf(sigma_exp[j] | 0, prior_sigma_sd) + eta[j];
    } else if (family[j] == 4) {
      target += lognormal_lpdf(mean_gam[j] | 0, prior_sigma_sd) + eta[j];
    } else {
      eta[j] ~ std_normal();
    }
  }

  // Likelihood
  for (i in 1:N) {
    for (t in 2:n_time) {
      row_vector[2] res = y[i][t] - y[i][t - 1] * Phi_T[i];
      vector[2] u;
      for (j in 1:2) {
        if (family[j] == 1) {
          target += normal_lpdf(res[j] | 0, sigma_eps[j]);
          u[j] = Phi(res[j] / sigma_eps[j]);
        } else if (family[j] == 2) {
          real x_shifted = sigma_exp[j] + skew_direction[j] * res[j];
          target += exponential_lpdf(x_shifted | rate_exp[j]);
          u[j] = exponential_cdf(x_shifted | rate_exp[j]);
          if (skew_direction[j] < 0) u[j] = 1.0 - u[j];
        } else if (family[j] == 3) {
          target += skew_normal_lpdf(res[j] | xi[j], omega[j], alpha[j]);
          u[j] = skew_normal_cdf(res[j] | xi[j], omega[j], alpha[j]);
        } else {
          real x_shifted = mean_gam[j] + skew_direction[j] * res[j];
          target += gamma_lpdf(x_shifted | shape_gam[j], rate_gam[j]);
          u[j] = gamma_cdf(x_shifted | shape_gam[j], rate_gam[j]);
          if (skew_direction[j] < 0) u[j] = 1.0 - u[j];
        }
      }
      target += gaussian_copula_uv_lpdf(u | rho);
    }
  }
}

generated quantities {
  array[N] matrix[2, 2] Phi;
  vector[N] spectral_radius;
  int<lower=0> n_nonstationary = 0;
  matrix[N, n_time_eff] log_lik;

  for (i in 1:N) {
    Phi[i] = Phi_T[i]';

    real tr = Phi[i][1, 1] + Phi[i][2, 2];
    real det_phi = Phi[i][1, 1] * Phi[i][2, 2] - Phi[i][1, 2] * Phi[i][2, 1];
    real disc = square(tr) - 4 * det_phi;
    if (disc >= 0) {
      real sqrt_disc = sqrt(disc);
      spectral_radius[i] = fmax(abs(0.5 * (tr + sqrt_disc)), abs(0.5 * (tr - sqrt_disc)));
    } else {
      spectral_radius[i] = sqrt(abs(det_phi));
    }
    if (spectral_radius[i] >= 1.0) n_nonstationary += 1;

    for (t in 2:n_time) {
      row_vector[2] res = y[i][t] - y[i][t - 1] * Phi_T[i];
      vector[2] u;
      log_lik[i, t - 1] = 0;
      for (j in 1:2) {
        if (family[j] == 1) {
          log_lik[i, t - 1] += normal_lpdf(res[j] | 0, sigma_eps[j]);
          u[j] = Phi(res[j] / sigma_eps[j]);
        } else if (family[j] == 2) {
          real x_shifted = sigma_exp[j] + skew_direction[j] * res[j];
          log_lik[i, t - 1] += exponential_lpdf(x_shifted | rate_exp[j]);
          u[j] = exponential_cdf(x_shifted | rate_exp[j]);
          if (skew_direction[j] < 0) u[j] = 1.0 - u[j];
        } else if (family[j] == 3) {
          log_lik[i, t - 1] += skew_normal_lpdf(res[j] | xi[j], omega[j], alpha[j]);
          u[j] = skew_normal_cdf(res[j] | xi[j], omega[j], alpha[j]);
        } else {
          real x_shifted = mean_gam[j] + skew_direction[j] * res[j];
          log_lik[i, t - 1] += gamma_lpdf(x_shifted | shape_gam[j], rate_gam[j]);
          u[j] = gamma_cdf(x_shifted | shape_gam[j], rate_gam[j]);
          if (skew_direction[j] < 0) u[j] = 1.0 - u[j];
        }
      }
      log_lik[i, t - 1] += gaussian_copula_uv_lpdf(u | rho);
    }
  }
}
