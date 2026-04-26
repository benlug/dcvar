// Multilevel Copula VAR(1) with random Phi_i, exponential margins, and global rho
// Person-mean centered data, no intercept.

functions {
#include functions/gaussian_copula_uv.stan
}

data {
  int<lower=1> N;
  int<lower=2> n_time;
  array[N] matrix[n_time, 2] y;
  vector[2] skew_direction;

  // Prior hyperparameters
  real<lower=0> prior_phi_bar_sd;
  real<lower=0> prior_tau_phi_scale;
  real<lower=0> prior_sigma_sd;
  real<lower=0> prior_rho_sd;
}

transformed data {
  int n_time_eff = n_time - 1;
}

parameters {
  vector[4] phi_bar;
  vector<lower=0>[4] tau_phi;
  matrix[N, 4] z_phi;
  vector[2] eta;
  real<lower=-1, upper=1> rho;
}

transformed parameters {
  matrix[N, 4] phi_unit;
  array[N] matrix[2, 2] Phi_T;

  for (i in 1:N) {
    for (k in 1:4) {
      phi_unit[i, k] = phi_bar[k] + tau_phi[k] * z_phi[i, k];
    }
    Phi_T[i][1, 1] = phi_unit[i, 1];
    Phi_T[i][2, 1] = phi_unit[i, 2];
    Phi_T[i][1, 2] = phi_unit[i, 3];
    Phi_T[i][2, 2] = phi_unit[i, 4];
  }
}

model {
  vector[2] sigma_lb;
  vector[2] sigma_exp_local;
  vector[2] rate_exp;

  phi_bar ~ normal(0, prior_phi_bar_sd);
  tau_phi ~ student_t(3, 0, prior_tau_phi_scale);
  to_vector(z_phi) ~ std_normal();
  rho ~ normal(0, prior_rho_sd);

  {
    row_vector[2] res0 = y[1][2] - y[1][1] * Phi_T[1];
    for (j in 1:2) {
      sigma_lb[j] = -skew_direction[j] * res0[j];
    }

    for (i in 1:N) {
      int t_start = (i == 1) ? 3 : 2;
      for (t in t_start:n_time) {
        row_vector[2] res_t = y[i][t] - y[i][t - 1] * Phi_T[i];
        for (j in 1:2) {
          sigma_lb[j] = fmax(sigma_lb[j], -skew_direction[j] * res_t[j]);
        }
      }
    }
  }

  for (j in 1:2) {
    sigma_lb[j] = fmax(sigma_lb[j], 0);
    sigma_exp_local[j] = sigma_lb[j] + exp(eta[j]) + 1e-9;
  }
  rate_exp = 1.0 ./ sigma_exp_local;

  for (j in 1:2) {
    target += lognormal_lpdf(sigma_exp_local[j] | 0, prior_sigma_sd) + eta[j];
  }

  for (i in 1:N) {
    for (t in 2:n_time) {
      row_vector[2] pred = y[i][t - 1] * Phi_T[i];
      row_vector[2] res = y[i][t] - pred;
      vector[2] u;

      for (j in 1:2) {
        real x_shifted = sigma_exp_local[j] + skew_direction[j] * res[j];
        target += exponential_lpdf(x_shifted | rate_exp[j]);
        u[j] = exponential_cdf(x_shifted | rate_exp[j]);
        if (skew_direction[j] < 0) {
          u[j] = 1.0 - u[j];
        }
      }

      target += gaussian_copula_uv_lpdf(u | rho);
    }
  }
}

generated quantities {
  vector[2] sigma_exp;
  vector[2] b_gq;
  vector[2] slack_gq;
  vector[2] rate_exp;
  array[N] matrix[2, 2] Phi;
  vector[N] spectral_radius;
  int<lower=0> n_nonstationary = 0;
  matrix[N, n_time_eff] log_lik;

  {
    row_vector[2] res0 = y[1][2] - y[1][1] * Phi_T[1];
    for (j in 1:2) {
      b_gq[j] = -skew_direction[j] * res0[j];
    }

    for (i in 1:N) {
      int t_start = (i == 1) ? 3 : 2;
      for (t in t_start:n_time) {
        row_vector[2] res_t = y[i][t] - y[i][t - 1] * Phi_T[i];
        for (j in 1:2) {
          b_gq[j] = fmax(b_gq[j], -skew_direction[j] * res_t[j]);
        }
      }
    }
  }

  for (j in 1:2) {
    b_gq[j] = fmax(b_gq[j], 0);
    slack_gq[j] = exp(eta[j]);
    sigma_exp[j] = b_gq[j] + slack_gq[j] + 1e-9;
    rate_exp[j] = 1.0 / sigma_exp[j];
  }

  for (i in 1:N) {
    Phi[i] = Phi_T[i]';

    {
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
    }

    if (spectral_radius[i] >= 1.0) {
      n_nonstationary += 1;
    }

    for (t in 2:n_time) {
      row_vector[2] pred = y[i][t - 1] * Phi_T[i];
      row_vector[2] res = y[i][t] - pred;
      vector[2] u;
      log_lik[i, t - 1] = 0;

      for (j in 1:2) {
        real x_shifted = sigma_exp[j] + skew_direction[j] * res[j];
        log_lik[i, t - 1] += exponential_lpdf(x_shifted | rate_exp[j]);
        u[j] = exponential_cdf(x_shifted | rate_exp[j]);
        if (skew_direction[j] < 0) {
          u[j] = 1.0 - u[j];
        }
      }

      log_lik[i, t - 1] += gaussian_copula_uv_lpdf(u | rho);
    }
  }
}
