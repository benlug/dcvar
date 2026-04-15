// DC-VAR Model with Gamma-Gaussian Margins (NCP)

functions {
#include functions/gaussian_copula_uv.stan
#include functions/var_residuals.stan
#include functions/ncp_random_walk.stan
}

data {
  int<lower=2> T;
  int<lower=2> D;
  matrix[T, D] Y;
  vector[D] skew_direction;

  real<lower=0> sigma_mu_prior;
  real<lower=0> sigma_phi_prior;
  real<lower=0> sigma_omega_prior;
  real<lower=0> rho_init_prior_sd;
}

transformed data {
  int T_eff = T - 1;
}

parameters {
  vector[D] mu;
  matrix[D, D] Phi;
  vector[D] eta;
  real<lower=0> shape_gam;

  real z_rho_init;
  real<lower=0.001> sigma_omega;
  vector[T_eff] omega_raw;
}

transformed parameters {
  matrix[T_eff, D] eps = compute_var_residuals(Y, mu, Phi, T_eff, D);
  vector[T_eff] z_rho = compute_z_rho_ncp(z_rho_init, sigma_omega, omega_raw, T_eff);
  vector[T_eff] rho;
  for (t in 1:T_eff) rho[t] = tanh(z_rho[t]);
}

model {
  vector[D] mean_lb;
  vector[D] mean_gam;
  vector[D] sigma_gam;
  vector[D] rate_gam;
  real sqrt_shape = sqrt(shape_gam);
  real sigma_eps = 1e-9;

  mu ~ normal(0, sigma_mu_prior);
  to_vector(Phi) ~ normal(0, sigma_phi_prior);
  z_rho_init ~ normal(0, rho_init_prior_sd);
  sigma_omega ~ exponential(1.0 / sigma_omega_prior);
  omega_raw ~ std_normal();
  shape_gam ~ lognormal(log(1), 0.5);

  for (i in 1:D) {
    real m = -skew_direction[i] * eps[1, i];
    for (t in 2:T_eff) m = fmax(m, -skew_direction[i] * eps[t, i]);
    mean_lb[i] = fmax(m, 0);
  }

  for (i in 1:D) mean_gam[i] = mean_lb[i] + exp(eta[i]) + sigma_eps;
  sigma_gam = mean_gam / sqrt_shape;
  rate_gam = shape_gam ./ mean_gam;

  for (i in 1:D) target += lognormal_lpdf(mean_gam[i] | 0, 0.5) + eta[i];

  for (t in 1:T_eff) {
    row_vector[D] res = eps[t];
    vector[2] u_vec;
    for (i in 1:D) {
      real x_shifted = mean_gam[i] + skew_direction[i] * res[i];
      target += gamma_lpdf(x_shifted | shape_gam, rate_gam[i]);
      u_vec[i] = gamma_cdf(x_shifted | shape_gam, rate_gam[i]);
      if (skew_direction[i] < 0) u_vec[i] = 1.0 - u_vec[i];
    }
    target += gaussian_copula_uv_lpdf(u_vec | rho[t]);
  }
}

generated quantities {
  vector[T_eff] log_lik;
  matrix[T_eff, D] eps_rep;
  vector[D] sigma_gam;
  vector[D] b_gq;
  vector[D] rate_gam;

  {
    real sqrt_shape = sqrt(shape_gam);
    real sigma_eps = 1e-9;
    vector[D] mean_gam;
    for (i in 1:D) {
      real m = -skew_direction[i] * eps[1, i];
      for (t in 2:T_eff) m = fmax(m, -skew_direction[i] * eps[t, i]);
      b_gq[i] = fmax(m, 0);
      mean_gam[i] = b_gq[i] + exp(eta[i]) + sigma_eps;
      sigma_gam[i] = mean_gam[i] / sqrt_shape;
      rate_gam[i] = shape_gam / mean_gam[i];
    }

    for (t in 1:T_eff) {
      log_lik[t] = 0;
      vector[2] u_vec;
      for (i in 1:D) {
        real x_shifted = mean_gam[i] + skew_direction[i] * eps[t, i];
        log_lik[t] += gamma_lpdf(x_shifted | shape_gam, rate_gam[i]);
        u_vec[i] = gamma_cdf(x_shifted | shape_gam, rate_gam[i]);
        if (skew_direction[i] < 0) u_vec[i] = 1.0 - u_vec[i];
      }
      log_lik[t] += gaussian_copula_uv_lpdf(u_vec | rho[t]);

      // NOTE: eps_rep contains copula-level z-scores, not gamma residuals,
      // because Stan lacks a gamma inverse CDF. plot_ppc() rejects gamma
      // fits until replicated residuals are available on the margin scale.
      {
        real z1_rep = std_normal_rng();
        real z2_rep = rho[t] * z1_rep + sqrt(1 - square(rho[t])) * std_normal_rng();
        eps_rep[t, 1] = z1_rep;
        eps_rep[t, 2] = z2_rep;
      }
    }
  }
}
