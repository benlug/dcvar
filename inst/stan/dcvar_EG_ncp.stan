// DC-VAR Model with Exponential-Gaussian Margins (NCP)
// Time-varying Gaussian copula with exponential innovation margins

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

  // Prior hyperparameters
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
  vector[D] eta;                     // unconstrained: sigma_exp = max(0, b) + exp(eta)

  // Time-varying copula (NCP)
  real z_rho_init;
  real<lower=0.001> sigma_omega;
  vector[T_eff] omega_raw;
}

transformed parameters {
  // VAR residuals
  matrix[T_eff, D] eps = compute_var_residuals(Y, mu, Phi, T_eff, D);

  // NCP random walk for rho
  vector[T_eff] z_rho = compute_z_rho_ncp(z_rho_init, sigma_omega, omega_raw, T_eff);
  vector[T_eff] rho;
  for (t in 1:T_eff) {
    rho[t] = tanh(z_rho[t]);
  }
}

model {
  vector[D] sigma_lb;
  vector[D] sigma_exp;
  vector[D] rate_exp;
  real sigma_eps = 1e-9;

  // Priors
  mu ~ normal(0, sigma_mu_prior);
  to_vector(Phi) ~ normal(0, sigma_phi_prior);
  z_rho_init ~ normal(0, rho_init_prior_sd);
  sigma_omega ~ exponential(1.0 / sigma_omega_prior);
  omega_raw ~ std_normal();

  // Feasibility bounds
  for (i in 1:D) {
    real m = -skew_direction[i] * eps[1, i];
    for (t in 2:T_eff) {
      m = fmax(m, -skew_direction[i] * eps[t, i]);
    }
    sigma_lb[i] = fmax(m, 0);
  }

  for (i in 1:D) {
    sigma_exp[i] = sigma_lb[i] + exp(eta[i]) + sigma_eps;
  }
  rate_exp = 1.0 ./ sigma_exp;

  // Induced prior on sigma_exp
  for (i in 1:D) {
    target += lognormal_lpdf(sigma_exp[i] | 0, 0.5) + eta[i];
  }

  // Likelihood
  for (t in 1:T_eff) {
    row_vector[D] res = eps[t];
    vector[2] u_vec;

    for (i in 1:D) {
      real x_shifted = sigma_exp[i] + skew_direction[i] * res[i];
      target += exponential_lpdf(x_shifted | rate_exp[i]);
      u_vec[i] = exponential_cdf(x_shifted | rate_exp[i]);
      if (skew_direction[i] < 0) u_vec[i] = 1.0 - u_vec[i];
    }

    target += gaussian_copula_uv_lpdf(u_vec | rho[t]);
  }
}

generated quantities {
  vector[T_eff] log_lik;
  matrix[T_eff, D] eps_rep;
  vector[D] sigma_exp;
  vector[D] b_gq;
  vector[D] rate_exp;

  {
    real sigma_eps = 1e-9;
    for (i in 1:D) {
      real m = -skew_direction[i] * eps[1, i];
      for (t in 2:T_eff) m = fmax(m, -skew_direction[i] * eps[t, i]);
      b_gq[i] = m;
      sigma_exp[i] = fmax(m, 0) + exp(eta[i]) + sigma_eps;
      rate_exp[i] = 1.0 / sigma_exp[i];
    }

    for (t in 1:T_eff) {
      log_lik[t] = 0;
      vector[2] u_vec;
      for (i in 1:D) {
        real x_shifted = sigma_exp[i] + skew_direction[i] * eps[t, i];
        log_lik[t] += exponential_lpdf(x_shifted | rate_exp[i]);
        u_vec[i] = exponential_cdf(x_shifted | rate_exp[i]);
        if (skew_direction[i] < 0) u_vec[i] = 1.0 - u_vec[i];
      }
      log_lik[t] += gaussian_copula_uv_lpdf(u_vec | rho[t]);

      // PPC: sample from bivariate normal copula, then invert through exponential margins
      {
        real z1_rep = std_normal_rng();
        real z2_rep = rho[t] * z1_rep + sqrt(1 - square(rho[t])) * std_normal_rng();
        real u1 = Phi(z1_rep);
        real u2 = Phi(z2_rep);
        eps_rep[t, 1] = skew_direction[1] * (-log1m(u1) / rate_exp[1] - sigma_exp[1]);
        eps_rep[t, 2] = skew_direction[2] * (-log1m(u2) / rate_exp[2] - sigma_exp[2]);
      }
    }
  }
}
