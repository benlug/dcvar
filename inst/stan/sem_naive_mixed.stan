// Naive SEM score VAR(1) with per-variable (mixed) margins and a Gaussian copula
// Observed y is the T x 2 matrix of row-mean factor scores; the VAR(1) residuals
// carry their own marginal family per dimension. rho = 0.97 * tanh(rho_raw).
//
// Family codes (per dimension): 1 = normal, 2 = exponential, 3 = skew_normal,
// 4 = gamma. All per-family parameters are declared unconditionally and a
// dimension that does not use a parameter samples it from its (proper) prior.

functions {
#include functions/gaussian_copula_uv.stan
}

data {
  int<lower=1> n_time;
  matrix[n_time, 2] y;                          // row-mean factor scores
  array[2] int<lower=1, upper=4> family;        // per-dimension margin family code
  vector[2] skew_direction;                     // consulted only for exp/gamma dims

  real<lower=0> prior_mu_sd;
  real<lower=0> prior_phi_sd;
  real<lower=0> prior_sigma_sd;
  real<lower=0> prior_rho_sd;
}

transformed data {
  real SQRT_2_OVER_PI = sqrt(2.0 / pi());
}

parameters {
  vector[2] mu;
  real<lower=-0.99, upper=0.99> phi11;
  real<lower=-0.99, upper=0.99> phi12;
  real<lower=-0.99, upper=0.99> phi21;
  real<lower=-0.99, upper=0.99> phi22;
  real rho_raw;                                 // rho = 0.97 * tanh(rho_raw)

  // Union of per-dimension marginal parameters
  vector<lower=0>[2] sigma_eps;
  vector[2] eta;
  vector<lower=0>[2] omega;
  vector<lower=-1, upper=1>[2] delta;
  vector<lower=0>[2] shape_gam;
}

transformed parameters {
  matrix[2, 2] Phi_T;
  matrix[n_time, 2] eps;
  real rho = 0.97 * tanh(rho_raw);

  vector[2] sigma_exp = rep_vector(0, 2);
  vector[2] rate_exp = rep_vector(0, 2);
  vector[2] mean_gam = rep_vector(0, 2);
  vector[2] rate_gam = rep_vector(0, 2);
  vector[2] sigma_gam = rep_vector(0, 2);
  vector[2] alpha = delta ./ sqrt(1 - square(delta));
  vector[2] xi = -omega .* (delta * SQRT_2_OVER_PI);

  Phi_T[1, 1] = phi11; Phi_T[2, 1] = phi12;
  Phi_T[1, 2] = phi21; Phi_T[2, 2] = phi22;

  {
    matrix[n_time, 2] ylag = rep_matrix(0.0, n_time, 2);
    if (n_time > 1) {
      ylag[2:n_time, ] = y[1:(n_time - 1), ];
    }
    eps = y - (rep_matrix(mu', n_time) + ylag * Phi_T);
  }

  for (i in 1:2) {
    if (family[i] == 2 || family[i] == 4) {
      real m = -skew_direction[i] * eps[1, i];
      for (t in 2:n_time) m = fmax(m, -skew_direction[i] * eps[t, i]);
      real lb = fmax(m, 0);
      if (family[i] == 2) {
        sigma_exp[i] = lb + exp(eta[i]) + 1e-9;
        rate_exp[i] = 1.0 / sigma_exp[i];
      } else {
        mean_gam[i] = lb + exp(eta[i]) + 1e-9;
        rate_gam[i] = shape_gam[i] / mean_gam[i];
        sigma_gam[i] = mean_gam[i] / sqrt(shape_gam[i]);
      }
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

  sigma_eps ~ lognormal(0, prior_sigma_sd);
  omega ~ normal(0, 1);
  delta ~ normal(0, 0.5);
  shape_gam ~ lognormal(log(1), 0.5);

  for (i in 1:2) {
    if (family[i] == 2) {
      target += lognormal_lpdf(sigma_exp[i] | 0, prior_sigma_sd) + eta[i];
    } else if (family[i] == 4) {
      target += lognormal_lpdf(mean_gam[i] | 0, prior_sigma_sd) + eta[i];
    } else {
      eta[i] ~ std_normal();
    }
  }

  for (t in 1:n_time) {
    vector[2] u;
    for (i in 1:2) {
      if (family[i] == 1) {
        target += normal_lpdf(eps[t, i] | 0, sigma_eps[i]);
        u[i] = Phi(eps[t, i] / sigma_eps[i]);
      } else if (family[i] == 2) {
        real x_shifted = sigma_exp[i] + skew_direction[i] * eps[t, i];
        target += exponential_lpdf(x_shifted | rate_exp[i]);
        u[i] = exponential_cdf(x_shifted | rate_exp[i]);
        if (skew_direction[i] < 0) u[i] = 1.0 - u[i];
      } else if (family[i] == 3) {
        target += skew_normal_lpdf(eps[t, i] | xi[i], omega[i], alpha[i]);
        u[i] = skew_normal_cdf(eps[t, i] | xi[i], omega[i], alpha[i]);
      } else {
        real x_shifted = mean_gam[i] + skew_direction[i] * eps[t, i];
        target += gamma_lpdf(x_shifted | shape_gam[i], rate_gam[i]);
        u[i] = gamma_cdf(x_shifted | shape_gam[i], rate_gam[i]);
        if (skew_direction[i] < 0) u[i] = 1.0 - u[i];
      }
    }
    target += gaussian_copula_uv_lpdf(u | rho);
  }
}

generated quantities {
  matrix[2, 2] Phi;
  vector[n_time] log_lik;

  Phi[1, 1] = phi11; Phi[1, 2] = phi12;
  Phi[2, 1] = phi21; Phi[2, 2] = phi22;

  for (t in 1:n_time) {
    vector[2] u;
    log_lik[t] = 0;
    for (i in 1:2) {
      if (family[i] == 1) {
        log_lik[t] += normal_lpdf(eps[t, i] | 0, sigma_eps[i]);
        u[i] = Phi(eps[t, i] / sigma_eps[i]);
      } else if (family[i] == 2) {
        real x_shifted = sigma_exp[i] + skew_direction[i] * eps[t, i];
        log_lik[t] += exponential_lpdf(x_shifted | rate_exp[i]);
        u[i] = exponential_cdf(x_shifted | rate_exp[i]);
        if (skew_direction[i] < 0) u[i] = 1.0 - u[i];
      } else if (family[i] == 3) {
        log_lik[t] += skew_normal_lpdf(eps[t, i] | xi[i], omega[i], alpha[i]);
        u[i] = skew_normal_cdf(eps[t, i] | xi[i], omega[i], alpha[i]);
      } else {
        real x_shifted = mean_gam[i] + skew_direction[i] * eps[t, i];
        log_lik[t] += gamma_lpdf(x_shifted | shape_gam[i], rate_gam[i]);
        u[i] = gamma_cdf(x_shifted | shape_gam[i], rate_gam[i]);
        if (skew_direction[i] < 0) u[i] = 1.0 - u[i];
      }
    }
    log_lik[t] += gaussian_copula_uv_lpdf(u | rho);
  }
}
