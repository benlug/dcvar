// SEM Copula VAR(1) with a fixed measurement model and per-variable (mixed) margins
// Latent innovations (zeta) carry their own marginal family per dimension; the
// latent state recursion and the Normal measurement model are unchanged from
// sem_copula_var.stan / sem_EG.stan. rho = 0.97 * tanh(rho_raw).
//
// Family codes (per dimension): 1 = normal, 2 = exponential, 3 = skew_normal,
// 4 = gamma. All per-family parameters are declared unconditionally and a
// dimension that does not use a parameter samples it from its (proper) prior.

functions {
#include functions/gaussian_copula_uv.stan
}

data {
  int<lower=1> n_time;
  int<lower=1> J;
  matrix[n_time, 2 * J] y;                     // indicators: y11..y1J y21..y2J
  vector[J] lambda;                            // fixed factor loadings
  real<lower=0> sigma_e;                        // fixed measurement error SD
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
  matrix[n_time, 2] zeta;                        // latent innovations

  // Union of per-dimension marginal parameters
  vector<lower=0>[2] sigma_eps;                 // normal innovation SDs
  vector[2] eta;                                // exponential / gamma scale
  vector<lower=0>[2] omega;                     // skew_normal scale
  vector<lower=-1, upper=1>[2] delta;           // skew_normal CP skewness
  vector<lower=0>[2] shape_gam;                 // gamma shape (per dimension)
}

transformed parameters {
  matrix[2, 2] B;
  matrix[n_time, 2] state;
  real rho = 0.97 * tanh(rho_raw);

  vector[2] sigma_exp = rep_vector(0, 2);
  vector[2] rate_exp = rep_vector(0, 2);
  vector[2] mean_gam = rep_vector(0, 2);
  vector[2] rate_gam = rep_vector(0, 2);
  vector[2] sigma_gam = rep_vector(0, 2);
  vector[2] alpha = delta ./ sqrt(1 - square(delta));
  vector[2] xi = -omega .* (delta * SQRT_2_OVER_PI);

  B[1, 1] = phi11; B[1, 2] = phi12;
  B[2, 1] = phi21; B[2, 2] = phi22;

  // Feasibility bounds and derived scales (per-dimension on the latent zeta)
  for (i in 1:2) {
    if (family[i] == 2 || family[i] == 4) {
      real m = -skew_direction[i] * zeta[1, i];
      for (t in 2:n_time) m = fmax(m, -skew_direction[i] * zeta[t, i]);
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

  // Latent state recursion with x_0 = 0
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

  // Union-parameter priors kept proper unconditionally
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

  // Latent innovation density: per-dimension margins + Gaussian copula
  for (t in 1:n_time) {
    vector[2] u_vec;
    for (i in 1:2) {
      if (family[i] == 1) {
        target += normal_lpdf(zeta[t, i] | 0, sigma_eps[i]);
        u_vec[i] = Phi(zeta[t, i] / sigma_eps[i]);
      } else if (family[i] == 2) {
        real x_shifted = sigma_exp[i] + skew_direction[i] * zeta[t, i];
        target += exponential_lpdf(x_shifted | rate_exp[i]);
        u_vec[i] = exponential_cdf(x_shifted | rate_exp[i]);
        if (skew_direction[i] < 0) u_vec[i] = 1.0 - u_vec[i];
      } else if (family[i] == 3) {
        target += skew_normal_lpdf(zeta[t, i] | xi[i], omega[i], alpha[i]);
        u_vec[i] = skew_normal_cdf(zeta[t, i] | xi[i], omega[i], alpha[i]);
      } else {
        real x_shifted = mean_gam[i] + skew_direction[i] * zeta[t, i];
        target += gamma_lpdf(x_shifted | shape_gam[i], rate_gam[i]);
        u_vec[i] = gamma_cdf(x_shifted | shape_gam[i], rate_gam[i]);
        if (skew_direction[i] < 0) u_vec[i] = 1.0 - u_vec[i];
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
      if (family[i] == 1) {
        log_lik[t] += normal_lpdf(zeta[t, i] | 0, sigma_eps[i]);
        u_vec[i] = Phi(zeta[t, i] / sigma_eps[i]);
      } else if (family[i] == 2) {
        real x_shifted = sigma_exp[i] + skew_direction[i] * zeta[t, i];
        log_lik[t] += exponential_lpdf(x_shifted | rate_exp[i]);
        u_vec[i] = exponential_cdf(x_shifted | rate_exp[i]);
        if (skew_direction[i] < 0) u_vec[i] = 1.0 - u_vec[i];
      } else if (family[i] == 3) {
        log_lik[t] += skew_normal_lpdf(zeta[t, i] | xi[i], omega[i], alpha[i]);
        u_vec[i] = skew_normal_cdf(zeta[t, i] | xi[i], omega[i], alpha[i]);
      } else {
        real x_shifted = mean_gam[i] + skew_direction[i] * zeta[t, i];
        log_lik[t] += gamma_lpdf(x_shifted | shape_gam[i], rate_gam[i]);
        u_vec[i] = gamma_cdf(x_shifted | shape_gam[i], rate_gam[i]);
        if (skew_direction[i] < 0) u_vec[i] = 1.0 - u_vec[i];
      }
    }
    log_lik[t] += gaussian_copula_uv_lpdf(u_vec | rho);

    for (j in 1:J) {
      log_lik[t] += normal_lpdf(y[t, j] | lambda[j] * state[t, 1], sigma_e);
      log_lik[t] += normal_lpdf(y[t, J + j] | lambda[j] * state[t, 2], sigma_e);
    }
  }
}
