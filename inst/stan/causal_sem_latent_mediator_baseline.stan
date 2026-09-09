// Study 3B measures the mediator baseline with three ordinal items.

functions {
#include functions/causal_measurement.stan
}

data {
  int<lower=1> N;
  int<lower=3, upper=3> J;
  array[N] int<lower=1, upper=2> group;
  array[3] int<lower=3, upper=5> K;
  array[N, 3] int<lower=1, upper=5> u;
  vector[N] y;
  vector[N] z;
  vector[N] m;
  vector[N] q1;
  vector[N] q2;
  real<lower=0> prior_beta_sd;
  real<lower=0> prior_intercept_sd;
  real<lower=0> prior_scale_sd;
  real prior_loading_mean;
  real<lower=0> prior_loading_sd;
  real<lower=0> prior_threshold_sd;
  real<lower=0> prior_item_log_sd;
  real<lower=0> prior_lkj_shape;
  int<lower=0, upper=1> prior_only;
}

transformed data {
  for (i in 1:N) {
    for (j in 1:3) {
      if (u[i, j] > K[j]) {
        reject("An item response exceeds its number of categories.");
      }
    }
  }
}

parameters {
  real mu_q1_treated;
  vector[2] mu_q2;
  matrix<lower=0>[2, 2] sigma_q;
  array[2] cholesky_factor_corr[2] L_q;
  vector[2] alpha_m;
  matrix[2, 2] B;
  vector[2] alpha_y;
  matrix[2, 2] D;
  vector[2] d;
  vector<lower=0>[2] sigma_m;
  vector<lower=0>[2] sigma_y;
  vector[2] lambda_free;
  ordered[K[1] - 1] threshold_1;
  ordered[K[2] - 1] threshold_2;
  ordered[K[3] - 1] threshold_3;
  vector<lower=0>[3] item_sd_treated;
  vector[N] latent_raw;
}

transformed parameters {
  matrix[2, 2] mu_q;
  vector[2] rho_q;
  vector[3] lambda;
  matrix[2, 3] item_sd;
  vector[N] latent;
  lambda[1] = 1;
  lambda[2:3] = lambda_free;
  item_sd[1] = rep_row_vector(1, 3);
  item_sd[2] = item_sd_treated';
  mu_q[1, 1] = 0;
  mu_q[2, 1] = mu_q1_treated;
  mu_q[, 2] = mu_q2;
  for (g in 1:2) {
    rho_q[g] = L_q[g][2, 1];
  }
  for (i in 1:N) {
    int g = group[i];
    // This conditional form retains the joint normal baseline distribution.
    real conditional_mean = mu_q[g, 1] + rho_q[g] * sigma_q[g, 1]
                            / sigma_q[g, 2] * (q2[i] - mu_q[g, 2]);
    real conditional_sd = sigma_q[g, 1] * L_q[g][2, 2];
    if (prior_only == 0) {
      // Condition the factor on M, then on Y, before the ordinal update.
      real m_mean = alpha_m[g] + B[g, 2] * q2[i] + B[g, 1] * conditional_mean;
      real m_sd = hypot(sigma_m[g], B[g, 1] * conditional_sd);
      real y_mean;
      real y_sd;
      conditional_mean += (conditional_sd / m_sd) * (B[g, 1] * conditional_sd / m_sd)
                          * (m[i] - m_mean);
      conditional_sd *= sigma_m[g] / m_sd;
      y_mean = alpha_y[g] + D[g, 2] * q2[i] + d[g] * m[i]
               + D[g, 1] * conditional_mean;
      y_sd = hypot(sigma_y[g], D[g, 1] * conditional_sd);
      conditional_mean += (conditional_sd / y_sd) * (D[g, 1] * conditional_sd / y_sd)
                          * (y[i] - y_mean);
      conditional_sd *= sigma_y[g] / y_sd;
    }
    latent[i] = conditional_mean + conditional_sd * latent_raw[i];
  }
}

model {
  mu_q1_treated ~ normal(0, prior_intercept_sd);
  mu_q2 ~ normal(0, prior_intercept_sd);
  to_vector(sigma_q) ~ normal(0, prior_scale_sd);
  for (g in 1:2) {
    L_q[g] ~ lkj_corr_cholesky(prior_lkj_shape);
  }
  alpha_m ~ normal(0, prior_intercept_sd);
  to_vector(B) ~ normal(0, prior_beta_sd);
  alpha_y ~ normal(0, prior_intercept_sd);
  to_vector(D) ~ normal(0, prior_beta_sd);
  d ~ normal(0, prior_beta_sd);
  sigma_m ~ normal(0, prior_scale_sd);
  sigma_y ~ normal(0, prior_scale_sd);
  lambda_free ~ normal(prior_loading_mean, prior_loading_sd);
  threshold_1 ~ normal(0, prior_threshold_sd);
  threshold_2 ~ normal(0, prior_threshold_sd);
  threshold_3 ~ normal(0, prior_threshold_sd);
  item_sd_treated ~ lognormal(0, prior_item_log_sd);
  latent_raw ~ std_normal();
  if (prior_only == 0) {
    for (i in 1:N) {
      int g = group[i];
      q2[i] ~ normal(mu_q[g, 2], sigma_q[g, 2]);
      {
        // These marginal densities retain the original joint normal model.
        real factor_mean = mu_q[g, 1] + rho_q[g] * sigma_q[g, 1]
                             / sigma_q[g, 2] * (q2[i] - mu_q[g, 2]);
        real factor_sd = sigma_q[g, 1] * L_q[g][2, 2];
        real m_mean = alpha_m[g] + B[g, 2] * q2[i] + B[g, 1] * factor_mean;
        real m_sd = hypot(sigma_m[g], B[g, 1] * factor_sd);
        real y_mean;
        real y_sd;
        m[i] ~ normal(m_mean, m_sd);
        factor_mean += (factor_sd / m_sd) * (B[g, 1] * factor_sd / m_sd)
                         * (m[i] - m_mean);
        factor_sd *= sigma_m[g] / m_sd;
        y_mean = alpha_y[g] + D[g, 2] * q2[i] + d[g] * m[i]
                   + D[g, 1] * factor_mean;
        y_sd = hypot(sigma_y[g], D[g, 1] * factor_sd);
        y[i] ~ normal(y_mean, y_sd);
      }
      target += causal_measurement_lpmf(u[i] | latent[i], lambda,
                    item_sd[g], threshold_1, threshold_2, threshold_3);
    }
  }
}

generated quantities {
  vector[N] latent_rep;
  vector[N] y_rep;
  array[N, 3] int u_rep;
  vector[N] m_rep;
  for (i in 1:N) {
    int g = group[i];
    // Prediction keeps Ypre and draws a new Mpre, mediator, and outcome.
    real conditional_mean = mu_q[g, 1] + rho_q[g] * sigma_q[g, 1]
                            / sigma_q[g, 2] * (q2[i] - mu_q[g, 2]);
    real conditional_sd = sigma_q[g, 1] * L_q[g][2, 2];
    latent_rep[i] = normal_rng(conditional_mean, conditional_sd);
    m_rep[i] = normal_rng(alpha_m[g] + B[g, 1] * latent_rep[i]
                           + B[g, 2] * q2[i], sigma_m[g]);
    y_rep[i] = normal_rng(alpha_y[g] + D[g, 1] * latent_rep[i]
                           + D[g, 2] * q2[i] + d[g] * m_rep[i], sigma_y[g]);
    u_rep[i] = causal_measurement_rng(latent_rep[i], lambda, item_sd[g],
                                      threshold_1, threshold_2, threshold_3);
  }
}
