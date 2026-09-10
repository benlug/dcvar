// Study 3A measures the mediator with three ordinal items.

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
  matrix[2, 2] mu_q;
  matrix<lower=0>[2, 2] sigma_q;
  array[2] cholesky_factor_corr[2] L_q;
  real alpha_m_treated;
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
  vector[2] rho_q;
  array[2] matrix[2, 2] L_cov_q;
  vector[2] alpha_m;
  vector[3] lambda;
  matrix[2, 3] item_sd;
  vector[N] latent;
  lambda[1] = 1;
  lambda[2:3] = lambda_free;
  item_sd[1] = rep_row_vector(1, 3);
  item_sd[2] = item_sd_treated';
  for (g in 1:2) {
    rho_q[g] = L_q[g][2, 1];
    L_cov_q[g] = diag_pre_multiply(to_vector(sigma_q[g]), L_q[g]);
  }
  // The marginal mediator mean is zero in the control group.
  alpha_m[1] = -dot_product(B[1], mu_q[1]);
  alpha_m[2] = alpha_m_treated;
  for (i in 1:N) {
    int g = group[i];
    real mean_m = alpha_m[g] + B[g, 1] * q1[i] + B[g, 2] * q2[i];
    if (prior_only == 1) {
      latent[i] = mean_m + sigma_m[g] * latent_raw[i];
    } else {
      real mean_y_base = alpha_y[g] + D[g, 1] * q1[i] + D[g, 2] * q2[i];
      real marginal_sd_y = hypot(sigma_y[g], d[g] * sigma_m[g]);
      real gain = (sigma_m[g] / marginal_sd_y)
                  * (d[g] * sigma_m[g] / marginal_sd_y);
      real conditional_mean = mean_m + gain * (y[i] - mean_y_base - d[g] * mean_m);
      real conditional_sd = sigma_m[g] * (sigma_y[g] / marginal_sd_y);
      latent[i] = conditional_mean + conditional_sd * latent_raw[i];
    }
  }
}

model {
  to_vector(mu_q) ~ normal(0, prior_intercept_sd);
  to_vector(sigma_q) ~ normal(0, prior_scale_sd);
  for (g in 1:2) {
    L_q[g] ~ lkj_corr_cholesky(prior_lkj_shape);
  }
  alpha_m_treated ~ normal(0, prior_intercept_sd);
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
      vector[2] q_i;
      q_i[1] = q1[i];
      q_i[2] = q2[i];
      q_i ~ multi_normal_cholesky(to_vector(mu_q[g]), L_cov_q[g]);
      y[i] ~ normal(alpha_y[g] + D[g, 1] * q1[i] + D[g, 2] * q2[i]
                    + d[g] * (alpha_m[g] + B[g, 1] * q1[i] + B[g, 2] * q2[i]),
                    hypot(sigma_y[g], d[g] * sigma_m[g]));
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
    // Prediction keeps both observed baselines and draws a new mediator.
    latent_rep[i] = normal_rng(alpha_m[g] + B[g, 1] * q1[i]
                                + B[g, 2] * q2[i], sigma_m[g]);
    m_rep[i] = latent_rep[i];
    y_rep[i] = normal_rng(alpha_y[g] + D[g, 1] * q1[i] + D[g, 2] * q2[i]
                           + d[g] * m_rep[i], sigma_y[g]);
    u_rep[i] = causal_measurement_rng(latent_rep[i], lambda, item_sd[g],
                                      threshold_1, threshold_2, threshold_3);
  }
}
