// Study 2 measures a latent outcome with three ordinal items.

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
  vector[2] mu_v;
  vector<lower=0>[2] sigma_v;
  real alpha_treated;
  vector[2] beta;
  vector<lower=0>[2] sigma_y;
  vector[2] lambda_free;
  ordered[K[1] - 1] threshold_1;
  ordered[K[2] - 1] threshold_2;
  ordered[K[3] - 1] threshold_3;
  vector<lower=0>[3] item_sd_treated;
  vector[N] latent_raw;
}

transformed parameters {
  vector[2] alpha;
  vector[3] lambda;
  matrix[2, 3] item_sd;
  vector[N] latent;
  lambda[1] = 1;
  lambda[2:3] = lambda_free;
  item_sd[1] = rep_row_vector(1, 3);
  item_sd[2] = item_sd_treated';
  // The marginal factor mean is zero in the control group.
  alpha[1] = -beta[1] * mu_v[1];
  alpha[2] = alpha_treated;
  for (i in 1:N) {
    int g = group[i];
    latent[i] = alpha[g] + beta[g] * z[i] + sigma_y[g] * latent_raw[i];
  }
}

model {
  mu_v ~ normal(0, prior_intercept_sd);
  sigma_v ~ normal(0, prior_scale_sd);
  alpha_treated ~ normal(0, prior_intercept_sd);
  beta ~ normal(0, prior_beta_sd);
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
      z[i] ~ normal(mu_v[g], sigma_v[g]);
      target += causal_measurement_lpmf(u[i] | latent[i], lambda,
                    item_sd[g], threshold_1, threshold_2, threshold_3);
    }
  }
}

generated quantities {
  vector[N] latent_rep;
  vector[N] y_rep;
  array[N, 3] int u_rep;
  for (i in 1:N) {
    int g = group[i];
    // Prediction keeps the observed covariate and draws a new outcome factor.
    latent_rep[i] = normal_rng(alpha[g] + beta[g] * z[i], sigma_y[g]);
    y_rep[i] = latent_rep[i];
    u_rep[i] = causal_measurement_rng(latent_rep[i], lambda, item_sd[g],
                                      threshold_1, threshold_2, threshold_3);
  }
}
