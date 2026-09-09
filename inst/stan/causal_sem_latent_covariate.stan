// Study 1 measures a latent covariate with three ordinal items.

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
  real mu_v_treated;
  vector<lower=0>[2] sigma_v;
  vector[2] alpha;
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
  vector[2] mu_v;
  vector[3] lambda;
  matrix[2, 3] item_sd;
  vector[N] latent;
  lambda[1] = 1;
  lambda[2:3] = lambda_free;
  item_sd[1] = rep_row_vector(1, 3);
  item_sd[2] = item_sd_treated';
  mu_v[1] = 0;
  mu_v[2] = mu_v_treated;
  for (i in 1:N) {
    int g = group[i];
    if (prior_only == 1) {
      latent[i] = mu_v[g] + sigma_v[g] * latent_raw[i];
    } else {
      real marginal_sd = hypot(beta[g] * sigma_v[g], sigma_y[g]);
      real gain = (sigma_v[g] / marginal_sd)
                    * (beta[g] * sigma_v[g] / marginal_sd);
      real conditional_sd = sigma_v[g] * (sigma_y[g] / marginal_sd);
      latent[i] = mu_v[g] + gain * (y[i] - alpha[g] - beta[g] * mu_v[g])
                    + conditional_sd * latent_raw[i];
    }
  }
}

model {
  mu_v_treated ~ normal(0, prior_intercept_sd);
  sigma_v ~ normal(0, prior_scale_sd);
  alpha ~ normal(0, prior_intercept_sd);
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
      // The marginal outcome and conditional factor give the same joint density.
      y[i] ~ normal(alpha[g] + beta[g] * mu_v[g],
                      hypot(beta[g] * sigma_v[g], sigma_y[g]));
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
    // A new person has a new covariate and a new outcome.
    latent_rep[i] = normal_rng(mu_v[g], sigma_v[g]);
    y_rep[i] = normal_rng(alpha[g] + beta[g] * latent_rep[i], sigma_y[g]);
    u_rep[i] = causal_measurement_rng(latent_rep[i], lambda, item_sd[g],
                                      threshold_1, threshold_2, threshold_3);
  }
}
