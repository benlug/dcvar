// Constant Copula VAR Model with Skew-Normal-Gaussian Margins

functions {
#include functions/gaussian_copula_uv.stan
#include functions/var_residuals.stan
}

data {
  int<lower=2> T;
  int<lower=2> D;
  matrix[T, D] Y;

  real<lower=0> sigma_mu_prior;
  real<lower=0> sigma_phi_prior;
  real<lower=0> z_rho_prior_sd;
}

transformed data {
  int T_eff = T - 1;
  real SQRT_2_OVER_PI = sqrt(2.0 / pi());
}

parameters {
  vector[D] mu;
  matrix[D, D] Phi;
  vector<lower=0>[D] omega;
  vector<lower=-1, upper=1>[D] delta;
  real z_rho;
}

transformed parameters {
  matrix[T_eff, D] eps = compute_var_residuals(Y, mu, Phi, T_eff, D);
  real rho = tanh(z_rho);
  vector[D] alpha;
  vector[D] xi;

  // CP -> DP transform
  alpha = delta ./ sqrt(1 - square(delta));
  xi = -omega .* (delta * SQRT_2_OVER_PI);
}

model {
  mu ~ normal(0, sigma_mu_prior);
  to_vector(Phi) ~ normal(0, sigma_phi_prior);
  z_rho ~ normal(0, z_rho_prior_sd);
  omega ~ normal(0, 1);
  delta ~ normal(0, 0.5);

  for (t in 1:T_eff) {
    row_vector[D] res = eps[t];
    vector[2] u_vec;

    for (i in 1:D) {
      target += skew_normal_lpdf(res[i] | xi[i], omega[i], alpha[i]);
      u_vec[i] = skew_normal_cdf(res[i] | xi[i], omega[i], alpha[i]);
    }
    target += gaussian_copula_uv_lpdf(u_vec | rho);
  }
}

generated quantities {
  vector[T_eff] log_lik;
  matrix[T_eff, D] eps_rep;

  for (t in 1:T_eff) {
    log_lik[t] = 0;
    vector[2] u_vec;
    for (i in 1:D) {
      log_lik[t] += skew_normal_lpdf(eps[t, i] | xi[i], omega[i], alpha[i]);
      u_vec[i] = skew_normal_cdf(eps[t, i] | xi[i], omega[i], alpha[i]);
    }
    log_lik[t] += gaussian_copula_uv_lpdf(u_vec | rho);

    // NOTE: eps_rep contains copula-level z-scores, not skew-normal residuals,
    // because Stan lacks a skew-normal inverse CDF. plot_ppc() rejects
    // skew-normal fits until replicated residuals are available on the
    // margin scale.
    {
      real z1_rep = std_normal_rng();
      real z2_rep = rho * z1_rep + sqrt(1 - square(rho)) * std_normal_rng();
      eps_rep[t, 1] = z1_rep;
      eps_rep[t, 2] = z2_rep;
    }
  }
}
