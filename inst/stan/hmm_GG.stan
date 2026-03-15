// HMM Copula VAR Model with Gamma-Gaussian Margins

functions {
#include functions/gaussian_copula_uv.stan
}

data {
  int<lower=2> T;
  int<lower=2> D;
  matrix[T, D] Y;
  int<lower=2> K;
  vector[D] skew_direction;

  real<lower=0> sigma_mu_prior;
  real<lower=0> sigma_phi_prior;
  real<lower=0> kappa;
  real<lower=0> alpha_off;
  real<lower=0> z_rho_prior_sd;
}

transformed data {
  int T_eff = T - 1;
  array[K] vector[K] dirichlet_prior;
  for (k in 1:K) {
    for (j in 1:K) {
      dirichlet_prior[k][j] = (k == j) ? kappa : alpha_off;
    }
  }
}

parameters {
  vector[D] mu;
  matrix[D, D] Phi;
  vector[D] eta;
  real<lower=0> shape_gam;
  ordered[K] z_rho;
  simplex[K] pi0;
  array[K] simplex[K] A;
}

transformed parameters {
  vector[K] rho_state;
  matrix[T_eff, D] eps;
  matrix[T_eff, K] obs_ll;
  matrix[K, K] log_A;
  vector[K] log_pi0;
  matrix[T_eff, K] log_alpha;

  for (k in 1:K) rho_state[k] = tanh(z_rho[k]);

  for (j in 1:K) {
    log_pi0[j] = log(pi0[j]);
    for (k in 1:K) log_A[j, k] = log(A[j][k]);
  }

  for (t in 1:T_eff) {
    vector[D] y_prev = to_vector(Y[t, ]);
    vector[D] y_curr = to_vector(Y[t + 1, ]);
    vector[D] y_hat = mu + Phi * (y_prev - mu);
    eps[t, ] = to_row_vector(y_curr - y_hat);
  }

  {
    real sqrt_shape = sqrt(shape_gam);
    real sigma_eps = 1e-9;
    vector[D] sigma_lb;
    vector[D] sigma_gam;
    vector[D] rate_gam;

    for (i in 1:D) {
      real m = -skew_direction[i] * eps[1, i];
      for (t in 2:T_eff) m = fmax(m, -skew_direction[i] * eps[t, i]);
      sigma_lb[i] = fmax(m / sqrt_shape, 0);
      sigma_gam[i] = sigma_lb[i] + exp(eta[i]) + sigma_eps;
    }
    rate_gam = sqrt_shape ./ sigma_gam;

    for (t in 1:T_eff) {
      real marginal_ll = 0;
      vector[2] u_vec;
      for (i in 1:D) {
        real mean_x = sqrt_shape * sigma_gam[i];
        real x_shifted = mean_x + skew_direction[i] * eps[t, i];
        marginal_ll += gamma_lpdf(x_shifted | shape_gam, rate_gam[i]);
        u_vec[i] = gamma_cdf(x_shifted | shape_gam, rate_gam[i]);
        if (skew_direction[i] < 0) u_vec[i] = 1.0 - u_vec[i];
      }
      for (k in 1:K) {
        obs_ll[t, k] = marginal_ll + gaussian_copula_uv_lpdf(u_vec | rho_state[k]);
      }
    }
  }

  for (k in 1:K) log_alpha[1, k] = log_pi0[k] + obs_ll[1, k];
  for (t in 2:T_eff) {
    for (k in 1:K) {
      vector[K] log_transition;
      for (j in 1:K) log_transition[j] = log_alpha[t - 1, j] + log_A[j, k];
      log_alpha[t, k] = log_sum_exp(log_transition) + obs_ll[t, k];
    }
  }
}

model {
  mu ~ normal(0, sigma_mu_prior);
  to_vector(Phi) ~ normal(0, sigma_phi_prior);
  z_rho ~ normal(0, z_rho_prior_sd);
  pi0 ~ dirichlet(rep_vector(1.0, K));
  for (k in 1:K) A[k] ~ dirichlet(dirichlet_prior[k]);
  shape_gam ~ lognormal(log(1), 0.5);

  // Induced prior on sigma_gam
  {
    real sqrt_shape = sqrt(shape_gam);
    real sigma_eps = 1e-9;
    vector[D] sigma_lb;
    vector[D] sigma_gam;
    for (i in 1:D) {
      real m = -skew_direction[i] * eps[1, i];
      for (t in 2:T_eff) m = fmax(m, -skew_direction[i] * eps[t, i]);
      sigma_lb[i] = fmax(m / sqrt_shape, 0);
      sigma_gam[i] = sigma_lb[i] + exp(eta[i]) + sigma_eps;
      target += lognormal_lpdf(sigma_gam[i] | 0, 0.5) + eta[i];
    }
  }

  target += log_sum_exp(to_vector(log_alpha[T_eff, ]));
}

generated quantities {
  matrix[T_eff, K] gamma;
  array[T_eff] int viterbi_state;
  vector[T_eff] rho_hmm;
  vector[T_eff] log_lik;
  matrix[T_eff, D] eps_rep;

  {
    matrix[T_eff, K] log_beta;
    for (k in 1:K) log_beta[T_eff, k] = 0;
    for (t_rev in 1:(T_eff - 1)) {
      int t = T_eff - t_rev;
      for (k in 1:K) {
        vector[K] log_terms;
        for (j in 1:K) log_terms[j] = log_A[k, j] + obs_ll[t + 1, j] + log_beta[t + 1, j];
        log_beta[t, k] = log_sum_exp(log_terms);
      }
    }
    for (t in 1:T_eff) {
      vector[K] log_gamma_t;
      for (k in 1:K) log_gamma_t[k] = log_alpha[t, k] + log_beta[t, k];
      real log_norm = log_sum_exp(log_gamma_t);
      for (k in 1:K) gamma[t, k] = exp(log_gamma_t[k] - log_norm);
    }
  }

  {
    matrix[T_eff, K] log_delta;
    array[T_eff, K] int psi;
    for (k in 1:K) { log_delta[1, k] = log_pi0[k] + obs_ll[1, k]; psi[1, k] = 0; }
    for (t in 2:T_eff) {
      for (k in 1:K) {
        real max_val = negative_infinity(); int max_idx = 1;
        for (j in 1:K) {
          real val = log_delta[t - 1, j] + log_A[j, k];
          if (val > max_val) { max_val = val; max_idx = j; }
        }
        log_delta[t, k] = max_val + obs_ll[t, k]; psi[t, k] = max_idx;
      }
    }
    { real max_val = negative_infinity();
      for (k in 1:K) if (log_delta[T_eff, k] > max_val) { max_val = log_delta[T_eff, k]; viterbi_state[T_eff] = k; }
    }
    for (t_rev in 1:(T_eff - 1)) { int t = T_eff - t_rev; viterbi_state[t] = psi[t + 1, viterbi_state[t + 1]]; }
  }

  for (t in 1:T_eff) {
    rho_hmm[t] = 0;
    for (k in 1:K) rho_hmm[t] += gamma[t, k] * rho_state[k];
  }

  log_lik[1] = log_sum_exp(to_vector(log_alpha[1, ]));
  for (t in 2:T_eff) log_lik[t] = log_sum_exp(to_vector(log_alpha[t, ])) - log_sum_exp(to_vector(log_alpha[t - 1, ]));

  // NOTE: eps_rep contains copula-level z-scores, not gamma residuals,
  // because Stan lacks a gamma inverse CDF. plot_ppc() rejects gamma
  // fits until replicated residuals are available on the margin scale.
  for (t in 1:T_eff) {
    real z1_rep = std_normal_rng();
    real z2_rep = rho_hmm[t] * z1_rep + sqrt(1 - square(rho_hmm[t])) * std_normal_rng();
    eps_rep[t, 1] = z1_rep;
    eps_rep[t, 2] = z2_rep;
  }
}
