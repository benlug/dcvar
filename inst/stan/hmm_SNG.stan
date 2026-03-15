// HMM Copula VAR Model with Skew-Normal-Gaussian Margins

functions {
#include functions/gaussian_copula_uv.stan
}

data {
  int<lower=2> T;
  int<lower=2> D;
  matrix[T, D] Y;
  int<lower=2> K;

  real<lower=0> sigma_mu_prior;
  real<lower=0> sigma_phi_prior;
  real<lower=0> kappa;
  real<lower=0> alpha_off;
  real<lower=0> z_rho_prior_sd;
}

transformed data {
  int T_eff = T - 1;
  real SQRT_2_OVER_PI = sqrt(2.0 / pi());
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
  vector<lower=0>[D] omega;
  vector<lower=-1, upper=1>[D] delta;
  ordered[K] z_rho;
  simplex[K] pi0;
  array[K] simplex[K] A;
}

transformed parameters {
  vector[K] rho_state;
  matrix[T_eff, D] eps;
  vector[D] sn_alpha;  // renamed to avoid conflict with alpha_off
  vector[D] xi;
  matrix[T_eff, K] obs_ll;
  matrix[K, K] log_A;
  vector[K] log_pi0;
  matrix[T_eff, K] log_alpha;

  sn_alpha = delta ./ sqrt(1 - square(delta));
  xi = -omega .* (delta * SQRT_2_OVER_PI);

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

  for (t in 1:T_eff) {
    real marginal_ll = 0;
    vector[2] u_vec;
    for (i in 1:D) {
      marginal_ll += skew_normal_lpdf(eps[t, i] | xi[i], omega[i], sn_alpha[i]);
      u_vec[i] = skew_normal_cdf(eps[t, i] | xi[i], omega[i], sn_alpha[i]);
    }
    for (k in 1:K) {
      obs_ll[t, k] = marginal_ll + gaussian_copula_uv_lpdf(u_vec | rho_state[k]);
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
  omega ~ normal(0, 1);
  delta ~ normal(0, 0.5);

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

  // NOTE: eps_rep contains copula-level z-scores, not skew-normal residuals,
  // because Stan lacks a skew-normal inverse CDF. plot_ppc() rejects
  // skew-normal fits until replicated residuals are available on the
  // margin scale.
  for (t in 1:T_eff) {
    real z1_rep = std_normal_rng();
    real z2_rep = rho_hmm[t] * z1_rep + sqrt(1 - square(rho_hmm[t])) * std_normal_rng();
    eps_rep[t, 1] = z1_rep;
    eps_rep[t, 2] = z2_rep;
  }
}
