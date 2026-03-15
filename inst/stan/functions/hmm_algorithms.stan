// HMM forward algorithm: compute log p(y_{1:t}, s_t = k | theta)
matrix hmm_forward(matrix obs_ll, matrix log_A, vector log_pi0, int T_eff, int K) {
  matrix[T_eff, K] log_alpha;
  for (k in 1:K)
    log_alpha[1, k] = log_pi0[k] + obs_ll[1, k];
  for (t in 2:T_eff) {
    for (k in 1:K) {
      vector[K] log_transition;
      for (j in 1:K)
        log_transition[j] = log_alpha[t - 1, j] + log_A[j, k];
      log_alpha[t, k] = log_sum_exp(log_transition) + obs_ll[t, k];
    }
  }
  return log_alpha;
}

// HMM forward-backward: compute state posteriors gamma[t,k] = p(s_t=k | y_{1:T})
matrix hmm_state_posteriors(matrix log_alpha, matrix obs_ll, matrix log_A, int T_eff, int K) {
  matrix[T_eff, K] gamma;
  matrix[T_eff, K] log_beta;
  for (k in 1:K)
    log_beta[T_eff, k] = 0;
  for (t_rev in 1:(T_eff - 1)) {
    int t = T_eff - t_rev;
    for (k in 1:K) {
      vector[K] log_terms;
      for (j in 1:K)
        log_terms[j] = log_A[k, j] + obs_ll[t + 1, j] + log_beta[t + 1, j];
      log_beta[t, k] = log_sum_exp(log_terms);
    }
  }
  for (t in 1:T_eff) {
    vector[K] log_gamma_t;
    for (k in 1:K)
      log_gamma_t[k] = log_alpha[t, k] + log_beta[t, k];
    real log_norm = log_sum_exp(log_gamma_t);
    for (k in 1:K)
      gamma[t, k] = exp(log_gamma_t[k] - log_norm);
  }
  return gamma;
}

// HMM Viterbi: MAP state sequence
array[] int hmm_viterbi(matrix obs_ll, matrix log_A, vector log_pi0, int T_eff, int K) {
  array[T_eff] int viterbi_state;
  matrix[T_eff, K] log_delta;
  array[T_eff, K] int psi;
  for (k in 1:K) {
    log_delta[1, k] = log_pi0[k] + obs_ll[1, k];
    psi[1, k] = 0;
  }
  for (t in 2:T_eff) {
    for (k in 1:K) {
      real max_val = negative_infinity();
      int max_idx = 1;
      for (j in 1:K) {
        real val = log_delta[t - 1, j] + log_A[j, k];
        if (val > max_val) {
          max_val = val;
          max_idx = j;
        }
      }
      log_delta[t, k] = max_val + obs_ll[t, k];
      psi[t, k] = max_idx;
    }
  }
  {
    real max_val = negative_infinity();
    for (k in 1:K) {
      if (log_delta[T_eff, k] > max_val) {
        max_val = log_delta[T_eff, k];
        viterbi_state[T_eff] = k;
      }
    }
  }
  for (t_rev in 1:(T_eff - 1)) {
    int t = T_eff - t_rev;
    viterbi_state[t] = psi[t + 1, viterbi_state[t + 1]];
  }
  return viterbi_state;
}

// HMM posterior-averaged rho
vector hmm_rho_average(matrix gamma_mat, vector rho_state, int T_eff, int K) {
  vector[T_eff] rho_hmm;
  for (t in 1:T_eff) {
    rho_hmm[t] = 0;
    for (k in 1:K)
      rho_hmm[t] += gamma_mat[t, k] * rho_state[k];
  }
  return rho_hmm;
}

// HMM log-likelihood from forward filter predictive density
vector hmm_log_lik(matrix log_alpha, int T_eff, int K) {
  vector[T_eff] log_lik;
  log_lik[1] = log_sum_exp(to_vector(log_alpha[1, ]));
  for (t in 2:T_eff)
    log_lik[t] = log_sum_exp(to_vector(log_alpha[t, ]))
                 - log_sum_exp(to_vector(log_alpha[t - 1, ]));
  return log_lik;
}
