// HMM Copula VAR Model - Regime-Switching Dependence
// K discrete hidden states with fixed correlation per state
// Forward algorithm to marginalize discrete states
// Same VAR(1) + Gaussian copula observation model as DC-VAR
//
// Author: Research Implementation
// Date: 2026

data {
  int<lower=2> T;                    // Number of time points
  int<lower=2> D;                    // Number of variables (typically 2)
  matrix[T, D] Y;                    // Observed data (T x D)
  int<lower=2> K;                    // Number of hidden states

  // Prior hyperparameters (VAR)
  real<lower=0> sigma_mu_prior;      // Prior SD for intercepts
  real<lower=0> sigma_phi_prior;     // Prior SD for VAR coefficients
  real<lower=0> sigma_eps_prior;     // Prior SD for innovation SDs

  // Prior hyperparameters (HMM)
  real<lower=0> kappa;               // Sticky Dirichlet: self-transition concentration
  real<lower=0> alpha_off;           // Sticky Dirichlet: off-diagonal concentration
  real<lower=0> z_rho_prior_sd;      // Prior SD for state-specific z_rho
}

transformed data {
  int T_eff = T - 1;

  // Build Dirichlet hyperparameter vectors for each row of transition matrix
  array[K] vector[K] dirichlet_prior;
  for (k in 1:K) {
    for (j in 1:K) {
      if (k == j) {
        dirichlet_prior[k][j] = kappa;
      } else {
        dirichlet_prior[k][j] = alpha_off;
      }
    }
  }
}

parameters {
  // VAR parameters (constant over time, same as DC-VAR)
  vector[D] mu;                      // Intercepts
  matrix[D, D] Phi;                  // VAR(1) coefficient matrix (unconstrained: stationarity
                                     // is not enforced, allowing unit roots / near-unit-root
                                     // dynamics. The normal prior provides soft regularization.)
  vector<lower=0.01>[D] sigma_eps;   // Marginal innovation SDs

  // HMM parameters
  // ordered[K] enforces z_rho[1] < z_rho[2] < ... < z_rho[K], which resolves
  // label switching: states are identified by monotonically increasing rho.
  // State 1 always has the lowest correlation, State K the highest.
  ordered[K] z_rho;                  // Fisher-z of state-specific rho (ordered for identification)
  simplex[K] pi0;                    // Initial state distribution
  array[K] simplex[K] A;              // Transition matrix: A[from] is a simplex over target states
}

transformed parameters {
  // State-specific correlations
  vector[K] rho_state;

  // VAR residuals
  matrix[T_eff, D] eps;

  // Standardized residuals as z-scores (used directly in copula density)
  matrix[T_eff, D] z_scores;

  // Observation log-likelihoods per state
  matrix[T_eff, K] obs_ll;

  // Log transition matrix and initial probs (precomputed once per iteration)
  matrix[K, K] log_A;
  vector[K] log_pi0;

  // Forward algorithm: log_alpha[t, k] = log p(y_{1:t}, s_t = k | theta)
  matrix[T_eff, K] log_alpha;

  // Precomputed state-specific quantities (for efficiency)
  vector[K] rho_sq_state;
  vector[K] log_one_m_rho_sq;

  // Transform z_rho to correlation scale via tanh
  // Maps to (-1, 1) smoothly; normal prior on z_rho provides soft regularization
  for (k in 1:K) {
    rho_state[k] = tanh(z_rho[k]);
  }

  // Precompute state-specific quantities (avoid recomputation in time loop)
  for (k in 1:K) {
    rho_sq_state[k] = square(rho_state[k]);
    log_one_m_rho_sq[k] = log1m(rho_sq_state[k]);
  }

  // Precompute log transition matrix and log initial probs
  for (j in 1:K) {
    log_pi0[j] = log(pi0[j]);
    for (k in 1:K) {
      log_A[j, k] = log(A[j][k]);
    }
  }

  // Compute VAR residuals (same as DC-VAR)
  for (t in 1:T_eff) {
    vector[D] y_prev = to_vector(Y[t, ]);
    vector[D] y_curr = to_vector(Y[t + 1, ]);
    vector[D] y_hat = mu + Phi * (y_prev - mu);
    eps[t, ] = to_row_vector(y_curr - y_hat);
  }

  // Compute z-scores directly from standardized residuals
  // (avoids the Phi -> clamp -> inv_Phi roundtrip that truncates tail information)
  for (t in 1:T_eff) {
    for (d in 1:D) {
      z_scores[t, d] = eps[t, d] / sigma_eps[d];
    }
  }

  // Compute observation log-likelihoods with inlined copula
  for (t in 1:T_eff) {
    // Marginal density (same for all states)
    real marginal_ll = 0;
    for (d in 1:D) {
      marginal_ll += normal_lpdf(eps[t, d] | 0, sigma_eps[d]);
    }

    // Precomputed z-scores and their products (shared across states)
    {
      real z1 = z_scores[t, 1];
      real z2 = z_scores[t, 2];
      real z1_sq = square(z1);
      real z2_sq = square(z2);
      real z1z2 = z1 * z2;
      real z_sq_sum = z1_sq + z2_sq;

      for (k in 1:K) {
        // Use precomputed quantities for efficiency
        real rho_k = rho_state[k];
        real rho_k_sq = rho_sq_state[k];
        real one_m_rho_sq = 1.0 - rho_k_sq;

        real log_cop = -0.5 * log_one_m_rho_sq[k]
                       - (rho_k_sq * z_sq_sum - 2.0 * rho_k * z1z2)
                         / (2.0 * one_m_rho_sq);

        obs_ll[t, k] = marginal_ll + log_cop;
      }
    }
  }

  // Forward algorithm using precomputed log_A and log_pi0
  for (k in 1:K) {
    log_alpha[1, k] = log_pi0[k] + obs_ll[1, k];
  }

  for (t in 2:T_eff) {
    for (k in 1:K) {
      vector[K] log_transition;
      for (j in 1:K) {
        log_transition[j] = log_alpha[t - 1, j] + log_A[j, k];
      }
      log_alpha[t, k] = log_sum_exp(log_transition) + obs_ll[t, k];
    }
  }
}

model {
  // Priors: VAR parameters (same as DC-VAR)
  mu ~ normal(0, sigma_mu_prior);
  to_vector(Phi) ~ normal(0, sigma_phi_prior);
  sigma_eps ~ exponential(1.0 / sigma_eps_prior);

  // Priors: HMM parameters
  z_rho ~ normal(0, z_rho_prior_sd);
  pi0 ~ dirichlet(rep_vector(1.0, K));
  for (k in 1:K) {
    A[k] ~ dirichlet(dirichlet_prior[k]);
  }

  // Marginal log-likelihood via forward algorithm
  target += log_sum_exp(to_vector(log_alpha[T_eff, ]));
}

generated quantities {
  // State posterior probabilities: gamma[t, k] = p(s_t = k | y_{1:T}, theta)
  matrix[T_eff, K] gamma;

  // Viterbi MAP state sequence
  array[T_eff] int viterbi_state;

  // Posterior-averaged time-varying rho
  vector[T_eff] rho_hmm;

  // Pointwise log-likelihood (for LOO-CV)
  vector[T_eff] log_lik;

  // Posterior predictive replicated residuals (for PPC)
  matrix[T_eff, D] eps_rep;

  // Forward-backward algorithm for state posteriors
  {
    matrix[T_eff, K] log_beta;

    for (k in 1:K) {
      log_beta[T_eff, k] = 0;
    }

    for (t_rev in 1:(T_eff - 1)) {
      int t = T_eff - t_rev;
      for (k in 1:K) {
        vector[K] log_terms;
        for (j in 1:K) {
          log_terms[j] = log_A[k, j] + obs_ll[t + 1, j] + log_beta[t + 1, j];
        }
        log_beta[t, k] = log_sum_exp(log_terms);
      }
    }

    for (t in 1:T_eff) {
      vector[K] log_gamma_t;
      for (k in 1:K) {
        log_gamma_t[k] = log_alpha[t, k] + log_beta[t, k];
      }
      {
        real log_norm = log_sum_exp(log_gamma_t);
        for (k in 1:K) {
          gamma[t, k] = exp(log_gamma_t[k] - log_norm);
        }
      }
    }
  }

  // Viterbi algorithm
  {
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
  }

  // Posterior-averaged rho
  for (t in 1:T_eff) {
    rho_hmm[t] = 0;
    for (k in 1:K) {
      rho_hmm[t] += gamma[t, k] * rho_state[k];
    }
  }

  // Pointwise log-likelihood via forward filter predictive density:
  // log p(y_t | y_{1:t-1}) = log p(y_{1:t}) - log p(y_{1:t-1})
  //                        = log_sum_exp(log_alpha[t,:]) - log_sum_exp(log_alpha[t-1,:])
  // This avoids the double-counting bug of using smoothed gamma (which conditions on y_t).
  log_lik[1] = log_sum_exp(to_vector(log_alpha[1, ]));
  for (t in 2:T_eff) {
    log_lik[t] = log_sum_exp(to_vector(log_alpha[t, ]))
                 - log_sum_exp(to_vector(log_alpha[t - 1, ]));
  }

  // Generate replicated residuals using posterior-averaged rho
  for (t in 1:T_eff) {
    real z1_rep = std_normal_rng();
    real z2_rep = rho_hmm[t] * z1_rep + sqrt(1 - square(rho_hmm[t])) * std_normal_rng();
    eps_rep[t, 1] = z1_rep * sigma_eps[1];
    eps_rep[t, 2] = z2_rep * sigma_eps[2];
  }
}
