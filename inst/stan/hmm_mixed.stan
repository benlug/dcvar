// HMM Copula VAR Model with per-variable (mixed) margins
// K discrete hidden states with state-specific correlation, where each
// dimension can use its own marginal family. This is the regime-switching
// counterpart of constant_mixed.stan: the per-dimension marginal dispatch is
// identical, only the copula correlation is state-dependent and marginalised
// over states with the forward algorithm. The R routing layer dispatches here
// whenever the requested per-dimension margins differ; same-family requests
// keep routing to the specialised models (hmm_copula_model / hmm_*).
//
// Family codes (per dimension): 1 = normal, 2 = exponential, 3 = skew_normal,
// 4 = gamma. All per-family parameters are declared unconditionally; on a
// dimension whose family does not use a parameter it samples from its (proper)
// prior. The Gaussian copula is applied via the CDF (u,v) form because the
// families can differ across dimensions.

functions {
#include functions/gaussian_copula_uv.stan
#include functions/var_residuals.stan
#include functions/hmm_algorithms.stan
}

data {
  int<lower=2> n_time;                         // Number of time points
  int<lower=2> D;                              // Number of variables (= 2)
  matrix[n_time, D] Y;                         // Observed data (n_time x D)
  int<lower=2> K;                              // Number of hidden states
  array[D] int<lower=1, upper=4> family;       // Per-dimension margin family code
  vector[D] skew_direction;                    // Consulted only for exp/gamma dims

  // Prior hyperparameters (VAR)
  real<lower=0> sigma_mu_prior;
  real<lower=0> sigma_phi_prior;
  real<lower=0> sigma_eps_prior;               // Prior mean for normal innovation SDs

  // Prior hyperparameters (HMM)
  real<lower=0> kappa;                         // Sticky Dirichlet self-transition concentration
  real<lower=0> alpha_off;                     // Sticky Dirichlet off-diagonal concentration
  real<lower=0> z_rho_prior_sd;                // Prior SD for state-specific z_rho
}

transformed data {
  int n_time_eff = n_time - 1;
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

  // HMM parameters (ordered z_rho resolves label switching)
  ordered[K] z_rho;
  simplex[K] pi0;
  array[K] simplex[K] A;

  // Union of per-dimension marginal parameters
  vector<lower=0.01>[D] sigma_eps;             // normal innovation SDs
  vector[D] eta;                               // exponential / gamma scale (unconstrained)
  vector<lower=0>[D] omega;                    // skew_normal scale
  vector<lower=-1, upper=1>[D] delta;          // skew_normal CP skewness
  vector<lower=0>[D] shape_gam;                // gamma shape (per dimension)
}

transformed parameters {
  vector[K] rho_state;
  matrix[n_time_eff, D] eps;
  matrix[n_time_eff, K] obs_ll;
  matrix[K, K] log_A;
  vector[K] log_pi0;
  matrix[n_time_eff, K] log_alpha;

  // Derived per-dimension scale parameters (filled only on their own family's
  // dimensions; 0 elsewhere). Declared as transformed parameters so the induced
  // priors (model block) and PPC (generated quantities) reuse them and so
  // sigma_exp / sigma_gam are available for extraction.
  vector[D] sigma_exp = rep_vector(0, D);
  vector[D] rate_exp = rep_vector(0, D);
  vector[D] mean_gam = rep_vector(0, D);
  vector[D] rate_gam = rep_vector(0, D);
  vector[D] sigma_gam = rep_vector(0, D);

  // Skew-normal CP -> DP transform (computed for all dims, used where family == 3)
  vector[D] alpha = delta ./ sqrt(1 - square(delta));
  vector[D] xi = -omega .* (delta * SQRT_2_OVER_PI);

  for (k in 1:K) rho_state[k] = tanh(z_rho[k]);
  for (j in 1:K) {
    log_pi0[j] = log(pi0[j]);
    for (k in 1:K) log_A[j, k] = log(A[j][k]);
  }

  eps = compute_var_residuals(Y, mu, Phi, n_time_eff, D);

  // Feasibility bounds and derived scales for exponential / gamma dimensions
  for (i in 1:D) {
    if (family[i] == 2 || family[i] == 4) {
      real m = -skew_direction[i] * eps[1, i];
      for (t in 2:n_time_eff) m = fmax(m, -skew_direction[i] * eps[t, i]);
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

  // Per-state observation log-likelihoods: shared marginal density + the
  // state-specific Gaussian copula on the per-dimension CDF values.
  for (t in 1:n_time_eff) {
    real marginal_ll = 0;
    vector[2] u_vec;
    for (i in 1:D) {
      if (family[i] == 1) {
        marginal_ll += normal_lpdf(eps[t, i] | 0, sigma_eps[i]);
        u_vec[i] = Phi(eps[t, i] / sigma_eps[i]);
      } else if (family[i] == 2) {
        real x_shifted = sigma_exp[i] + skew_direction[i] * eps[t, i];
        marginal_ll += exponential_lpdf(x_shifted | rate_exp[i]);
        u_vec[i] = exponential_cdf(x_shifted | rate_exp[i]);
        if (skew_direction[i] < 0) u_vec[i] = 1.0 - u_vec[i];
      } else if (family[i] == 3) {
        marginal_ll += skew_normal_lpdf(eps[t, i] | xi[i], omega[i], alpha[i]);
        u_vec[i] = skew_normal_cdf(eps[t, i] | xi[i], omega[i], alpha[i]);
      } else {
        real x_shifted = mean_gam[i] + skew_direction[i] * eps[t, i];
        marginal_ll += gamma_lpdf(x_shifted | shape_gam[i], rate_gam[i]);
        u_vec[i] = gamma_cdf(x_shifted | shape_gam[i], rate_gam[i]);
        if (skew_direction[i] < 0) u_vec[i] = 1.0 - u_vec[i];
      }
    }
    for (k in 1:K) {
      obs_ll[t, k] = marginal_ll + gaussian_copula_uv_lpdf(u_vec | rho_state[k]);
    }
  }

  log_alpha = hmm_forward(obs_ll, log_A, log_pi0, n_time_eff, K);
}

model {
  // Shared priors
  mu ~ normal(0, sigma_mu_prior);
  to_vector(Phi) ~ normal(0, sigma_phi_prior);

  // HMM priors
  z_rho ~ normal(0, z_rho_prior_sd);
  pi0 ~ dirichlet(rep_vector(1.0, K));
  for (k in 1:K) A[k] ~ dirichlet(dirichlet_prior[k]);

  // Union-parameter priors kept proper unconditionally
  sigma_eps ~ exponential(1.0 / sigma_eps_prior);
  omega ~ normal(0, 1);
  delta ~ normal(0, 0.5);
  shape_gam ~ lognormal(log(1), 0.5);

  // Per-dimension induced scale priors / eta priors
  for (i in 1:D) {
    if (family[i] == 2) {
      target += lognormal_lpdf(sigma_exp[i] | 0, 0.5) + eta[i];
    } else if (family[i] == 4) {
      target += lognormal_lpdf(mean_gam[i] | 0, 0.5) + eta[i];
    } else {
      eta[i] ~ std_normal();
    }
  }

  // Marginal log-likelihood via the forward algorithm
  target += log_sum_exp(to_vector(log_alpha[n_time_eff, ]));
}

generated quantities {
  matrix[n_time_eff, K] gamma;
  array[n_time_eff] int viterbi_state;
  vector[n_time_eff] rho_hmm;
  vector[n_time_eff] log_lik;
  matrix[n_time_eff, D] eps_rep;

  gamma = hmm_state_posteriors(log_alpha, obs_ll, log_A, n_time_eff, K);
  viterbi_state = hmm_viterbi(obs_ll, log_A, log_pi0, n_time_eff, K);
  rho_hmm = hmm_rho_average(gamma, rho_state, n_time_eff, K);
  log_lik = hmm_log_lik(log_alpha, n_time_eff, K);

  // Posterior predictive using the posterior-averaged rho; invert per dimension
  // where an inverse CDF exists (normal, exponential). skew_normal / gamma dims
  // carry copula-level z-scores (plot_ppc guards fits containing such a dim).
  for (t in 1:n_time_eff) {
    real z1_rep = std_normal_rng();
    real z2_rep = rho_hmm[t] * z1_rep + sqrt(1 - square(rho_hmm[t])) * std_normal_rng();
    array[2] real z_rep = {z1_rep, z2_rep};
    for (i in 1:D) {
      real u_i = Phi(z_rep[i]);
      if (family[i] == 1) {
        eps_rep[t, i] = z_rep[i] * sigma_eps[i];
      } else if (family[i] == 2) {
        // The likelihood uses u = 1 - F(x_shifted) on left-skewed dims,
        // so invert at the flipped uniform.
        real u_eff = skew_direction[i] < 0 ? 1 - u_i : u_i;
        eps_rep[t, i] = skew_direction[i] * (-log1m(u_eff) / rate_exp[i] - sigma_exp[i]);
      } else {
        eps_rep[t, i] = z_rep[i];
      }
    }
  }
}
