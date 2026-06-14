// Full Markov-switching copula VAR(1) engine.
//
// Generalises hmm_mixed.stan: in addition to the always state-specific copula
// correlation, the VAR intercepts mu, the VAR(1) coefficients Phi (per-coefficient
// mask), the residual scales, and the marginal FAMILY itself can be regime-specific.
// Which components switch is selected by data flags; the per-state family assignment
// is supplied as DATA (family[k, d], not estimated). Backward compatibility: the
// default dcvar_hmm() (only rho switches, global margins) keeps routing to the
// specialised legacy files; this generic engine is selected internally only when
// extra switching or per-state families are requested.
//
// Family codes (per state, per dimension): 1 = normal, 2 = exponential,
// 3 = skew_normal, 4 = gamma. The Gaussian copula is applied via the CDF (u,v)
// form because families can differ across states and dimensions.
//
// Identifiability: rho is ALWAYS state-specific and `ordered[K] z_rho` is the sole
// label-switching anchor (states ordered by increasing correlation). mu/Phi/margin
// switching are additive on top. The R layer forbids dropping rho while other
// components switch, and orders any per-state `margins` list in increasing-rho order.
//
// Switching idioms (conditional sizing, masked-coefficient baseline+deviation,
// zero-extent containers) mirror dcvar_dynamic.stan.

functions {
#include functions/gaussian_copula_uv.stan
#include functions/var_residuals.stan
#include functions/hmm_algorithms.stan
}

data {
  int<lower=2> n_time;                          // Number of time points
  int<lower=2, upper=2> D;  // copula code is hard-wired bivariate
  matrix[n_time, D] Y;                          // Observed data (n_time x D)
  int<lower=2> K;                               // Number of hidden states

  // Per-state family assignment (DATA, not estimated). Rows are identical when
  // margins are global; they differ for a per-state margins specification.
  array[K, D] int<lower=1, upper=4> family;
  array[K, D] int<lower=-1, upper=1> skew_direction; // consulted only on exp/gamma dims

  // Switch flags (rho always switches; not a flag). phi_switch_mask is row-major
  // (phi11, phi12, phi21, phi22), as in dcvar_dynamic.stan.
  int<lower=0, upper=1> switch_mu;
  int<lower=0, upper=1> switch_phi;
  array[4] int<lower=0, upper=1> phi_switch_mask;
  int<lower=0, upper=1> switch_margins;

  // Prior hyperparameters (scalars; broadcast across states)
  real<lower=0> sigma_mu_prior;
  real<lower=0> sigma_phi_prior;
  real<lower=0> prior_phi_dev_sd;               // SD of the per-state Phi deviations
  real<lower=0> sigma_eps_prior;
  real<lower=0> kappa;                          // Sticky Dirichlet self-transition
  real<lower=0> alpha_off;                      // Sticky Dirichlet off-diagonal
  real<lower=0> z_rho_prior_sd;
}

transformed data {
  int n_time_eff = n_time - 1;
  real SQRT_2_OVER_PI = sqrt(2.0 / pi());
  int Mu_K  = switch_mu ? K : 1;                // location extent
  int Mrg_K = switch_margins ? K : 1;           // margin-config extent
  array[K] vector[K] dirichlet_prior;
  for (k in 1:K) {
    for (j in 1:K) {
      dirichlet_prior[k][j] = (k == j) ? kappa : alpha_off;
    }
  }
}

parameters {
  // VAR location: K-sized iff switching, else shared (length 1)
  array[Mu_K] vector[D] mu;

  // VAR dynamics: shared baseline + per-state deviations on masked coefficients
  matrix[D, D] Phi_base;
  array[switch_phi ? K : 0] matrix[D, D] Phi_dev;

  // Copula correlation: always state-specific, ordered for label identification
  ordered[K] z_rho;
  simplex[K] pi0;
  array[K] simplex[K] A;

  // Per-margin-config marginal union (config == state when switch_margins, else shared)
  array[Mrg_K] vector<lower=0.01>[D] sigma_eps; // normal innovation SDs
  array[Mrg_K] vector[D] eta;                   // exponential / gamma log-scale
  array[Mrg_K] vector<lower=0>[D] omega;        // skew_normal scale
  array[Mrg_K] vector<lower=-1, upper=1>[D] delta; // skew_normal CP skewness
  array[Mrg_K] vector<lower=0>[D] shape_gam;    // gamma shape
}

transformed parameters {
  vector[K] rho_state;
  matrix[K, K] log_A;
  vector[K] log_pi0;
  matrix[n_time_eff, K] obs_ll;
  matrix[n_time_eff, K] log_alpha;

  // Per-state effective VAR parameters and residuals
  array[K] vector[D] mu_eff;
  array[K] matrix[D, D] Phi_eff;
  array[K] matrix[n_time_eff, D] eps_k;

  // Derived per-config scales (filled only on their own family's dimensions; 0 else)
  array[Mrg_K] vector[D] sigma_exp = rep_array(rep_vector(0, D), Mrg_K);
  array[Mrg_K] vector[D] rate_exp  = rep_array(rep_vector(0, D), Mrg_K);
  array[Mrg_K] vector[D] mean_gam  = rep_array(rep_vector(0, D), Mrg_K);
  array[Mrg_K] vector[D] rate_gam  = rep_array(rep_vector(0, D), Mrg_K);
  array[Mrg_K] vector[D] sigma_gam = rep_array(rep_vector(0, D), Mrg_K);
  array[Mrg_K] vector[D] sn_alpha;
  array[Mrg_K] vector[D] sn_xi;

  for (k in 1:K) rho_state[k] = tanh(z_rho[k]);
  for (j in 1:K) {
    log_pi0[j] = log(pi0[j]);
    for (k in 1:K) log_A[j, k] = log(A[j][k]);
  }

  // Effective mu / Phi per state and the resulting residual series.
  for (k in 1:K) {
    mu_eff[k] = mu[switch_mu ? k : 1];
    matrix[D, D] P = Phi_base;
    if (switch_phi == 1) {
      if (phi_switch_mask[1] == 1) P[1, 1] += Phi_dev[k][1, 1];
      if (phi_switch_mask[2] == 1) P[1, 2] += Phi_dev[k][1, 2];
      if (phi_switch_mask[3] == 1) P[2, 1] += Phi_dev[k][2, 1];
      if (phi_switch_mask[4] == 1) P[2, 2] += Phi_dev[k][2, 2];
    }
    Phi_eff[k] = P;
    eps_k[k] = compute_var_residuals(Y, mu_eff[k], Phi_eff[k], n_time_eff, D);
  }

  // Skew-normal CP -> DP transform and exp/gamma feasibility-bounded scales, per
  // margin config m (residuals of state m; when margins are shared, m = 1 and all
  // residual series coincide because mu/Phi are then global as well).
  for (m in 1:Mrg_K) {
    sn_alpha[m] = delta[m] ./ sqrt(1 - square(delta[m]));
    sn_xi[m] = -omega[m] .* (delta[m] * SQRT_2_OVER_PI);
    for (d in 1:D) {
      int fam = family[m, d];
      if (fam == 2 || fam == 4) {
        int sd = skew_direction[m, d];
        real mx = -sd * eps_k[m][1, d];
        for (t in 2:n_time_eff) mx = fmax(mx, -sd * eps_k[m][t, d]);
        real lb = fmax(mx, 0);
        if (fam == 2) {
          sigma_exp[m][d] = lb + exp(eta[m][d]) + 1e-9;
          rate_exp[m][d] = 1.0 / sigma_exp[m][d];
        } else {
          mean_gam[m][d] = lb + exp(eta[m][d]) + 1e-9;
          rate_gam[m][d] = shape_gam[m][d] / mean_gam[m][d];
          sigma_gam[m][d] = mean_gam[m][d] / sqrt(shape_gam[m][d]);
        }
      }
    }
  }

  // Per-state observation log-likelihoods: state-specific marginal density on the
  // state's residuals + the state-specific Gaussian copula.
  for (k in 1:K) {
    int mk = switch_margins ? k : 1;
    for (t in 1:n_time_eff) {
      real marginal_ll = 0;
      vector[2] u_vec;
      for (d in 1:D) {
        int fam = family[k, d];
        int sd = skew_direction[k, d];
        real e = eps_k[k][t, d];
        if (fam == 1) {
          marginal_ll += normal_lpdf(e | 0, sigma_eps[mk][d]);
          u_vec[d] = Phi(e / sigma_eps[mk][d]);
        } else if (fam == 2) {
          real x = sigma_exp[mk][d] + sd * e;
          marginal_ll += exponential_lpdf(x | rate_exp[mk][d]);
          u_vec[d] = exponential_cdf(x | rate_exp[mk][d]);
          if (sd < 0) u_vec[d] = 1.0 - u_vec[d];
        } else if (fam == 3) {
          marginal_ll += skew_normal_lpdf(e | sn_xi[mk][d], omega[mk][d], sn_alpha[mk][d]);
          u_vec[d] = skew_normal_cdf(e | sn_xi[mk][d], omega[mk][d], sn_alpha[mk][d]);
        } else {
          real x = mean_gam[mk][d] + sd * e;
          marginal_ll += gamma_lpdf(x | shape_gam[mk][d], rate_gam[mk][d]);
          u_vec[d] = gamma_cdf(x | shape_gam[mk][d], rate_gam[mk][d]);
          if (sd < 0) u_vec[d] = 1.0 - u_vec[d];
        }
      }
      obs_ll[t, k] = marginal_ll + gaussian_copula_uv_lpdf(u_vec | rho_state[k]);
    }
  }

  log_alpha = hmm_forward(obs_ll, log_A, log_pi0, n_time_eff, K);
}

model {
  // Location / dynamics priors (loop over the actual extent so length-1 and
  // K-sized groups each get exactly the right number of proper priors).
  for (i in 1:Mu_K) mu[i] ~ normal(0, sigma_mu_prior);
  to_vector(Phi_base) ~ normal(0, sigma_phi_prior);
  for (i in 1:size(Phi_dev)) to_vector(Phi_dev[i]) ~ normal(0, prior_phi_dev_sd);

  // HMM priors
  z_rho ~ normal(0, z_rho_prior_sd);
  pi0 ~ dirichlet(rep_vector(1.0, K));
  for (k in 1:K) A[k] ~ dirichlet(dirichlet_prior[k]);

  // Margin-union priors per config, kept proper unconditionally (mirror hmm_mixed)
  for (m in 1:Mrg_K) {
    sigma_eps[m] ~ exponential(1.0 / sigma_eps_prior);
    omega[m] ~ normal(0, 1);
    delta[m] ~ normal(0, 0.5);
    shape_gam[m] ~ lognormal(log(1), 0.5);
  }

  // Induced lognormal scale priors per (config, dim) for exp/gamma; eta std_normal
  // otherwise. The + eta Jacobian is the change of variables for lb + exp(eta).
  for (m in 1:Mrg_K) {
    for (d in 1:D) {
      int fam = family[m, d];
      if (fam == 2) {
        target += lognormal_lpdf(sigma_exp[m][d] | 0, 0.5) + eta[m][d];
      } else if (fam == 4) {
        target += lognormal_lpdf(mean_gam[m][d] | 0, 0.5) + eta[m][d];
      } else {
        eta[m][d] ~ std_normal();
      }
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

  // Posterior predictive from the regime mixture: draw a state from the smoothed
  // gamma, then use THAT state's family/params. Invert per dimension where an
  // inverse CDF exists (normal, exponential); skew_normal/gamma dims carry
  // copula-level z-scores (plot_ppc guards fits containing such a dim).
  for (t in 1:n_time_eff) {
    vector[K] gt = to_vector(gamma[t, ]);
    int k = categorical_rng(gt / sum(gt));
    int mk = switch_margins ? k : 1;
    real rho_rep = rho_state[k];
    real z1_rep = std_normal_rng();
    real z2_rep = rho_rep * z1_rep + sqrt(1 - square(rho_rep)) * std_normal_rng();
    array[2] real z_rep = {z1_rep, z2_rep};
    for (d in 1:D) {
      real u_i = Phi(z_rep[d]);
      int fam = family[k, d];
      int sd = skew_direction[k, d];
      if (fam == 1) {
        eps_rep[t, d] = z_rep[d] * sigma_eps[mk][d];
      } else if (fam == 2) {
        real u_eff = sd < 0 ? 1 - u_i : u_i;
        eps_rep[t, d] = sd * (-log1m(u_eff) / rate_exp[mk][d] - sigma_exp[mk][d]);
      } else {
        eps_rep[t, d] = z_rep[d];
      }
    }
  }
}
