// Constant Copula VAR Model with per-variable (mixed) margins and a Clayton copula
// VAR(1) with a time-invariant Clayton copula (single dependence parameter
// theta) where each dimension can use its own marginal family. This is the
// Clayton counterpart of constant_mixed.stan: the per-dimension marginal
// dispatch is identical, only the copula is swapped. The R routing layer
// dispatches here for mixed margins with copula = "clayton"; same-family normal
// requests keep routing to constant_NCl.stan.
//
// Family codes (per dimension): 1 = normal, 2 = exponential, 3 = skew_normal,
// 4 = gamma. The Clayton copula naturally takes CDF values, so mixed margins
// work directly; all per-family parameters are declared unconditionally and a
// dimension that does not use a parameter samples it from its (proper) prior.

functions {
#include functions/clayton_copula.stan
#include functions/var_residuals.stan
}

data {
  int<lower=2> n_time;                         // Number of time points
  int<lower=2> D;                              // Number of variables (= 2)
  matrix[n_time, D] Y;                         // Observed data (n_time x D)
  array[D] int<lower=1, upper=4> family;       // Per-dimension margin family code
  vector[D] skew_direction;                    // Consulted only for exp/gamma dims

  // Prior hyperparameters
  real<lower=0> sigma_mu_prior;                // Prior SD for intercepts
  real<lower=0> sigma_phi_prior;               // Prior SD for VAR coefficients
  real<lower=0> sigma_eps_prior;               // Prior mean for normal innovation SDs (exponential prior scale)
}

transformed data {
  int n_time_eff = n_time - 1;
  real SQRT_2_OVER_PI = sqrt(2.0 / pi());
}

parameters {
  vector[D] mu;                                // Intercepts
  matrix[D, D] Phi;                            // VAR(1) coefficient matrix
  real<lower=1e-6> theta;                      // Clayton dependence parameter

  // Union of per-dimension marginal parameters
  vector<lower=0.01>[D] sigma_eps;             // normal innovation SDs
  vector[D] eta;                               // exponential / gamma scale (unconstrained)
  vector<lower=0>[D] omega;                    // skew_normal scale
  vector<lower=-1, upper=1>[D] delta;          // skew_normal CP skewness
  vector<lower=0>[D] shape_gam;                // gamma shape (per dimension)
}

transformed parameters {
  matrix[n_time_eff, D] eps = compute_var_residuals(Y, mu, Phi, n_time_eff, D);

  // Skew-normal CP -> DP transform (computed for all dims, used where family == 3)
  vector[D] alpha = delta ./ sqrt(1 - square(delta));
  vector[D] xi = -omega .* (delta * SQRT_2_OVER_PI);
}

model {
  real tiny = 1e-9;
  vector[D] sigma_exp = rep_vector(0, D);
  vector[D] rate_exp = rep_vector(0, D);
  vector[D] mean_gam = rep_vector(0, D);
  vector[D] rate_gam = rep_vector(0, D);

  // Shared priors
  mu ~ normal(0, sigma_mu_prior);
  to_vector(Phi) ~ normal(0, sigma_phi_prior);
  theta ~ lognormal(0, 1);

  // Union-parameter priors kept proper unconditionally
  sigma_eps ~ exponential(1.0 / sigma_eps_prior);
  omega ~ normal(0, 1);
  delta ~ normal(0, 0.5);
  shape_gam ~ lognormal(log(1), 0.5);

  for (i in 1:D) {
    if (family[i] == 2 || family[i] == 4) {
      real m = -skew_direction[i] * eps[1, i];
      for (t in 2:n_time_eff) m = fmax(m, -skew_direction[i] * eps[t, i]);
      real lb = fmax(m, 0);
      if (family[i] == 2) {
        sigma_exp[i] = lb + exp(eta[i]) + tiny;
        rate_exp[i] = 1.0 / sigma_exp[i];
        target += lognormal_lpdf(sigma_exp[i] | 0, 0.5) + eta[i];
      } else {
        mean_gam[i] = lb + exp(eta[i]) + tiny;
        rate_gam[i] = shape_gam[i] / mean_gam[i];
        target += lognormal_lpdf(mean_gam[i] | 0, 0.5) + eta[i];
      }
    } else {
      eta[i] ~ std_normal();
    }
  }

  // Likelihood: per-dimension margins + Clayton copula on the CDF values
  for (t in 1:n_time_eff) {
    row_vector[D] res = eps[t];
    vector[2] u_vec;

    for (i in 1:D) {
      if (family[i] == 1) {
        target += normal_lpdf(res[i] | 0, sigma_eps[i]);
        u_vec[i] = Phi(res[i] / sigma_eps[i]);
      } else if (family[i] == 2) {
        real x_shifted = sigma_exp[i] + skew_direction[i] * res[i];
        target += exponential_lpdf(x_shifted | rate_exp[i]);
        u_vec[i] = exponential_cdf(x_shifted | rate_exp[i]);
        if (skew_direction[i] < 0) u_vec[i] = 1.0 - u_vec[i];
      } else if (family[i] == 3) {
        target += skew_normal_lpdf(res[i] | xi[i], omega[i], alpha[i]);
        u_vec[i] = skew_normal_cdf(res[i] | xi[i], omega[i], alpha[i]);
      } else {
        real x_shifted = mean_gam[i] + skew_direction[i] * res[i];
        target += gamma_lpdf(x_shifted | shape_gam[i], rate_gam[i]);
        u_vec[i] = gamma_cdf(x_shifted | shape_gam[i], rate_gam[i]);
        if (skew_direction[i] < 0) u_vec[i] = 1.0 - u_vec[i];
      }
    }

    target += clayton_copula_ld(u_vec[1], u_vec[2], theta);
  }
}

generated quantities {
  vector[n_time_eff] log_lik;
  matrix[n_time_eff, D] eps_rep;

  // Reported per-dimension scale parameters (0 on dims of a different family).
  vector[D] sigma_exp = rep_vector(0, D);
  vector[D] sigma_gam = rep_vector(0, D);

  {
    real tiny = 1e-9;
    vector[D] rate_exp = rep_vector(0, D);
    vector[D] mean_gam = rep_vector(0, D);
    vector[D] rate_gam = rep_vector(0, D);

    for (i in 1:D) {
      if (family[i] == 2 || family[i] == 4) {
        real m = -skew_direction[i] * eps[1, i];
        for (t in 2:n_time_eff) m = fmax(m, -skew_direction[i] * eps[t, i]);
        real lb = fmax(m, 0);
        if (family[i] == 2) {
          sigma_exp[i] = lb + exp(eta[i]) + tiny;
          rate_exp[i] = 1.0 / sigma_exp[i];
        } else {
          mean_gam[i] = lb + exp(eta[i]) + tiny;
          sigma_gam[i] = mean_gam[i] / sqrt(shape_gam[i]);
          rate_gam[i] = shape_gam[i] / mean_gam[i];
        }
      }
    }

    for (t in 1:n_time_eff) {
      log_lik[t] = 0;
      vector[2] u_vec;
      for (i in 1:D) {
        if (family[i] == 1) {
          log_lik[t] += normal_lpdf(eps[t, i] | 0, sigma_eps[i]);
          u_vec[i] = Phi(eps[t, i] / sigma_eps[i]);
        } else if (family[i] == 2) {
          real x_shifted = sigma_exp[i] + skew_direction[i] * eps[t, i];
          log_lik[t] += exponential_lpdf(x_shifted | rate_exp[i]);
          u_vec[i] = exponential_cdf(x_shifted | rate_exp[i]);
          if (skew_direction[i] < 0) u_vec[i] = 1.0 - u_vec[i];
        } else if (family[i] == 3) {
          log_lik[t] += skew_normal_lpdf(eps[t, i] | xi[i], omega[i], alpha[i]);
          u_vec[i] = skew_normal_cdf(eps[t, i] | xi[i], omega[i], alpha[i]);
        } else {
          real x_shifted = mean_gam[i] + skew_direction[i] * eps[t, i];
          log_lik[t] += gamma_lpdf(x_shifted | shape_gam[i], rate_gam[i]);
          u_vec[i] = gamma_cdf(x_shifted | shape_gam[i], rate_gam[i]);
          if (skew_direction[i] < 0) u_vec[i] = 1.0 - u_vec[i];
        }
      }
      log_lik[t] += clayton_copula_ld(u_vec[1], u_vec[2], theta);

      // Posterior predictive: sample (u, v) from the Clayton copula via the
      // conditional inverse, then invert per dimension where possible. normal ->
      // inv_Phi(u) * sigma; exponential -> inverse-CDF; skew_normal / gamma lack
      // a Stan inverse CDF, so those dims carry copula-level z-scores.
      {
        real u1 = uniform_rng(0, 1);
        real q = uniform_rng(0, 1);
        real u2 = pow(1 + pow(u1, -theta) * (pow(q, -theta / (1 + theta)) - 1), -1 / theta);
        array[2] real u_rep = {u1, u2};
        for (i in 1:D) {
          real ui = fmax(1e-9, fmin(1 - 1e-9, u_rep[i]));
          if (family[i] == 1) {
            eps_rep[t, i] = inv_Phi(ui) * sigma_eps[i];
          } else if (family[i] == 2) {
            eps_rep[t, i] = skew_direction[i] * (-log1m(ui) / rate_exp[i] - sigma_exp[i]);
          } else {
            eps_rep[t, i] = inv_Phi(ui);  // skew_normal / gamma: copula-level z-score
          }
        }
      }
    }
  }
}
