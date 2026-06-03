// Constant Copula VAR Model with per-variable (mixed) margins
// VAR(1) with a time-invariant Gaussian copula (single fixed rho) where each
// dimension can use its own marginal family. This is the generic model that
// the R routing layer dispatches to whenever the requested per-dimension
// margins differ. Same-family requests keep routing to the specialised models
// (constant_copula_var / constant_EG / constant_SNG / constant_GG), which
// preserves their exact results and the gamma shared-scalar shape semantics.
//
// Family codes (per dimension): 1 = normal, 2 = exponential, 3 = skew_normal,
// 4 = gamma. The marginal math for each family is identical to the specialised
// model of the same name; here the per-dimension blocks are gated on family[i].
//
// All per-family parameters are declared unconditionally (the union of the
// specialised models). On a dimension whose family does not use a given
// parameter, that parameter is unused in the likelihood and simply samples
// from its (proper) prior. This keeps the posterior proper and is negligible
// at D = 2. The Gaussian copula is always applied via the CDF (u,v) form
// because the families can differ across dimensions.

functions {
#include functions/gaussian_copula_uv.stan
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
  real<lower=0> sigma_eps_prior;               // Prior mean for normal innovation SDs
  real<lower=0> z_rho_prior_sd;                // Prior SD for rho (Fisher-z scale)
}

transformed data {
  int n_time_eff = n_time - 1;
  real SQRT_2_OVER_PI = sqrt(2.0 / pi());
}

parameters {
  // VAR parameters (constant over time)
  vector[D] mu;                                // Intercepts
  matrix[D, D] Phi;                            // VAR(1) coefficient matrix
  real z_rho;                                  // Copula correlation (Fisher-z scale)

  // Union of per-dimension marginal parameters
  vector<lower=0.01>[D] sigma_eps;             // normal innovation SDs
  vector[D] eta;                               // exponential / gamma scale (unconstrained)
  vector<lower=0>[D] omega;                    // skew_normal scale
  vector<lower=-1, upper=1>[D] delta;          // skew_normal CP skewness
  vector<lower=0>[D] shape_gam;                // gamma shape (per dimension)
}

transformed parameters {
  matrix[n_time_eff, D] eps = compute_var_residuals(Y, mu, Phi, n_time_eff, D);
  real rho = tanh(z_rho);

  // Skew-normal CP -> DP transform (computed for all dims, used where family == 3)
  vector[D] alpha = delta ./ sqrt(1 - square(delta));
  vector[D] xi = -omega .* (delta * SQRT_2_OVER_PI);
}

model {
  real tiny = 1e-9;
  vector[D] sigma_exp = rep_vector(0, D);      // exponential scale (filled for exp dims)
  vector[D] rate_exp = rep_vector(0, D);
  vector[D] mean_gam = rep_vector(0, D);       // gamma mean (filled for gamma dims)
  vector[D] rate_gam = rep_vector(0, D);

  // Shared priors
  mu ~ normal(0, sigma_mu_prior);
  to_vector(Phi) ~ normal(0, sigma_phi_prior);
  z_rho ~ normal(0, z_rho_prior_sd);

  // Union-parameter priors kept proper unconditionally (so unused params on a
  // given dimension remain proper draws). eta is handled per-dimension below.
  sigma_eps ~ exponential(1.0 / sigma_eps_prior);
  omega ~ normal(0, 1);
  delta ~ normal(0, 0.5);
  shape_gam ~ lognormal(log(1), 0.5);

  // Per-dimension feasibility bounds, induced scale priors, and eta priors
  for (i in 1:D) {
    if (family[i] == 2 || family[i] == 4) {
      // Shifted-positive support bound shared by exponential and gamma margins
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
      // eta unused on this dimension; keep it proper
      eta[i] ~ std_normal();
    }
  }

  // Likelihood: per-dimension margins + Gaussian copula on the CDF values
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

    target += gaussian_copula_uv_lpdf(u_vec | rho);
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
      log_lik[t] += gaussian_copula_uv_lpdf(u_vec | rho);

      // Posterior predictive: draw from the bivariate-normal copula, then invert
      // per dimension where an inverse CDF exists. normal -> z * sigma;
      // exponential -> inverse-CDF on the shifted exponential. skew_normal and
      // gamma lack a Stan inverse CDF, so those dims carry copula-level z-scores
      // (plot_ppc() guards fits that contain such a dimension).
      {
        real z1_rep = std_normal_rng();
        real z2_rep = rho * z1_rep + sqrt(1 - square(rho)) * std_normal_rng();
        array[2] real z_rep = {z1_rep, z2_rep};
        for (i in 1:D) {
          real u_i = Phi(z_rep[i]);
          if (family[i] == 1) {
            eps_rep[t, i] = z_rep[i] * sigma_eps[i];
          } else if (family[i] == 2) {
            eps_rep[t, i] = skew_direction[i] * (-log1m(u_i) / rate_exp[i] - sigma_exp[i]);
          } else {
            eps_rep[t, i] = z_rep[i];  // skew_normal / gamma: copula-level z-score
          }
        }
      }
    }
  }
}
