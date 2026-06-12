// TV-VAR Copula Model - fully time-varying parameters, per-variable margins
// Non-centered parameterization throughout.
//
// Extends dcvar_mixed_ncp.stan with optional time-varying VAR coefficients
// (tv_phi: 4 independent random walks on phi11, phi12, phi21, phi22) and
// optional time-varying residual scales (tv_sigma: log-scale random walks on
// the normal / skew-normal dimensions). The copula correlation rho(t) is
// always time-varying (the Fisher-z random walk of the dcvar family).
//
// Every time-varying component is parameterized as a constant baseline plus a
// non-centered random-walk deviation, and the deviation containers are
// conditionally sized to zero when the corresponding flag is off. With
// tv_phi = tv_sigma = 0 the target density is identical term-by-term to
// dcvar_mixed_ncp.stan, so this model exactly nests the plain DC-VAR.
//
// Family codes (per dimension): 1 = normal, 2 = exponential, 3 = skew_normal,
// 4 = gamma. All per-family parameters are declared unconditionally; on a
// dimension whose family does not use a parameter it samples from its (proper)
// prior. The Gaussian copula is applied via the CDF (u,v) form because the
// families can differ across dimensions.
//
// Exponential / gamma dimensions keep a CONSTANT scale even when tv_sigma = 1
// (deviation forced to zero): with the shifted margins, a pointwise
// feasibility bound makes the likelihood insensitive to all below-mean
// residuals, and a global bound floors the scale path. The R layer warns
// about this restriction; a soft-barrier construction is the designated
// future extension.

functions {
#include functions/gaussian_copula_uv.stan
#include functions/var_residuals.stan
#include functions/ncp_random_walk.stan
}

data {
  int<lower=2> n_time;                         // Number of time points
  int<lower=2, upper=2> D;  // copula code is hard-wired bivariate
  matrix[n_time, D] Y;                         // Observed data (n_time x D)
  array[D] int<lower=1, upper=4> family;       // Per-dimension margin family code
  vector[D] skew_direction;                    // Consulted only for exp/gamma dims

  int<lower=0, upper=1> tv_phi;                // time-varying Phi(t)
  int<lower=0, upper=1> tv_sigma;              // time-varying scales (normal/skew-normal dims)

  // Prior hyperparameters
  real<lower=0> sigma_mu_prior;                // Prior SD for intercepts
  real<lower=0> sigma_phi_prior;               // Prior SD for baseline VAR coefficients
  real<lower=0> sigma_eps_prior;               // Prior mean for normal innovation SDs
  real<lower=0> sigma_omega_prior;             // Prior mean for rho random-walk SD
  real<lower=0> rho_init_prior_sd;             // Prior SD for initial rho (Fisher-z)
  real<lower=0> tau_phi_prior;                 // Prior mean for Phi random-walk SDs
  real<lower=0> tau_sigma_prior;               // Prior mean for log-scale random-walk SDs
}

transformed data {
  int n_time_eff = n_time - 1;
  real SQRT_2_OVER_PI = sqrt(2.0 / pi());
  int n_phi = tv_phi ? n_time_eff : 0;
  int n_sig = tv_sigma ? n_time_eff : 0;
}

parameters {
  // ---- constant baselines: identical to dcvar_mixed_ncp ----
  vector[D] mu;                                // Intercepts
  matrix[D, D] Phi;                            // Baseline VAR(1) coefficient matrix

  // Time-varying copula parameter (non-centred random walk)
  real z_rho_init;                             // Initial Fisher-z rho
  real<lower=0.001> sigma_omega;               // Random-walk innovation SD
  vector[n_time_eff] omega_raw;                // Raw innovations (std_normal)

  // Union of per-dimension marginal parameters
  vector<lower=0.01>[D] sigma_eps;             // normal baseline innovation SDs
  vector[D] eta;                               // exponential / gamma scale (unconstrained)
  vector<lower=0>[D] omega;                    // skew_normal baseline scale
  vector<lower=-1, upper=1>[D] delta;          // skew_normal CP skewness (constant in time)
  vector<lower=0>[D] shape_gam;                // gamma shape (per dimension)

  // ---- time-varying deviations (zero-sized when the flag is off) ----
  vector<lower=0>[tv_phi ? 4 : 0] tau_phi;     // per-coefficient walk innovation SDs
  matrix[n_phi, 4] phi_raw;                    // std-normal innovations, cols row-major (11,12,21,22)
  vector<lower=0>[tv_sigma ? D : 0] tau_sigma; // per-dimension log-scale walk SDs
  matrix[n_sig, D] sigma_raw;                  // std-normal innovations (inert on exp/gamma dims)
}

transformed parameters {
  // Non-centred random walk for the Fisher-z correlation
  vector[n_time_eff] z_rho = compute_rw_ncp(z_rho_init, sigma_omega, omega_raw, n_time_eff);
  vector[n_time_eff] rho;

  // Deviation walks (first state carries an innovation, package convention)
  matrix[n_phi, 4] phi_dev;
  matrix[n_sig, D] sigma_dev;                  // log-scale deviations; forced 0 on exp/gamma dims

  matrix[n_time_eff, D] eps;

  // Skew-normal CP -> DP transform (shape part is constant in time)
  vector[D] alpha = delta ./ sqrt(1 - square(delta));
  vector[D] xi = -omega .* (delta * SQRT_2_OVER_PI);

  for (t in 1:n_time_eff) rho[t] = tanh(z_rho[t]);

  if (tv_phi == 1) {
    matrix[n_time_eff, 4] phi_t_local;
    for (k in 1:4) {
      phi_dev[, k] = compute_rw_ncp(0, tau_phi[k], phi_raw[, k], n_time_eff);
    }
    phi_t_local[, 1] = phi_dev[, 1] + Phi[1, 1];
    phi_t_local[, 2] = phi_dev[, 2] + Phi[1, 2];
    phi_t_local[, 3] = phi_dev[, 3] + Phi[2, 1];
    phi_t_local[, 4] = phi_dev[, 4] + Phi[2, 2];
    eps = compute_var_residuals_tv(Y, mu, phi_t_local, n_time_eff, D);
  } else {
    eps = compute_var_residuals(Y, mu, Phi, n_time_eff, D);
  }

  if (tv_sigma == 1) {
    for (i in 1:D) {
      if (family[i] == 1 || family[i] == 3) {
        sigma_dev[, i] = compute_rw_ncp(0, tau_sigma[i], sigma_raw[, i], n_time_eff);
      } else {
        // Shifted exp/gamma margins keep a constant scale (see file header).
        sigma_dev[, i] = rep_vector(0, n_time_eff);
      }
    }
  }
}

model {
  real tiny = 1e-9;
  vector[D] sigma_exp = rep_vector(0, D);
  vector[D] rate_exp = rep_vector(0, D);
  vector[D] mean_gam = rep_vector(0, D);
  vector[D] rate_gam = rep_vector(0, D);

  // Shared priors (verbatim from dcvar_mixed_ncp)
  mu ~ normal(0, sigma_mu_prior);
  to_vector(Phi) ~ normal(0, sigma_phi_prior);
  z_rho_init ~ normal(0, rho_init_prior_sd);
  sigma_omega ~ exponential(1.0 / sigma_omega_prior);
  omega_raw ~ std_normal();

  // Union-parameter priors kept proper unconditionally
  sigma_eps ~ exponential(1.0 / sigma_eps_prior);
  omega ~ normal(0, 1);
  delta ~ normal(0, 0.5);
  shape_gam ~ lognormal(log(1), 0.5);

  // Time-varying deviation priors (no-ops on zero-sized containers)
  tau_phi ~ exponential(1.0 / tau_phi_prior);
  to_vector(phi_raw) ~ std_normal();
  tau_sigma ~ exponential(1.0 / tau_sigma_prior);
  to_vector(sigma_raw) ~ std_normal();

  // Per-dimension feasibility bounds, induced scale priors, and eta priors
  // (constant scales; eps already reflects Phi(t) when tv_phi = 1)
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

  // Likelihood: per-dimension margins (per-t scales where enabled) +
  // time-varying Gaussian copula
  for (t in 1:n_time_eff) {
    row_vector[D] res = eps[t];
    vector[2] u_vec;

    for (i in 1:D) {
      real dev = (tv_sigma == 1) ? sigma_dev[t, i] : 0;
      if (family[i] == 1) {
        real s_t = sigma_eps[i] * exp(dev);
        target += normal_lpdf(res[i] | 0, s_t);
        u_vec[i] = Phi(res[i] / s_t);
      } else if (family[i] == 2) {
        real x_shifted = sigma_exp[i] + skew_direction[i] * res[i];
        target += exponential_lpdf(x_shifted | rate_exp[i]);
        u_vec[i] = exponential_cdf(x_shifted | rate_exp[i]);
        if (skew_direction[i] < 0) u_vec[i] = 1.0 - u_vec[i];
      } else if (family[i] == 3) {
        real omega_t = omega[i] * exp(dev);
        real xi_t = -omega_t * (delta[i] * SQRT_2_OVER_PI);
        target += skew_normal_lpdf(res[i] | xi_t, omega_t, alpha[i]);
        u_vec[i] = skew_normal_cdf(res[i] | xi_t, omega_t, alpha[i]);
      } else {
        real x_shifted = mean_gam[i] + skew_direction[i] * res[i];
        target += gamma_lpdf(x_shifted | shape_gam[i], rate_gam[i]);
        u_vec[i] = gamma_cdf(x_shifted | shape_gam[i], rate_gam[i]);
        if (skew_direction[i] < 0) u_vec[i] = 1.0 - u_vec[i];
      }
    }

    target += gaussian_copula_uv_lpdf(u_vec | rho[t]);
  }
}

generated quantities {
  vector[n_time_eff] log_lik;
  matrix[n_time_eff, D] eps_rep;

  // Reported per-dimension scale parameters (0 on dims of a different family).
  vector[D] sigma_exp = rep_vector(0, D);
  vector[D] sigma_gam = rep_vector(0, D);

  // Trajectories (zero-sized when the component is disabled; the R layer
  // tiles the constant baselines instead).
  matrix[n_phi, 4] phi_t;                      // (phi11, phi12, phi21, phi22)(t)
  matrix[n_sig, D] sigma_t;                    // per-dim scale path (family's natural scale)

  // Stationarity monitoring: per t when tv_phi, single value otherwise.
  vector[tv_phi ? n_time_eff : 1] spectral_radius;
  int<lower=0> n_nonstationary = 0;

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

    // Phi trajectory + per-t spectral radius (2x2 closed form)
    if (tv_phi == 1) {
      for (t in 1:n_time_eff) {
        phi_t[t, 1] = Phi[1, 1] + phi_dev[t, 1];
        phi_t[t, 2] = Phi[1, 2] + phi_dev[t, 2];
        phi_t[t, 3] = Phi[2, 1] + phi_dev[t, 3];
        phi_t[t, 4] = Phi[2, 2] + phi_dev[t, 4];
        real tr = phi_t[t, 1] + phi_t[t, 4];
        real det_phi = phi_t[t, 1] * phi_t[t, 4] - phi_t[t, 2] * phi_t[t, 3];
        real disc = square(tr) - 4 * det_phi;
        if (disc >= 0) {
          real sqrt_disc = sqrt(disc);
          spectral_radius[t] = fmax(abs(0.5 * (tr + sqrt_disc)), abs(0.5 * (tr - sqrt_disc)));
        } else {
          spectral_radius[t] = sqrt(abs(det_phi));
        }
        if (spectral_radius[t] >= 1.0) n_nonstationary += 1;
      }
    } else {
      real tr = Phi[1, 1] + Phi[2, 2];
      real det_phi = Phi[1, 1] * Phi[2, 2] - Phi[1, 2] * Phi[2, 1];
      real disc = square(tr) - 4 * det_phi;
      if (disc >= 0) {
        real sqrt_disc = sqrt(disc);
        spectral_radius[1] = fmax(abs(0.5 * (tr + sqrt_disc)), abs(0.5 * (tr - sqrt_disc)));
      } else {
        spectral_radius[1] = sqrt(abs(det_phi));
      }
      if (spectral_radius[1] >= 1.0) n_nonstationary = 1;
    }

    // Scale trajectories on each family's natural reporting scale
    if (tv_sigma == 1) {
      for (i in 1:D) {
        for (t in 1:n_time_eff) {
          if (family[i] == 1) {
            sigma_t[t, i] = sigma_eps[i] * exp(sigma_dev[t, i]);
          } else if (family[i] == 3) {
            sigma_t[t, i] = omega[i] * exp(sigma_dev[t, i]);
          } else if (family[i] == 2) {
            sigma_t[t, i] = sigma_exp[i];
          } else {
            sigma_t[t, i] = sigma_gam[i];
          }
        }
      }
    }

    // log_lik + posterior predictive (mirrors the model block, incl. per-t
    // scales and the u = 1 - F flip on left-skewed dims)
    for (t in 1:n_time_eff) {
      log_lik[t] = 0;
      vector[2] u_vec;
      array[D] real s_norm;                    // per-t normal SDs for eps_rep
      for (i in 1:D) {
        real dev = (tv_sigma == 1) ? sigma_dev[t, i] : 0;
        if (family[i] == 1) {
          s_norm[i] = sigma_eps[i] * exp(dev);
          log_lik[t] += normal_lpdf(eps[t, i] | 0, s_norm[i]);
          u_vec[i] = Phi(eps[t, i] / s_norm[i]);
        } else if (family[i] == 2) {
          real x_shifted = sigma_exp[i] + skew_direction[i] * eps[t, i];
          log_lik[t] += exponential_lpdf(x_shifted | rate_exp[i]);
          u_vec[i] = exponential_cdf(x_shifted | rate_exp[i]);
          if (skew_direction[i] < 0) u_vec[i] = 1.0 - u_vec[i];
        } else if (family[i] == 3) {
          real omega_ti = omega[i] * exp(dev);
          real xi_ti = -omega_ti * (delta[i] * SQRT_2_OVER_PI);
          log_lik[t] += skew_normal_lpdf(eps[t, i] | xi_ti, omega_ti, alpha[i]);
          u_vec[i] = skew_normal_cdf(eps[t, i] | xi_ti, omega_ti, alpha[i]);
        } else {
          real x_shifted = mean_gam[i] + skew_direction[i] * eps[t, i];
          log_lik[t] += gamma_lpdf(x_shifted | shape_gam[i], rate_gam[i]);
          u_vec[i] = gamma_cdf(x_shifted | shape_gam[i], rate_gam[i]);
          if (skew_direction[i] < 0) u_vec[i] = 1.0 - u_vec[i];
        }
      }
      log_lik[t] += gaussian_copula_uv_lpdf(u_vec | rho[t]);

      // Posterior predictive: draw from the bivariate-normal copula at rho[t],
      // invert per dimension where an inverse CDF exists (normal, exponential);
      // skew_normal/gamma dims carry copula-level z-scores (plot_ppc guards them).
      {
        real z1_rep = std_normal_rng();
        real z2_rep = rho[t] * z1_rep + sqrt(1 - square(rho[t])) * std_normal_rng();
        array[2] real z_rep = {z1_rep, z2_rep};
        for (i in 1:D) {
          real u_i = Phi(z_rep[i]);
          if (family[i] == 1) {
            eps_rep[t, i] = z_rep[i] * s_norm[i];
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
  }
}
