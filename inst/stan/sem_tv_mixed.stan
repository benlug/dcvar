// Time-varying SEM Copula VAR(1) with fixed measurement model and mixed margins.
//
// The measurement model matches sem_mixed.stan. The latent VAR can carry
// optional time-varying Phi coefficients and optional time-varying innovation
// scales. TV SEM fits also use a Fisher-z random walk for the Gaussian copula
// correlation on the transition time axis.

functions {
#include functions/gaussian_copula_uv.stan
#include functions/ncp_random_walk.stan
#include functions/softplus.stan
}

data {
  int<lower=2> n_time;
  int<lower=1> J;
  matrix[n_time, 2 * J] y;
  vector[J] lambda;
  real<lower=0> sigma_e;
  array[2] int<lower=1, upper=4> family;
  vector[2] skew_direction;

  int<lower=0, upper=1> tv_phi;
  array[4] int<lower=0, upper=1> phi_tv_mask;
  int<lower=0, upper=1> tv_sigma;

  real<lower=0> prior_mu_sd;
  real<lower=0> prior_phi_sd;
  real<lower=0> prior_sigma_sd;
  real<lower=0> prior_rho_sd;
  real<lower=0> sigma_omega_prior;
  real<lower=0> tau_phi_prior;
  real<lower=0> tau_sigma_prior;
  real<lower=0> barrier_k;
}

transformed data {
  int n_time_eff = n_time - 1;
  real SQRT_2_OVER_PI = sqrt(2.0 / pi());
  int n_phi_tv = 0;
  for (k in 1:4) n_phi_tv += phi_tv_mask[k];
  int n_phi = tv_phi ? n_time_eff : 0;
  int n_sig = tv_sigma ? n_time_eff : 0;
}

parameters {
  vector[2] mu;
  matrix<lower=-0.99, upper=0.99>[2, 2] Phi;

  real z_rho_init;
  real<lower=0.001> sigma_omega;
  vector[n_time_eff] omega_raw;

  matrix[n_time, 2] zeta;

  vector<lower=0>[2] sigma_eps;
  vector[2] eta;
  vector<lower=0>[2] omega;
  vector<lower=-1, upper=1>[2] delta;
  vector<lower=0>[2] shape_gam;

  vector<lower=0>[n_phi_tv] tau_phi;
  matrix[n_phi, n_phi_tv] phi_raw;
  vector<lower=0>[tv_sigma ? 2 : 0] tau_sigma;
  matrix[n_sig, 2] sigma_raw;
}

transformed parameters {
  vector[n_time_eff] z_rho = compute_rw_ncp(z_rho_init, sigma_omega, omega_raw, n_time_eff);
  vector[n_time_eff] rho;
  matrix[n_phi, 4] phi_dev;
  matrix[n_sig, 2] sigma_dev;
  matrix[n_time, 2] state;

  vector[2] alpha = delta ./ sqrt(1 - square(delta));

  for (t in 1:n_time_eff) rho[t] = 0.97 * tanh(z_rho[t]);

  if (tv_phi == 1) {
    int idx = 0;
    for (k in 1:4) {
      if (phi_tv_mask[k] == 1) {
        idx += 1;
        phi_dev[, k] = compute_rw_ncp(0, tau_phi[idx], phi_raw[, idx], n_time_eff);
      } else {
        phi_dev[, k] = rep_vector(0, n_time_eff);
      }
    }
  }

  if (tv_sigma == 1) {
    for (i in 1:2) {
      sigma_dev[, i] = compute_rw_ncp(0, tau_sigma[i], sigma_raw[, i], n_time_eff);
    }
  }

  {
    vector[2] s;
    for (t in 1:n_time) {
      vector[2] zt = to_vector(zeta[t]);
      if (t == 1) {
        s = mu + zt;
      } else {
        matrix[2, 2] B = Phi;
        if (tv_phi == 1) {
          B[1, 1] += phi_dev[t - 1, 1];
          B[1, 2] += phi_dev[t - 1, 2];
          B[2, 1] += phi_dev[t - 1, 3];
          B[2, 2] += phi_dev[t - 1, 4];
        }
        s = mu + B * s + zt;
      }
      state[t, 1] = s[1];
      state[t, 2] = s[2];
    }
  }
}

model {
  real tiny = 1e-9;
  vector[2] sigma_exp = rep_vector(0, 2);
  vector[2] rate_exp = rep_vector(0, 2);
  vector[2] mean_gam = rep_vector(0, 2);
  vector[2] rate_gam = rep_vector(0, 2);

  mu ~ normal(0, prior_mu_sd);
  to_vector(Phi) ~ normal(0, prior_phi_sd);
  z_rho_init ~ normal(0, prior_rho_sd);
  sigma_omega ~ exponential(1.0 / sigma_omega_prior);
  omega_raw ~ std_normal();

  sigma_eps ~ lognormal(0, prior_sigma_sd);
  omega ~ normal(0, 1);
  delta ~ normal(0, 0.5);
  shape_gam ~ lognormal(log(1), 0.5);

  tau_phi ~ exponential(1.0 / tau_phi_prior);
  to_vector(phi_raw) ~ std_normal();
  tau_sigma ~ exponential(1.0 / tau_sigma_prior);
  to_vector(sigma_raw) ~ std_normal();

  for (i in 1:2) {
    if (family[i] == 2 || family[i] == 4) {
      if (tv_sigma == 1) {
        eta[i] ~ normal(0, prior_sigma_sd);
      } else {
        real m = -skew_direction[i] * zeta[1, i];
        for (t in 2:n_time) m = fmax(m, -skew_direction[i] * zeta[t, i]);
        real lb = fmax(m, 0);
        if (family[i] == 2) {
          sigma_exp[i] = lb + exp(eta[i]) + tiny;
          rate_exp[i] = 1.0 / sigma_exp[i];
          target += lognormal_lpdf(sigma_exp[i] | 0, prior_sigma_sd) + eta[i];
        } else {
          mean_gam[i] = lb + exp(eta[i]) + tiny;
          rate_gam[i] = shape_gam[i] / mean_gam[i];
          target += lognormal_lpdf(mean_gam[i] | 0, prior_sigma_sd) + eta[i];
        }
      }
    } else {
      eta[i] ~ std_normal();
    }
  }

  for (t in 1:n_time) {
    int tt = t == 1 ? 1 : t - 1;
    vector[2] u_vec;
    for (i in 1:2) {
      real dev = (tv_sigma == 1) ? sigma_dev[tt, i] : 0;
      if (family[i] == 1) {
        real s_t = sigma_eps[i] * exp(dev);
        target += normal_lpdf(zeta[t, i] | 0, s_t);
        u_vec[i] = Phi(zeta[t, i] / s_t);
      } else if (family[i] == 2) {
        if (tv_sigma == 1) {
          real m_t = exp(eta[i] + dev);
          real arg = m_t + skew_direction[i] * zeta[t, i];
          real x_shifted = softplus_k(arg, barrier_k);
          target += exponential_lpdf(x_shifted | 1.0 / m_t)
                    + log_softplus_k_jac(arg, barrier_k);
          u_vec[i] = exponential_cdf(x_shifted | 1.0 / m_t);
        } else {
          real x_shifted = sigma_exp[i] + skew_direction[i] * zeta[t, i];
          target += exponential_lpdf(x_shifted | rate_exp[i]);
          u_vec[i] = exponential_cdf(x_shifted | rate_exp[i]);
        }
        if (skew_direction[i] < 0) u_vec[i] = 1.0 - u_vec[i];
      } else if (family[i] == 3) {
        real omega_t = omega[i] * exp(dev);
        real xi_t = -omega_t * (delta[i] * SQRT_2_OVER_PI);
        target += skew_normal_lpdf(zeta[t, i] | xi_t, omega_t, alpha[i]);
        u_vec[i] = skew_normal_cdf(zeta[t, i] | xi_t, omega_t, alpha[i]);
      } else {
        if (tv_sigma == 1) {
          real m_t = exp(eta[i] + dev);
          real arg = m_t + skew_direction[i] * zeta[t, i];
          real x_shifted = softplus_k(arg, barrier_k);
          target += gamma_lpdf(x_shifted | shape_gam[i], shape_gam[i] / m_t)
                    + log_softplus_k_jac(arg, barrier_k);
          u_vec[i] = gamma_cdf(x_shifted | shape_gam[i], shape_gam[i] / m_t);
        } else {
          real x_shifted = mean_gam[i] + skew_direction[i] * zeta[t, i];
          target += gamma_lpdf(x_shifted | shape_gam[i], rate_gam[i]);
          u_vec[i] = gamma_cdf(x_shifted | shape_gam[i], rate_gam[i]);
        }
        if (skew_direction[i] < 0) u_vec[i] = 1.0 - u_vec[i];
      }
    }
    target += gaussian_copula_uv_lpdf(u_vec | rho[tt]);
  }

  for (t in 1:n_time) {
    for (j in 1:J) {
      y[t, j]     ~ normal(lambda[j] * state[t, 1], sigma_e);
      y[t, J + j] ~ normal(lambda[j] * state[t, 2], sigma_e);
    }
  }
}

generated quantities {
  vector[n_time] log_lik;
  vector[2] sigma_exp = rep_vector(0, 2);
  vector[2] sigma_gam = rep_vector(0, 2);
  matrix[n_phi, 4] phi_t;
  matrix[n_sig, 2] sigma_t;
  vector[tv_phi ? n_time_eff : 1] spectral_radius;
  int<lower=0> n_nonstationary = 0;

  {
    real tiny = 1e-9;
    vector[2] rate_exp = rep_vector(0, 2);
    vector[2] mean_gam = rep_vector(0, 2);
    vector[2] rate_gam = rep_vector(0, 2);

    for (i in 1:2) {
      if (family[i] == 2 || family[i] == 4) {
        if (tv_sigma == 1) {
          real m_base = exp(eta[i]);
          if (family[i] == 2) {
            sigma_exp[i] = m_base;
          } else {
            sigma_gam[i] = m_base / sqrt(shape_gam[i]);
          }
        } else {
          real m = -skew_direction[i] * zeta[1, i];
          for (t in 2:n_time) m = fmax(m, -skew_direction[i] * zeta[t, i]);
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
    }

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

    if (tv_sigma == 1) {
      for (i in 1:2) {
        for (t in 1:n_time_eff) {
          if (family[i] == 1) {
            sigma_t[t, i] = sigma_eps[i] * exp(sigma_dev[t, i]);
          } else if (family[i] == 3) {
            sigma_t[t, i] = omega[i] * exp(sigma_dev[t, i]);
          } else if (family[i] == 2) {
            sigma_t[t, i] = exp(eta[i] + sigma_dev[t, i]);
          } else {
            sigma_t[t, i] = exp(eta[i] + sigma_dev[t, i]) / sqrt(shape_gam[i]);
          }
        }
      }
    }

    for (t in 1:n_time) {
      int tt = t == 1 ? 1 : t - 1;
      vector[2] u_vec;
      log_lik[t] = 0;
      for (i in 1:2) {
        real dev = (tv_sigma == 1) ? sigma_dev[tt, i] : 0;
        if (family[i] == 1) {
          real s_t = sigma_eps[i] * exp(dev);
          log_lik[t] += normal_lpdf(zeta[t, i] | 0, s_t);
          u_vec[i] = Phi(zeta[t, i] / s_t);
        } else if (family[i] == 2) {
          if (tv_sigma == 1) {
            real m_t = exp(eta[i] + dev);
            real arg = m_t + skew_direction[i] * zeta[t, i];
            real x_shifted = softplus_k(arg, barrier_k);
            log_lik[t] += exponential_lpdf(x_shifted | 1.0 / m_t)
                          + log_softplus_k_jac(arg, barrier_k);
            u_vec[i] = exponential_cdf(x_shifted | 1.0 / m_t);
          } else {
            real x_shifted = sigma_exp[i] + skew_direction[i] * zeta[t, i];
            log_lik[t] += exponential_lpdf(x_shifted | rate_exp[i]);
            u_vec[i] = exponential_cdf(x_shifted | rate_exp[i]);
          }
          if (skew_direction[i] < 0) u_vec[i] = 1.0 - u_vec[i];
        } else if (family[i] == 3) {
          real omega_ti = omega[i] * exp(dev);
          real xi_ti = -omega_ti * (delta[i] * SQRT_2_OVER_PI);
          log_lik[t] += skew_normal_lpdf(zeta[t, i] | xi_ti, omega_ti, alpha[i]);
          u_vec[i] = skew_normal_cdf(zeta[t, i] | xi_ti, omega_ti, alpha[i]);
        } else {
          if (tv_sigma == 1) {
            real m_t = exp(eta[i] + dev);
            real arg = m_t + skew_direction[i] * zeta[t, i];
            real x_shifted = softplus_k(arg, barrier_k);
            log_lik[t] += gamma_lpdf(x_shifted | shape_gam[i], shape_gam[i] / m_t)
                          + log_softplus_k_jac(arg, barrier_k);
            u_vec[i] = gamma_cdf(x_shifted | shape_gam[i], shape_gam[i] / m_t);
          } else {
            real x_shifted = mean_gam[i] + skew_direction[i] * zeta[t, i];
            log_lik[t] += gamma_lpdf(x_shifted | shape_gam[i], rate_gam[i]);
            u_vec[i] = gamma_cdf(x_shifted | shape_gam[i], rate_gam[i]);
          }
          if (skew_direction[i] < 0) u_vec[i] = 1.0 - u_vec[i];
        }
      }
      log_lik[t] += gaussian_copula_uv_lpdf(u_vec | rho[tt]);
      for (j in 1:J) {
        log_lik[t] += normal_lpdf(y[t, j] | lambda[j] * state[t, 1], sigma_e);
        log_lik[t] += normal_lpdf(y[t, J + j] | lambda[j] * state[t, 2], sigma_e);
      }
    }
  }
}
