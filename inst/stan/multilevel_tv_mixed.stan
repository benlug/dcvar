// Time-varying multilevel Copula VAR(1) with mixed margins.
//
// Unit-specific Phi baselines are retained from multilevel_mixed.stan. Optional
// time variation is a shared population drift around those unit baselines, not a
// separate per-unit random walk, which avoids stacking two weakly identified
// sources of unit-level Phi variation.

functions {
#include functions/gaussian_copula_uv.stan
#include functions/ncp_random_walk.stan
#include functions/softplus.stan
}

data {
  int<lower=1> N;
  int<lower=2> n_time;
  array[N] matrix[n_time, 2] y;
  array[2] int<lower=1, upper=4> family;
  vector[2] skew_direction;

  int<lower=0, upper=1> tv_phi;
  array[4] int<lower=0, upper=1> phi_tv_mask;
  int<lower=0, upper=1> tv_sigma;

  real<lower=0> prior_phi_bar_sd;
  real<lower=0> prior_tau_phi_scale;
  real<lower=0> prior_sigma_sd;
  real<lower=0> prior_rho_sd;
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
  vector[4] phi_bar;
  vector<lower=0>[4] tau_phi;
  matrix[N, 4] z_phi;
  real<lower=-1, upper=1> rho;

  vector<lower=0.01>[2] sigma_eps;
  vector[2] eta;
  vector<lower=0>[2] omega;
  vector<lower=-1, upper=1>[2] delta;
  vector<lower=0>[2] shape_gam;

  vector<lower=0>[n_phi_tv] tau_phi_tv;
  matrix[n_phi, n_phi_tv] phi_raw;
  vector<lower=0>[tv_sigma ? 2 : 0] tau_sigma;
  matrix[n_sig, 2] sigma_raw;
}

transformed parameters {
  matrix[N, 4] phi_unit;
  matrix[n_phi, 4] phi_dev;
  matrix[n_sig, 2] sigma_dev;

  vector[2] alpha = delta ./ sqrt(1 - square(delta));

  for (i in 1:N) {
    for (k in 1:4) phi_unit[i, k] = phi_bar[k] + tau_phi[k] * z_phi[i, k];
  }

  if (tv_phi == 1) {
    int idx = 0;
    for (k in 1:4) {
      if (phi_tv_mask[k] == 1) {
        idx += 1;
        phi_dev[, k] = compute_rw_ncp(0, tau_phi_tv[idx], phi_raw[, idx], n_time_eff);
      } else {
        phi_dev[, k] = rep_vector(0, n_time_eff);
      }
    }
  }

  if (tv_sigma == 1) {
    for (j in 1:2) {
      sigma_dev[, j] = compute_rw_ncp(0, tau_sigma[j], sigma_raw[, j], n_time_eff);
    }
  }
}

model {
  vector[2] sigma_exp = rep_vector(0, 2);
  vector[2] rate_exp = rep_vector(0, 2);
  vector[2] mean_gam = rep_vector(0, 2);
  vector[2] rate_gam = rep_vector(0, 2);
  real tiny = 1e-9;

  phi_bar ~ normal(0, prior_phi_bar_sd);
  tau_phi ~ student_t(3, 0, prior_tau_phi_scale);
  to_vector(z_phi) ~ std_normal();
  rho ~ normal(0, prior_rho_sd);

  sigma_eps ~ normal(0, prior_sigma_sd);
  omega ~ normal(0, 1);
  delta ~ normal(0, 0.5);
  shape_gam ~ lognormal(log(1), 0.5);

  tau_phi_tv ~ exponential(1.0 / tau_phi_prior);
  to_vector(phi_raw) ~ std_normal();
  tau_sigma ~ exponential(1.0 / tau_sigma_prior);
  to_vector(sigma_raw) ~ std_normal();

  for (j in 1:2) {
    if (family[j] == 2 || family[j] == 4) {
      if (tv_sigma == 1) {
        eta[j] ~ normal(0, prior_sigma_sd);
      } else {
        real m = negative_infinity();
        for (i in 1:N) {
          for (t in 2:n_time) {
            matrix[2, 2] Phi_T;
            row_vector[2] res;
            int tt = t - 1;
            Phi_T[1, 1] = phi_unit[i, 1] + (tv_phi == 1 ? phi_dev[tt, 1] : 0);
            Phi_T[2, 1] = phi_unit[i, 2] + (tv_phi == 1 ? phi_dev[tt, 2] : 0);
            Phi_T[1, 2] = phi_unit[i, 3] + (tv_phi == 1 ? phi_dev[tt, 3] : 0);
            Phi_T[2, 2] = phi_unit[i, 4] + (tv_phi == 1 ? phi_dev[tt, 4] : 0);
            res = y[i][t] - y[i][t - 1] * Phi_T;
            m = fmax(m, -skew_direction[j] * res[j]);
          }
        }
        real lb = fmax(m, 0);
        if (family[j] == 2) {
          sigma_exp[j] = lb + exp(eta[j]) + tiny;
          rate_exp[j] = 1.0 / sigma_exp[j];
          target += lognormal_lpdf(sigma_exp[j] | 0, prior_sigma_sd) + eta[j];
        } else {
          mean_gam[j] = lb + exp(eta[j]) + tiny;
          rate_gam[j] = shape_gam[j] / mean_gam[j];
          target += lognormal_lpdf(mean_gam[j] | 0, prior_sigma_sd) + eta[j];
        }
      }
    } else {
      eta[j] ~ std_normal();
    }
  }

  for (i in 1:N) {
    for (t in 2:n_time) {
      int tt = t - 1;
      matrix[2, 2] Phi_T;
      row_vector[2] res;
      vector[2] u;
      Phi_T[1, 1] = phi_unit[i, 1] + (tv_phi == 1 ? phi_dev[tt, 1] : 0);
      Phi_T[2, 1] = phi_unit[i, 2] + (tv_phi == 1 ? phi_dev[tt, 2] : 0);
      Phi_T[1, 2] = phi_unit[i, 3] + (tv_phi == 1 ? phi_dev[tt, 3] : 0);
      Phi_T[2, 2] = phi_unit[i, 4] + (tv_phi == 1 ? phi_dev[tt, 4] : 0);
      res = y[i][t] - y[i][t - 1] * Phi_T;

      for (j in 1:2) {
        real dev = (tv_sigma == 1) ? sigma_dev[tt, j] : 0;
        if (family[j] == 1) {
          real s_t = sigma_eps[j] * exp(dev);
          target += normal_lpdf(res[j] | 0, s_t);
          u[j] = Phi(res[j] / s_t);
        } else if (family[j] == 2) {
          if (tv_sigma == 1) {
            real m_t = exp(eta[j] + dev);
            real arg = m_t + skew_direction[j] * res[j];
            real x_shifted = softplus_k(arg, barrier_k);
            target += exponential_lpdf(x_shifted | 1.0 / m_t)
                      + log_softplus_k_jac(arg, barrier_k);
            u[j] = exponential_cdf(x_shifted | 1.0 / m_t);
          } else {
            real x_shifted = sigma_exp[j] + skew_direction[j] * res[j];
            target += exponential_lpdf(x_shifted | rate_exp[j]);
            u[j] = exponential_cdf(x_shifted | rate_exp[j]);
          }
          if (skew_direction[j] < 0) u[j] = 1.0 - u[j];
        } else if (family[j] == 3) {
          real omega_t = omega[j] * exp(dev);
          real xi_t = -omega_t * (delta[j] * SQRT_2_OVER_PI);
          target += skew_normal_lpdf(res[j] | xi_t, omega_t, alpha[j]);
          u[j] = skew_normal_cdf(res[j] | xi_t, omega_t, alpha[j]);
        } else {
          if (tv_sigma == 1) {
            real m_t = exp(eta[j] + dev);
            real arg = m_t + skew_direction[j] * res[j];
            real x_shifted = softplus_k(arg, barrier_k);
            target += gamma_lpdf(x_shifted | shape_gam[j], shape_gam[j] / m_t)
                      + log_softplus_k_jac(arg, barrier_k);
            u[j] = gamma_cdf(x_shifted | shape_gam[j], shape_gam[j] / m_t);
          } else {
            real x_shifted = mean_gam[j] + skew_direction[j] * res[j];
            target += gamma_lpdf(x_shifted | shape_gam[j], rate_gam[j]);
            u[j] = gamma_cdf(x_shifted | shape_gam[j], rate_gam[j]);
          }
          if (skew_direction[j] < 0) u[j] = 1.0 - u[j];
        }
      }
      target += gaussian_copula_uv_lpdf(u | rho);
    }
  }
}

generated quantities {
  array[N] matrix[2, 2] Phi;
  vector[2] sigma_exp = rep_vector(0, 2);
  vector[2] sigma_gam = rep_vector(0, 2);
  matrix[n_phi, 4] phi_t;
  array[N] matrix[n_phi, 4] phi_unit_t;
  matrix[n_sig, 2] sigma_t;
  matrix[N, tv_phi ? n_time_eff : 1] spectral_radius;
  int<lower=0> n_nonstationary = 0;
  matrix[N, n_time_eff] log_lik;

  {
    vector[2] rate_exp = rep_vector(0, 2);
    vector[2] mean_gam = rep_vector(0, 2);
    vector[2] rate_gam = rep_vector(0, 2);
    real tiny = 1e-9;

    for (j in 1:2) {
      if (family[j] == 2 || family[j] == 4) {
        if (tv_sigma == 1) {
          real m_base = exp(eta[j]);
          if (family[j] == 2) {
            sigma_exp[j] = m_base;
          } else {
            sigma_gam[j] = m_base / sqrt(shape_gam[j]);
          }
        } else {
          real m = negative_infinity();
          for (i in 1:N) {
            for (t in 2:n_time) {
              int tt = t - 1;
              matrix[2, 2] Phi_T;
              row_vector[2] res;
              Phi_T[1, 1] = phi_unit[i, 1] + (tv_phi == 1 ? phi_dev[tt, 1] : 0);
              Phi_T[2, 1] = phi_unit[i, 2] + (tv_phi == 1 ? phi_dev[tt, 2] : 0);
              Phi_T[1, 2] = phi_unit[i, 3] + (tv_phi == 1 ? phi_dev[tt, 3] : 0);
              Phi_T[2, 2] = phi_unit[i, 4] + (tv_phi == 1 ? phi_dev[tt, 4] : 0);
              res = y[i][t] - y[i][t - 1] * Phi_T;
              m = fmax(m, -skew_direction[j] * res[j]);
            }
          }
          real lb = fmax(m, 0);
          if (family[j] == 2) {
            sigma_exp[j] = lb + exp(eta[j]) + tiny;
            rate_exp[j] = 1.0 / sigma_exp[j];
          } else {
            mean_gam[j] = lb + exp(eta[j]) + tiny;
            sigma_gam[j] = mean_gam[j] / sqrt(shape_gam[j]);
            rate_gam[j] = shape_gam[j] / mean_gam[j];
          }
        }
      }
    }

    for (i in 1:N) {
      Phi[i][1, 1] = phi_unit[i, 1];
      Phi[i][1, 2] = phi_unit[i, 2];
      Phi[i][2, 1] = phi_unit[i, 3];
      Phi[i][2, 2] = phi_unit[i, 4];
    }

    if (tv_phi == 1) {
      for (t in 1:n_time_eff) {
        phi_t[t, 1] = phi_bar[1] + phi_dev[t, 1];
        phi_t[t, 2] = phi_bar[2] + phi_dev[t, 2];
        phi_t[t, 3] = phi_bar[3] + phi_dev[t, 3];
        phi_t[t, 4] = phi_bar[4] + phi_dev[t, 4];
        for (i in 1:N) {
          for (k in 1:4) phi_unit_t[i][t, k] = phi_unit[i, k] + phi_dev[t, k];
        }
      }
    }

    if (tv_sigma == 1) {
      for (j in 1:2) {
        for (t in 1:n_time_eff) {
          if (family[j] == 1) {
            sigma_t[t, j] = sigma_eps[j] * exp(sigma_dev[t, j]);
          } else if (family[j] == 3) {
            sigma_t[t, j] = omega[j] * exp(sigma_dev[t, j]) *
              sqrt(1 - 2 * square(delta[j]) / pi());
          } else if (family[j] == 2) {
            sigma_t[t, j] = exp(eta[j] + sigma_dev[t, j]);
          } else {
            sigma_t[t, j] = exp(eta[j] + sigma_dev[t, j]) / sqrt(shape_gam[j]);
          }
        }
      }
    }

    // n_nonstationary counts nonstationary UNITS (max N), invariant to tv_phi:
    // a unit is flagged if its spectral radius reaches 1 at any time point.
    for (i in 1:N) {
      int sr_len = tv_phi ? n_time_eff : 1;
      int unit_nonstationary = 0;
      for (s in 1:sr_len) {
        int tt = tv_phi ? s : 1;
        real p11 = phi_unit[i, 1] + (tv_phi == 1 ? phi_dev[tt, 1] : 0);
        real p12 = phi_unit[i, 2] + (tv_phi == 1 ? phi_dev[tt, 2] : 0);
        real p21 = phi_unit[i, 3] + (tv_phi == 1 ? phi_dev[tt, 3] : 0);
        real p22 = phi_unit[i, 4] + (tv_phi == 1 ? phi_dev[tt, 4] : 0);
        real tr = p11 + p22;
        real det_phi = p11 * p22 - p12 * p21;
        real disc = square(tr) - 4 * det_phi;
        if (disc >= 0) {
          real sqrt_disc = sqrt(disc);
          spectral_radius[i, s] = fmax(abs(0.5 * (tr + sqrt_disc)), abs(0.5 * (tr - sqrt_disc)));
        } else {
          spectral_radius[i, s] = sqrt(abs(det_phi));
        }
        if (spectral_radius[i, s] >= 1.0) unit_nonstationary = 1;
      }
      n_nonstationary += unit_nonstationary;
    }

    for (i in 1:N) {
      for (t in 2:n_time) {
        int tt = t - 1;
        matrix[2, 2] Phi_T;
        row_vector[2] res;
        vector[2] u;
        Phi_T[1, 1] = phi_unit[i, 1] + (tv_phi == 1 ? phi_dev[tt, 1] : 0);
        Phi_T[2, 1] = phi_unit[i, 2] + (tv_phi == 1 ? phi_dev[tt, 2] : 0);
        Phi_T[1, 2] = phi_unit[i, 3] + (tv_phi == 1 ? phi_dev[tt, 3] : 0);
        Phi_T[2, 2] = phi_unit[i, 4] + (tv_phi == 1 ? phi_dev[tt, 4] : 0);
        res = y[i][t] - y[i][t - 1] * Phi_T;
        log_lik[i, tt] = 0;
        for (j in 1:2) {
          real dev = (tv_sigma == 1) ? sigma_dev[tt, j] : 0;
          if (family[j] == 1) {
            real s_t = sigma_eps[j] * exp(dev);
            log_lik[i, tt] += normal_lpdf(res[j] | 0, s_t);
            u[j] = Phi(res[j] / s_t);
          } else if (family[j] == 2) {
            if (tv_sigma == 1) {
              real m_t = exp(eta[j] + dev);
              real arg = m_t + skew_direction[j] * res[j];
              real x_shifted = softplus_k(arg, barrier_k);
              log_lik[i, tt] += exponential_lpdf(x_shifted | 1.0 / m_t)
                                + log_softplus_k_jac(arg, barrier_k);
              u[j] = exponential_cdf(x_shifted | 1.0 / m_t);
            } else {
              real x_shifted = sigma_exp[j] + skew_direction[j] * res[j];
              log_lik[i, tt] += exponential_lpdf(x_shifted | rate_exp[j]);
              u[j] = exponential_cdf(x_shifted | rate_exp[j]);
            }
            if (skew_direction[j] < 0) u[j] = 1.0 - u[j];
          } else if (family[j] == 3) {
            real omega_ti = omega[j] * exp(dev);
            real xi_ti = -omega_ti * (delta[j] * SQRT_2_OVER_PI);
            log_lik[i, tt] += skew_normal_lpdf(res[j] | xi_ti, omega_ti, alpha[j]);
            u[j] = skew_normal_cdf(res[j] | xi_ti, omega_ti, alpha[j]);
          } else {
            if (tv_sigma == 1) {
              real m_t = exp(eta[j] + dev);
              real arg = m_t + skew_direction[j] * res[j];
              real x_shifted = softplus_k(arg, barrier_k);
              log_lik[i, tt] += gamma_lpdf(x_shifted | shape_gam[j], shape_gam[j] / m_t)
                                + log_softplus_k_jac(arg, barrier_k);
              u[j] = gamma_cdf(x_shifted | shape_gam[j], shape_gam[j] / m_t);
            } else {
              real x_shifted = mean_gam[j] + skew_direction[j] * res[j];
              log_lik[i, tt] += gamma_lpdf(x_shifted | shape_gam[j], rate_gam[j]);
              u[j] = gamma_cdf(x_shifted | shape_gam[j], rate_gam[j]);
            }
            if (skew_direction[j] < 0) u[j] = 1.0 - u[j];
          }
        }
        log_lik[i, tt] += gaussian_copula_uv_lpdf(u | rho);
      }
    }
  }
}
