// Generic scalar non-centered random walk:
// z[t] = init + tau * cumsum(raw[1:t])
// Package-wide convention: the first state already carries an innovation.
vector compute_rw_ncp(real init, real tau, vector raw, int n) {
  vector[n] z;
  real cs = 0;
  for (t in 1:n) {
    cs += raw[t];
    z[t] = init + tau * cs;
  }
  return z;
}

// Non-centered parameterization for random walk on Fisher-z scale
// (backward-compatible wrapper used by the existing dcvar_* models)
vector compute_z_rho_ncp(real z_rho_init, real sigma_omega, vector omega_raw, int n_time_eff) {
  return compute_rw_ncp(z_rho_init, sigma_omega, omega_raw, n_time_eff);
}
