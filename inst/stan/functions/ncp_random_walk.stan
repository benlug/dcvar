// Non-centered parameterization for random walk on Fisher-z scale
// z_rho[t] = z_rho_init + sigma_omega * cumsum(omega_raw[1:t])
vector compute_z_rho_ncp(real z_rho_init, real sigma_omega, vector omega_raw, int T_eff) {
  vector[T_eff] z_rho;
  real cumsum_omega = 0;
  for (t in 1:T_eff) {
    cumsum_omega += omega_raw[t];
    z_rho[t] = z_rho_init + sigma_omega * cumsum_omega;
  }
  return z_rho;
}
