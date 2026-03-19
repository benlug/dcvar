// Bivariate Gaussian copula log-density (z-score parameterization)
// z1, z2 are standardized residuals (equivalent to inv_Phi(u) but computed directly,
// avoiding the Phi -> clamp -> inv_Phi roundtrip that truncates tail information)
// rho is the copula correlation parameter
real gaussian_copula_z_lpdf(vector z, real rho) {
  real rho_sq = square(rho);
  real one_m_rho_sq = 1.0 - rho_sq;
  real numerator = rho_sq * (square(z[1]) + square(z[2])) - 2.0 * rho * z[1] * z[2];
  return -0.5 * log1m(rho_sq) - numerator / (2.0 * one_m_rho_sq);
}
