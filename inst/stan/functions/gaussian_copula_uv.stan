// Bivariate Gaussian copula log-density (u,v parameterization)
// For non-Gaussian margins where we have CDF values u, v in (0,1)
// rather than z-scores directly.
// Clamps u, v away from boundaries to avoid numerical issues with inv_Phi.
real gaussian_copula_uv_lpdf(vector uv, real rho) {
  real eps = 1e-9;
  real uu = fmax(eps, fmin(1 - eps, uv[1]));
  real vv = fmax(eps, fmin(1 - eps, uv[2]));
  real z1 = inv_Phi(uu);
  real z2 = inv_Phi(vv);
  real rho_sq = square(rho);
  real one_m_rho_sq = 1.0 - rho_sq;
  return -0.5 * log1m(rho_sq)
         - 0.5 / one_m_rho_sq * (square(z1) - 2.0 * rho * z1 * z2 + square(z2))
         + 0.5 * (square(z1) + square(z2));
}
