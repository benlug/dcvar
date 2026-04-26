// Bivariate Clayton copula log-density.
// Inputs are clamped away from 0 and 1 for numerical stability.
real clayton_copula_ld(real u, real v, real theta) {
  real eps = 1e-9;
  real uu = fmax(eps, fmin(1 - eps, u));
  real vv = fmax(eps, fmin(1 - eps, v));
  real up = pow(uu, -theta);
  real vp = pow(vv, -theta);
  real s = up + vp - 1;
  return log1p(theta)
         - (1 + theta) * (log(uu) + log(vv))
         - (2 + 1 / theta) * log(s);
}

// z-score parameterization for normal margins.
real clayton_copula_z_lpdf(vector z, real theta) {
  return clayton_copula_ld(Phi(z[1]), Phi(z[2]), theta);
}
