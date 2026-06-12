// Bivariate Clayton copula log-density.
// Inputs are clamped away from 0 and 1 for numerical stability.
real clayton_copula_ld(real u, real v, real theta) {
  real eps = 1e-9;
  real uu = fmax(eps, fmin(1 - eps, u));
  real vv = fmax(eps, fmin(1 - eps, v));
  real log_u = log(uu);
  real log_v = log(vv);
  // log(u^-theta + v^-theta - 1) in log space: pow(uu, -theta) overflows to
  // +inf (with NaN gradients) for theta > ~34 at the lower clamp boundary,
  // which theta's unbounded prior can reach during warmup. The sum exceeds 2
  // because u, v < 1, so log_diff_exp against 0 is always well defined.
  real log_sum = log_sum_exp(-theta * log_u, -theta * log_v);
  real log_s = log_diff_exp(log_sum, 0);
  return log1p(theta)
         - (1 + theta) * (log_u + log_v)
         - (2 + 1 / theta) * log_s;
}

// z-score parameterization for normal margins.
real clayton_copula_z_lpdf(vector z, real theta) {
  return clayton_copula_ld(Phi(z[1]), Phi(z[2]), theta);
}
