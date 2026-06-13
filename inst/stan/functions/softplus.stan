// Soft-barrier helpers for time-varying scales of the shifted exponential /
// gamma margins. The shifted variate x = scale + skew * eps has a hard
// support boundary at 0 that, with a time-varying scale, would couple the
// scale pointwise to the residual. Replacing the affine shift with a softplus
// transform, x = softplus_k(scale + skew * eps), keeps x strictly positive,
// matches the affine shift in the interior (argument >> 0), and rounds the
// boundary smoothly so the residual has unbounded support and the scale can
// vary freely. k controls sharpness (k -> infinity recovers the hard ReLU
// boundary). See the dcvar_tv_mixed.stan header.

// softplus_k(a) = log(1 + exp(k a)) / k  (numerically stable via log1p_exp)
real softplus_k(real a, real k) {
  return log1p_exp(k * a) / k;
}

// Inverse, defined for y > 0: a such that softplus_k(a) = y.
real inv_softplus_k(real y, real k) {
  return log_diff_exp(k * y, 0) / k;
}

// log of the transform Jacobian |d softplus_k(a)/da| = inv_logit(k a)
real log_softplus_k_jac(real a, real k) {
  return log_inv_logit(k * a);
}
