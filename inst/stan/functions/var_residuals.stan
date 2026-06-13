// Compute VAR(1) residuals: eps[t] = Y[t+1] - (mu + Phi * (Y[t] - mu))
matrix compute_var_residuals(matrix Y, vector mu, matrix Phi, int n_time_eff, int D) {
  matrix[n_time_eff, D] eps;
  for (t in 1:n_time_eff) {
    vector[D] y_prev = to_vector(Y[t, ]);
    vector[D] y_curr = to_vector(Y[t + 1, ]);
    vector[D] y_hat = mu + Phi * (y_prev - mu);
    eps[t, ] = to_row_vector(y_curr - y_hat);
  }
  return eps;
}

// Time-varying-Phi residuals: eps[t] = Y[t+1] - (mu + Phi(t) * (Y[t] - mu)).
// phi_t columns are row-major (phi11, phi12, phi21, phi22); bivariate only.
matrix compute_var_residuals_tv(matrix Y, vector mu, matrix phi_t, int n_time_eff, int D) {
  matrix[n_time_eff, D] eps;
  for (t in 1:n_time_eff) {
    matrix[2, 2] P;
    P[1, 1] = phi_t[t, 1];
    P[1, 2] = phi_t[t, 2];
    P[2, 1] = phi_t[t, 3];
    P[2, 2] = phi_t[t, 4];
    vector[D] y_prev = to_vector(Y[t, ]);
    vector[D] y_curr = to_vector(Y[t + 1, ]);
    vector[D] y_hat = mu + P * (y_prev - mu);
    eps[t, ] = to_row_vector(y_curr - y_hat);
  }
  return eps;
}
