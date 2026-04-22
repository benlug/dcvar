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
