// The three items share one factor. Each item has its own ordered thresholds.
real causal_measurement_lpmf(array[] int response, real factor,
                            vector lambda, row_vector item_sd,
                            vector threshold_1, vector threshold_2,
                            vector threshold_3) {
  return ordered_probit_lpmf(response[1] | lambda[1] * factor / item_sd[1],
                            threshold_1 / item_sd[1])
       + ordered_probit_lpmf(response[2] | lambda[2] * factor / item_sd[2],
                            threshold_2 / item_sd[2])
       + ordered_probit_lpmf(response[3] | lambda[3] * factor / item_sd[3],
                            threshold_3 / item_sd[3]);
}

// A common factor draw preserves the dependence between replicated items.
array[] int causal_measurement_rng(real factor, vector lambda,
                                  row_vector item_sd, vector threshold_1,
                                  vector threshold_2, vector threshold_3) {
  array[3] int response;
  response[1] = ordered_probit_rng(lambda[1] * factor / item_sd[1],
                                  threshold_1 / item_sd[1]);
  response[2] = ordered_probit_rng(lambda[2] * factor / item_sd[2],
                                  threshold_2 / item_sd[2]);
  response[3] = ordered_probit_rng(lambda[3] * factor / item_sd[3],
                                  threshold_3 / item_sd[3]);
  return response;
}
