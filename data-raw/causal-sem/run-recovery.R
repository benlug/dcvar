source("data-raw/causal-sem/validation-helpers.R")
validation_setup()
config <- validation_config("recovery")
results <- list()
index <- 0L
cases <- strsplit(Sys.getenv("DCVAR_RECOVERY_CASES", "nonzero,null"), ",", fixed = TRUE)[[1]]
if (any(!cases %in% c("nonzero", "null"))) stop("Use nonzero or null recovery cases.")
for (model in config$models) {
  for (effect_case in cases) {
    for (replication in seq_len(config$reps)) {
      index <- index + 1L
      parameters <- list()
      if (effect_case == "null") {
        baseline <- dcvar::simulate_dcvar_causal_sem(n = config$n, model = model,
                                                     seed = config$seed)$parameters
        if (model %in% c("latent_covariate", "latent_outcome")) {
          baseline$mu_v[2] <- baseline$mu_v[1]
          baseline$alpha[2] <- baseline$alpha[1]
          baseline$beta[2] <- baseline$beta[1]
        } else {
          baseline$mu_q[2, ] <- baseline$mu_q[1, ]
          baseline$alpha_m[2] <- baseline$alpha_m[1]
          baseline$alpha_y[2] <- baseline$alpha_y[1]
          baseline$B[2, ] <- baseline$B[1, ]
          baseline$D[2, ] <- baseline$D[1, ]
          baseline$d[2] <- baseline$d[1]
        }
        parameters <- baseline
      }
      sim <- dcvar::simulate_dcvar_causal_sem(n = config$n, model = model,
                    parameters = parameters, seed = config$seed + index)
      result <- validation_fit(sim, config, seed = config$seed + 100000L + index)
      id <- paste0(model, "-", effect_case, "-rep", replication)
      results[[index]] <- validation_record(result, id, config)
    }
  }
}
validation_report(results, config)
