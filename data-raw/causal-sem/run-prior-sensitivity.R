source("data-raw/causal-sem/validation-helpers.R")
validation_setup()
config <- validation_config("prior-sensitivity", default_reps = 5L)
prior_sets <- list(default = list(), regularized = list(beta_sd = 0.5),
                    wide = list(beta_sd = 2))
results <- list()
index <- 0L
for (model in config$models) {
  for (replication in seq_len(config$reps)) {
    sim <- dcvar::simulate_dcvar_causal_sem(n = config$n, model = model,
                                           seed = config$seed + replication)
    for (prior_name in names(prior_sets)) {
      index <- index + 1L
      result <- validation_fit(sim, config, seed = config$seed + 100000L + index,
                                priors = prior_sets[[prior_name]])
      result$prior <- prior_sets[[prior_name]]
      id <- paste0(model, "-", prior_name, "-rep", replication)
      results[[index]] <- validation_record(result, id, config)
    }
  }
}
validation_report(results, config)
