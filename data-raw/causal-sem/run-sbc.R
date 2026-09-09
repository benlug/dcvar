source("data-raw/causal-sem/validation-helpers.R")
validation_setup()
config <- validation_config("sbc")
results <- list()
index <- 0L
for (model in config$models) {
  for (replication in seq_len(config$reps)) {
    index <- index + 1L
    set.seed(config$seed + index)
    # The generator uses exactly the fit priors on the fixed internal scale.
    sim <- dcvar::simulate_dcvar_causal_sem(n = config$n, model = model,
                      parameter_prior = TRUE, seed = config$seed + index)
    result <- validation_fit(sim, config, seed = config$seed + 100000L + index)
    if (identical(result$status, "ok")) {
      parameter_matrix <- posterior::as_draws_matrix(result$parameter_draws)
      effect_matrix <- posterior::as_draws_matrix(result$effect_draws)
      combined <- cbind(parameter_matrix, effect_matrix)
      truth <- c(validation_flatten(sim$parameters), sim$effects)
      wanted <- min(config$ranks, nrow(combined))
      # Spread retained draws over each chain. Review ESS before reading ranks.
      keep <- unique(round(seq(1, nrow(combined), length.out = wanted)))
      result$ranks <- data.frame(variable = colnames(combined),
        rank = vapply(colnames(combined), function(variable) {
          sum(combined[keep, variable] < truth[[variable]])
        }, integer(1)), posterior_draws = length(keep),
        truth = unname(truth[colnames(combined)]))
      result$ranks <- result$ranks[is.finite(result$ranks$truth), ]
    } else {
      # Do not replace failed or invalid samples. They define the denominator.
      result$constant_indicators <- names(sim$data)[vapply(sim$data, function(value) {
        is.ordered(value) && length(unique(value)) < 2L
      }, logical(1))]
    }
    result$generating_parameters <- sim$parameters
    id <- paste0(model, "-prior-rep", replication)
    results[[index]] <- validation_record(result, id, config)
  }
}
validation_report(results, config)
validation_rank_report(results, config)
