source("data-raw/causal-sem/validation-helpers.R")
validation_setup()
config <- validation_config("priors", default_reps = 100L)
results <- list()
index <- 0L
for (model in config$models) {
  for (replication in seq_len(config$reps)) {
    index <- index + 1L
    started <- proc.time()[["elapsed"]]
    sim <- dcvar::simulate_dcvar_causal_sem(n = config$n, model = model,
                        parameter_prior = TRUE, seed = config$seed + index)
    indicators <- unlist(sim$args$indicators, use.names = FALSE)
    frequencies <- lapply(sim$data[indicators], function(x) prop.table(table(x)))
    continuous <- names(sim$data)[vapply(sim$data, is.numeric, logical(1))]
    continuous <- setdiff(continuous, sim$args$treatment)
    result <- list(status = "ok", effects = sim$effects, frequencies = frequencies,
                   constant = vapply(sim$data[indicators], function(x) length(unique(x)) < 2L,
                                     logical(1)),
                   quantiles = lapply(sim$data[continuous], stats::quantile,
                                      probs = c(0.01, 0.5, 0.99)),
                   parameters = sim$parameters,
                   seconds = proc.time()[["elapsed"]] - started)
    results[[index]] <- validation_record(result, paste0(model, "-rep", replication), config)
  }
}
summary <- do.call(rbind, lapply(config$models, function(model) {
  rows <- results[startsWith(vapply(results, `[[`, character(1), "id"), paste0(model, "-rep"))]
  data.frame(model = model, simulations = length(rows),
               fraction_constant = mean(vapply(rows, function(x) any(x$constant), logical(1))))
}))
utils::write.csv(summary, file.path(config$output, "summary.csv"), row.names = FALSE)
saveRDS(list(config = config, summary = summary, session = utils::sessionInfo()),
        file.path(config$output, "run.rds"))
print(summary)
