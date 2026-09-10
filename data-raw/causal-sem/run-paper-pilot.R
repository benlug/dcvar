source("data-raw/causal-sem/validation-helpers.R")
source("data-raw/causal-sem/paper-vita.R")
validation_setup()
config <- validation_config("paper-pilot", default_reps = 1L)
if (config$n %% 2L != 0L) stop("The balanced paper pilot requires an even total sample size.")
results <- list()
index <- 0L
generate_only <- identical(Sys.getenv("DCVAR_PAPER_GENERATE_ONLY"), "true")
Nmax <- as.integer(Sys.getenv("DCVAR_PAPER_NMAX", "100000"))
paper_root <- Sys.getenv("DCVAR_PAPER_ROOT", "../sem_causal_effects")
families <- strsplit(Sys.getenv("DCVAR_PAPER_COPULAS", "gauss,clayton,joe"), ",", fixed = TRUE)[[1]]
for (model in config$models) {
  for (family in families) {
    for (replication in seq_len(config$reps)) {
      index <- index + 1L
      condition <- paper_condition(model, n = as.integer(config$n / 2), copula = rep(family, 2))
      started <- proc.time()[["elapsed"]]
      result <- tryCatch({
        sim <- paper_vita_sample(condition, seed = config$seed + index,
                                  paper_root = paper_root, Nmax = Nmax)
        if (generate_only) {
          fitted <- list(status = "generated", effects = NULL,
                          seconds = proc.time()[["elapsed"]] - started)
        } else {
          fitted <- validation_fit(sim, config, seed = config$seed + 100000L + index)
        }
        fitted$generation_seconds <- proc.time()[["elapsed"]] - started -
          if (generate_only) 0 else fitted$seconds
        fitted$paper <- sim$metadata
        fitted$data <- sim$data
        fitted$covariance <- sim$covariance
        fitted$calibration <- sim$calibration
        fitted
      }, error = function(error) list(status = "error", error = conditionMessage(error),
                                       seconds = proc.time()[["elapsed"]] - started))
      id <- paste0(model, "-", family, "-rep", replication)
      results[[index]] <- validation_record(result, id, config)
    }
  }
}
validation_report(results, config)
