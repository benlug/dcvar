source("data-raw/causal-sem/validation-helpers.R")
source("data-raw/causal-sem/paper-vita.R")
validation_setup()
config <- validation_config("paper-fits", default_reps = 1L)
input <- Sys.getenv("DCVAR_PAPER_INPUT", "")
if (!nzchar(input)) stop("Set DCVAR_PAPER_INPUT to a generation result directory.")
paths <- list.files(input, pattern = "-rep[0-9]+[.]rds$", full.names = TRUE)
if (!length(paths)) stop("The input has no paper datasets.")
results <- list()
resume <- identical(Sys.getenv("DCVAR_PAPER_RESUME"), "true")
for (index in seq_along(paths)) {
  original <- readRDS(paths[index])
  model <- sub("-(gauss|clayton|joe)-rep[0-9]+$", "", original$id)
  if (!model %in% config$models) next
  existing <- file.path(config$output, paste0(original$id, ".rds"))
  if (resume && file.exists(existing)) {
    results[[length(results) + 1L]] <- readRDS(existing)
    message(original$id, ": use the saved fit")
    next
  }
  if (is.null(original$data)) {
    result <- list(status = "generation_error", error = original$error,
                    seconds = 0, generation_seconds = original$seconds,
                    source_result = paths[index])
  } else {
    stopifnot(identical(model, original$paper$condition$model))
    sim <- list(data = original$data, parameters = list(),
                  args = paper_fit_args(original$data, model),
                  effects = original$paper$identified_reference_effects)
    result <- validation_fit(sim, config, seed = config$seed + 100000L + index)
    result$paper <- original$paper
    result$source_result <- paths[index]
    result$generation_seconds <- original$generation_seconds
  }
  results[[length(results) + 1L]] <- validation_record(result, original$id, config)
}
validation_report(results, config)
