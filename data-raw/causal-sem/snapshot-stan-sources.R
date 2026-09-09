# Save the source held by each fitted object. The package source can change later.
source("data-raw/causal-sem/validation-helpers.R")
validation_setup()
arguments <- commandArgs(trailingOnly = TRUE)
if (length(arguments) != 1L) stop("Supply one result directory.")
input <- arguments[1]
output <- file.path(input, "stan-sources")
dir.create(file.path(output, "functions"), recursive = TRUE, showWarnings = FALSE)
includes <- list.files("inst/stan/functions", pattern = "causal", full.names = TRUE)
file.copy(includes, file.path(output, "functions"), overwrite = FALSE)
paths <- list.files(input, pattern = "-rep[0-9]+[.]rds$", full.names = TRUE)
rows <- list()
for (path in paths) {
  result <- readRDS(path)
  if (!identical(result$status, "ok")) next
  fit_file <- file.path(input, paste0("fit-seed-", result$seed, ".rds"))
  if (!file.exists(fit_file)) stop("A complete saved fit is missing: ", fit_file)
  fit <- readRDS(fit_file)
  code <- if (fit$backend == "cmdstanr") fit$fit$code() else fit$fit@stanmodel@model_code
  code <- strsplit(paste(code, collapse = "\n"), "\n", fixed = TRUE)[[1]]
  marker <- switch(fit$model, latent_covariate = "marginal_sd",
                   latent_mediator = "marginal_sd_y",
                   latent_mediator_baseline = "hypot(sigma_m", "unused_marker")
  parameterization <- if (any(grepl(marker, code, fixed = TRUE))) {
    "conditioned_normal"
  } else "original_noncentered"
  source_file <- file.path(output, paste0(result$id, ".stan"))
  if (file.exists(source_file) && !identical(readLines(source_file), code)) {
    stop("A different source snapshot exists: ", source_file)
  }
  if (!file.exists(source_file)) writeLines(code, source_file, useBytes = TRUE)
  data_file <- file.path(output, paste0(result$id, "-data.rds"))
  if (file.exists(data_file) && !identical(readRDS(data_file), fit$stan_data)) {
    stop("A different data snapshot exists: ", data_file)
  }
  if (!file.exists(data_file)) saveRDS(fit$stan_data, data_file, version = 3)
  rows[[length(rows) + 1L]] <- data.frame(id = result$id, model = fit$model,
    parameterization = parameterization, normalized_md5 = unname(tools::md5sum(source_file)),
    source = file.path("stan-sources", basename(source_file)),
    data_md5 = unname(tools::md5sum(data_file)),
    data = file.path("stan-sources", basename(data_file)))
}
report <- do.call(rbind, rows)
utils::write.csv(report, file.path(input, "stan-sources.csv"), row.names = FALSE)
print(report[c("id", "parameterization", "normalized_md5")])
