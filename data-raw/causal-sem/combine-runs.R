# Combine disjoint fit batches. Keep each original batch and its full fits.
source("data-raw/causal-sem/validation-helpers.R")
arguments <- commandArgs(trailingOnly = TRUE)
if (length(arguments) < 3L) stop("Supply a new output directory and at least two run directories.")
output <- arguments[1]
inputs <- arguments[-1]
if (dir.exists(output)) stop("Use a new output directory.")
manifests <- lapply(file.path(inputs, "run.rds"), readRDS)
records <- list()
for (input in inputs) {
  files <- list.files(input, pattern = "-rep[0-9]+[.]rds$", full.names = TRUE)
  for (file in files) {
    value <- readRDS(file)
    value$source_result <- normalizePath(file)
    records[[length(records) + 1L]] <- value
  }
}
ids <- vapply(records, `[[`, character(1), "id")
if (anyDuplicated(ids)) stop("The batches contain duplicate result IDs.")
dir.create(output, recursive = TRUE)
config <- manifests[[1]]$config
config$mode <- "combined-fits"
config$output <- output
config$models <- unique(unlist(lapply(manifests, function(x) x$config$models)))
config$batches <- stats::setNames(lapply(manifests, `[[`, "config"), normalizePath(inputs))
for (record in records) saveRDS(record, file.path(output, paste0(record$id, ".rds")))
validation_report(records, config)
for (name in c("predictive-frequencies.csv", "predictive-correlations.csv",
               "stan-sources.csv")) {
  files <- file.path(inputs, name)
  if (all(file.exists(files))) {
    rows <- do.call(rbind, lapply(files, utils::read.csv))
    utils::write.csv(rows, file.path(output, name), row.names = FALSE)
  }
}
for (input in inputs) {
  source_root <- file.path(input, "stan-sources")
  sources <- list.files(source_root, recursive = TRUE)
  for (source in sources) {
    source_file <- file.path(source_root, source)
    output_file <- file.path(output, "stan-sources", source)
    dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
    if (file.exists(output_file)) {
      if (unname(tools::md5sum(source_file)) != unname(tools::md5sum(output_file))) {
        stop("The batches contain different source files: ", source)
      }
    } else if (!file.copy(source_file, output_file)) {
      stop("A source snapshot cannot be copied: ", source)
    }
  }
}
