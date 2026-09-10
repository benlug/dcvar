# Rebuild aggregate reports from stored records. Do not run a new fit.
source("data-raw/causal-sem/validation-helpers.R")
arguments <- commandArgs(trailingOnly = TRUE)
if (length(arguments) != 1L) stop("Supply one completed run directory.")
input <- arguments[1]
manifest <- readRDS(file.path(input, "run.rds"))
if (identical(manifest$config$mode, "priors")) stop("The prior draw run has a separate report.")
paths <- list.files(input, pattern = "-rep[0-9]+[.]rds$", full.names = TRUE)
records <- lapply(paths, readRDS)
manifest$config$output <- input
validation_report(records, manifest$config)
if (identical(manifest$config$mode, "sbc")) validation_rank_report(records, manifest$config)
