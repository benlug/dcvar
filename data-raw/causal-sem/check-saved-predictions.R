# Check saved predictive draws. This script does not run a new sampler.
source("data-raw/causal-sem/validation-helpers.R")
validation_setup()
arguments <- commandArgs(trailingOnly = TRUE)
if (length(arguments) != 1L) stop("Supply one result directory.")
input <- arguments[1]
paths <- list.files(input, pattern = "-rep[0-9]+[.]rds$", full.names = TRUE)
categories <- correlations <- list()
for (path in paths) {
  result <- readRDS(path)
  if (!identical(result$status, "ok")) next
  fit_path <- file.path(input, paste0("fit-seed-", result$seed, ".rds"))
  if (!file.exists(fit_path)) stop("A complete saved fit is missing: ", fit_path)
  fit <- readRDS(fit_path)
  data <- data.frame(x = fit$stan_data$group - 1L)
  for (j in 1:3) data[[paste0("u", j)]] <- ordered(fit$stan_data$u[, j], levels = seq_len(fit$stan_data$K[j]))
  sim <- list(data = data, args = list(treatment = "x", control = 0, treated = 1,
                                        indicators = list(L = paste0("u", 1:3))))
  categories[[length(categories) + 1L]] <- data.frame(id = result$id,
                                                       validation_category_check(fit, sim))
  correlations[[length(correlations) + 1L]] <- data.frame(id = result$id,
                                                          validation_dependence_check(fit))
}
write_report <- function(rows, name) {
  if (!length(rows)) return(invisible(NULL))
  data <- do.call(rbind, rows)
  data$outside <- data$observed < data$lower | data$observed > data$upper
  utils::write.csv(data, file.path(input, name), row.names = FALSE)
  cat(name, ":", sum(data$outside, na.rm = TRUE), "of", sum(!is.na(data$outside)),
      "observed values lie outside the intervals.\n")
}
write_report(categories, "predictive-frequencies.csv")
write_report(correlations, "predictive-correlations.csv")
