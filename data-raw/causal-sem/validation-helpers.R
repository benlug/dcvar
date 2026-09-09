# Use these functions from the repository root.

validation_models <- c("latent_covariate", "latent_outcome", "latent_mediator",
                       "latent_mediator_baseline")

validation_setup <- function() {
  extra_library <- Sys.getenv("DCVAR_VALIDATION_LIB", "")
  if (nzchar(extra_library)) .libPaths(c(extra_library, .libPaths()))
  if (file.exists("DESCRIPTION") && requireNamespace("pkgload", quietly = TRUE)) {
    pkgload::load_all(".", quiet = TRUE)
  } else {
    library(dcvar)
  }
  invisible(NULL)
}

validation_config <- function(mode, default_reps = 20L) {
  integer_setting <- function(name, default) {
    value <- suppressWarnings(as.integer(Sys.getenv(name, as.character(default))))
    if (length(value) != 1L || is.na(value) || value < 1L) stop(name, " must be positive.")
    value
  }
  output <- Sys.getenv("DCVAR_VALIDATION_OUTPUT",
                       file.path(tempdir(), paste0("dcvar-causal-", mode)))
  dir.create(output, recursive = TRUE, showWarnings = FALSE)
  models <- strsplit(Sys.getenv("DCVAR_VALIDATION_MODELS",
                                paste(validation_models, collapse = ",")), ",", fixed = TRUE)[[1]]
  if (any(!models %in% validation_models)) stop("Use documented model names.")
  source_files <- c(list.files("R", pattern = "causal", full.names = TRUE),
                    list.files("inst/stan", pattern = "causal", full.names = TRUE),
                    list.files("inst/stan/functions", pattern = "causal", full.names = TRUE))
  list(mode = mode, output = output, models = models,
       source_md5 = tools::md5sum(source_files),
       reps = integer_setting("DCVAR_VALIDATION_REPS", default_reps),
       n = integer_setting("DCVAR_VALIDATION_N", 400L),
       seed = integer_setting("DCVAR_VALIDATION_SEED", 90421L),
       ranks = integer_setting("DCVAR_VALIDATION_RANK_DRAWS", 100L),
       sampler = list(backend = Sys.getenv("DCVAR_VALIDATION_BACKEND", "cmdstanr"),
                      chains = integer_setting("DCVAR_VALIDATION_CHAINS", 4L),
                      cores = integer_setting("DCVAR_VALIDATION_CORES", 2L),
                      iter_warmup = integer_setting("DCVAR_VALIDATION_WARMUP", 1000L),
                      iter_sampling = integer_setting("DCVAR_VALIDATION_SAMPLING", 1000L),
                      adapt_delta = 0.99, max_treedepth = 13L, refresh = 0L))
}

validation_flatten <- function(parameters) {
  unlist(lapply(names(parameters), function(name) {
    value <- parameters[[name]]
    if (is.matrix(value)) {
      index <- arrayInd(seq_along(value), dim(value))
      names_out <- paste0(name, "[", index[, 1], ",", index[, 2], "]")
    } else {
      names_out <- paste0(name, "[", seq_along(value), "]")
    }
    stats::setNames(as.numeric(value), names_out)
  }), use.names = TRUE)
}

validation_effect_array <- function(fit) {
  value <- dcvar::causal_effects(fit, summary = FALSE)
  posterior::as_draws_array(value)
}

validation_effect_truth <- function(parameters, model, weights = c(0.5, 0.5)) {
  p <- parameters
  if (model %in% c("latent_covariate", "latent_outcome")) {
    target_mean <- sum(weights * p$mu_v)
    adjusted <- p$alpha + p$beta * target_mean
    return(c(ATE = adjusted[2] - adjusted[1], interaction = p$beta[2] - p$beta[1]))
  }
  target_q <- colSums(p$mu_q * weights)
  group_m <- p$alpha_m + rowSums(p$B * p$mu_q)
  target_m <- sum(weights * group_m)
  total <- direct <- numeric(2)
  for (g in 1:2) {
    adjusted_m <- p$alpha_m[g] + sum(p$B[g, ] * target_q)
    total[g] <- p$alpha_y[g] + sum(p$D[g, ] * target_q) + p$d[g] * adjusted_m
    direct[g] <- p$alpha_y[g] + sum(p$D[g, ] * target_q) + p$d[g] * target_m
  }
  c(ATE = total[2] - total[1], ADE = direct[2] - direct[1],
     AIE = (total[2] - direct[2]) - (total[1] - direct[1]))
}

validation_category_check <- function(fit, sim) {
  replicated <- posterior::as_draws_matrix(stats::predict(fit, type = "indicators"))
  items <- unlist(sim$args$indicators, use.names = FALSE)
  output <- list()
  row <- 0L
  for (group in c(sim$args$control, sim$args$treated)) {
    persons <- which(sim$data[[sim$args$treatment]] == group)
    for (item in seq_along(items)) {
      observed <- sim$data[[items[item]]][persons]
      selected <- replicated[, paste0("u_rep[", persons, ",", item, "]"), drop = FALSE]
      for (category in seq_along(levels(observed))) {
        frequencies <- rowMeans(selected == category)
        row <- row + 1L
        output[[row]] <- data.frame(group = group, item = items[item], category = category,
          observed = mean(as.integer(observed) == category),
          predicted = mean(frequencies),
          lower = unname(stats::quantile(frequencies, 0.025)),
          upper = unname(stats::quantile(frequencies, 0.975)))
      }
    }
  }
  do.call(rbind, output)
}

validation_dependence_check <- function(fit) {
  replicated <- posterior::as_draws_matrix(stats::predict(fit, type = "indicators"))
  observed <- fit$stan_data$u
  groups <- fit$stan_data$group
  pairs <- utils::combn(1:3, 2)
  output <- list()
  row <- 0L
  for (g in 1:2) {
    persons <- which(groups == g)
    for (pair in seq_len(ncol(pairs))) {
      a <- pairs[1, pair]
      b <- pairs[2, pair]
      left <- replicated[, paste0("u_rep[", persons, ",", a, "]"), drop = FALSE]
      right <- replicated[, paste0("u_rep[", persons, ",", b, "]"), drop = FALSE]
      left_mean <- rowMeans(left)
      right_mean <- rowMeans(right)
      covariance <- rowMeans(left * right) - left_mean * right_mean
      product_sd <- sqrt(pmax(0, rowMeans(left^2) - left_mean^2) *
                           pmax(0, rowMeans(right^2) - right_mean^2))
      correlation <- covariance / product_sd
      finite <- is.finite(correlation)
      observed_r <- if (stats::sd(observed[persons, a]) > 0 &&
                        stats::sd(observed[persons, b]) > 0) {
        stats::cor(observed[persons, a], observed[persons, b])
      } else NA_real_
      interval <- if (any(finite)) stats::quantile(correlation[finite], c(0.025, 0.975)) else c(NA, NA)
      row <- row + 1L
      output[[row]] <- data.frame(group = g, item_a = a, item_b = b,
        observed = observed_r, predicted = if (any(finite)) mean(correlation[finite]) else NA_real_,
        lower = interval[1], upper = interval[2], finite_fraction = mean(finite))
    }
  }
  do.call(rbind, output)
}

validation_fit <- function(sim, config, seed, priors = list()) {
  messages <- character()
  started <- proc.time()[["elapsed"]]
  result <- tryCatch(withCallingHandlers({
    fit_path <- file.path(config$output, paste0("fit-seed-", seed, ".rds"))
    if (file.exists(fit_path)) stop("A saved fit for this seed already exists: ", fit_path)
    stan_files <- c(system.file("stan", paste0("causal_sem_", sim$args$model, ".stan"),
                                package = "dcvar", mustWork = TRUE),
                    list.files(system.file("stan", "functions", package = "dcvar"),
                               pattern = "causal", full.names = TRUE))
    source_snapshot <- tools::md5sum(stan_files)
    if (length(sim$parameters)) {
      independent_truth <- validation_effect_truth(sim$parameters, sim$args$model,
                                                    sim$args$target_weights)
      stopifnot(isTRUE(all.equal(sim$effects, independent_truth, tolerance = 1e-10)))
      sim$effects <- independent_truth
    }
    arguments <- utils::modifyList(sim$args, config$sampler)
    arguments$seed <- seed
    arguments$priors <- priors
    arguments$standardize <- FALSE
    fit <- do.call(dcvar::dcvar_causal_sem, arguments)
    if (file.exists(fit_path)) stop("Another fit uses this output seed: ", fit_path)
    saveRDS(fit, fit_path)
    effects <- validation_effect_array(fit)
    all_draws <- posterior::as_draws_array(dcvar::draws(fit))
    truth <- validation_flatten(sim$parameters)
    available <- intersect(names(truth), posterior::variables(all_draws))
    parameter_sd <- apply(posterior::as_draws_matrix(all_draws)[, available, drop = FALSE],
                           2, stats::sd)
    available <- available[is.finite(parameter_sd) & parameter_sd > 0]
    selected <- parameter_summary <- NULL
    if (length(available)) {
      selected <- posterior::subset_draws(all_draws, variable = available)
      parameter_summary <- posterior::summarise_draws(selected, "mean", "sd", "rhat",
                                                      "ess_bulk", "ess_tail", "mcse_mean")
      parameter_summary$truth <- unname(truth[parameter_summary$variable])
    }
    effect_summary <- posterior::summarise_draws(effects, "mean", "sd", "rhat",
                                                 "ess_bulk", "ess_tail", "mcse_mean")
    matrix_effects <- posterior::as_draws_matrix(effects)
    intervals <- apply(matrix_effects, 2, stats::quantile, c(0.025, 0.975))
    effect_summary$lower <- intervals[1, effect_summary$variable]
    effect_summary$upper <- intervals[2, effect_summary$variable]
    effect_summary$truth <- unname(sim$effects[effect_summary$variable])
    diagnostic <- dcvar::dcvar_diagnostics(fit)
    list(status = "ok", effects = effect_summary, parameters = parameter_summary,
         effect_draws = effects, parameter_draws = selected, diagnostics = diagnostic,
         category_check = validation_category_check(fit, sim),
         correlation_check = validation_dependence_check(fit),
         stan_source_md5 = source_snapshot,
         sample_args = sim$args[names(sim$args) != "data"])
  }, warning = function(warning) {
    messages <<- c(messages, conditionMessage(warning))
    invokeRestart("muffleWarning")
  }), error = function(error) list(status = "error", error = conditionMessage(error)))
  result$warnings <- unique(messages)
  result$seconds <- proc.time()[["elapsed"]] - started
  result$seed <- seed
  result
}

validation_record <- function(result, id, config) {
  result$id <- id
  path <- file.path(config$output, paste0(id, ".rds"))
  if (file.exists(path)) stop("The result file exists: ", path)
  saveRDS(result, path)
  message(id, ": ", result$status, " (", round(result$seconds, 1), " s)")
  result
}

validation_report <- function(results, config) {
  status <- do.call(rbind, lapply(results, function(x) {
    diagnostic <- x$diagnostics
    metric <- function(name) if (is.null(diagnostic[[name]])) NA_real_ else diagnostic[[name]]
    pass <- if (is.null(diagnostic) || is.null(x$effects)) NA else
      all(is.finite(c(metric("max_rhat"), x$effects$rhat, x$effects$ess_bulk,
                        x$effects$ess_tail, x$effects$mcse_mean / x$effects$sd))) &&
      metric("n_divergent") == 0 && metric("n_max_treedepth") == 0 &&
      !isTRUE(diagnostic$incomplete_diagnostics) &&
      is.data.frame(diagnostic$parameters) && nrow(diagnostic$parameters) > 0L &&
      all(is.finite(diagnostic$parameters$rhat)) &&
      all(diagnostic$parameters$rhat < 1.01) &&
      metric("max_rhat") < 1.01 && all(x$effects$rhat < 1.01) &&
      all(x$effects$ess_bulk >= 400) && all(x$effects$ess_tail >= 400) &&
      all(x$effects$mcse_mean / x$effects$sd <= 0.05) &&
      length(diagnostic$ebfmi) > 0L && all(is.finite(diagnostic$ebfmi)) &&
      all(diagnostic$ebfmi > 0.3)
    data.frame(id = x$id, status = x$status, seconds = x$seconds,
               divergences = metric("n_divergent"),
               max_treedepth = metric("n_max_treedepth"),
               max_rhat = metric("max_rhat"),
               min_ebfmi = if (is.null(diagnostic$ebfmi)) NA_real_ else min(diagnostic$ebfmi),
               incomplete_diagnostics = if (is.null(diagnostic$incomplete_diagnostics)) NA else diagnostic$incomplete_diagnostics,
               diagnostic_pass = pass, warnings = paste(x$warnings, collapse = " | "),
               error = if (is.null(x$error)) "" else x$error)
  }))
  estimates <- lapply(results, function(x) {
    if (!identical(x$status, "ok")) return(NULL)
    value <- as.data.frame(x$effects)
    value$id <- x$id
    value
  })
  estimates <- do.call(rbind, estimates)
  summary <- NULL
  if (!is.null(estimates) && nrow(estimates)) {
    split_rows <- split(estimates, sub("-rep[0-9]+$", "", estimates$id))
    summary <- do.call(rbind, lapply(names(split_rows), function(condition) {
      rows <- split_rows[[condition]]
      do.call(rbind, lapply(split(rows, rows$variable), function(values) {
        deviation <- values$mean - values$truth
        coverage <- mean(values$lower <= values$truth & values$upper >= values$truth)
        n <- nrow(values)
        z <- stats::qnorm(0.975)
        denominator <- 1 + z^2 / n
        coverage_center <- (coverage + z^2 / (2 * n)) / denominator
        coverage_half <- z * sqrt(coverage * (1 - coverage) / n +
                                    z^2 / (4 * n^2)) / denominator
        data.frame(condition = condition, effect = values$variable[1],
                   successful = n, bias = mean(deviation), absolute_bias = abs(mean(deviation)),
                   rmse = sqrt(mean(deviation^2)), coverage = coverage,
                   coverage_mcse = sqrt(coverage * (1 - coverage) / n),
                   nominal_coverage_mcse = sqrt(0.95 * 0.05 / n),
                   coverage_wilson_lower = max(0, coverage_center - coverage_half),
                   coverage_wilson_upper = min(1, coverage_center + coverage_half),
                   interval_width = mean(values$upper - values$lower),
                   max_rhat = max(values$rhat), min_ess_bulk = min(values$ess_bulk),
                   min_ess_tail = min(values$ess_tail),
                   max_mcse_fraction = max(values$mcse_mean / values$sd))
      }))
    }))
  }
  utils::write.csv(status, file.path(config$output, "status.csv"), row.names = FALSE)
  if (!is.null(estimates)) utils::write.csv(estimates, file.path(config$output, "effects.csv"), row.names = FALSE)
  if (!is.null(summary)) utils::write.csv(summary, file.path(config$output, "summary.csv"), row.names = FALSE)
  saveRDS(list(config = config, session = utils::sessionInfo(), status = status,
               summary = summary), file.path(config$output, "run.rds"))
  print(status)
  if (!is.null(summary)) print(summary)
  invisible(summary)
}

validation_rank_report <- function(results, config) {
  rows <- do.call(rbind, lapply(results, function(result) {
    if (is.null(result$ranks)) return(NULL)
    data.frame(id = result$id, result$ranks)
  }))
  if (!is.null(rows)) {
    status <- utils::read.csv(file.path(config$output, "status.csv"))
    rows$diagnostic_pass <- status$diagnostic_pass[match(rows$id, status$id)]
    utils::write.csv(rows, file.path(config$output, "ranks.csv"), row.names = FALSE)
  }
  invisible(rows)
}
