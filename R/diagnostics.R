# ============================================================================
# MCMC Diagnostics
# ============================================================================

#' Extract MCMC diagnostics
#'
#' Returns a summary of sampling diagnostics including divergences,
#' tree depth warnings, Rhat, and effective sample size.
#'
#' @param object A fitted model object.
#' @param ... Additional arguments (unused).
#'
#' @return A named list with:
#'   - `n_divergent`: total number of divergent transitions
#'   - `n_max_treedepth`: transitions hitting max tree depth
#'   - `max_rhat`: worst (highest) Rhat across all parameters
#'   - `min_ess_bulk`: smallest bulk ESS
#'   - `min_ess_tail`: smallest tail ESS
#'   - `mean_accept_prob`: mean acceptance probability
#' @export
dcVar_diagnostics <- function(object, ...) {
  UseMethod("dcVar_diagnostics")
}

#' @rdname dcVar_diagnostics
#' @export
dcVar_diagnostics.default <- function(object, ...) {
  cli_abort("{.fun dcVar_diagnostics} is not defined for objects of class {.cls {class(object)[[1]]}}.")
}

#' Internal: extract common diagnostics from a CmdStanMCMC fit
#' @noRd
.sampling_diagnostics_from_fit <- function(fit) {
  diag_summ <- fit$diagnostic_summary(quiet = TRUE)
  summ <- suppressWarnings(fit$summary())
  sampler_diags <- fit$sampler_diagnostics()
  accept_stat <- as.numeric(sampler_diags[, , "accept_stat__", drop = TRUE])
  accept_stat <- accept_stat[is.finite(accept_stat)]

  rhat <- summ$rhat[is.finite(summ$rhat)]
  ess_bulk <- summ$ess_bulk[is.finite(summ$ess_bulk)]
  ess_tail <- summ$ess_tail[is.finite(summ$ess_tail)]

  list(
    n_divergent = sum(diag_summ$num_divergent),
    n_max_treedepth = sum(diag_summ$num_max_treedepth),
    max_rhat = if (length(rhat) > 0) max(rhat) else NA_real_,
    min_ess_bulk = if (length(ess_bulk) > 0) min(ess_bulk) else NA_real_,
    min_ess_tail = if (length(ess_tail) > 0) min(ess_tail) else NA_real_,
    mean_accept_prob = if (length(accept_stat) > 0) mean(accept_stat) else NA_real_
  )
}

#' Internal: report post-sampling diagnostics to the user
#' @noRd
.report_sampling_outcome <- function(fit, model_label, chains = NA_integer_,
                                     rhat_threshold = 1.10) {
  diag <- .sampling_diagnostics_from_fit(fit)
  issues <- character()

  if (diag$n_divergent > 0) {
    issues <- c(
      issues,
      sprintf("%d divergent transition%s.", diag$n_divergent,
              if (diag$n_divergent == 1) "" else "s")
    )
  }
  if (diag$n_max_treedepth > 0) {
    issues <- c(
      issues,
      sprintf("%d transition%s hit the maximum treedepth.", diag$n_max_treedepth,
              if (diag$n_max_treedepth == 1) "" else "s")
    )
  }
  if (isTRUE(chains > 1) && is.finite(diag$max_rhat) && diag$max_rhat > rhat_threshold) {
    issues <- c(
      issues,
      sprintf("Max R-hat is %.3f.", diag$max_rhat)
    )
  }

  if (length(issues) == 0) {
    cli_alert_success("{model_label} sampling complete.")
  } else {
    cli_warn(c(
      "{model_label} sampling finished with diagnostic issues.",
      setNames(issues, rep("!", length(issues))),
      "i" = "Inspect {.fun dcVar_diagnostics} before using this fit for inference."
    ))
  }

  invisible(diag)
}

#' @rdname dcVar_diagnostics
#' @export
dcVar_diagnostics.dcVar_model_fit <- function(object, ...) {
  .sampling_diagnostics_from_fit(object$fit)
}
