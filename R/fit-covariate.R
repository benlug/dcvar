#' Prepare data for the covariate DC-VAR model
#'
#' Transforms a data frame into a list suitable for the Gaussian covariate
#' DC-VAR Stan models. Outcome rows are sorted by `time_var`; rows with missing
#' outcomes are removed with the same adjacency rules used by [prepare_dcvar_data()],
#' and covariates are filtered in the same order so `X[t + 1, ]` aligns with the
#' outcome occasion of transition `Y[t, ] -> Y[t + 1, ]`.
#'
#' @inheritParams prepare_dcvar_data
#' @param covariates Character vector of covariate column names.
#' @param standardize_covariates Logical; whether to z-score covariates
#'   (default: `FALSE`). Binary phase indicators should usually be left on
#'   their original scale.
#' @param prior_beta_sd Prior SD for the covariate effects:
#'   `beta ~ normal(0, prior_beta_sd)`.
#' @param zero_init_eta Logical; if `TRUE` (default), fixes the first residual
#'   drift state at zero (`eta[1] = 0`). If `FALSE`, the first transition can
#'   receive an immediate residual random-walk shock.
#'
#' @return A named list suitable as Stan data input.
#' @export
prepare_dcvar_covariate_data <- function(data, vars, covariates,
                                         time_var = "time",
                                         standardize = TRUE,
                                         standardize_covariates = FALSE,
                                         allow_gaps = FALSE,
                                         prior_mu_sd = 2,
                                         prior_phi_sd = 0.5,
                                         prior_sigma_eps_rate = 1,
                                         prior_sigma_omega_rate = 0.1,
                                         prior_rho_init_sd = 1,
                                         prior_beta_sd = 1,
                                         zero_init_eta = TRUE) {
  if (!is.character(covariates) || length(covariates) == 0L || anyNA(covariates)) {
    cli_abort("{.arg covariates} must be a non-empty character vector.")
  }
  .prep_validate_unique_vars(covariates, "covariates")
  .prep_validate_scalar_logical(standardize_covariates, "standardize_covariates")
  .prep_validate_scalar_logical(zero_init_eta, "zero_init_eta")
  .prep_validate_positive_scalar(prior_mu_sd, "prior_mu_sd")
  .prep_validate_positive_scalar(prior_phi_sd, "prior_phi_sd")
  .prep_validate_positive_scalar(prior_sigma_eps_rate, "prior_sigma_eps_rate")
  .prep_validate_positive_scalar(prior_sigma_omega_rate, "prior_sigma_omega_rate")
  .prep_validate_positive_scalar(prior_rho_init_sd, "prior_rho_init_sd")
  .prep_validate_positive_scalar(prior_beta_sd, "prior_beta_sd")

  missing_covariates <- setdiff(covariates, names(data))
  if (length(missing_covariates) > 0L) {
    cli_abort("Covariate column{?s} not found in data: {.val {missing_covariates}}")
  }

  prep <- .prepare_var_data(data, vars, time_var, standardize, allow_gaps)
  aligned_data <- prep$data

  X_df <- aligned_data[, covariates, drop = FALSE]
  for (v in covariates) {
    if (!is.numeric(X_df[[v]])) {
      cli_abort("Covariate {.val {v}} must be numeric, got {.cls {class(X_df[[v]])[1]}}.")
    }
  }

  X <- as.matrix(X_df)
  if (any(is.na(X) | is.nan(X) | is.infinite(X))) {
    cli_abort("{.arg covariates} must be complete finite numeric values after outcome filtering.")
  }

  X_means <- NULL
  X_sds <- NULL
  if (standardize_covariates) {
    X_means <- colMeans(X)
    X_sds <- apply(X, 2, stats::sd)
    if (any(X_sds < 1e-10)) {
      cli_abort("Cannot standardize covariates: one or more covariates have zero (or near-zero) variance.")
    }
    X <- scale(X)
  }

  stan_data <- list(
    n_time = prep$T_obs,
    D = prep$D,
    Y = prep$Y,
    sigma_mu_prior = prior_mu_sd,
    sigma_phi_prior = prior_phi_sd,
    sigma_eps_prior = prior_sigma_eps_rate,
    sigma_omega_prior = prior_sigma_omega_rate,
    rho_init_prior_sd = prior_rho_init_sd,
    P = ncol(X),
    X = X,
    sigma_beta_prior = prior_beta_sd,
    zero_init_eta = as.integer(zero_init_eta)
  )

  attr(stan_data, "vars") <- prep$vars
  attr(stan_data, "covariates") <- covariates
  attr(stan_data, "standardized") <- prep$standardized
  attr(stan_data, "standardized_covariates") <- standardize_covariates
  attr(stan_data, "time_values") <- prep$time_values
  attr(stan_data, "zero_init_eta") <- zero_init_eta
  if (prep$standardized) {
    attr(stan_data, "Y_means") <- prep$Y_means
    attr(stan_data, "Y_sds") <- prep$Y_sds
  }
  if (standardize_covariates) {
    attr(stan_data, "X_means") <- X_means
    attr(stan_data, "X_sds") <- X_sds
  }

  stan_data
}


#' Fit the covariate DC-VAR model
#'
#' Fits a Gaussian Dynamic Copula VAR(1) model in which contemporaneous
#' innovation dependence varies on the Fisher-z scale as a function of observed
#' covariates. With `drift = TRUE`, the model adds residual random-walk drift:
#' `zeta_i = beta_0 + x_{i+1}' beta + eta_i`. With `drift = FALSE`, the
#' dependence trajectory is explained entirely by the covariates.
#'
#' @inheritParams dcvar
#' @param covariates Character vector of covariate column names.
#' @param standardize_covariates Logical; whether to z-score covariates
#'   (default: `FALSE`).
#' @param drift Logical; if `TRUE` (default), include residual random-walk drift
#'   after the covariate effect. If `FALSE`, fit the no-drift sensitivity model.
#' @param zero_init_eta Logical; if `TRUE` (default), fixes `eta[1] = 0` in the
#'   residual-drift model.
#' @param prior_beta_sd Prior SD for covariate effects.
#'
#' @return A `dcvar_covariate_fit` object.
#'
#' @seealso [prepare_dcvar_covariate_data()], [covariate_effects()],
#'   [rho_trajectory()], and [dcvar()] for the covariate-free random-walk model.
#' @export
#'
#' @examples
#' \donttest{
#' sim <- simulate_dcvar(
#'   n_time = 12,
#'   rho_trajectory = rho_step(12),
#'   seed = 1
#' )
#' sim$Y_df$phase <- as.numeric(sim$Y_df$time > 6)
#' fit <- dcvar_covariate(
#'   sim$Y_df,
#'   vars = c("y1", "y2"),
#'   covariates = "phase",
#'   chains = 1,
#'   iter_warmup = 10,
#'   iter_sampling = 10,
#'   refresh = 0,
#'   seed = 1
#' )
#' covariate_effects(fit)
#'
#' fit_nodrift <- dcvar_covariate(
#'   sim$Y_df,
#'   vars = c("y1", "y2"),
#'   covariates = "phase",
#'   drift = FALSE,
#'   chains = 1,
#'   iter_warmup = 10,
#'   iter_sampling = 10,
#'   refresh = 0,
#'   seed = 1
#' )
#' rho_trajectory(fit_nodrift)
#' }
dcvar_covariate <- function(data, vars, covariates, time_var = "time",
                            standardize = TRUE,
                            standardize_covariates = FALSE,
                            drift = TRUE,
                            zero_init_eta = TRUE,
                            allow_gaps = FALSE,
                            prior_mu_sd = 2,
                            prior_phi_sd = 0.5,
                            prior_sigma_eps_rate = 1,
                            prior_sigma_omega_rate = 0.1,
                            prior_rho_init_sd = 1,
                            prior_beta_sd = 1,
                            chains = 4,
                            iter_warmup = 2000,
                            iter_sampling = 4000,
                            adapt_delta = 0.99,
                            max_treedepth = 12,
                            seed = NULL,
                            cores = NULL,
                            refresh = 500,
                            init = NULL,
                            stan_file = NULL,
                            backend = getOption("dcvar.backend", "auto"),
                            ...) {
  .prep_validate_scalar_logical(drift, "drift")
  backend <- .resolve_backend(backend)
  .validate_sampling_args(chains, iter_warmup, iter_sampling,
                          adapt_delta, max_treedepth)

  stan_data <- prepare_dcvar_covariate_data(
    data = data,
    vars = vars,
    covariates = covariates,
    time_var = time_var,
    standardize = standardize,
    standardize_covariates = standardize_covariates,
    allow_gaps = allow_gaps,
    prior_mu_sd = prior_mu_sd,
    prior_phi_sd = prior_phi_sd,
    prior_sigma_eps_rate = prior_sigma_eps_rate,
    prior_sigma_omega_rate = prior_sigma_omega_rate,
    prior_rho_init_sd = prior_rho_init_sd,
    prior_beta_sd = prior_beta_sd,
    zero_init_eta = zero_init_eta
  )

  model_type <- if (drift) "dcvar_covariate" else "dcvar_covariate_nodrift"
  drift_label <- if (drift) "with residual drift" else "no residual drift"
  cli_inform("Fitting covariate DC-VAR model ({drift_label}; n_time = {stan_data$n_time}, D = {stan_data$D}, P = {stan_data$P})...")

  model <- .compile_model(model_type, margins = "normal", stan_file = stan_file,
                          backend = backend)

  if (!drift) {
    stan_data$sigma_omega_prior <- NULL
    stan_data$zero_init_eta <- NULL
  }

  if (is.null(init)) {
    D <- stan_data$D
    n_time_obs <- stan_data$n_time
    P <- stan_data$P
    init <- function() .init_dcvar_covariate_params(
      D = D,
      T_obs = n_time_obs,
      P = P,
      drift = drift,
      zero_init_eta = zero_init_eta
    )
  }

  cores <- .normalize_cores(cores, chains)

  fit <- .sample_model(
    compiled_model = model,
    stan_data = stan_data,
    backend = backend,
    chains = chains,
    iter_warmup = iter_warmup,
    iter_sampling = iter_sampling,
    adapt_delta = adapt_delta,
    max_treedepth = max_treedepth,
    seed = seed,
    cores = cores,
    init = init,
    refresh = refresh,
    ...
  )

  .report_sampling_outcome(fit, "Covariate DC-VAR", chains = chains, backend = backend)

  new_dcvar_covariate_fit(
    fit = fit,
    stan_data = stan_data,
    vars = vars,
    covariates = covariates,
    standardized = standardize,
    standardized_covariates = standardize_covariates,
    drift = drift,
    zero_init_eta = zero_init_eta,
    backend = backend,
    priors = list(
      mu_sd = prior_mu_sd,
      phi_sd = prior_phi_sd,
      sigma_eps_rate = prior_sigma_eps_rate,
      sigma_omega_rate = prior_sigma_omega_rate,
      rho_init_sd = prior_rho_init_sd,
      beta_sd = prior_beta_sd
    ),
    meta = list(
      chains = chains,
      iter_warmup = iter_warmup,
      iter_sampling = iter_sampling,
      adapt_delta = adapt_delta,
      max_treedepth = max_treedepth,
      seed = seed
    )
  )
}
