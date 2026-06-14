#' Fit the HMM copula model
#'
#' Fits a Hidden Markov Model copula VAR(1) with K discrete states and
#' state-specific correlations. Uses ordered z_rho constraint to prevent
#' label switching and a sticky Dirichlet prior to encourage state persistence.
#'
#' @section State-specific parameters and per-state margins:
#' By default only the copula correlation switches by state. The `switch`
#' argument turns a full Markov-switching VAR on: the intercepts, the VAR
#' coefficients, and the residual scales can each be made state-specific. The
#' marginal **family** can also switch by state by passing a length-`K`
#' `margins` list. Because the states are identified by an ordered correlation
#' (`z_rho`), `rho` is the mandatory anchor and a per-state `margins` list is
#' consumed in increasing-correlation order (`margins[[1]]` is the lowest-rho
#' state). Richly-parameterized configurations on a single short bivariate series
#' are weakly identified; fit a ladder (`switch = "rho"` ->
#' `switch = c("rho", "ar")` -> full) and widen `prior_phi_sd`. See
#' [hmm_state_params()] for the per-state parameter view.
#'
#' @inheritParams dcvar
#' @param margins Marginal distribution specification. Either a single string
#'   applied to both variables, a length-2 character vector giving a per-variable
#'   (mixed) margin, e.g. `c("normal", "exponential")`, or -- for a state-specific
#'   family -- a length-`K` list of such specifications, one per hidden state,
#'   e.g. `list(c("normal", "normal"), c("exponential", "gamma"))`. Each entry is
#'   one of `"normal"` (default), `"exponential"`, `"skew_normal"`, or `"gamma"`.
#'   A per-state list is consumed in increasing-`rho` order.
#' @param skew_direction Skew direction for asymmetric margins: each element `1`
#'   (right-skewed) or `-1` (left-skewed). A length-2 vector (recycled across
#'   states) or, for per-state margins, a length-`K` list of length-2 vectors.
#'   Required whenever any (state, dimension) uses an `"exponential"` or
#'   `"gamma"` margin; only those entries consult it.
#' @param switch Which parameters become state-specific. `"rho"` (default) keeps
#'   only the copula correlation switching (the classic HMM). A character vector
#'   drawn from `c("rho", "mu", "phi", "ar", "cross", "sigma")` (with `"ar"` /
#'   `"cross"` / coefficient names like `"phi11"` selecting a subset of the VAR
#'   coefficients), or `TRUE` for everything. `rho` must be included whenever any
#'   other component switches (it is the label-switching anchor).
#' @param K Number of hidden states (default: 2).
#' @param prior_kappa Sticky Dirichlet self-transition concentration (default: 10).
#' @param prior_alpha_off Sticky Dirichlet off-diagonal concentration (default: 1).
#' @param prior_z_rho_sd Prior SD for state-specific z_rho (default: 1.0).
#'
#' @param backend Character: `"auto"` (default, uses rstan), `"rstan"`, or
#'   `"cmdstanr"`. Can also be set globally via
#'   `options(dcvar.backend = "cmdstanr")`.
#' @param ... Additional backend-specific sampling arguments.
#'
#' @return A `dcvar_hmm_fit` object.
#'
#' @seealso [dcvar()] for the smooth time-varying model,
#'   [dcvar_constant()] for the time-invariant baseline,
#'   [hmm_states()] for state extraction, [plot_hmm_states()] for visualisation,
#'   [dcvar_compare()] for LOO-CV model comparison.
#' @export
#'
#' @examples
#' \donttest{
#' sim <- simulate_dcvar(
#'   n_time = 12,
#'   rho_trajectory = rho_step(12),
#'   seed = 1
#' )
#' fit <- dcvar_hmm(
#'   sim$Y_df,
#'   vars = c("y1", "y2"),
#'   K = 2,
#'   chains = 1,
#'   iter_warmup = 10,
#'   iter_sampling = 10,
#'   refresh = 0,
#'   seed = 1
#' )
#' print(fit)
#' hmm_states(fit)
#' }
dcvar_hmm <- function(data, vars, K = 2,
                      time_var = "time",
                      standardize = TRUE,
                      margins = "normal",
                      skew_direction = NULL,
                      switch = "rho",
                      allow_gaps = FALSE,
                      prior_mu_sd = 2,
                      prior_phi_sd = 0.5,
                      prior_sigma_eps_rate = 1,
                      prior_kappa = 10,
                      prior_alpha_off = 1,
                      prior_z_rho_sd = 1.0,
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
  backend <- .resolve_backend(backend)
  .validate_sampling_args(chains, iter_warmup, iter_sampling,
                          adapt_delta, max_treedepth)

  if (!is.numeric(K) || length(K) != 1 || K < 2 || K != as.integer(K)) {
    cli_abort("{.arg K} must be an integer >= 2, got {.val {K}}.")
  }
  K <- as.integer(K)

  # Resolve the switching spec and the (global or per-state) margin configuration.
  sw <- .resolve_switch_spec(switch)
  cfg <- .hmm_margin_config(margins, skew_direction, K)

  # Force state-specific margin parameters when residuals are state-specific and
  # any exp/gamma family is present (the feasibility-bounded scale needs one
  # bijection per state), when families differ per state, or on explicit request.
  any_exp_gamma <- any(cfg$family %in% c(2L, 4L))
  if (cfg$per_state_differ ||
      ((sw$mu == 1L || any(sw$phi_mask > 0L)) && any_exp_gamma)) {
    sw$margins <- 1L
  }

  needs_engine <- cfg$per_state || sw$mu == 1L || any(sw$phi_mask > 0L) || sw$margins == 1L
  uses_engine <- .uses_hmm_switching_stan_file(stan_file) ||
    (is.null(stan_file) && needs_engine)

  if (needs_engine && !is.null(stan_file) && !uses_engine) {
    cli_abort(c(
      "{.arg switch} / per-state {.arg margins} require the bundled switching engine.",
      "i" = "Omit {.arg stan_file}, or pass {.code stan_file = dcvar_stan_path(\"hmm_switching\")}."
    ))
  }

  # The engine overrides family/skew via .as_hmm_switching_stan_data(), so a
  # normal placeholder suffices for prepare_hmm_data(); the legacy path uses the
  # real global margins.
  base_margins <- if (uses_engine) "normal" else cfg$margins_global
  base_skew <- if (uses_engine) NULL else cfg$skew_global

  stan_data <- prepare_hmm_data(
    data, vars, K, time_var, standardize, base_margins, base_skew,
    prior_mu_sd, prior_phi_sd, prior_sigma_eps_rate,
    prior_kappa, prior_alpha_off, prior_z_rho_sd,
    allow_gaps = allow_gaps
  )

  if (uses_engine) {
    comps <- c(
      "rho",
      if (sw$mu == 1L) "mu",
      if (any(sw$phi_mask > 0L)) {
        if (all(sw$phi_mask == 1L)) "Phi" else paste0("Phi:", paste(names(sw$phi_mask)[sw$phi_mask == 1L], collapse = "/"))
      },
      if (sw$margins == 1L) "margins"
    )
    cli_inform("Fitting Markov-switching HMM [state-specific: {paste(comps, collapse = ', ')}] (n_time = {stan_data$n_time}, D = {stan_data$D}, K = {K})...")
  } else {
    margins_label <- if (all(cfg$margins_char == "normal")) {
      ""
    } else {
      paste0(" [", paste(cfg$margins_char[1, ], collapse = ", "), "]")
    }
    cli_inform("Fitting HMM copula model{margins_label} (n_time = {stan_data$n_time}, D = {stan_data$D}, K = {K})...")
  }

  # Compile model (the engine compiles the bundled hmm_switching.stan)
  model_type <- if (uses_engine) "hmm_switching" else "hmm"
  compile_stan_file <- if (uses_engine) dcvar_stan_path("hmm_switching") else stan_file
  model <- .compile_model(model_type, margins = base_margins,
                          stan_file = compile_stan_file, backend = backend)

  # Augment the prepared data into the engine's superset block.
  if (uses_engine) {
    stan_data <- .as_hmm_switching_stan_data(stan_data, sw, cfg$family, cfg$skew,
                                             prior_phi_dev_sd = prior_phi_sd)
  }

  # Weak-identification warnings for richly-parameterized switching configs.
  if (uses_engine) {
    if (sw$mu == 1L && all(sw$phi_mask == 1L) && stan_data$n_time < 150) {
      cli_warn(c(
        "All-switching (mu + Phi + rho) on a short series (n_time = {stan_data$n_time}) is weakly identified.",
        "i" = "Consider a switching ladder ({.code \"rho\"} -> {.code c(\"rho\", \"ar\")} -> full) and widen {.arg prior_phi_sd}."
      ))
    }
    if (K >= 4L) {
      cli_warn("Recovering {K} regimes from a single bivariate series is hard; prefer K = 2-3 unless the design is strong.")
    }
    if (cfg$per_state_differ) {
      cli_warn(c(
        "Per-state {.arg margins} or {.arg skew_direction} differ across states.",
        "i" = "States are ordered by increasing correlation, so the first state specification is the lowest-rho state."
      ))
    }
  }

  # Default init
  if (is.null(init)) {
    D <- stan_data$D
    init <- if (uses_engine) {
      function() .init_hmm_switching_params(D, K, sw)
    } else {
      function() .init_hmm_params(D, K, base_margins)
    }
  }

  cores <- .normalize_cores(cores, chains)

  # Fit
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

  .report_sampling_outcome(fit, if (uses_engine) "Markov-switching HMM" else "HMM copula",
                           chains = chains, backend = backend)

  # Wrap in S3 class
  new_dcvar_hmm_fit(
    fit = fit,
    stan_data = stan_data,
    K = K,
    vars = vars,
    standardized = standardize,
    margins = cfg$margins_global %||% "normal",
    skew_direction = skew_direction,
    switch = sw,
    switching = uses_engine,
    margins_matrix = if (uses_engine) cfg$margins_char else NULL,
    backend = backend,
    priors = list(
      mu_sd = prior_mu_sd,
      phi_sd = prior_phi_sd,
      sigma_eps_rate = prior_sigma_eps_rate,
      kappa = prior_kappa,
      alpha_off = prior_alpha_off,
      z_rho_sd = prior_z_rho_sd
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
