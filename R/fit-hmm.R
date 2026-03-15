#' Fit the HMM copula model
#'
#' Fits a Hidden Markov Model copula VAR(1) with K discrete states and
#' state-specific correlations. Uses ordered z_rho constraint to prevent
#' label switching and a sticky Dirichlet prior to encourage state persistence.
#'
#' @inheritParams dcvar
#' @param margins Character string specifying the marginal distribution.
#'   One of `"normal"` (default), `"exponential"`, `"skew_normal"`, or `"gamma"`.
#' @param skew_direction Integer vector of length D indicating skew direction
#'   for asymmetric margins. Each element must be `1` (right-skewed) or `-1`
#'   (left-skewed). Required for `"exponential"` and `"gamma"` margins.
#' @param K Number of hidden states (default: 2).
#' @param prior_kappa Sticky Dirichlet self-transition concentration (default: 10).
#' @param prior_alpha_off Sticky Dirichlet off-diagonal concentration (default: 1).
#' @param prior_z_rho_sd Prior SD for state-specific z_rho (default: 1.0).
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
#' \dontrun{
#' sim <- simulate_dcvar(T = 100, rho_trajectory = rho_step(100))
#' fit <- dcvar_hmm(sim$Y_df, vars = c("y1", "y2"), K = 2)
#' print(fit)
#' hmm_states(fit)
#' }
dcvar_hmm <- function(data, vars, K = 2,
                      time_var = "time",
                      standardize = TRUE,
                      margins = "normal",
                      skew_direction = NULL,
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
                      ...) {
  .check_cmdstanr()
  .validate_sampling_args(chains, iter_warmup, iter_sampling,
                          adapt_delta, max_treedepth)
  .validate_margins(margins, skew_direction)

  if (!is.numeric(K) || length(K) != 1 || K < 2 || K != as.integer(K)) {
    cli_abort("{.arg K} must be an integer >= 2, got {.val {K}}.")
  }

  # Prepare data
  stan_data <- prepare_hmm_data(
    data, vars, K, time_var, standardize, margins, skew_direction,
    prior_mu_sd, prior_phi_sd, prior_sigma_eps_rate,
    prior_kappa, prior_alpha_off, prior_z_rho_sd,
    allow_gaps = allow_gaps
  )

  margins_label <- if (margins == "normal") "" else paste0(" [", margins, "]")
  cli_inform("Fitting HMM copula model{margins_label} (T = {stan_data$T}, D = {stan_data$D}, K = {K})...")

  # Compile model
  model <- .compile_model("hmm", margins = margins, stan_file = stan_file)

  # Default init
  if (is.null(init)) {
    D <- stan_data$D
    init <- function() .init_hmm_params(D, K, margins)
  }

  if (is.null(cores)) cores <- parallel::detectCores(logical = FALSE)

  # Fit
  fit <- model$sample(
    data = stan_data,
    chains = chains,
    iter_warmup = iter_warmup,
    iter_sampling = iter_sampling,
    adapt_delta = adapt_delta,
    max_treedepth = max_treedepth,
    seed = seed,
    parallel_chains = cores,
    init = init,
    refresh = refresh,
    ...
  )

  .report_sampling_outcome(fit, "HMM copula", chains = chains)

  # Wrap in S3 class
  new_dcvar_hmm_fit(
    fit = fit,
    stan_data = stan_data,
    K = K,
    vars = vars,
    standardized = standardize,
    margins = margins,
    skew_direction = skew_direction,
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
