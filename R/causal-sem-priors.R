#' Internal: validate the causal SEM model name
#' @noRd
.causal_sem_validate_model <- function(model) {
  match.arg(model, c(
    "latent_covariate", "latent_outcome", "latent_mediator",
    "latent_mediator_baseline"
  ))
}

#' Internal: validate and complete causal SEM priors
#' @noRd
.causal_sem_priors <- function(priors = list()) {
  defaults <- list(
    beta_sd = 1, intercept_sd = 2.5, scale_sd = 1.5,
    loading_mean = 1, loading_sd = 1, threshold_sd = 2.5,
    item_log_sd = 0.5, lkj_shape = 2
  )
  if (!is.list(priors) || (length(priors) &&
      (is.null(names(priors)) || anyNA(names(priors)) ||
       any(!nzchar(names(priors))) || anyDuplicated(names(priors))))) {
    cli_abort("{.arg priors} must be a list with unique non-empty names.")
  }
  unknown <- setdiff(names(priors), names(defaults))
  if (length(unknown)) {
    cli_abort("Unknown prior name{?s}: {.val {unknown}}.")
  }
  for (name in names(priors)) defaults[name] <- priors[name]
  for (name in names(defaults)) {
    value <- defaults[[name]]
    if (!is.numeric(value) || !is.null(dim(value)) ||
        length(value) != 1L || !is.finite(value) ||
        (name != "loading_mean" && value <= 0)) {
      requirement <- if (name == "loading_mean") "finite" else "positive finite"
      cli_abort("Prior {.val {name}} must be one {requirement} number.")
    }
    defaults[[name]] <- as.numeric(value)
  }
  defaults
}

#' Internal: make the Stan prior fields
#' @noRd
.causal_sem_stan_priors <- function(priors = list()) {
  out <- .causal_sem_priors(priors)
  names(out) <- paste0("prior_", names(out))
  out
}

#' Internal: draw parameters from the causal SEM prior
#' @noRd
.causal_sem_draw_prior <- function(model, priors = list(), K = c(3L, 3L, 3L)) {
  model <- .causal_sem_validate_model(model)
  priors <- .causal_sem_priors(priors)
  if (!is.numeric(K) || !is.null(dim(K)) || length(K) != 3L || anyNA(K) ||
      any(!is.finite(K)) || any(K != floor(K)) || any(K < 3L | K > 5L)) {
    cli_abort("{.arg K} must contain three integers from 3 to 5.")
  }
  mean_draw <- function(n) stats::rnorm(n, 0, priors$intercept_sd)
  beta_draw <- function(n) stats::rnorm(n, 0, priors$beta_sd)
  scale_draw <- function(n) abs(stats::rnorm(n, 0, priors$scale_sd))
  out <- list(
    lambda = c(1, stats::rnorm(2, priors$loading_mean, priors$loading_sd)),
    item_sd = rbind(rep(1, 3), stats::rlnorm(3, 0, priors$item_log_sd))
  )
  for (j in seq_len(3L)) {
    out[[paste0("threshold_", j)]] <- sort(stats::rnorm(K[j] - 1L, 0, priors$threshold_sd))
  }
  if (model %in% c("latent_covariate", "latent_outcome")) {
    out$mu_v <- mean_draw(2)
    out$sigma_v <- scale_draw(2)
    out$alpha <- mean_draw(2)
    out$beta <- beta_draw(2)
    out$sigma_y <- scale_draw(2)
    if (model == "latent_covariate") out$mu_v[1] <- 0
    if (model == "latent_outcome") out$alpha[1] <- -out$beta[1] * out$mu_v[1]
  } else {
    out$mu_q <- matrix(mean_draw(4), 2, 2)
    out$sigma_q <- matrix(scale_draw(4), 2, 2)
    # For a 2 by 2 matrix, LKJ(eta) gives (rho + 1) / 2 ~ beta(eta, eta).
    out$rho_q <- 2 * stats::rbeta(2, priors$lkj_shape, priors$lkj_shape) - 1
    out$alpha_m <- mean_draw(2)
    out$B <- matrix(beta_draw(4), 2, 2)
    out$alpha_y <- mean_draw(2)
    out$D <- matrix(beta_draw(4), 2, 2)
    out$d <- beta_draw(2)
    out$sigma_m <- scale_draw(2)
    out$sigma_y <- scale_draw(2)
    if (model == "latent_mediator") {
      out$alpha_m[1] <- -sum(out$B[1, ] * out$mu_q[1, ])
    } else {
      out$mu_q[1, 1] <- 0
    }
  }
  out
}
