source("data-raw/causal-sem/paper-vita.R")

paper_covariance_formula <- function(condition, group) {
  p <- condition
  loading <- p$loading
  if (p$model %in% c("latent_covariate", "latent_outcome")) {
    result <- matrix(loading^2, 4, 4)
    result[4, 1:3] <- result[1:3, 4] <- loading * p$beta[group]
    diag(result) <- 1
    return(result)
  }
  q <- matrix(c(1, p$rho, p$rho, 1), 2)
  B <- c(p$B[group, 1], 0)
  D <- c(0, p$D[group, 2])
  d <- p$d[group]
  residual_m <- if (p$model == "latent_mediator") 1 - B[1]^2 else 1
  covariance_qm <- as.vector(q %*% B)
  variance_m <- sum(B * covariance_qm) + residual_m
  covariance_qy <- as.vector(q %*% D) + d * covariance_qm
  covariance_my <- sum(D * covariance_qm) + d * variance_m
  variance_y <- as.numeric(t(D) %*% q %*% D) +
    2 * d * sum(D * covariance_qm) + d^2 * variance_m + 1
  joint <- matrix(0, 4, 4)
  joint[1:2, 1:2] <- q
  joint[1:2, 3] <- joint[3, 1:2] <- covariance_qm
  joint[1:2, 4] <- joint[4, 1:2] <- covariance_qy
  joint[3, 3] <- variance_m
  joint[4, 4] <- variance_y
  joint[3, 4] <- joint[4, 3] <- covariance_my
  measurement <- matrix(0, 6, 4)
  if (p$model == "latent_mediator") {
    measurement[1:3, 3] <- loading
    measurement[4, 4] <- measurement[5, 1] <- measurement[6, 2] <- 1
  } else {
    measurement[1:3, 1] <- loading
    measurement[4, 3] <- measurement[5, 4] <- measurement[6, 2] <- 1
  }
  measurement %*% joint %*% t(measurement) + diag(c(rep(1 - loading^2, 3), 0, 0, 0))
}

check_paper_covariance <- function() {
  extra_library <- Sys.getenv("DCVAR_VALIDATION_LIB", "")
  if (nzchar(extra_library)) .libPaths(c(extra_library, .libPaths()))
  models <- c("latent_covariate", "latent_outcome", "latent_mediator", "latent_mediator_baseline")
  errors <- stats::setNames(numeric(length(models)), models)
  for (model in models) {
    condition <- paper_condition(model, n = 10L)
    generated <- paper_vita_sample(condition, seed = 2026, gaussian_reference = TRUE)
    for (g in 1:2) {
      errors[model] <- max(errors[model], abs(generated$covariance[[g]] -
                            paper_covariance_formula(condition, g)))
    }
    if (model %in% c("latent_mediator", "latent_mediator_baseline")) {
      changed <- condition
      changed$B[, 2] <- 0.9
      changed$D[, 1] <- -0.5
      ignored <- paper_vita_sample(changed, seed = 2026, gaussian_reference = TRUE)
      stopifnot(identical(generated$covariance, ignored$covariance))
    }
  }
  stopifnot(all(errors < 1e-12))
  list(passed = TRUE, max_error = errors, session = utils::sessionInfo())
}

if (sys.nframe() == 0L) {
  result <- check_paper_covariance()
  print(result[c("passed", "max_error")])
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args)) saveRDS(result, args[1])
}
