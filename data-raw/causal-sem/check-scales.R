# Check the item scale without using the fitted model or its simulator.

category_prob <- function(factor, loading, thresholds, residual_sd) {
  diff(stats::pnorm((c(-Inf, thresholds, Inf) - loading * factor) /
                    residual_sd))
}

ordinal_moments <- function(mu, sd, loading, thresholds, residual_sd) {
  conditional <- function(z, j) {
    p <- category_prob(mu + sd * z, loading[j], thresholds[[j]], residual_sd[j])
    sum(seq_along(p) * p)
  }
  expected <- vapply(seq_along(loading), function(j) {
    stats::integrate(Vectorize(function(z) conditional(z, j) * stats::dnorm(z)),
                     -9, 9, rel.tol = 1e-10)$value
  }, numeric(1))
  pairs <- utils::combn(seq_along(loading), 2)
  covariance <- apply(pairs, 2, function(pair) {
    cross <- stats::integrate(Vectorize(function(z) {
      conditional(z, pair[1]) * conditional(z, pair[2]) * stats::dnorm(z)
    }), -9, 9, rel.tol = 1e-10)$value
    cross - prod(expected[pair])
  })
  c(expected, covariance)
}

check_scales <- function(tolerance = 1e-8) {
  loading <- c(0.8, 0.9, 0.7)
  threshold <- list(c(-0.7, 0.6), c(-1.1, -0.2, 0.9), c(-1.2, -0.3, 0.4, 1.1))
  mu <- c(0.4, 0.9)
  sd <- c(0.8, 1.3)
  residual <- rbind(c(0.6, 0.5, 0.7), c(0.9, 0.4, 1.1))
  factor_scale <- loading[1] / residual[1, 1]
  lambda_new <- loading / residual[1, ] / factor_scale
  threshold_new <- lapply(seq_along(loading), function(j) {
    (threshold[[j]] - loading[j] * mu[1]) / residual[1, j]
  })
  residual_new <- sweep(residual, 2, residual[1, ], "/")
  error <- c(conditional = 0, marginal = 0, covariance = 0, ordinal = 0)
  for (g in 1:2) {
    total_sd <- sqrt(loading^2 * sd[g]^2 + residual[g, ]^2)
    delta <- 1 / total_sd
    for (j in seq_along(loading)) {
      for (l in seq(-2, 2, length.out = 11)) {
        original <- category_prob(l, loading[j], threshold[[j]], residual[g, j])
        theta <- category_prob(factor_scale * (l - mu[1]), lambda_new[j],
                               threshold_new[[j]], residual_new[g, j])
        delta_prob <- category_prob(l, loading[j] * delta[j],
                                    threshold[[j]] * delta[j],
                                    residual[g, j] * delta[j])
        error["conditional"] <- max(error["conditional"], abs(original - theta),
                                      abs(original - delta_prob))
      }
      original <- diff(stats::pnorm((c(-Inf, threshold[[j]], Inf) -
                                    loading[j] * mu[g]) / total_sd[j]))
      transformed <- diff(stats::pnorm((c(-Inf, threshold_new[[j]], Inf) -
                      lambda_new[j] * factor_scale * (mu[g] - mu[1])) /
                      sqrt((lambda_new[j] * factor_scale * sd[g])^2 +
                           residual_new[g, j]^2)))
      error["marginal"] <- max(error["marginal"], abs(original - transformed))
    }
    original_cov <- tcrossprod(loading) * sd[g]^2 + diag(residual[g, ]^2)
    theta_cov <- tcrossprod(lambda_new) * (factor_scale * sd[g])^2 +
      diag(residual_new[g, ]^2)
    expected_cov <- original_cov / outer(residual[1, ], residual[1, ])
    delta_cov <- original_cov * outer(delta, delta)
    error["covariance"] <- max(error["covariance"], abs(theta_cov - expected_cov),
                               abs(diag(delta_cov) - 1))
    original_moments <- ordinal_moments(mu[g], sd[g], loading, threshold,
                                        residual[g, ])
    theta_moments <- ordinal_moments(factor_scale * (mu[g] - mu[1]),
                                     factor_scale * sd[g], lambda_new,
                                     threshold_new, residual_new[g, ])
    error["ordinal"] <- max(error["ordinal"], abs(original_moments - theta_moments))
  }
  stopifnot(lambda_new[1] == 1, all(residual_new[1, ] == 1), all(error < tolerance))

  # Transform a latent covariate and a latent outcome in separate checks.
  alpha <- c(-0.2, 0.7)
  beta <- c(0.3, 0.8)
  weights <- c(0.3, 0.7)
  ate <- diff(alpha) + diff(beta) * sum(weights * mu)
  new_alpha <- alpha + beta * mu[1]
  new_beta <- beta / factor_scale
  new_mu <- factor_scale * (mu - mu[1])
  stopifnot(abs(diff(new_alpha) + diff(new_beta) * sum(weights * new_mu) - ate) < tolerance)
  outcome_alpha <- factor_scale * (alpha - mu[1])
  outcome_beta <- factor_scale * beta
  stopifnot(abs(diff(outcome_alpha) + diff(outcome_beta) * sum(weights * mu) -
                 factor_scale * ate) < tolerance)

  # Swap labels algebraically on one fixed scale.
  reversed <- diff(rev(alpha)) + diff(rev(beta)) * sum(rev(weights) * rev(mu))
  stopifnot(abs(reversed + ate) < tolerance)

  mediation <- function(p, weights) {
    mean_q <- colSums(p$mu_q * weights)
    mean_m <- sum(weights * (p$a_m + rowSums(p$B * p$mu_q)))
    base_y <- p$a_y + as.vector(p$D %*% mean_q)
    total_y <- base_y + p$d * (p$a_m + as.vector(p$B %*% mean_q))
    direct_y <- base_y + p$d * mean_m
    c(ATE = diff(total_y), ADE = diff(direct_y), AIE = diff(total_y - direct_y))
  }
  p <- list(mu_q = matrix(c(0.3, -0.2, 0.7, 0.4), 2, byrow = TRUE),
            a_m = c(-0.2, 0.6), B = matrix(c(0.5, 0.2, 0.8, -0.1), 2, byrow = TRUE),
            a_y = c(0.1, 0.9), D = matrix(c(0.2, 0.4, -0.1, 0.6), 2, byrow = TRUE),
            d = c(0.3, 0.7))
  effects <- mediation(p, weights)
  transformed_m <- p
  transformed_m$a_m <- factor_scale * (p$a_m - 0.4)
  transformed_m$B <- factor_scale * p$B
  transformed_m$a_y <- p$a_y + p$d * 0.4
  transformed_m$d <- p$d / factor_scale
  stopifnot(max(abs(mediation(transformed_m, weights) - effects)) < tolerance)
  transformed_q <- p
  transformed_q$mu_q[, 1] <- factor_scale * (p$mu_q[, 1] - p$mu_q[1, 1])
  transformed_q$a_m <- p$a_m + p$B[, 1] * p$mu_q[1, 1]
  transformed_q$B[, 1] <- p$B[, 1] / factor_scale
  transformed_q$a_y <- p$a_y + p$D[, 1] * p$mu_q[1, 1]
  transformed_q$D[, 1] <- p$D[, 1] / factor_scale
  stopifnot(max(abs(mediation(transformed_q, weights) - effects)) < tolerance)
  swapped <- lapply(p, function(value) {
    if (is.matrix(value)) value[2:1, , drop = FALSE] else rev(value)
  })
  stopifnot(max(abs(mediation(swapped, rev(weights)) + effects)) < tolerance,
            abs(effects["ATE"] - effects["ADE"] - effects["AIE"]) < tolerance)

  # Retain the change in the manifest mediator mean in paper study 3B.
  loading_3b <- 0.8
  nominal_mean_m <- 0.3
  mediator_path <- 0.5
  generated_aie <- mediator_path * loading_3b * nominal_mean_m
  nominal_aie <- mediator_path * nominal_mean_m
  stopifnot(abs(generated_aie - 0.12) < tolerance,
            abs(nominal_aie - 0.15) < tolerance,
            abs(generated_aie - nominal_aie + 0.03) < tolerance)
  list(passed = TRUE, max_error = error, factor_scale = factor_scale,
       paper_3b = c(nominal_AIE = nominal_aie, generated_AIE = generated_aie),
       checked_at = Sys.time(), session = utils::sessionInfo())
}

if (sys.nframe() == 0L) {
  result <- check_scales()
  print(result[c("passed", "max_error", "factor_scale", "paper_3b")])
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args)) saveRDS(result, args[1])
}
