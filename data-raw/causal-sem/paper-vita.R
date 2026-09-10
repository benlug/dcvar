# Load only the two generator functions. Do not run the paper simulation grid.

paper_condition <- function(model, n = 250L, loading = 0.8,
                            copula = c("gauss", "gauss"), categories = 3L,
                            threshold_shape = "symmetric") {
  stopifnot(length(model) == 1L, length(n) == 1L, is.finite(n), n == floor(n),
            length(loading) == 1L, length(categories) == 1L,
            model %in% c("latent_covariate", "latent_outcome", "latent_mediator",
                         "latent_mediator_baseline"), n > 0,
            loading %in% c(0.8, 0.9), categories %in% 3:5,
            all(copula %in% c("gauss", "clayton", "joe")), length(copula) == 2L)
  bound <- if (model == "latent_mediator") 1.8 else 3
  if (threshold_shape == "symmetric") {
    thresholds <- seq(-bound, bound, length.out = categories + 1L)
  } else {
    probabilities <- seq(0, 1, length.out = categories + 1L)^1.5
    thresholds <- stats::qnorm(probabilities)
    if (threshold_shape == "right") thresholds <- -rev(thresholds)
    if (!threshold_shape %in% c("left", "right")) stop("Use symmetric, left, or right.")
  }
  thresholds[c(1, length(thresholds))] <- c(-Inf, Inf)
  list(model = model, n = as.integer(n), loading = loading, copula = copula,
       beta = c(0.3, 0.5), mean_v = c(0, 0), mean_y = c(0, 1),
       mean_q = matrix(0, 2, 2),
       mean_m = c(0, if (model == "latent_mediator_baseline") 0.3 else 1),
       B = matrix(c(0.8, 0, 0.8, 0), 2, 2, byrow = TRUE),
       D = matrix(c(0, 0.75, 0, 0.75), 2, 2, byrow = TRUE),
       d = c(0.5, 0.5), rho = 0.3, thresholds = thresholds,
       threshold_shape = threshold_shape, margins = "normal")
}

paper_reference_effects <- function(condition) {
  p <- condition
  loading <- p$loading
  residual <- sqrt(1 - loading^2)
  factor_scale <- loading / residual
  if (p$model == "latent_covariate") {
    mean_v <- p$mean_v / loading
    alpha <- p$mean_y - p$beta * mean_v
    nominal_alpha <- p$mean_y - p$beta * p$mean_v
    return(list(nominal = c(ATE = diff(nominal_alpha) + diff(p$beta) * mean(p$mean_v),
                             interaction = diff(p$beta)),
                 identified = c(ATE = diff(alpha) + diff(p$beta) * mean(mean_v),
                                  interaction = diff(p$beta) / factor_scale)))
  }
  if (p$model == "latent_outcome") {
    mean_factor <- c(p$mean_y[1] / loading, p$mean_y[2])
    alpha <- mean_factor - p$beta * p$mean_v
    nominal_alpha <- p$mean_y - p$beta * p$mean_v
    return(list(nominal = c(ATE = diff(nominal_alpha) + diff(p$beta) * mean(p$mean_v),
                             interaction = diff(p$beta)),
                 identified = c(ATE = (diff(alpha) + diff(p$beta) * mean(p$mean_v)) * factor_scale,
                                  interaction = diff(p$beta) * factor_scale)))
  }
  mediation_effect <- function(mean_q, mean_m, B, D) {
    alpha_m <- mean_m - rowSums(B * mean_q)
    alpha_y <- p$mean_y - rowSums(D * mean_q) - p$d * mean_m
    target_q <- colMeans(mean_q)
    adjusted_m <- alpha_m + as.vector(B %*% target_q)
    total <- alpha_y + as.vector(D %*% target_q) + p$d * adjusted_m
    direct <- alpha_y + as.vector(D %*% target_q) + p$d * mean(mean_m)
    c(ATE = diff(total), ADE = diff(direct), AIE = diff(total) - diff(direct))
  }
  nominal <- mediation_effect(p$mean_q, p$mean_m, p$B, p$D)
  actual_q <- p$mean_q
  if (p$model == "latent_mediator") {
    # The source passes this mean to each continuous indicator.
    actual_m <- p$mean_m / loading
  } else {
    # The source scales the treated manifest mediator mean only.
    actual_m <- p$mean_m
    actual_m[2] <- actual_m[2] * loading
    actual_q[, 1] <- actual_q[, 1] / loading
  }
  actual_B <- p$B
  actual_D <- p$D
  actual_B[, 2] <- 0
  actual_D[, 1] <- 0
  list(nominal = nominal,
       identified = mediation_effect(actual_q, actual_m, actual_B, actual_D))
}

paper_fit_args <- function(data, model) {
  prefix <- switch(model, latent_covariate = "z", latent_outcome = "y", "m")
  indicators <- paste0(prefix, 1:3)
  roles <- switch(model,
    latent_covariate = list(outcome = "y", covariate = "xi", indicators = list(xi = indicators)),
    latent_outcome = list(outcome = "eta", covariate = "z", indicators = list(eta = indicators)),
    latent_mediator = list(outcome = "Y", mediator = "M", mediator_baseline = "Mpre",
      outcome_baseline = "Ypre", indicators = list(M = indicators)),
    latent_mediator_baseline = list(outcome = "Y", mediator = "M", mediator_baseline = "Mpre",
      outcome_baseline = "Ypre", indicators = list(Mpre = indicators)))
  c(list(data = data, model = model, treatment = "x", control = 0, treated = 1,
           standardize = FALSE, target_weights = c(0.5, 0.5)), roles)
}

paper_vita_sample <- function(condition, seed,
                               paper_root = "../sem_causal_effects",
                               Nmax = 100000L, gaussian_reference = FALSE) {
  packages <- c("lavaan", if (!gaussian_reference) c("covsim", "rvinecopulib"))
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop("Install these optional packages: ", paste(missing, collapse = ", "))
  if (gaussian_reference && any(condition$copula != "gauss")) {
    stop("The Gaussian reference requires Gaussian copulas in both groups.")
  }
  study <- switch(condition$model, latent_covariate = "study_1",
                   latent_outcome = "study_2", latent_mediator = "study_path_A",
                   latent_mediator_baseline = "study_path_B")
  paths <- file.path(paper_root, "code", study, "main", "functions",
                    c("generate_model.R", "generate_data.R"))
  if (!all(file.exists(paths))) stop("The paper generator files are missing.")
  environment <- new.env(parent = baseenv())
  environment$sem <- lavaan::sem
  environment$fitted <- lavaan::fitted
  environment$as_tibble <- function(x) as.data.frame(x)
  if (!gaussian_reference) {
    for (name in c("bicop_dist", "dvine_structure", "vinecop_dist", "rvine")) {
      environment[[name]] <- getExportedValue("rvinecopulib", name)
    }
    environment$vita <- covsim::vita
  }
  for (path in paths) sys.source(path, envir = environment)
  set.seed(seed)
  p <- condition
  groups <- covariance <- calibration <- vector("list", 2)
  for (g in 1:2) {
    simple <- p$model %in% c("latent_covariate", "latent_outcome")
    model_args <- if (simple) {
      list(beta_z = p$beta[g], load_z = p$loading)
    } else {
      list(sMZ = p$B[g, 1], sMW = p$B[g, 2], sYW = p$D[g, 2],
           sYM = p$d[g], sYZ = p$D[g, 1], MpreYpre_corr = p$rho, loading = p$loading)
    }
    covariance[[g]] <- do.call(environment$generate_model, model_args)
    if (simple) {
      mean_y <- p$mean_y[g]
      if (p$model == "latent_outcome" && g == 2L) mean_y <- p$loading * mean_y
      generator_args <- list(n = p$n, mean_y = mean_y, mean_z = p$mean_v[g],
        Sigma = covariance[[g]], margin_dist_y = "norm", margin_dist_z = "norm",
        copula_type_ZY = p$copula[g], Nmax = Nmax)
      generator_args[[if (p$model == "latent_covariate") "num_z_vars" else "num_y_vars"]] <- 3L
      means <- if (p$model == "latent_covariate") c(rep(p$mean_v[g], 3), mean_y) else c(rep(mean_y, 3), p$mean_v[g])
    } else {
      mean_m <- p$mean_m[g]
      if (p$model == "latent_mediator_baseline" && g == 2L) mean_m <- mean_m * p$loading
      generator_args <- list(n = p$n, mean_Ypre = p$mean_q[g, 2],
        mean_Y = p$mean_y[g], mean_Mpre = p$mean_q[g, 1], mean_M = mean_m,
        Sigma = covariance[[g]], margin_dist_M = "norm", margin_dist_Ypre = "norm",
        copula = p$copula[g], Nmax = Nmax, num_m_vars = 3L)
      means <- if (p$model == "latent_mediator") {
        c(rep(mean_m, 3), p$mean_y[g], p$mean_q[g, ])
      } else {
        c(rep(p$mean_q[g, 1], 3), mean_m, p$mean_y[g], p$mean_q[g, 2])
      }
    }
    if (gaussian_reference) {
      sample <- matrix(stats::rnorm(p$n * length(means)), p$n) %*% chol(covariance[[g]])
      sample <- sweep(sample, 2, means, "+")
      colnames(sample) <- colnames(covariance[[g]])
      groups[[g]] <- as.data.frame(sample)
    } else {
      generated <- do.call(environment$generate_data, generator_args)
      groups[[g]] <- as.data.frame(generated[[1]])
      calibration[[g]] <- generated[[2]]
    }
    groups[[g]]$x <- g - 1L
  }
  continuous <- do.call(rbind, groups)
  data <- continuous
  prefix <- switch(p$model, latent_covariate = "z", latent_outcome = "y", "m")
  indicators <- paste0(prefix, 1:3)
  for (name in indicators) {
    data[[name]] <- ordered(cut(data[[name]], p$thresholds, labels = FALSE),
                            levels = seq_len(length(p$thresholds) - 1L))
  }
  effects <- paper_reference_effects(p)
  list(data = data, continuous = continuous, parameters = list(), effects = effects$identified,
       args = paper_fit_args(data, p$model),
       metadata = list(condition = condition, seed = seed, Nmax = Nmax,
                         engine = if (gaussian_reference) "exact Gaussian reference" else "paper VITA",
                         source_paths = normalizePath(paths), source_md5 = tools::md5sum(paths),
                         package_versions = vapply(packages, function(name) as.character(utils::packageVersion(name)), character(1)),
                         nominal_effects = effects$nominal,
                         identified_reference_effects = effects$identified,
                         effect_scope = "Gaussian SEM reference; no observed person factors from VITA"),
       covariance = covariance, calibration = calibration)
}
