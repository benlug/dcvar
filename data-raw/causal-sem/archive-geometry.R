# Archive the completed geometry experiments. Do not run a new fit.
source("data-raw/causal-sem/validation-helpers.R")
validation_setup()
arguments <- commandArgs(trailingOnly = TRUE)
if (length(arguments) != 1L) stop("Supply a new archive directory.")
output <- arguments[1]
if (dir.exists(output)) stop("Use a new archive directory.")
dir.create(output, recursive = TRUE)

copy_artifacts <- function(model, root, files) {
  for (name in files) {
    source <- file.path(root, name)
    destination <- file.path(output, model, name)
    if (!file.exists(source)) stop("An experiment artifact is missing: ", source)
    dir.create(dirname(destination), recursive = TRUE, showWarnings = FALSE)
    if (!file.copy(source, destination)) stop("An artifact cannot be copied.")
    if (grepl("[.]log$", destination)) {
      lines <- readLines(destination, warn = FALSE)
      writeLines(sub("[[:blank:]]+$", "", lines), destination, useBytes = TRUE)
    }
  }
}
copy_artifacts("latent_covariate", "/tmp/dcvar-conditional-covariate",
  c("summary.md", "derivation.md", "density-equivalence.csv", "diagnostics-all.csv",
    "comparison-matched.csv", "reference/comparison.csv", "prior-only-check.txt",
    "run.R", "run-original.R", "run-reference.R", "prior-only.R", "report.R"))
copy_artifacts("latent_mediator", "/tmp/dcvar-causal-centered-check",
  c("geometry-report.md", "geometry-comparison.csv", "stress-conditioned-effects.csv",
    "conditioned-normal-full-comparison.csv", "divergence-parameter-summary.csv",
    "check-density-and-prior.R", "summarise-geometry-checks.R",
    "density-and-prior-check.log"))
copy_artifacts("latent_mediator", "/tmp",
  paste0("run-dcvar-causal-", c("centered-check", "conditioned-check",
    "conditioned-long", "conditioned-normal"), ".R"))
copy_artifacts("latent_mediator_baseline", "/tmp",
  c("dcvar-causal-lmb-density.R", "dcvar-causal-lmb-geometry.R",
    "dcvar-causal-lmb-conditioned.R", "dcvar-causal-lmb-conditioned-prior.R",
    "dcvar-causal-lmb-conditioned-reference.R", "dcvar-causal-lmb-conditioned-prior.log",
    "dcvar-causal-lmb-conditioned-reference.csv",
    "dcvar-causal-lmb-conditioned-reference-parameters.csv"))

fits <- list(
  latent_covariate = c(
    original_short = "/tmp/dcvar-causal-sbc-run/fit-seed-190423.rds",
    original_long = "/tmp/dcvar-conditional-covariate/original/conditional-fit.rds",
    conditioned_long = "/tmp/dcvar-conditional-covariate/conditional-fit.rds",
    conditioned_normal = "/tmp/dcvar-conditional-covariate/reference/conditional-fit.rds"),
  latent_mediator = c(
    original_short = "/tmp/dcvar-causal-sbc-run/fit-seed-190427.rds",
    centered_short = "/tmp/dcvar-causal-centered-check/fit-centered-500.rds",
    conditioned_short = "/tmp/dcvar-causal-centered-check/fit-conditioned-500.rds",
    conditioned_long = "/tmp/dcvar-causal-centered-check/fit-conditioned-2000.rds",
    conditioned_normal = "/tmp/dcvar-causal-centered-check/fit-conditioned-normal-2000.rds"),
  latent_mediator_baseline = c(
    original_short = "/tmp/dcvar-causal-sbc-run/fit-seed-190428.rds",
    original_long = "/tmp/dcvar-causal-lmb-geometry.rds",
    conditioned_long = "/tmp/dcvar-causal-lmb-conditioned.rds",
    conditioned_normal = "/tmp/dcvar-causal-lmb-conditioned-reference.rds")
)
rows <- list()
for (model in names(fits)) {
  directory <- file.path(output, model)
  dir.create(file.path(directory, "functions"), showWarnings = FALSE)
  file.copy("inst/stan/functions/causal_measurement.stan",
    file.path(directory, "functions"))
  file.copy(file.path("inst/stan", paste0("causal_sem_", model, ".stan")),
    file.path(directory, "production.stan"))
  for (case in names(fits[[model]])) {
    fit_path <- fits[[model]][[case]]
    saved <- readRDS(fit_path)
    fit <- if (inherits(saved, "dcvar_causal_sem_fit")) saved else saved$fit
    diagnostic <- dcvar_diagnostics(fit)
    metadata <- fit$fit$metadata()
    settings <- metadata[intersect(c("seed", "iter_warmup", "iter_sampling",
      "max_treedepth", "adapt_delta", "num_chains", "id", "thin"), names(metadata))]
    source_path <- file.path(directory, paste0(case, ".stan"))
    writeLines(strsplit(paste(fit$fit$code(), collapse = "\n"), "\n",
      fixed = TRUE)[[1]], source_path, useBytes = TRUE)
    input_path <- file.path(directory, paste0(case, "-input.rds"))
    saveRDS(list(stan_data = fit$stan_data, meta = fit$meta, priors = fit$priors,
      sampler = settings), input_path, version = 3)
    utils::write.csv(causal_effects(fit),
      file.path(directory, paste0(case, "-effects.csv")), row.names = FALSE)
    rows[[length(rows) + 1L]] <- data.frame(model = model, case = case,
      source_fit = fit_path, seed = metadata$seed[1],
      warmup = metadata$iter_warmup[1], sampling = metadata$iter_sampling[1],
      chains = fit$fit$num_chains(),
      divergences = diagnostic$n_divergent,
      max_treedepth = diagnostic$n_max_treedepth,
      max_rhat = diagnostic$max_rhat, min_ess_bulk = diagnostic$min_ess_bulk,
      min_ess_tail = diagnostic$min_ess_tail,
      min_ebfmi = min(diagnostic$ebfmi), max_ebfmi = max(diagnostic$ebfmi),
      incomplete_diagnostics = diagnostic$incomplete_diagnostics,
      source = file.path(model, basename(source_path)),
      source_md5 = unname(tools::md5sum(source_path)),
      input = file.path(model, basename(input_path)),
      input_md5 = unname(tools::md5sum(input_path)))
  }
}
utils::write.csv(do.call(rbind, rows), file.path(output, "diagnostics.csv"), row.names = FALSE)
files <- list.files(output, recursive = TRUE, full.names = TRUE)
utils::write.csv(data.frame(file = substring(files, nchar(output) + 2L),
  md5 = unname(tools::md5sum(files))), file.path(output, "files.csv"), row.names = FALSE)
manifest_lines <- c("Geometry experiment archive", "",
  "The CSV contains diagnostics from complete saved fits.",
  "Each Stan source comes from its saved fit object.",
  "Each input file contains Stan data, priors, role metadata, and sampler settings.",
  "production.stan is the final package source at archive time.",
  "The original experiment scripts retain their original paths.",
  "Some source builders require the original package source.",
  "Use the stored Stan snapshots to repeat a fit after a source change.",
  "The complete fits stay at the source paths in diagnostics.csv.", "",
  capture.output(sessionInfo()))
writeLines(sub("[[:blank:]]+$", "", manifest_lines), file.path(output, "manifest.txt"))
print(do.call(rbind, rows)[c("model", "case", "divergences", "max_rhat", "min_ebfmi")])
