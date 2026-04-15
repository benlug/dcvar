## Test environments

- Local Ubuntu development branch, R 4.5.2
- GitHub Actions:
  - ubuntu-latest (release, devel, oldrel-1)
  - macos-latest (release)
  - windows-latest (release)

## R CMD check results

Local `R CMD check --as-cran` on `dcvar_0.1.0.tar.gz`:

- 0 errors
- 0 package-caused warnings
- 0 package-caused notes

The local sandboxed machine reported:

- 1 warning: `qpdf` not available
- 3 notes:
  - no network access for CRAN incoming / URL validation
  - unable to verify current time
  - `tidy` / `V8` unavailable for HTML validation helpers

These are local environment limitations rather than package issues.

## Notes for the reviewer

1. `cmdstanr` is in `Suggests`, not `Imports`

`dcvar` uses `rstan` as the default and fully supported backend in `Imports`.
`cmdstanr` is optional and only enables an additional backend path for users
who have CmdStan installed. The package installs, loads, runs examples, builds
vignettes, and passes checks without `cmdstanr`.

`cmdstanr` is declared in `Suggests` and resolved through:

- `Additional_repositories: https://stan-dev.r-universe.dev`

2. Stan backend behavior

Bundled models are compiled at use time, not at install time. Examples and
tests are written to keep CRAN-time execution lightweight, while vignettes and
the check-visible surface still exercise real model-fitting paths.

3. Scope of the package

The package implements Gaussian-copula VAR workflows. The SEM extension now
supports normal and exponential latent innovation margins; multilevel remains
normal-only.
