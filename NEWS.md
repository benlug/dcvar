# dcvar 0.1.0

## Scope and documentation

- Clarified that the package currently implements Gaussian-copula models only.
- Marked the multilevel and SEM variants as experimental extensions with a
  narrower post-estimation interface than the core single-level models.
- Documented PIT diagnostics as posterior-mean plug-in diagnostics and made the
  unsupported multilevel and SEM paths explicit in the help pages.
- Documented the current scope of `fitted()`, `predict()`, and `loo()`.

## Build and submission hygiene

- Excluded `data-raw/`, `cran-comments.md`, and the source-only Quarto
  walkthrough from source builds.
- Added package citation metadata and a draft `cran-comments.md`.
- Ignored local `*-test-local.log` artifacts in Git.

## Testing

- Added gamma and skew-normal fit coverage for `dcvar()` and `dcvar_hmm()`.
- Added PIT smoke coverage for gamma and skew-normal margins.
