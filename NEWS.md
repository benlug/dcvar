# dcvar 0.1.0

## Scope and documentation

- Clarified that the package currently implements Gaussian-copula models only.
- Marked the multilevel and SEM variants as experimental extensions with
  narrower diagnostic support than the core single-level models.
- Documented PIT diagnostics as posterior-mean plug-in diagnostics and made the
  unsupported multilevel and SEM paths explicit in the help pages.
- Clarified the current scope of `loo()` support across model classes.

## Build and submission hygiene

- Excluded local `*-test-local.log` artifacts from source builds.
- Added package citation metadata.

## Testing

- Added skew-normal fit coverage for `dcvar()` and `dcvar_hmm()`.
- Added PIT smoke coverage for gamma and skew-normal margins.
