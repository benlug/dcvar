# dcvar 0.2.0

## Simulation model parity

- Added a constant Clayton-copula DC-VAR for normal margins via
  `dcvar_constant(copula = "clayton")`.
- Added exponential-margin support for `dcvar_multilevel()`.
- Added naive SEM score models via `dcvar_sem(method = "naive")` for normal
  and exponential margins.
- Added `dependence_summary()` for Kendall's tau summaries across Gaussian and
  Clayton copula fits.

## Infrastructure

- Added copula-family dispatch alongside the existing margin dispatch.
- Added bundled Stan models for the new Clayton, multilevel exponential, and
  naive SEM variants.
- Updated extractors, summaries, diagnostics, LOO support, and tests for the
  new model variants.

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
