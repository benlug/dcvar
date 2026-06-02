# dcvar 0.3.0

## Per-variable (mixed) margins

- `dcvar_constant()`, `dcvar()`, and `dcvar_hmm()` now accept a length-2
  `margins` vector so each variable can use its own marginal family, e.g.
  `margins = c("normal", "exponential")`. This exploits the copula's separation
  of the marginal distributions from the dependence structure. Mixed margins
  currently require the Gaussian copula.
- Added generic mixed-margins Stan models (`constant_mixed.stan`,
  `dcvar_mixed_ncp.stan`, `hmm_mixed.stan`) that dispatch each dimension to its
  own marginal family and apply the Gaussian copula on the CDF scale. The mixed
  models cover the constant (single rho), dynamic (time-varying rho random
  walk), and regime-switching (state-specific rho) dependence structures.
- `simulate_dcvar()` accepts the same length-2 `margins` vector so mixed-family
  data can be generated (for example for parameter-recovery studies).
- `coef()`, `var_params()`, `pit_values()`, and the diagnostics/plots report
  each dimension under its own family for mixed fits.

## Backward compatibility

- A single `margins` string is unchanged, and an all-identical `margins` vector
  (such as `c("normal", "normal")`) routes to the existing specialised
  single-family model, so prior results, tests, and the gamma shared-shape
  parameterisation are preserved exactly.
- Mixed margins are implemented for the single-level time-series models
  (`dcvar_constant()`, `dcvar()`, `dcvar_hmm()`); the multilevel and SEM models
  still require a single margin family and report a clear error if given a
  per-variable vector.

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
