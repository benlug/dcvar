# dcvar 0.3.1

## Bug fixes

- `summary()` printouts for HMM and multilevel fits now include the marginal
  scale/shape parameters (e.g. `sigma_exp`, `omega`, `delta`, `sigma_gam`,
  `shape_gam`) for non-normal and per-variable (mixed) margins. Previously the
  HMM summary omitted all margin scale parameters, and the multilevel summary
  dropped them for mixed margins.
- `simulate_breakpoint_data()` now records the breakpoint specification
  (`type`, plus `breakpoint` or `breakpoints`) in `true_params`, so the
  documented access no longer returns `NULL`.

## Validation

- `rho_decreasing()` and `rho_increasing()` now validate that `rho_start` and
  `rho_end` are single finite values in `[-1, 1]`, matching the other rho
  trajectory generators.

## Documentation

- Extensive roxygen, vignette, README, and CITATION updates so the
  documentation reflects per-variable (mixed) margin support across all model
  families (including the Clayton-copula constant model) and the correct
  `coef()` / `fitted()` / data-preparation contracts.

## Internal

- Removed unused helper code, made the exponential-margin diagnostic generated
  quantity (`b_gq`) report the clamped lower bound consistently across the
  exponential models (inference-neutral), and corrected Stan prior-hyperparameter
  comments. The Clayton-normal constant model's `sigma_eps` prior is left
  unchanged and its intentional divergence from the other constant models is now
  documented.

# dcvar 0.3.0

## Per-variable (mixed) margins

- All copula VAR model families now accept a length-2 `margins` vector so each
  variable can use its own marginal family, e.g.
  `margins = c("normal", "exponential")`. This exploits the copula's separation
  of the marginal distributions from the dependence structure. Supported across
  `dcvar_constant()`, `dcvar()`, `dcvar_hmm()`, `dcvar_multilevel()`, and
  `dcvar_sem()` (both the indicator and naive methods).
- Added generic mixed-margins Stan models that dispatch each dimension to its
  own marginal family and apply the copula on the CDF scale:
  `constant_mixed.stan`, `dcvar_mixed_ncp.stan`, `hmm_mixed.stan`,
  `multilevel_mixed.stan`, `sem_mixed.stan`, and `sem_naive_mixed.stan`.
- Mixed margins are also available with the **Clayton** copula for the constant
  model (`constant_mixed_clayton.stan`), previously limited to normal margins.
- `simulate_dcvar()`, `simulate_dcvar_multilevel()`, and `simulate_dcvar_sem()`
  accept the same per-variable `margins` vector so mixed-family data can be
  generated (for example for parameter-recovery studies).
- `coef()`, `var_params()`, `pit_values()`, and the diagnostics/plots report
  each dimension under its own family for mixed fits across all model families.

## Backward compatibility

- A single `margins` string is unchanged, and an all-identical `margins` vector
  (such as `c("normal", "normal")`) routes to the existing specialised
  single-family model, so prior results, tests, and the gamma shared-shape
  parameterisation are preserved exactly.
- The mixed multilevel and SEM models support all four families per dimension.
  Single-family multilevel and SEM fits keep their existing
  normal/exponential-only restriction (there is no specialised gamma or
  skew-normal model for those structures); request the other families through a
  per-variable `margins` vector instead.

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
