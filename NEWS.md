# dcvar 0.4.0

## Breaking / scientifically relevant fixes

- **Simulators no longer negate the dependence for left-skewed margins.**
  The Stan likelihoods treat the copula uniform on a `skew_direction = -1`
  exponential/gamma dimension as `u = 1 - F(x_shifted)`, but `simulate_dcvar()`
  (and the shared mixed-margin helper used by `simulate_dcvar_multilevel()` and
  `simulate_dcvar_sem()`, plus the homogeneous exponential SEM path) negated
  the quantile-transformed variable *after* applying the copula. Data simulated
  with exactly one negative `skew_direction` therefore implied copula
  correlation `-rho` while `true_params` recorded `+rho`, so simulate-then-fit
  studies silently recovered the negated dependence trajectory. Simulations
  that used `c(1, 1)` or `c(-1, -1)` are unaffected (the flips cancel); any
  results generated with asymmetric skew directions should be re-run.
- The `eps_rep` posterior-predictive blocks in seven Stan models had the
  mirror-image inversion bug and replicated negated dependence on left-skewed
  dimensions; they now invert at the flipped uniform.
- `dependence_summary()` for HMM fits now averages Kendall's tau over states
  per draw (`sum_k gamma[t,k] * (2/pi) asin(rho_state[k])`) instead of applying
  `asin` to the gamma-weighted mean rho, which understated tau during regime
  transitions.
- `hmm_states()$viterbi` is now a genuine Viterbi (MAP) decoding on
  posterior-mean emission/transition log-probabilities. The previous
  most-frequent-sampled-path estimator degenerated to an arbitrary
  (lexicographically first) single draw's path whenever parameter uncertainty
  made nearly all sampled paths unique.
- HMM posterior predictive replicates (`eps_rep`) are now drawn from the
  regime mixture (a state sampled from the smoothed probabilities, then that
  state's rho) instead of the gamma-weighted mean rho.
- `.relative_bias()` (used by `compute_rho_metrics()` /
  `compute_param_metrics()`) now normalizes the mean error by the mean
  absolute true value. The previous pointwise form exploded to ~1e12 % when
  any true value was near zero (e.g. null-dependence conditions), silently
  corrupting aggregated summaries. It returns `NA` with a warning when all
  true values are near zero.

## Bug fixes

- `hmm_EG.stan` now exposes `sigma_exp` (plus `b_gq` and `rate_exp`) as
  saved generated quantities. Previously they were brace-local, so
  `summary()`, `var_params()`, `coef()`, and `plot_diagnostics()` errored and
  `pit_values()` silently returned all-NA for exponential-margin HMM fits.
- Tibble input no longer breaks multilevel fits: per-unit time-grid
  validation was silently skipped and the stored time axis corrupted, crashing
  `rho_trajectory()`, `fitted()`, and `predict()` downstream. All preparation
  functions now coerce to a plain data frame up front.
- Character and unordered-factor time columns are now rejected: they sort
  lexicographically ("T10" before "T2"), silently scrambling the VAR(1)
  ordering while making the spacing checks vacuous.
- Leading/trailing runs of two or more missing rows are now treated as edge
  trims rather than adjacency-breaking interior gaps, so they no longer abort
  (or warn misleadingly) when removing them preserves adjacency.
- `loo()` and `dcvar_compare()` now support all multilevel fits:
  `multilevel_copula_var.stan` stores per-observation `log_lik` like its EG
  and mixed siblings, and the margin-based whitelist was removed.
  `dcvar_compare()` warns when an HMM fit (state-marginalized predictive
  density) is compared against a dynamic fit (conditioned on the smoothed
  latent rho trajectory), which can systematically favor the dynamic model.
- The `seed` argument now makes fits reproducible: default per-chain initial
  values are generated under a deterministic RNG (previously they came from
  the unseeded global R RNG, so two fits with the same seed differed).
- cmdstanr fits now survive `saveRDS()` and an R restart: posterior draws are
  eagerly read into the fit object after sampling instead of being lazily
  re-read from CSVs in the session tempdir.
- `interpret_rho_trajectory()` no longer errors on Clayton constant fits; it
  interprets the dependence via Kendall's tau.
- The Clayton copula log-density is computed in log space, avoiding overflow
  to `-Inf` (with NaN gradients) for `theta` greater than about 34 during
  warmup.
- PIT helpers abort with a clear message when a margin parameter group is
  missing from the Stan output instead of silently returning all-NA values.
- `simulate_dcvar_sem(n_time = 1)` no longer crashes with a subscript error;
  `prepare_multilevel_data()` requires at least 3 occasions per unit;
  `.safe_cor()` handles length-1 inputs; `Y[cc, ]` row removal keeps the
  matrix shape with a single remaining row.
- `print()` for constant-fit summaries no longer swaps the CI bounds when
  `probs` are not passed in ascending order.
- `draws()` validates its `format` argument instead of silently returning a
  `draws_array`.
- Divergence/treedepth counts are reported as `NA` (unknown) rather than 0
  for fits without sampler diagnostics.
- `plot_hmm_states()` pins the state factor levels so Viterbi point colors
  match the fill colors when the MAP path skips a state;
  `plot_trajectories()` forwards `...` only to the scenario generators that
  accept each argument.
- Warnings are emitted for prior arguments that a configuration ignores
  (`prior_sigma_eps_rate` with no normal margin; `prior_z_rho_sd` with the
  Clayton copula), and the recorded `priors` list omits unused entries.

## Validation

- All bivariate Stan models now declare `int<lower=2, upper=2> D`, rejecting
  D > 2 data that the hard-wired bivariate copula code would silently
  mishandle.

## Documentation

- Documented the plug-in nature of `predict()` intervals, the raw-data scale
  of simulator `true_params` (fit with `standardize = FALSE` for round-trip
  comparisons), the lognormal scale prior for exponential/gamma multilevel
  margins, and `prior_rho_init_sd` as the covariate-model dependence
  intercept (`beta_0`) prior.
- Added the covariate model family to the README table and the pkgdown
  reference index (which previously failed to build due to missing topics).

## Tests

- New regression tests pin the simulator copula orientation for every
  `skew_direction` combination, seed reproducibility, tibble input, time
  column validation, and edge-gap handling. The Clayton-normal constant model
  and `dcvar_covariate()` (drift and no-drift) now get real MCMC smoke fits
  in CI.

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
