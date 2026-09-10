# Causal SEM validation

Run these scripts from the repository root.
They test all four models through the public R API.
They use a fixed scale with `standardize = FALSE`.
They store each result before they start the next fit.
They keep failed fits in `status.csv`.
They do not replace rejected datasets.
Use a new output directory for each run.

Git excludes the `results/` directory.
It contains local output and optional run archives.
The execution record below describes those local runs.
The package does not load files from `data-raw`.

## Fast scale checks

```sh
TZ=Europe/Berlin Rscript data-raw/causal-sem/check-scales.R /tmp/causal-scales.rds
```

This check does not use the package simulator.
It compares conditional category probabilities and marginal category probabilities.
It also compares continuous covariance matrices and ordinal moments.
It checks affine factor changes and the sign after a group swap.
The group swap uses one fixed scale.
It does not fit a model with a new reference group.
`check-paper-covariance.R` compares independent structural covariance formulas
with the original lavaan generator for all four models.
It also checks the two unused mediation arguments.
This check needs the optional lavaan package.

Let the paper response be `U_j = lambda_j * L + error_j`.
Let `s_j0` be its residual SD in the control group.
Let `mu_0` be the factor mean in that group.
The new factor is `L_new = c * (L - mu_0)`.
Here, `c = lambda_1 / s_10`.
The new item response is `(U_j - lambda_j * mu_0) / s_j0`.
Thus:

- The new loading is `lambda_j / (s_j0 * c)`.
- The new threshold is `(tau_jk - lambda_j * mu_0) / s_j0`.
- The new residual SD is `s_jg / s_j0`.
- The Delta scale is `1 / sqrt(lambda_j^2 * var(L_g) + s_jg^2)`.

Delta scales the total item response.
It is not an item residual SD.
An outcome effect on the factor scale changes by `c`.
A slope for a latent covariate changes by `1 / c`.
A manifest outcome effect keeps its units.

`check-paper-measurement.R` also calls the original measurement generator.
Use an installed `EffectLiteR` package or a source file from its archive:

```sh
DCVAR_EFFECTLITER_SOURCE=/tmp/EffectLiteR/R/elr_generate_measurement_model.R \
  TZ=Europe/Berlin Rscript data-raw/causal-sem/check-paper-measurement.R \
  /tmp/causal-measurement.rds
```

The check records the source hash or package version.
The generated syntax fixes each control Delta scale at one.
It leaves each treatment Delta scale free.
It shares loadings and thresholds across groups.
The paper then replaces the first loading and first threshold with known values.
The package instead anchors the factor mean and loading.

## Prior draws, recovery, and calibration

```sh
DCVAR_VALIDATION_OUTPUT=/tmp/causal-priors \
  Rscript data-raw/causal-sem/run-priors.R
DCVAR_VALIDATION_OUTPUT=/tmp/causal-recovery \
  Rscript data-raw/causal-sem/run-recovery.R
DCVAR_VALIDATION_OUTPUT=/tmp/causal-sbc \
  Rscript data-raw/causal-sem/run-sbc.R
DCVAR_VALIDATION_OUTPUT=/tmp/causal-prior-sensitivity \
  Rscript data-raw/causal-sem/run-prior-sensitivity.R
```

The prior script needs no Stan fit.
It stores category frequencies, constant items, continuous quantiles, and effects.
The recovery script uses null effects and nonzero effects.
It checks the simulator effects against a separate formula before each fit.
Set `DCVAR_RECOVERY_CASES=nonzero` to run only the nonzero cases.
The sensitivity script fits the same dataset with slope prior SDs of 0.5, 1, and 2.

Simulation-based calibration, or SBC, draws parameters from the fit prior.
It then draws data from the model and fits those data.
It records ranks for structural parameters, measurement parameters, and effects.
It excludes fixed parameters from these ranks.
The rank table includes the diagnostic pass flag for each fit.
It retains at most 100 posterior draws for each rank by default.
Check effective sample size before you read the ranks.
Autocorrelation can distort a rank histogram.
Two replications check execution only.
They do not establish calibration or interval coverage.
Use many replications for a calibration claim.
The [Stan SBC guide](https://mc-stan.org/docs/stan-users-guide/simulation-based-calibration.html)
describes the rank check.

SBC can produce globally constant indicators.
The data contract rejects these datasets.
The result records the rejection and the constant item names.
Calibration for accepted datasets does not establish unconditional calibration.

Each fit stores effect draws, parameter draws, truth, warnings, and diagnostics.
Each new fit also stores its complete fit object before output extraction.
It also compares observed and predicted category frequencies in each group.
The CSV report includes bias, RMSE, coverage, coverage Monte Carlo error, and interval width.
Bias is the mean estimate minus the true effect.
Its unit is the effect unit.
The report also includes its absolute value.
Coverage uses successful fits only.
Read its denominator together with `status.csv`.
The plug-in coverage Monte Carlo error is zero when every fit covers the truth.
That value does not imply precise coverage.
The report also gives a Wilson interval and the Monte Carlo error at nominal 95% coverage.
For two replications, nominal coverage Monte Carlo error is 0.154.
Two covered values give a Wilson interval of about `[0.342, 1]`.

Reference fits require no divergences and no maximum tree depth events.
Relevant parameters and effects require `R-hat < 1.01`.
Effects require bulk ESS and tail ESS of at least 400.
Their mean Monte Carlo error must be at most 5% of their SD.
Inspect chain energy and predictions as well.
The scripts report these values.
They do not label a short run as a passed reference run.
The `diagnostic_pass` column applies these numeric limits.
It rejects missing diagnostics when the fit reports them.
Older records can have `NA` in `incomplete_diagnostics`.
Their pass flag still checks the stored numeric diagnostics.
It does not replace checks of identification or calibration.

Set these environment variables to change a run:

| Variable | Default |
| --- | --- |
| `DCVAR_VALIDATION_MODELS` | All four names, separated by commas |
| `DCVAR_VALIDATION_REPS` | 20; prior draws use 100; sensitivity uses 5 |
| `DCVAR_VALIDATION_N` | 400 persons in total |
| `DCVAR_VALIDATION_SEED` | 90421 |
| `DCVAR_VALIDATION_BACKEND` | `cmdstanr` |
| `DCVAR_VALIDATION_CHAINS` | 4 |
| `DCVAR_VALIDATION_CORES` | 2 |
| `DCVAR_VALIDATION_WARMUP` | 1000 per chain |
| `DCVAR_VALIDATION_SAMPLING` | 1000 per chain |
| `DCVAR_VALIDATION_RANK_DRAWS` | 100 |
| `DCVAR_VALIDATION_LIB` | No extra package library |
| `DCVAR_VALIDATION_OUTPUT` | A new path below the R temporary directory |

## Paper pilot

`paper-vita.R` loads only `generate_model.R` and `generate_data.R`.
It does not source `one_rep.R` or a main script.
It does not install packages, change the working directory, or start a cluster.
The adapter needs optional `lavaan`, `covsim`, and `rvinecopulib` packages.
They are not fit dependencies.
The [covsim manual](https://cran.r-project.org/web/packages/covsim/covsim.pdf)
defines the VITA calibration options.

```sh
DCVAR_VALIDATION_LIB=/tmp/dcvar-causal-validation-lib \
  DCVAR_VALIDATION_OUTPUT=/tmp/causal-paper \
  DCVAR_VALIDATION_N=500 \
  Rscript data-raw/causal-sem/run-paper-pilot.R
```

The pilot uses 250 persons per group, loading 0.8, and three symmetric categories.
It runs all four models with Gaussian, Clayton, and Joe copulas in both groups.
It uses normal margins from the main paper scripts.
It preserves the original vine structures.
Study 2 lets VITA select from one supplied family.
The other studies supply a vine object with mixed Gaussian pairs.
The metadata stores the seed, source hashes, package versions, condition, and calibration size.

Set `DCVAR_PAPER_ROOT` to change the paper path.
Its default is `../sem_causal_effects`.
Set `DCVAR_PAPER_COPULAS` to a comma-separated family list.
Set `DCVAR_PAPER_GENERATE_ONLY=true` to check generation without a fit.
Set `DCVAR_PAPER_NMAX` to change the calibration size.
Its default is 100000, as in the original code.
A smaller value reduces calibration precision.
`gaussian_reference = TRUE` in `paper_vita_sample()` draws exact normal data from
the original target covariance matrix.
This mode checks the reference moments.
It does not call VITA.

To fit saved pilot datasets, use `fit-paper-pilot.R`:

```sh
DCVAR_PAPER_INPUT=/tmp/causal-paper-generation \
  DCVAR_VALIDATION_OUTPUT=/tmp/causal-paper-fits \
  Rscript data-raw/causal-sem/fit-paper-pilot.R
```

This script keeps generation errors in its status table.
It fits the stored data without a new calibration draw.
Set `DCVAR_PAPER_RESUME=true` to reuse complete records in the output directory.
Stop the previous R process and its Stan children before a resume.
An interrupted shell can leave its child processes active.
Use separate output directories for concurrent workers.
The model filter also applies to generation errors.
Sampler seeds use positions in the complete source list.
They remain the same when you split models across workers.
`combine-runs.R NEW_OUTPUT RUN_1 RUN_2 ...` combines disjoint batches.
It rejects duplicate result IDs and keeps the source batch settings.
Use `archive-run.R RUN_DIRECTORY NEW_ARCHIVE_DIRECTORY` to retain small CSV reports
and a text manifest.
The complete fits stay at the run location.
`snapshot-stan-sources.R RUN_DIRECTORY` stores the actual source from each saved fit.
It also stores the Stan input data.
Its CSV records a hash for each source and data file.
This check works when the package source changes after a fit.
The archive script also copies these snapshots.
The archive removes trailing spaces from logs and manifests.

The local source has these details:

- Study 3A uses only loading 0.8 in its active main grid.
- Study 3B uses sample sizes 1000 and 250 in its active main grid.
- Both mediation generators ignore the `sMW` and `sYZ` arguments.
  The fit still estimates those paths.
- Study 3A passes `mean_M` directly to each indicator mean.
  With loading 0.8, a treatment indicator mean of 1 implies a factor mean of 1.25.
  The Gaussian reference AIE is then 0.625.
  The nominal AIE is 0.5.
- Study 3B multiplies the treated manifest mediator mean by the loading.
  With loading 0.8, its mean changes from 0.3 to 0.24.
  The Gaussian reference AIE is then 0.12.
  The nominal AIE is 0.15.

The adapter keeps these details.
It stores nominal effects and Gaussian reference effects separately.
The package simulator uses every documented structural path.
It simulates the identified Bayes model directly.
VITA provides no observed true person factors.
The VITA pilot tests sensitivity to distribution assumptions.
Its effect reports use the Gaussian reference targets.
For copula data, the reported bias is displacement from that reference.
It does not establish a true causal effect for the copula generator.
It does not replace SBC.
The full paper grid and exponential margins require a separate run.

## Local execution record

Date: 9 September 2026.
R 4.5.2, CmdStan 2.38.0, lavaan 0.7.2, covsim 1.1.0, and rvinecopulib 0.7.3.1.0 are available.
The three paper packages use `/tmp/dcvar-causal-validation-lib`.
The installed global libraries do not change.
The EffectLiteR 0.5-1 source supplies the measurement function.
Its full installation needs an unavailable CMake dependency.

The independent scale check passes.
The largest numerical difference is `1.8e-15`.
The EffectLiteR restriction check passes.
The independent paper covariance check passes with a largest difference of `2.3e-16`.
These checks take less than one second each.
One hundred prior datasets per model at `n = 400` produce constant items in
6% of latent covariate datasets and 3% of latent mediator datasets.
The other two models produce no constant items in this run.
These rates are estimates from this finite run.

The recovery pilot completes 16 fits.
It uses two replications of each null and nonzero case at `n = 300`.
It uses four chains with 500 warmup steps and 500 retained steps per chain.
Fit times range from 29 to 66 seconds.
This pilot detects no divergence.
Some parameter R-hat values exceed the reference limit.
The coverage Wilson intervals show the large uncertainty from two replications.
This early run retains draws, diagnostics, and settings in its full result records.
It precedes the complete fit snapshots and source hash records.

The longer reference run uses one nonzero dataset per model at `n = 300`.
It uses four chains with 1000 warmup steps and 2000 retained steps per chain.
All four fits meet the stated diagnostic limits.
The recovery, reference, SBC, and sensitivity pilots use the original factor representation.
The package later adopts an exact conditional representation for three models.
The separate geometry checks below assess that final representation.

| Model | Maximum R-hat | Minimum E-BFMI | Fit time in seconds |
| --- | ---: | ---: | ---: |
| Latent covariate | 1.00437 | 0.594 | 65 |
| Latent outcome | 1.00433 | 0.450 | 53 |
| Latent mediator | 1.00315 | 0.567 | 87 |
| Latent mediator baseline | 1.00264 | 0.547 | 87 |

These runs check model execution and diagnostics at specified parameter values.
They do not establish performance across the full parameter space.
All 72 category frequencies and all 24 item correlations from the four reference
datasets lie inside their respective 95% predictive intervals.
The reports are in `results/reference` and `results/recovery-pilot`.

The separate backend comparison uses one latent outcome dataset with 200 persons.
All 23 nonconstant parameter and effect differences are within three combined
Monte Carlo standard errors.
The largest absolute ratio is 1.727.
The RStan fit has maximum R-hat 1.00243 and no divergences or tree depth events.
Its category and response pattern checks place all observed frequencies inside
the corresponding 95% predictive intervals.
See `results/backend-comparison`.
Separate CmdStan checks for the latent outcome and latent mediator models place
all 36 category frequencies and all 12 item correlations inside those intervals.
See `results/predictive-reference`.

The final package checks report 1323 passed fast assertions and 217 skipped test blocks.
They report no failed assertion and no test warning.
The final real workflows with both Stan backends pass another 62 assertions.
These workflow tests have no failed assertion, warning, or skip.
`R CMD check` completes with no errors, warnings, or notes.

The VITA generation pilot completes 11 of 12 cases at `Nmax = 10000`.
The latent outcome case with Joe copulas returns an invalid calibration object.
That failure stays in `results/paper-generation-pilot`.
A separate run with the same seed, 90427, and `Nmax = 100000` succeeds.
See `results/paper-generation-joe-refinement`.
Generation takes 1 to 49 seconds per initial case.
The refined Joe case takes 10 seconds.
All 11 generated pilot datasets complete their fits.
Each fit uses four chains with 1000 warmup steps and 1000 retained steps.
All 11 fits meet the diagnostic limits.
They have no divergence or maximum tree depth event.
Maximum R-hat is 1.00982 across these fits.
Fit times range from 159 to 287 seconds under concurrent execution.
All 198 category frequencies and all 66 item correlations lie inside their
respective 95% predictive intervals.
Three of 28 Gaussian reference effect targets lie outside their effect intervals.
They are the latent outcome interaction and the latent mediator ADE and AIE
in the Clayton cases.
One dataset per condition cannot estimate interval coverage.
See `results/paper-fits`.

Ten pilot fits use the original factor representation.
The latent mediator Joe fit uses the final conditional representation.
The archived source from each fit records this distinction.
The separate latent outcome Joe fit at `Nmax = 100000` also meets the diagnostic limits.
It takes 101 seconds and has maximum R-hat 1.00929.
Its 18 category frequencies and six item correlations lie inside their intervals.
See `results/paper-joe-refined-fit`.
This separate result does not replace the original generation failure.

The SBC pilot completes two prior datasets per model at `n = 200`.
It uses four chains with 500 warmup steps and 500 retained steps per chain.
All eight fits violate at least one diagnostic limit.
The second latent mediator fit has 43 divergences.
The second latent mediator baseline fit has 14 divergences.
Some fits also have low E-BFMI and high R-hat.
These results detect difficult prior regions.
They do not establish calibration.
All ranks and failed diagnostic checks remain in `results/sbc-pilot`.
Separate geometry checks use the same difficult datasets.
The longer standard fit for the first latent mediator baseline dataset uses
four chains with 2000 warmup steps and 2000 retained steps per chain.
It still has 63 divergences, maximum R-hat 1.024, and E-BFMI from 0.111 to 0.206.
Its effect diagnostics alone look better.
They do not resolve the structural parameter problem.

## Final geometry checks

The latent covariate, latent mediator, and latent mediator baseline models now
use an exact conditional normal representation.
The representation conditions the factor on the continuous observations.
It uses the same joint density, priors, and causal effects.
The latent outcome model keeps its original representation.
The density checks include the Jacobian.
Their largest absolute log density difference is below `7e-12`.
Separate `prior_only` checks change the outcomes and item values.
All prior draws remain identical at the same seed.

The final stress fits use the same difficult prior datasets as the SBC pilot.
They use four chains with 1000 warmup steps and 2000 retained steps.
They have no divergence or maximum tree depth event.
They also report complete diagnostics.

| Model | Data seed | Maximum R-hat | Minimum bulk ESS | Minimum tail ESS | Minimum E-BFMI |
| --- | ---: | ---: | ---: | ---: | ---: |
| Latent covariate | 90423 | 1.00290 | 1198 | 1527 | 0.503 |
| Latent mediator | 90427 | 1.00436 | 1080 | 1451 | 0.696 |
| Latent mediator baseline | 90428 | 1.00328 | 1629 | 1449 | 0.722 |

Longer runs with the original representation do not resolve all problems.
The latent covariate long run has 42 divergences.
The latent mediator baseline long run has 63 divergences.
A centered mediator trial also retains divergences and poor mixing.
These failed comparisons stay in the archive.

Separate normal cases compare the original and final representations.
All 22 latent covariate, 42 latent mediator, and 41 latent mediator baseline
parameter and effect differences are within three combined Monte Carlo errors.
The largest absolute ratios are 2.238, 2.478, and 2.528, respectively.
These checks support the exact transformation and its use in these cases.
They do not establish calibration across the full prior.
The original SBC pilot remains unchanged.
A large SBC study and the full paper grid remain separate research runs.

`results/geometry` contains diagnostics, effects, source snapshots, and input data.
It also contains the original experiment scripts and their reports.
Its manifest records file hashes and the locations of the complete fits.
The original scripts retain their original paths and source assumptions.
Use the stored Stan source and input snapshots when you repeat an old experiment.
`archive-geometry.R NEW_ARCHIVE_DIRECTORY` rebuilds this archive from the complete fits.

## Prior sensitivity result

The prior sensitivity pilot uses one latent mediator dataset with 300 persons.
It fits slope prior SDs of 0.5, 1, and 2.
Each fit uses four chains with 1000 warmup steps and 2000 retained steps per chain.
All three fits meet the diagnostic limits.
Maximum R-hat is at most 1.00442.
ATE means range from 1.00385 to 1.00613.
ADE means range from 0.73277 to 0.73831.
AIE means range from 0.26554 to 0.27335.
This result applies to that dataset and those priors.
See `results/prior-sensitivity`.
