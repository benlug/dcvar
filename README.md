# dcvar: Dynamic Copula VAR Models for Time-Varying Dependence

[<img src="man/figures/logo.png" align="right" width="15%" height="15%" alt="dcvar logo"/>](https://github.com/benlug/dcvar)

<!-- badges: start -->
[![R-CMD-check](https://github.com/benlug/dcvar/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/benlug/dcvar/actions/workflows/R-CMD-check.yaml)
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
<!-- badges: end -->

`dcvar` is an R package for fitting Bayesian copula VAR(1) models to bivariate
time series with time-varying dependence. It supports continuous random-walk,
regime-switching, and constant-copula specifications, with extensions for
multilevel and SEM settings, all estimated through
[CmdStan](https://mc-stan.org/users/interfaces/cmdstan).

## Installation

`dcvar` depends on [`cmdstanr`](https://mc-stan.org/cmdstanr/) and a working
CmdStan installation.

Install `cmdstanr` and CmdStan:

```r
install.packages(
  "cmdstanr",
  repos = c("https://stan-dev.r-universe.dev", getOption("repos"))
)
cmdstanr::install_cmdstan()
```

Install `dcvar` from GitHub:

```r
install.packages("remotes")
remotes::install_github("benlug/dcvar")
```

For skew-normal margins, install `sn`:

```r
install.packages("sn")
```

## Example

The example below simulates a bivariate time series with decreasing
dependence, fits the baseline DC-VAR model, and compares it to HMM and
constant-copula alternatives.

```r
library(dcvar)

# simulate data with decreasing coupling
sim <- simulate_dcvar(
  T = 150,
  rho_trajectory = rho_decreasing(150, rho_start = 0.7, rho_end = 0.3)
)

# fit the DC-VAR model
fit <- dcvar(sim$Y_df, vars = c("y1", "y2"))

# inspect results
summary(fit)
plot_rho(fit, true_rho = sim$true_params$rho)

# compare models via LOO-CV
fit_hmm <- dcvar_hmm(sim$Y_df, vars = c("y1", "y2"), K = 2)
fit_con <- dcvar_constant(sim$Y_df, vars = c("y1", "y2"))
dcvar_compare(dcvar = fit, hmm = fit_hmm, constant = fit_con)
```

## Supported Models

| Model | Function | Dependence Structure |
| --- | --- | --- |
| **DC-VAR** | `dcvar()` | Continuous random-walk on Fisher-z scale |
| **HMM Copula** | `dcvar_hmm()` | Discrete regime-switching with K states |
| **Constant Copula** | `dcvar_constant()` | Time-invariant baseline |
| **Multilevel** | `dcvar_multilevel()` | Random VAR coefficients for panel data |
| **SEM** | `dcvar_sem()` | Fixed measurement model for latent processes |

All models use Gaussian copulas. The core three time-series models
(`dcvar()`, `dcvar_hmm()`, and `dcvar_constant()`) support four marginal
distributions: **normal**, **exponential**, **skew-normal**, and **gamma**.
The multilevel and SEM variants currently support normal margins only.

`fitted()` and `predict()` are currently implemented for the three core
time-series models. `plot_ppc()` is available for normal and exponential
margins; gamma and skew-normal fits do not yet have replicated residuals on the
observed margin scale.

## Documentation

- Getting started vignette: [vignettes/getting-started.Rmd](vignettes/getting-started.Rmd)
- Model comparison vignette: [vignettes/model-comparison.Rmd](vignettes/model-comparison.Rmd)
- Simulation tools vignette: [vignettes/simulation-tools.Rmd](vignettes/simulation-tools.Rmd)
- Source code and issue tracker: <https://github.com/benlug/dcvar>

## Citation

If you use `dcvar` in your work, cite it with:

```r
citation("dcvar")
```

## Getting Help

- Report bugs or request features at <https://github.com/benlug/dcvar/issues>
- For usage questions, include a minimal reproducible example when possible
