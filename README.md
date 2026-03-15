# dcVar: Dynamic Copula VAR Models for Time-Varying Dependence

<!-- badges: start -->
[![R-CMD-check](https://github.com/benlug/dcvar/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/benlug/dcvar/actions/workflows/R-CMD-check.yaml)
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
<!-- badges: end -->

**dcVar** fits Bayesian copula VAR(1) models for bivariate time series, estimating how dependence between two processes evolves over time. Estimation is via Hamiltonian Monte Carlo through [CmdStan](https://mc-stan.org/users/interfaces/cmdstan).

## Models

| Model | Function | Dependence Structure |
|---|---|---|
| **DC-VAR** | `dcVar()` | Continuous random-walk on Fisher-z scale |
| **HMM Copula** | `dcVar_hmm()` | Discrete regime-switching with K states |
| **Constant Copula** | `dcVar_constant()` | Time-invariant baseline |
| **Multilevel** | `dcVar_multilevel()` | Random VAR coefficients for panel data |
| **SEM** | `dcVar_sem()` | Fixed measurement model for latent processes |

All models use Gaussian copulas. The core three time-series models (`dcVar()`,
`dcVar_hmm()`, and `dcVar_constant()`) support four marginal distributions:
**normal**, **exponential**, **skew-normal**, and **gamma**. The multilevel and
SEM variants currently support normal margins only.

`fitted()` and `predict()` are currently implemented for the three core
time-series models. `plot_ppc()` is available for normal and exponential
margins; gamma and skew-normal fits do not yet have replicated residuals on the
observed margin scale.

## Installation

```r
# install cmdstanr (not on CRAN, use R-universe)
install.packages("cmdstanr", repos = c("https://stan-dev.r-universe.dev", getOption("repos")))

# install CmdStan
cmdstanr::install_cmdstan()

# install remotes for GitHub installation
install.packages("remotes")

# install dcVar from GitHub
remotes::install_github("benlug/dcvar")

# optional: for skew-normal margins
install.packages("sn")
```

## Quick Start

```r
library(dcVar)

# simulate data with decreasing coupling
sim <- simulate_dcVar(
  T = 150,
  rho_trajectory = rho_decreasing(150, rho_start = 0.7, rho_end = 0.3)
)

# fit the DC-VAR model
fit <- dcVar(sim$Y_df, vars = c("y1", "y2"))

# inspect results
summary(fit)
plot_rho(fit, true_rho = sim$true_params$rho)

# compare models via LOO-CV
fit_hmm <- dcVar_hmm(sim$Y_df, vars = c("y1", "y2"), K = 2)
fit_con <- dcVar_constant(sim$Y_df, vars = c("y1", "y2"))
dcVar_compare(dcVar = fit, hmm = fit_hmm, constant = fit_con)
```

## Learning More

- **[Getting Started](vignettes/getting-started.Rmd)** - Installation, basic workflow, and interpretation
- **[Model Comparison](vignettes/model-comparison.Rmd)** - Comparing DC-VAR, HMM, and Constant models via LOO-CV
- **[Simulation Tools](vignettes/simulation-tools.Rmd)** - Trajectory generators, parameter recovery, and metrics

## Citation

If you use dcVar in your research, please cite it:

```r
citation("dcVar")
```

## Getting Help

- For **usage questions** or **bug reports**, open an [issue](https://github.com/benlug/dcvar/issues) with a minimal reproducible example.
