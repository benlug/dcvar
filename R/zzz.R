.onAttach <- function(libname, pkgname) {
  packageStartupMessage(
    "dcvar: Dynamic Copula VAR Models for Time-Varying Dependence\n",
    "- Use dcvar() to fit the DC-VAR model (random-walk rho)\n",
    "- Use dcvar_hmm() to fit the HMM copula model (regime-switching)\n",
    "- Use dcvar_constant() to fit the constant copula baseline\n",
    "- Use dcvar_multilevel() for panel data with random VAR coefficients\n",
    "- Use dcvar_sem() for latent processes with measurement models\n",
    "- Requires cmdstanr + CmdStan for model fitting"
  )
}
