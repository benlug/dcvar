.onAttach <- function(libname, pkgname) {
  packageStartupMessage(
    "dcVar: Dynamic Copula VAR Models for Time-Varying Dependence\n",
    "- Use dcVar() to fit the DC-VAR model (random-walk rho)\n",
    "- Use dcVar_hmm() to fit the HMM copula model (regime-switching)\n",
    "- Use dcVar_constant() to fit the constant copula baseline\n",
    "- Use dcVar_multilevel() for panel data with random VAR coefficients\n",
    "- Use dcVar_sem() for latent processes with measurement models\n",
    "- Requires cmdstanr + CmdStan for model fitting"
  )
}
