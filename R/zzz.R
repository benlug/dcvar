.onAttach <- function(libname, pkgname) {
  if (interactive()) {
    packageStartupMessage(
      "dcvar: Gaussian-copula VAR models for bivariate time series. See ?dcvar for the core workflow."
    )
  }
}
