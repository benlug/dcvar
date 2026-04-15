# CRAN submission comments

## Resubmission

- This is the first submission of `dcvar`.

## Test environments

- Local Ubuntu 24.04, R 4.5.x

## R CMD check results

- Local `R CMD check --as-cran` completed with no errors.
- A local warning about `qpdf` may appear when the utility is not installed.
- URL-check failures on restricted hosts are expected when outbound network
  access is blocked.

## Notes for reviewers

- The package uses `cmdstanr` from the Stan r-universe and requires an external
  CmdStan installation (`SystemRequirements: CmdStan`).
- Stan-backed fitting examples and tests skip gracefully when `cmdstanr` or
  CmdStan is unavailable. The non-Stan parts of the package continue to check
  normally in those environments.
- The source tarball excludes development-only files such as `data-raw/` and
  the Quarto walkthrough source.
