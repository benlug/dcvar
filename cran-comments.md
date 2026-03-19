## Test environments

- Local macOS Tahoe 26.3.1, R 4.5.2
- GitHub Actions:
  - ubuntu-latest (release, devel, oldrel-1)
  - macos-latest (release)
  - windows-latest (release)

## R CMD check results

- 0 errors
- 0 warnings
- 2 notes

## Notes

1. `Suggests or Enhances not in mainstream repositories: cmdstanr`

`cmdstanr` is an optional backend only. The default backend remains `rstan`,
which is in `Imports`. The package installs, checks, and runs examples without
`cmdstanr`. When available, `cmdstanr` enables an additional backend path that
is exercised by dedicated regression tests. It is declared in `Suggests` and
resolved via `Additional_repositories: https://stan-dev.r-universe.dev`.

2. `unable to verify current time`

This note was produced by the local macOS check environment and is not caused
by package code or metadata.
