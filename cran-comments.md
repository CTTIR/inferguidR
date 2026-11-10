## Submission

This is a new submission: inferguidR 0.1.0.

It imports infercalcR, which is submitted to CRAN first; this package can
only be submitted once infercalcR is available on CRAN.

## Test environments

* Local: macOS 27 on Apple silicon, R 4.6.1 (aarch64-apple-darwin),
  `R CMD check --as-cran --no-manual` on the built tarball with `NOT_CRAN`
  unset. The PDF manual was built locally with TinyTeX (TeX Live 2026) via
  `R CMD Rd2pdf` without errors; the check run itself used `--no-manual`.
* GitHub Actions (r-lib/actions v2, `rcmdcheck` with `--as-cran` on one
  shared source archive): macos-latest (R release), windows-latest
  (R release), ubuntu-latest (R devel, R release, R oldrel-1), plus an
  R 4.1.0 core-compatibility job.

## R CMD check results

0 errors | 0 warnings | 1 note

* New submission.

Until infercalcR is on CRAN, the incoming feasibility check also reports
it as a strong dependency not in a mainstream repository.

The installed size is about 7.8 Mb (R 3.0 Mb, help 1.7 Mb, advisor 1.3 Mb)
and may be reported on some platforms. The package ships a registry of
power-analysis procedures, versioned JSON Schemas and the method-advisor
rule set as structured data in `inst/`; these are needed at run time.

## Notes

* Optional packages in Suggests ('shiny', 'ggplot2', 'jsonvalidate',
  'Exact', 'pROC', 'yardstick' and others) are used conditionally, and the
  tests that need them skip when they are absent.
* Examples that build or export a full report take a few seconds each and
  are wrapped in `\donttest{}`; the interactive workbench example runs only
  in interactive sessions.
* Tests use at most two worker processes. The slowest validation test
  files skip on CRAN, so the tests finish in a few minutes; they run in
  full in continuous integration.
