## Submission

This is a new submission: inferguidR 0.1.0.

It imports infercalcR, which is submitted to CRAN first; this package can
only be submitted once infercalcR is available on CRAN.

## Test environments

* Local: macOS (aarch64-apple-darwin), R 4.6.1, 14 November 2026:
  full `R CMD check --as-cran` on the built source tarball, including
  PDF manual compilation with TinyTeX (TeX Live 2026) and HTML manual
  validation. This supersedes the earlier local check environment.
* GitHub Actions (r-lib/actions v2, `rcmdcheck` with `--as-cran` on one
  shared source archive): macos-latest (R release), windows-latest
  (R release), ubuntu-latest (R devel, R release, R oldrel-1), plus an
  R 4.1.0 core-compatibility job.

## R CMD check results

0 errors | 1 warning | 0 notes (local check before infercalcR is on CRAN)

* The incoming feasibility WARNING includes the new-submission notice and
  the strong dependency infercalcR, which is not yet in the CRAN or
  Bioconductor software repositories (the check confirms its availability
  from the `Additional_repositories` entry, https://cttir.r-universe.dev).
  Submission remains blocked until infercalcR is available on CRAN.
* The same WARNING lists the reserved Zenodo DOI in the README as an invalid
  URL (404); the DOI resolves once the archived release is published.

The installed size is about 6.8 Mb (R 2.0 Mb, help 1.7 Mb, advisor 1.2 Mb)
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
