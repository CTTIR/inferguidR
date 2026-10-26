# Continuous checks

`R-CMD-check.yaml` checks Linux release/devel, Windows release and macOS release
with all declared suggested packages required. A source-build job creates one
archive with a SHA-256 manifest and source commit. Every platform downloads and
verifies that same archive before checking it. Before distribution,
`tools/release/audit-source.py` checks the archive's inventory and distribution
hygiene; its receipt is retained even on failure. Each full-check job retains
its logs, installs the checked archive, and runs `tools/validation/check-public.R`,
`tools/validation/check-numerics.R`, `tools/validation/check-coverage.R`,
`tools/product/check-bootstrap.R` and `tools/release/check-package.R` against the
installed package. The method campaign samples every registered numerical
method; it does not exhaust every input or option combination. The bootstrap
campaign checks bounded classifier resamples against independent calculations.
A failure in any checker fails the job. Their result directories remain under
`check/` for the always-run evidence upload. Numerical libraries
use one thread and package tests use at most two workers.

The checkers receive the verified archive SHA-256 and commit through
`INFERGUIDR_SOURCE_SHA256` and `INFERGUIDR_SOURCE_COMMIT`. The public checker's
exact-rational oracle uses Python 3.12 installed
by the pinned setup action; `PYTHON` points to that interpreter explicitly,
including on Windows.

Both platform profiles use `install-candidate.R` to install into a newly created
library outside the dependency cache. It calls `R CMD INSTALL`, records the
exit status and log, and exports the library and source hash only after a
successful installation. A fresh R process verifies that `find.package()`
resolves to this candidate library before running further checks. Installation
failure cannot fall back to a cached copy of inferguidR. The core smoke script
also enforces `INFERGUIDR_EXPECTED_LIBRARY` when that variable is supplied.

The checks use `--as-cran --no-manual`, with remote incoming checks disabled.
They are regression and packaging evidence, not CRAN acceptance or complete
G*Power agreement. Platform jobs share one immutable source archive; their
interpreter and dependency environments are recorded separately in check logs.

The separate minimum-R job uses exactly R 4.1.0. Its dependency installer creates
a fresh library containing only `digest` and `jsonlite`, and hides ordinary site
and user libraries from subsequent R processes. It verifies and installs the
same source archive as the full-check jobs, exercises the core API, then runs
the installed release preflight with `--minimal`. That preflight requires its
selected optional packages to be absent and checks portable report exports and
typed refusals for unavailable display functions. Recommended packages shipped
with R remain available. For a local source installation:

```sh
R CMD INSTALL .
Rscript tools/ci/minimum-r-smoke.R /path/to/new-or-empty/results
```

The script checks registry loading, a standard two-group design, agreement with
`stats::power.t.test()`, invalid probability rejection, classification arithmetic,
advisor startup and JSON hash round-trips. It writes a receipt and session
information. This is **core compatibility**, not the full test suite or a check
of optional UI/report packages on R 4.1.0. Use an isolated library appropriate to
the interpreter; do not reuse compiled packages from a newer R version.

`pkgdown.yaml` builds the website for pull requests and main. The build job has
read-only repository permissions. Deployment requires a successful build, runs
only for main outside pull requests, and uses the `github-pages` environment
with Pages and identity-token permissions. The `website` dependency profile in
DESCRIPTION installs the shared CTTIR theme. The repository must have GitHub
Pages configured to use GitHub Actions.

Workflow files do not establish hosted success by themselves. Record actual
workflow results, deployment URL and HTTP readback before claiming that CI or
the documentation site is operational.

Workstation commands must also follow the local shared CPU and temporary-file
policy; those host-specific settings are not needed on GitHub-hosted runners.
