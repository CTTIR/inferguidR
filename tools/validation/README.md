# Public independent checks

This directory provides a small, reproducible mathematical check of an
**installed** inferguidR candidate. Inputs are synthetic. No reference
application, captured output, private project file or development namespace is
used. Passing these checks is not full numerical certification or evidence of
agreement with G*Power.

## Run

Install the candidate with `R CMD INSTALL` into a suitable library, then make
that library visible through `R_LIBS`. This prepends the candidate while keeping
the usual compatible dependency libraries. Run from any directory:

```sh
R_LIBS=/path/to/candidate-library \
  INFERGUIDR_EXPECTED_LIBRARY=/path/to/candidate-library \
  Rscript --vanilla tools/validation/check-public.R /path/to/new-results-directory
```

The first argument is a required output directory. It must be new or empty:
the driver refuses to overwrite earlier evidence. The checked package is the
one found on R's library path; the driver does not call `pkgload::load_all()`.
Use the same R executable that installed the package. Dependencies are the
installed package and its normal imports, plus Python 3 (standard library only).
`PYTHON` can name the Python executable if it is not found as `python3` or
`python` on `PATH`.

Set `INFERGUIDR_EXPECTED_LIBRARY` to require that the loaded package comes from
that candidate library. A different installation then fails before numerical
checks. CI installs into a fresh library, checks the install exit status and
uses this guard to prevent falling back to an older cached installation.

Set `INFERGUIDR_SOURCE_SHA256` to the tested source archive's SHA-256 when
available. The summary labels it as caller-supplied: this value is not proof
that the installed package came from that archive. Retain installation logs
and the archive alongside the results. Runtime/package paths, versions,
reference-source hashes and session information are also recorded.

Locally, follow the workstation's shared CPU policy. CI should use one thread
for OpenBLAS, OpenMP, MKL and BLIS; the checker itself does not spawn R workers.

## Output and failure behavior

- `results.csv`: each comparison, observed and reference values, absolute
  error, fixed tolerances, pass/fail result and oracle description.
- `summary.json`: case/comparison counts, run times, package/runtime/source
  identity, checker hashes, negative-control marker and any fatal error.
- `session-info.txt`: the actual R environment.
- `boschloo-ranks/`: exact rational ranks and gap diagnostics produced from
  the specified synthetic designs; `python-output.txt` retains Python output.

The command exits zero only when all comparisons pass and no fatal error
occurs. Missing dependencies, non-`ok` results on supported cases and ambiguous
reference decisions fail the run. Partial CSV/JSON results are retained when
an error occurs after setup. Installation/setup errors before JSON support is
available can terminate without a summary.

To exercise the same failure path without modifying the package:

```sh
R_LIBS=/path/to/candidate-library \
  Rscript --vanilla tools/validation/check-public.R /path/to/new-negative-results --negative-control
```

This deliberately adds 0.01 to the first observed t-test power. The expected
outcome is exit status **1**, `status: "FAIL"`, `negative_control_injected:
true`, and a failed `t-forward-01` / `power` row. These results must not be
presented as an ordinary validation failure or as a successful candidate run.

## Cases, oracles and tolerances

| Cases | Independent calculation | Acceptance |
|---|---|---|
| 3 two-sample t forward designs, both tail variants, unequal allocation and negative effect | Integral of the normal CDF over chi-square density; critical t value found by root search on that integral | Power: absolute `1e-8` plus relative `1e-6`; critical value: absolute `1e-8`; df exact |
| 2 equal-allocation t sample-size searches | Full reference-integral scan for group sizes 2 through 60 | First total N exact; every earlier design below target; achieved power under the forward tolerance |
| 1 t sensitivity search at fixed unequal group sizes | Reference power evaluated at the returned effect | Target residual at most `1e-8` |
| 3 generic chi-square designs | Convolution of squared shifted normal with central chi-square, including the df=1 special case | Null critical tail, power and beta: absolute `1e-8` plus relative `1e-6` |
| 3 one-sided binomial designs with positive and negative effects | Pascal coefficients and direct polynomial probability sums; null rejection tail enumerated independently | Attained alpha, power and beta: absolute `1e-12` |
| 8 Boschloo design/variant combinations: designs (5,7), (4,6); less, greater, minlike, central | Exact rational conditional orderings; direct polynomial nuisance probabilities; grid plus golden-section maximisation | All observed-table p-values and attained size: absolute `1e-9`; engine upper bounds contain reference within `1e-12`; region membership exact; power at 4 alternatives per combination: absolute `1e-12` |
| 1 exact-unconditional compromise request | Declared unresolved-mode contract | Status remains `reference_pending`, with no returned alpha estimate |

The driver records **21 cases**. It records scalar comparisons separately:
all tables in each Boschloo design are checked, not a selected subset that
agrees. A reference p-value within `1e-8` of the tested alpha makes its region
comparison inconclusive and fails the run. The tolerances are fixed before
comparison; they are not fitted to each observed discrepancy.

The continuous references use different mathematical representations from
the production routines, but share central base-R probability primitives and
QUADPACK integration. They therefore provide algorithmic cross-checks, not
independence from every numerical library. The Boschloo maximiser supplies a
numerical supremum estimate rather than a rigorous reference interval; the
package's reported upper bounds are compared separately. The grid is fixed at
2001 points and every near-top grid local maximum is refined.

## Provenance and boundaries

Copyright (c) 2026 Raban Heller. The reference source is covered by the
repository's [MIT licence](../../LICENSE.md).

- `reference/continuous.R` contains the relevant pure reference algorithms
  from `tests/testthat/helper-reference.R`: `ref_log_integrate`,
  `ref_chisq_const`, `ref_pt`, `ref_qt_upper`, `ref_power_t2` and
  `ref_pchisq`. Registry-mutating test helpers are excluded.
- `reference/boschloo.R` and `reference/exact_classes.py` preserve the
  independent checker and exact-rational rank generator from the existing
  Boschloo validation work. The public driver replaces the earlier exhaustive
  driver; no historical result CSVs, binaries, manual material or captures are
  included.

The scope does not cover every method/mode/option, large-design numerical
limits, complete inverse searches, classifier bootstrap coverage, advisor
semantics or cross-platform acceptance. The eleven package-test skips for
excluded development provenance are separate and are not closed by these
selected checks. See the repository's `VALIDATION.md` and `ROADMAP.md` for the
broader evidence and qualification requirements.
