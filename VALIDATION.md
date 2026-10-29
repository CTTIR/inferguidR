# Validation and supported scope

inferguidR is development software. A method being registered, a calculation
returning `ok`, an independent numerical comparison and a complete external
reference certification are different statements.

| Evidence | What it establishes | What it does not establish |
|---|---|---|
| `list_methods()` implementation status | An engine and generic solvers are registered | Every parameter branch is available or qualified |
| Method `independently_checked` status | Independent checks exist for the method | Exhaustive coverage of all modes, options or numerical boundaries |
| Public package regression tests | The tested behaviors agree with their regression expectations | Complete G*Power parity; some reference tests intentionally use agreeing cases |
| Public independent checker | The documented synthetic cases agree with separate reference calculations | Certification of untested methods or a proof of all numerical error bounds |
| G*Power `reference_pending` status | Complete external comparison is not certified | Absence of all captured reference evidence |
| Passing platform checks | The identified package candidate passes on the recorded environments | Scientific approval, hosted deployment or CRAN acceptance |

## Run the checks

Use a clean checkout and a separate candidate library, keeping compatible
dependency libraries available. From the checkout root:
Choose a new library directory whose parent already exists; a failed build or
installation must stop the sequence before any checker can load stale code.

```sh
set -eu
candidate_library=/absolute/path/to/library
mkdir "$candidate_library"
R CMD build .
R CMD check --as-cran inferguidR_0.0.0.9000.tar.gz
R CMD INSTALL --library="$candidate_library" inferguidR_0.0.0.9000.tar.gz
export INFERGUIDR_EXPECTED_LIBRARY="$candidate_library"
export R_LIBS="$candidate_library${R_LIBS:+:$R_LIBS}"
Rscript --vanilla tools/validation/check-public.R /absolute/path/to/results
Rscript --vanilla tools/validation/check-numerics.R /absolute/path/to/numerical-results
Rscript --vanilla tools/validation/check-coverage.R /absolute/path/to/coverage-results
Rscript --vanilla tools/validation/check-inverses.R /absolute/path/to/inverse-results
```

The independent checker uses an installed package and public synthetic inputs.
Its [instructions](https://github.com/CTTIR/inferguidR/blob/main/tools/validation/README.md) specify the reference algorithms,
tolerances, dependencies, output and failure behavior. Retain the source archive
and its SHA-256 alongside the result files. The checker is deliberately scoped;
it does not read private captures or certify the full G*Power inventory.

Repository automation runs full checks on Linux release/devel, Windows and
macOS, plus a separate minimum-R core profile. The R 4.1.0 core profile uses
hard dependencies and functional examples; it is not the full optional-package
test suite. Actual workflow runs are required before claiming hosted coverage;
all five profiles passed for `c83e08c` on 2026-10-29.

The expanded numerical campaign streams individual results and a complete
method/mode/option declaration inventory. Declared combinations are not
automatically valid input designs or covered numerical cases. Uncovered,
inactive and pending branches remain explicit in the inventory.
The three numerical campaigns jointly contain 728 cases and 2,655 comparisons:
241 registered method/mode pairs have numerical samples and two exact
unconditional compromise pairs have intentional refusal checks. This covers
all 243 registered pairs at that limited level, while 825 marginal option
values remain unsampled. The separate public smoke campaign contains 21 cases
and 761 comparisons. See the
[inverse coverage contract](https://github.com/CTTIR/inferguidR/blob/main/tools/validation/INVERSE-COVERAGE.md)
for integer-design scans and discrete-boundary conventions.

Installed API/help, startup and optional-dependency absence checks and source
archive audits are described in the [release preflight instructions](https://github.com/CTTIR/inferguidR/blob/main/tools/release/README.md).

## Why reference checks can be skipped

The initial public package check had eleven skips. Their purposes are:

| Group | Count | Files or reason |
|---|---:|---|
| Development protocol specifications | 1 | `test-protocol-export.R` |
| MANOVA/G10 reference fixture files | 2 | `test-methods-g4-manova.R`, `test-methods-g10-fixtures.R` |
| Registry source inventory | 2 | `test-core-registry.R` |
| G7 development corpus | 1 | `test-methods-g7-fixtures.R` |
| Advisor source specifications | 1 | `test-advisor-rules.R` |
| G1/G6/G11 captured fixture provenance | 3 | Corresponding `test-methods-*-fixtures.R` files |
| Validation ledger | 1 | `test-slice-t-two-independent.R` |

These materials are not shipped in the source package. A skip is not a pass.
The embedded numerical regression cases still run; their provenance checks
and complete comparison summaries are separate. G*Power executables, manuals
and raw captures must not be inferred to be redistributable from the package's
MIT licence. The public checker contains separate code and synthetic cases.

## Current explicit limits

- Exact unconditional compromise retains `reference_pending`: the generic
  solver's nominal alpha and the exact test's attained size need a separately
  validated contract.
- Tetrachoric approximation-mode confidence intervals retain
  `reference_pending`.
- Classifier bootstrap validation remains `development`.
- Standalone Boschloo calculations have independent numerical checks. The
  explicit `advisor_power_spec()` bridge preserves central, minlike and
  prespecified directional variants for mapped, completed recommendations.
  It requires supplied event probabilities, design parameters and allocation
  policy; it does not infer these from the answers or certify every numerical
  branch. Unmapped recommendations remain unavailable for automatic planning.
- `unsupported`, `infeasible`, `invalid_input`, `numerical_failure` and
  `reference_pending` outcomes must be handled explicitly, not coerced into
  usable design estimates.

The [roadmap](ROADMAP.md) gives the remaining numerical, advisor, workbench,
report, platform and release gates. No supported-scope expansion should be
made solely by changing a status label.

## Advisor rules revision

Rules `2026-10-26.1` corrected validation explanations. Rules `2026-10-26.2`
add explicit Boschloo planning mappings and the additive `available` planning
state. Questions and estimands are unchanged. Saved sessions retain their
original rules version. Restoration across versions requires explicit
migration/replay; do not rewrite old JSON or silently relabel its provenance.
Existing recommendation documents remain readable; an older reader may refuse
the new availability value.

## Reporting a discrepancy

Include the method and mode, complete synthetic specification, returned typed
status, warnings, package revision/version, platform and `sessionInfo()`. For a
reference comparison, identify the reference algorithm or application version,
display precision, tolerance, input mapping and which values disagree. Retain
failed cases alongside passing ones.
