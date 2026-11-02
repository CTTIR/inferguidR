# Inference suite validation, 31 October 2026

The three development packages have completed the required local package and integration checks on R 4.6.1, arm64 macOS 27.0.1. Package version is 0.0.0.9000 and shared contract revision is 1.0.0. This records bounded development validation, not CRAN release readiness or universal statistical certification.

## Tested source artifacts

| Package | Implementation commit | Source archive SHA-256 |
|---|---|---|
| infercalcR | `70a3ecbd626df69da4b9b9ec101da2eb68b94d12` | `da124ad9bdce13156d0ea701568eed6a8fbf45b0571c8653ad3c494772df0635` |
| infereportR | `a6ddf365b99f3c445d073de91f69bdf58b056794` | `9ad0279021c4c6f6dbc87570898e9984434c14039d37653493019e1440fbf0e7` |
| inferguidR | `afe4323d72df5f0a456b1fec52be17b8f0416939` | `0f0d2aca513c3f59bd911f15e22a24bb5cd4f78f7c4fbf6d2eeb6458a891fe74` |

The later shared-contract copy and this evidence summary are excluded from the source package by `.Rbuildignore`; they do not change the tested runtime, tests, schemas or packaged documentation. Mutable command records, logs, archived failed attempts and the publication manifest are local under `runs/20261031-001/`. `candidate-manifest.csv` binds the tested commits, trees and archive digests; remote publication is recorded separately.

## Executed checks

| Lane | Observed result |
|---|---|
| Calculation full package check | Exit 0; 7,177 retained assertions, zero failures/warnings/skips; all 71 mandatory capabilities execute, with 68 primary mappings and three helpers. One CRAN incoming feasibility NOTE only. |
| Reporting full package check | Exit 0; actual pinned DescrTab2 import/render and 12 ggplot behavior groups pass; six plot examples rendered. One CRAN incoming feasibility WARNING only. |
| Guidance local package check | Status OK, exit 0; 38,838 retained assertions, zero failures/warnings, 17 explicitly classified skips below. Examples, vignettes, installed resources and new suite contract pass. |
| Independently isolated minimal core | Fresh calculation/reporting source installs, absence-of-optional-backend probes and both package checks pass with contributed packages hidden. |
| Independent integration | 11 public behavior groups, six attained-coverage groups, and 20 historical whole-result objects/identities pass. |
| Saved documents and planning | Five historical JSON documents remain byte-identical over two round trips; seven complete planning comparisons pass; valid-hex identity tampering is rejected. |
| Documented workflows and public API | Three suite workflows and 26 installed package checks pass; all 99 guidance exports have help, including all 98 historical exports. |
| Targeted rank regression | Five existing guidance groups, 183 assertions, zero failures/warnings/skips. Requested confidence remains in (0,1); valid exact intervals retain attained coverage 1 and unbounded endpoints. |
| Bounded numerical campaigns | 761 public comparisons, 990 numerical comparisons over 313 cases, 1,067 extended comparisons over 271 cases, 598 inverse comparisons over 144 cases, and 57 bootstrap comparisons pass. |

The bounded numerical campaigns used final-02 candidates. Exact source diffs establish that their legacy numerical kernels, test drivers and resources remain unchanged in final-05; the subsequent calculation runtime repair changes only the common finalized confidence-level validator. Current full package and integration checks were rerun after that repair. Bootstrap coverage is limited to the documented synthetic settings, six metrics, 250 draws and one seed; it is not general nominal-coverage certification.

The independently identified design-provenance, prediction-estimand, omnibus-bracket, adjustment-family reuse and attained-confidence issues were repaired and rechecked. Native method settings, estimates, intervals and historical hashes were preserved. Tests were not relaxed to close these findings.

## Explicit exceptions

The final guidance check has 17 skips:

- Six optional comparison tests in `test-classif-compare.R`: three require pROC and three require yardstick. These are unexecuted comparator checks, not supported runtime backends. Required native classifier behavior is tested independently.
- Eleven development-only provenance checks: protocol source specification (one), registry builder/inventory (two), advisor source specifications (one), original G1/G4/G6/G7/G10/G11 fixture records (six), and the original documented-example transcription ledger (one). The corresponding development inputs are absent from the fetched baseline and original checkout. Embedded numerical fixture assertions and installed runtime/schema checks still execute.

A strict guidance CRAN submission preflight failed on absent optional MESS, pROC, pwr, pwrss, quarto, WebPower and yardstick. That attempt remains failed. The separate required same-archive local check explicitly sets `_R_CHECK_FORCE_SUGGESTS_=false`, `_R_CHECK_CRAN_INCOMING_=false` and `INFER_SUITE_FULL=true`; a positive preflight verifies all 19 required suite/production/test dependencies and exact installation paths. Missing optional reference comparators are not counted as passes.

Calculation/reporting incoming notes and warnings concern development versions, private or not-yet-created repository URLs, Remotes and dependencies unavailable through CRAN/Bioconductor, including the pinned archived DescrTab2 source. All actual supported-backend lanes executed. Hosted platform/minimum-R checks and private-provider authentication must be observed after repository publication; workflow YAML alone is not a passing result.

## Reproduction and installation

Use the exact archives in `runs/20261031-001/final-05/`, installing infercalcR before infereportR and inferguidR. Set both R_LIBS and R_LIBS_USER to the selected library. Plain R CMD INSTALL does not fetch dependencies; guidance and legacy identity-bearing APIs retain digest/jsonlite. The full suite also deliberately installs optional ggplot2, DescrTab2 and specialist backends; plain calculation/reporting works without them. See [INFER-SUITE.md](../../INFER-SUITE.md) for runnable workflows and the base-R local check runner.

Detailed evidence is recorded by the command JSON/log pairs in `final-05/`, the independent review under `independent-review/`, and the dated `test-results.csv`, `phase-status.csv`, `dispositions.json` and final publication handoff. Earlier failed or superseded runs are retained and are not represented as successful gates.
