# Inference suite portability verification, 2 November 2026

The three development packages use version 0.0.0.9000 and common contract 1.0.0. This supplement records the repairs prompted by the first hosted matrix, while retaining the original [31 October validation](VALIDATION-20261031.md) and its failed/superseded evidence.

## Frozen source artifacts

| Package | Tested implementation commit | Source archive SHA-256 |
|---|---|---|
| infercalcR | `d89081f5664c69930e2ee6291ae4b7d9670ad6a8` | `69b132efc871921c805350139b10cf7ecce77f0f893c7fe41da67bc264aae01e` |
| infereportR | `32f3be6f668fd18b9eb345705a6d7fc23cc8a2db` | `33b4da7b22856f5c48c5587d901c39540708f013fa3d3789761f063e31e7acfe` |
| inferguidR | `be242e23499a43477d206efa4f3fc156dfc154d5` | `fd95522c1a6f87ea5174a640e094b1786c5d809777db7f8c12422587a4d90a37` |

Exact final source archives and fresh command logs are under `runs/20261031-001/final-07/`. Compared with final-06, calculation changes only two test files and their test-policy README; DESCRIPTION also carries the generated build timestamp. Both consumer archives differ from final-06 only in DESCRIPTION, which pins the new calculation commit. The parity ledger checks every archive member. Prior evidence remains in its original final-05/final-06 directories; it is not relabelled as freshly executed final-07 work. Later changes to this excluded evidence file do not change package content.

## Repairs and their scientific scope

Older R can materialize shared deferred factor levels or names while a model is summarized or validated. Calculation now isolates the affected model conversion and ordinal/multinomial inputs before those operations. Tests retain strict original-input serialization checks. Extracted numerical modules, original golden fixtures, estimates, intervals and native diagnostics are unchanged.

Wilcoxon results now label the confidence level as native-reported. R 4.1 can retain a nominal or rounded attribute: the enumerated rank-sum example covers 33/35 while reporting .95, and a one-sided signed-rank example covers 7/8 while reporting .87. R 4.6 uses different conservative endpoint/attribute policies for these examples. The legacy field and every native number are preserved; reporting exposes the qualification instead of substituting a newly calculated interval or coverage value.

Frozen backend comparisons now distinguish representation from numerical regression. Whole factor columns and their correlations are aligned together; Ward clades, heights, cophenetic distances and the requested partition are checked despite reordered tied merges. EFA and multinomial optimizer comparisons use narrow field-specific bounds with independent covariance, likelihood, score and information checks. GLMM checks additionally integrate the random-intercept likelihood independently and retain strict current-runtime refit and Wald-algebra comparisons. Deliberately corrupted results must still fail. These fixture-specific checks do not establish universal statistical accuracy or nominal coverage.

The final ordinal check addresses a platform difference of approximately 1e-14 in a small Wald p-value. It requires the current value to equal the native normal-tail formula, retains the original coefficient/standard-error/z tolerances, and bounds the frozen p-value drift by the intervening normal density and z drift plus a small floating-point allowance. It rejects even a one-representable-step modification of the current p-value. Actual captures from all four hosted lanes and deliberately corrupted outputs were checked independently. This changes test policy only; production results and golden fixtures are unchanged.

## Fresh final-07 evidence

The full local calculation package check (`R CMD check --no-manual --no-vignettes`) is Status OK. Its checked source is byte-identical to the final source archive except the generated DESCRIPTION Packaged timestamp; a fresh exact-archive installation and independent installed comparator check follow that source-parity proof. This final check is not a new `--as-cran` run. All 71 mandatory capabilities and 7,177 legacy assertions execute. A clean three-package installation verifies exact library paths and provider pins. The 15 real rank-result reporting cases, shared planning contract, three documented workflows and all 26 installed release checks pass; 99 exports have documentation. The ordinal policy passes the actual four-lane captures and negative corruption tests.

## Retained final-06 evidence

- Calculation R 4.6.1 full package check: Status OK; 7,177 legacy assertions, all 71 mandatory capabilities, zero failures/warnings/skips in those lanes.
- Reporting full package check: Status OK; 15 calculated/imported rank results preserve confidence provenance and numbers under a no-refit guard; 12 ggplot behavior groups, six rendered examples and pinned DescrTab2 import/table/text/HTML checks pass.
- Guidance affected installed regression: 131 groups and 10,447 assertions pass without failures, warnings or skips, including advisor templates/scenarios, rank-version policies, shared planning contract, core JSON/results and numeric JSON precision.
- Independent exact-archive probes on R 4.1.3 and R 4.6.2 pass model/factor-input byte preservation, native parity, all nine Wilcoxon design/direction cases and rendering under a no-refit guard.
- Guidance documentation, examples, vignette rebuild and package check with `--no-tests`: Status OK. Historical JSON identities, seven planning baselines, three suite workflows and all 26 installed release checks pass; 99 exports have documentation.
- Independent constrained-library calculation/reporting installs and package checks: Status OK, with optional contributed backends confirmed absent.
- Calculation Linux R 4.6.2 full package check: Status OK; all 7,177 legacy assertions and 71 mandatory capabilities execute with zero failures/warnings/skips in those lanes, including current specialist backends.

Retained final-06 integration package checks explicitly disable CRAN incoming feasibility and force-Suggests checks, while requiring the full supported integration lane. This is development validation, not a CRAN submission result. The previous complete guidance run (38,838 assertions; 17 classified optional/development skips) remains evidence for unchanged source. The retained final-06 guidance package check deliberately omits the entire test suite and is supplemented by the affected installed checks; it is not reported as another complete 38,838-assertion run.

## Remote and credential status

The private `CTTIR/infercalcR` and `CTTIR/infereportR` repositories exist. Guidance remains reviewable in [pull request 2](https://github.com/CTTIR/inferguidR/pull/2); source main is not merged. Earlier calculation matrices failed portability comparisons and led to the repairs described above. The final candidates have been pushed and their remote commits verified. At 16:46 UTC on 2 November, the final [calculation matrix](https://github.com/CTTIR/infercalcR/actions/runs/41483183470) has passed the Windows and Ubuntu current-R full lanes, both with Status OK, and the Ubuntu R 4.1 minimal lane. The minimal lane has the expected note for absent optional jsonlite/lme4/geepack/testthat packages. Its macOS job remains queued with no runner assigned; an overall hosted pass is not claimed. The workflow link provides subsequent status.

The final [reporting run](https://github.com/CTTIR/infereportR/actions/runs/41483191637), [guidance checks](https://github.com/CTTIR/inferguidR/actions/runs/41483211708) and [documentation](https://github.com/CTTIR/inferguidR/actions/runs/41483212067) stop at the missing-provider-credential guard before package execution (the reporting macOS job is also queued). These are hosted authentication blockers, not package-validation passes.

Consumer checks require `INFER_SUITE_TOKEN` with read access to the private calculation repository. Credential creation was explicitly authorized, but browser control remained unavailable. A tested alternative using separate read-only deploy keys was rejected by GitHub with HTTP 422, “Deploy keys are disabled for this repository.” No GitHub deploy key or repository secret was created, no broad local credential was copied, and no organization/repository policy was changed. The unpublished deploy-key workflow prototype is retained only in local evidence; shipped workflows retain the token guard.
