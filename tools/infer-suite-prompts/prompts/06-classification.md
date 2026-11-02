# P06 Reusable classifier calculations

Role: classifier computation owner. Prerequisite: G03. Work from the actual classif-* ownership map.

## Objective

Move reusable diagnostic/classification calculations into infercalcR while keeping interactive guidance in guid and publication formatting in report. Preserve the current distinction between theoretical inputs, observed estimates and prevalence scenarios.

## Actions

1. Extract the existing supported metric calculations, interval procedures, operating-point evaluation, ROC/PR summaries and AUC comparisons, calibration/net-benefit computations, cutoff selection, theoretical models and related numerical data preparation as assigned by P01.
2. Keep UI controls, teaching explanations, user decisions and session state in guid. Keep drawing/layout in the owning UI or reporting layer. A function generating numeric plot data is computational; its labels/layout need separate ownership.
3. Replace internal cross-file/namespace assumptions with public calc accessors. Base data.frames/lists are sufficient. Keep existing exported guid calls functional through P07 compatibility wrappers.
4. Preserve row-level status/reason, denominators, metric meaning, interval type/level, event and score orientation, comparator, threshold origin, selection/evaluation split, conditioning and scenario provenance.
5. Preserve data attachment/fingerprint and report privacy rules. Aggregate-only input does not identify a score distribution, AUC or cutoff curve. Do not fabricate those outputs from four counts.
6. Retain explicit held-out/tuning/external/same-data labels. A selected cutoff evaluated on its training/tuning data must not acquire a validation claim during conversion.

## Dependency treatment

Use base-R arithmetic and existing native algorithms. Specialist reference packages may remain optional comparators or optional engines. Do not replace a validated method with an improvised formula to claim zero dependencies. Record necessary packages by feature and separate their role as production backend versus validation reference.

## Verification

Create installed-package base-R tests for synthetic count tables, zero denominators, no positives/negatives, perfect/reversed/random classifiers, tied scores, missing scores, prevalence scenarios, one-class and boundary cases, cutoff tie rules, denominator and row-status retention, and seeded uncertainty procedures actually supported.

For inferential comparisons, verify pairing/covariance assumptions and sample identity are retained. Compare to baseline fixtures in a clean process; retain optional independent reference comparisons where already present. Do not weaken scoped validation claims.

Confirm round trips/accessors do not accidentally retain raw observations in portable report fields. Do not equate serialization success with privacy or statistical validity.

## Gate G06

All assigned existing classifier computations have a canonical core owner, documented public interfaces and evidence of preserved numerical/semantic behavior. UI/report consumers have the numeric contracts they need; no dependency back from calc to guid or report exists.

