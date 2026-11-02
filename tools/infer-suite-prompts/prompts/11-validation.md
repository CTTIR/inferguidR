# P11 Independent validation and integration gate

Role: independent reviewer/validator. Prerequisite: G10. Read `90-independent-review.md` as well.

## Objective

Test the actual integrated candidates, identify regressions and determine readiness. Do not accept another agent's summary as proof. Validation is scoped to the implemented/claimed capabilities; it does not imply universal statistical certification.

## Candidate freeze

Require clean, committed implementation candidates before testing; the coordinator integrates reviewed worker commits first. Record exact repositories, branches, commits, source tree IDs and contract revision. Mutable evidence lives in the ignored run directory defined in P02 and must not alter checked content. Build source tarballs. Install calc first, then guid/report consumers into isolated libraries. Do not load baseline and candidate versions of the same package in one R session. If candidate code or packaged content changes after testing, invalidate the affected gate, rerun its checks and refresh the candidate manifest.

Use the R executable belonging to each selected runtime. Run R CMD build on each package source directory, then R CMD INSTALL and R CMD check on the resulting source tarball. Run an appropriate --as-cran lane when its environment is available; distinguish infrastructure notes from package defects. Inspect check logs, not just the shell exit code. A package can exit with warnings that still need disposition.

## Required review matrix

### Architecture and dependencies

- No calc dependency on guid/report, no report requirement for guid, and no hard cycles.
- No cross-package private API access, sourcing of implementation files, runtime package installation or hidden conditional package dependency.
- Minimal core works in an actually constrained library. Base/recommended/contributed dependencies are accurately classified; every exception has a functional reason.
- Shared schemas/accessors have one owner; metadata and registries are consistent across packages.

### Numerical and method behavior

- Boschloo variants, ties, bounds, decisions and rejection-region/planning consistency.
- Standard-R adapter correctness with explicit matched settings and nontrivial edge cases.
- Pairing, event/reference orientation, contrast sign, missingness and denominators.
- Estimates and intervals target the same estimand; valid unbounded intervals remain distinct from numeric failure.
- Adjustment is performed on the declared family and preserved across filtering.
- Supported power planning refers to the intended analysis; approximations and unsupported mappings are explicit.
- Classifier rows retain status, denominator, scenario and selection/evaluation context.

### Compatibility

- Legacy public APIs/classes/conditions have their declared behavior.
- Historical JSON fixtures retain content and recorded identity without silent hash rewriting.
- New provenance is honest and numerical equivalence is tested separately from historical hash equality.
- Source resources/schema files are included in installed packages and resolved from the correct owner.
- Existing meaningful tests were retained or ported with recorded executed coverage.

### Reporting

- Tables, text and plot labels derive from the same stored values/formatter.
- Rendering performs no inference, refit, interval recomputation or p-value adjustment.
- Inputs remain unchanged; result-only reporting works without raw observations.
- Requested ggplot lane actually ran: transformations, facets, order, inherited aesthetics, horizontal/dodged explicit positions and missing results.
- Optional backend absent/available paths have truthful status; no archived dependency silently becomes mandatory.

## Evidence rules

New software tests are base-R scripts with assertions, expected condition classes and documented numerical tolerances. Preserve existing framework coverage until equivalence was executed. A direct stats parity test demonstrates wrapper behavior; independent arithmetic/property/reference tests support the underlying extraction where feasible.

Report each check PASS, FAIL, SKIP, BLOCKED or NOT_RUN, with commands, environment, exit status and log path. Required test exceptions must fail the lane. Require a positive count of expected executed checks so an empty test selection cannot pass. Compare actual run identifiers with the required inventory; unused tests and hidden skips are review findings.

## Repair loop

Return findings by severity with exact reproducible cases and owning phase. The coordinator assigns repairs to file owners. Preserve acceptance criteria; do not relax assertions or redefine existing behavior simply to close findings. Reinspect the fixes and rerun affected integrated checks. Stop optional expansion when the concrete remaining risks and required gates are resolved.

## Gate G11

G11 is the LOCAL candidate gate. Approve it only if required local package checks, legacy compatibility, ggplot integration, the requested pinned DescrTab2 import/render check and every other claimed/mandatory production capability lane actually passed. Missing mandatory local execution is BLOCKED and cannot satisfy G11. Track hosted checks separately as NOT_RUN/BLOCKED/PASS/FAIL; hosted checks awaiting new repositories do not prevent a local G11 PASS or the P12 create/push actions. A remaining hosted failure/blocker prevents an unqualified full-completion claim. Do not claim a fully green suite, supported optional backend, or release readiness from skipped execution.

