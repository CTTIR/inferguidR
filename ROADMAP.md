# inferguidR improvement roadmap

Baseline: public revision `d03e87b`, audited 26 October 2026. Version
`0.0.0.9000` is experimental. The baseline Linux package check passed 43,887
assertions with 11 reference/specification checks skipped. This roadmap does
not turn that result into numerical certification or release approval.

## Strategy

First make validation states accurate and the public package reproducible.
Then establish an automated platform matrix and documentation build. Numerical
and product qualification follow explicit contracts, with unresolved outcomes
retained rather than hidden by broader success labels. A release candidate is
one immutable source archive with evidence from every required environment.

## Constraints

- Preserve R >= 4.1.0 and the existing public API and JSON formats.
- Keep `reference_pending` and other typed refusals until their mathematical
  and interface contracts are independently verified.
- Keep G*Power executables, manuals, raw captures and private project records
  outside the distributable repository and source package.
- Do not alter accepted historical evidence. New results identify their source
  archive, runtime, command, scope, skips and exit status.
- Run local computation through the shared workstation resource wrapper,
  with at most two test workers and numerical-library threads set to one.
- A passing local build does not establish hosted CI, site deployment, Windows
  compatibility, external-reference parity or CRAN acceptance.

## Phase 1 — accurate status and public reproducibility

### 1. Correct validation explanations

Files: `R/core-registry.R`, `inst/advisor/rules.json`, relevant registry and
advisor tests, and `man/list_methods.Rd` / `man/validation_status.Rd` if their
documented interpretation changes.

- [x] Explain that G*Power certification is pending without claiming that no
  captured evidence exists.
- [x] Explain that `implemented` means a registered engine and mode solver;
  individual parameter branches can still return a typed refusal.
- [x] Describe Boschloo's independently checked standalone engine separately
  from its unfinished advisor-to-planner bridge. Keep the bridge's conservative
  compatibility state until it has a validated mapping.
- [x] Remove obsolete phase-number language from the two Boschloo planning notes.
- [x] Read back the public explanations in both languages and run existing
  regression coverage for unchanged status codes, rule-version migration,
  and saved recommendation round-trips.

Acceptance: focused registry/advisor tests pass; rendered/API readback confirms
that the false explanations are removed; no claim of G*Power parity or automatic
Boschloo planning is introduced. Changes to rules content carry a deliberate
version decision.

### 2. Provide reproducible public validation

Files: `tools/validation/`, `VALIDATION.md`, `.Rbuildignore`, README validation
links, and any required vignette/website references.

- [x] Add an installed-package validation entry point using synthetic inputs
  and independent reference calculations, with explicit case counts/tolerances.
- [x] Include representative forward and inverse cases, allocation and tail
  variants, integer sample-size minimality and safe refusal behavior.
- [x] Document the distinction between public regression tests, independent
  mathematical checks, captured-reference comparison and full certification.
- [x] Explain all eleven baseline skips and the missing private provenance
  material; do not relabel the skipped checks as passed.
- [x] Make the checker fail on a numerical mismatch and write an auditable
  machine-readable result to a caller-selected output directory.

Acceptance: the checker runs using an installed package without local private
files or G*Power; every claimed case has an identified independent oracle;
negative-control mismatches are detected. Public files contain no private
paths, credentials, binaries or raw application material.

## Phase 2 — repeatable checks and documentation

### 3. Add repository automation

Files: `.github/workflows/R-CMD-check.yaml`, `.github/workflows/pkgdown.yaml`,
`_pkgdown.yml`, and supporting validation instructions.

- [x] Configure package checks for Linux release/devel, Windows release and
  macOS release; include a deliberate minimum-R check with compatible dependencies.
- [x] Run the public independent checker against the installed candidate.
- [x] Bound test/numerical thread concurrency and retain logs on failure.
- [x] Pin external actions; give read-only jobs minimal permissions and restrict
  deployment permissions to the Pages deployment job.
- [x] Build documentation on pull requests; deploy the reviewed main branch
  using the already configured GitHub Pages workflow mechanism.
- [x] Set the canonical pkgdown URL and verify the local website output.

Acceptance: workflow syntax/security lint and local equivalents pass. Hosted
acceptance requires actual successful workflow runs after publication. Site
acceptance requires a successful deployment and HTTP readback; writing a YAML
file alone does not close the observed 404.

### 4. Validate a frozen local candidate

- [x] Build one source archive after code changes and record its SHA-256.
- [x] Run the complete installed Linux check and regression suite.
- [x] Run the public checker on that installed candidate.
- [x] Verify installation and core examples on the locally available R 4.1.0
  toolchain with isolated, compatible libraries.
- [x] Run the same archive on the available R-devel toolchain after Linux,
  recording the actual devel build date and dependency profile.
- [x] Render the website and inspect its home, validation and reference pages.
- [x] Review the final changes independently and resolve material findings.

Acceptance: logs and receipts identify the same archive and exact limitations.
Windows/macOS, hosted workflows, browser interaction and manual generation
remain separate unless actually exercised.

## Phase 3 — numerical qualification

These tasks extend the preserved readiness ledger. Each task produces a new
evidence set; none is satisfied merely by repeating the ordinary test suite.

| Order | Work package | Concrete deliverable | Acceptance gate |
|---|---|---|---|
| N1 | Exact unconditional compromise | Written attained-size versus nominal-alpha contract, independently enumerated small-design cases, integer-boundary tests, solver implementation only after the contract is fixed | Every returned design satisfies the chosen error-ratio contract; discontinuities and no-solution cases are explicit; current guard remains until then |
| N2 | Exact search and inverse modes | Case grid covering effect representations, direction, boundaries, allocation, a priori/sensitivity/criterion and solver warnings | Independent enumeration/bracketing agrees within declared tolerances; sample-size minimality holds; unresolved reference disagreements have adjudications |
| N3 | Rank and GLM families | Streaming/reduced-memory G5/G6/G8/G12 comparisons with per-case receipts | All included cases scored; no selection only among agreeing cases; failed/time-limited cases retained; memory budget respected |
| N4 | Multivariate and utilities | Independent checks for RM-ANOVA/MANOVA, Hotelling, correlations, distribution/effect utilities and option interactions | Branch-specific numerical and domain contracts pass; approximation limitations are exposed |
| N5 | Tetrachoric and bootstrap limits | Separate mathematical specifications and simulation/enumeration studies | Coverage/error tolerances and failure behavior justified before upgrading pending/development states |
| N6 | Certification inventory | Method × mode × option coverage and evidence links | Every supported branch accounted for; no global certified label inferred from method registration |

## Phase 4 — advisor, workbench and reports

| Order | Work package | Concrete deliverable | Acceptance gate |
|---|---|---|---|
| P1 | Advisor semantics | Review of all 69 analysis variants, assumptions, estimands and source links | Route-specific scientific and implementation review, including intentionally unavailable routes |
| P2 | Planner bridges | Explicit same-test mappings, starting with Boschloo variants | Analysis-test variant, tails, allocation, nuisance control and planning calculation agree; overrides do not silently inherit unrelated mappings |
| P3 | Classifier workspaces | Finish calibration/net-benefit panels and uncertainty contracts in an isolated prototype | Independent arithmetic, data-role separation, no scenario fabrication, failed-input/stale-result and serialization tests |
| P4 | Planner and history | Mode/option/helper/sweep/history coverage | Stale results cannot be exported as current; independent session state and reproducible exports |
| P5 | Teaching and reports | Verified bilingual sources/glossary, complete classifier/XY/lattice report coverage and visual references | Real sources and anchors; checked pagination/fonts/labels; adversarial and interrupted export cases |
| P6 | Installed browser review | Offline desktop/mobile widths, keyboard access, two sessions, downloads and error recovery | Recorded browser behavior against the installed candidate; automated server tests alone do not close this gate |

## Phase 5 — release qualification

1. Audit API/help/dependencies/licences, startup effects and runtime, preserving
   meaningful regression coverage while reducing duplicated expensive cases.
2. Freeze one release candidate; run Linux, Windows, macOS, R-devel and minimum
   R profiles, including optional-dependency absence and installed/offline use.
3. Review all remaining numerical and product gates. Unfinished functionality
   must be explicitly excluded or retain its typed refusal and documentation.
4. Only then set a stable version, write release notes and prepare the release
   dossier. Publication, release tagging and CRAN submission require their
   own recorded authorization and acceptance evidence.

## Progress and evidence

Baseline audit: 26 October 2026, `R CMD check --no-manual` OK, 43,887 passing
assertions, 11 documented skips. Phase checkboxes refer to the current
remediation, while numerical and product tables are acceptance requirements.
Update checkboxes only after inspecting the corresponding evidence.


### Executed local milestone — 26 October 2026

- Phase 1: status correction, version-aware advisor notes and public independent
  checks implemented and verified.
- Phase 2 local work: pinned shared-archive workflows and Pages configuration
  prepared; documentation rendered and reviewed; workflow lint passed.
- Frozen archive SHA-256:
  `697447a5e4a98ed80af3cd49ffe4231505c01a5128f91a2b8c01588f57ad8265`.
- Linux R 4.6.1 and the available R-devel snapshot (2026-10-10 r90659):
  `R CMD check --as-cran --no-manual` returned `Status: OK`; each passed
  43,887 assertions with 11 documented skips and no failures or warnings.
  Incoming CRAN checks were disabled.
- R 4.1.0: isolated installation and 19 core compatibility checks passed.
  This is not a full minimum-R optional-dependency regression run.
- Public independent checks: 21 synthetic cases / 761 comparisons passed
  against that archive on all three runtimes. Intentional mismatch, wrong
  installed-package identity, non-OK result and missing-reference probes
  behaved as required.
- Local site: 116 HTML pages, 5,143 local links/assets/anchors checked with
  no missing targets; validation article rendered without horizontal overflow.
- Remaining Phase 2 acceptance: publish the reviewed changes, run hosted
  Windows/macOS/release/devel CI and verify the deployed Pages endpoint.
- Phases 3–5 remain open. This milestone does not certify complete numerical
  parity, complete workbench behavior or release readiness.
