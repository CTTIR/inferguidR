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

### Follow-on execution — 26 October 2026

The reviewed local milestone was merged to `main` and published as `0595042`.
The follow-on work addresses every work package below. Its final verification
record is completed after the source archive is frozen; a row describing an
implemented check does not imply that all scientific qualification gates pass.

| Work package | Implemented work and evidence entry points | Remaining acceptance boundary |
|---|---|---|
| N1 | Enumerated every deterministic rejection state for an attained-size compromise example; documented why the requested equality can be unattainable in `tools/validation/NUMERICAL-CONTRACTS.md` | Exact attained-size compromise stays excluded with `reference_pending`; no nearest-ratio or randomized test is substituted |
| N2 | Corrected split rational ties in pooled, corrected pooled and Mantel-Haenszel exact statistics; independent enumeration, integer-boundary and inverse checks, including the 144-case inverse campaign | Historical native-reference disagreements remain visible; representative tests do not certify every continuous input domain |
| N3 | Added bounded, streaming GLM, rank and regression checks; replayed all scheduled historical G5/G6/G8/G12 mathematical campaigns against the unchanged numerical source, retaining per-case receipts | G12 reuses retained high-precision references; two historical G8 refusals remain explicit. These results do not certify the entire native capture corpus or erase failed/OOM runs |
| N4 | Added design-derived multivariate parameter mappings, independent distribution integrals and utility checks | Approximation choices and finite candidate option grids are distinguished from exhaustive continuous-domain qualification |
| N5 | Independently checked tetrachoric quadrature/variance and bootstrap weighted arithmetic; enumerated finite-sample percentile-interval coverage diagnostics | Bootstrap remains development; undercoverage and finite-sample limitations are retained |
| N6 | Published all registered method/mode/option declarations and scoped sample coverage | Uncovered branches and unresolved applicability remain explicit; no global certification label is introduced |
| P1 | Reviewed variant identities, actual executable templates, estimands, assumptions, overrides and source relationships; corrected stale interval and data-handling promises | Technical consistency and source provenance do not establish external scientific approval or source equivalence |
| P2 | Added explicit `advisor_power_spec()` mappings for the four Boschloo variants, including direction, allocation and supported numerical controls | Other method bridges remain unavailable unless explicitly mapped; older saved sessions require replay |
| P3 | Added observed calibration/net-benefit panels and settings controls, independent arithmetic and adversarial aggregate-report validation | No prevalence-scenario fabrication; frozen probability summaries are opt-in and bootstrap coverage remains qualified |
| P4 | Checked setting-aware staleness, history, individual sweep results, reproducible exports and every advertised helper layout | Failed or changed settings must not export an old result as current; UI editors cover 20 layouts for 11 methods, while other API drawers remain explicitly unavailable in the UI |
| P5 | Added an original bilingual glossary, live resource/anchor audit and frozen aggregate probability and lattice figures alongside existing XY sweeps; replaced destructive export copying with a checked sibling-file rename | Full scored-row ROC/PR/density tables remain intentionally excluded from portable reports; filesystem rename behavior does not establish crash durability or concurrent-writer isolation |
| P6 | Installed-browser review covers offline use, desktop/mobile widths, keyboard operation, separate sessions, downloads and recovery | Acceptance requires actual browser receipts and visual inspection for the final source candidate |
| Release | Added source archive inspection, exported-help/startup/report preflight, isolated optional-dependency absence, licence review and runtime measurement | Final local checks, hosted platform results and deployed HTTP readback are recorded separately; no stable version or submission is inferred |

The portable-report extension is versioned as `inferguidR.report_document/0.4.0`.
Existing report formats remain readable. Probability summaries require explicit
inclusion; lattice-only reports do not imply permission to include scored-data
aggregates. The numerical repair can change criterion results where an earlier
floating-point calculation split a mathematically tied rejection class.

Release policy: retain the experimental version while documented exclusions,
external-reference adjudications or required platform acceptance remain open.
The roadmap's stable-release step is not satisfied by passing local regression
tests or by merely enumerating the registered methods.

The three numerical campaigns now contain 728 cases and 2,655 comparisons,
with numerical representatives for all 50 methods. Of 243 registered
method/mode pairs, 241 have numerical samples and two have only an intentional
exact-compromise refusal check. No method/mode pair is unsampled, but option
interactions and continuous input domains are not exhaustively certified. See
[`EXTENDED-COVERAGE.md`](https://github.com/CTTIR/inferguidR/blob/main/tools/validation/EXTENDED-COVERAGE.md) for the branch
contracts and the distinction between sampled evidence and certification, and
[`INVERSE-COVERAGE.md`](https://github.com/CTTIR/inferguidR/blob/main/tools/validation/INVERSE-COVERAGE.md)
for independent inverse residuals, admissible sample-size scans and discrete
jump conventions. The separate compact public campaign retains its 21 cases
and 761 comparisons.

The [advisor review tables](https://github.com/CTTIR/inferguidR/blob/main/tools/product/reviews/2026-10-26/README.md) make the
69-variant review, 30-reference bibliographic verification and 19 lesson-gap
dispositions inspectable. External candidates are not silently registered as
exact lessons, and bibliographic verification is separate from full-content
or scientific acceptance. Historical source labels remain historical records.

### Final local candidate — 26 October 2026

The locally checked package source archive has SHA-256
`48242513e561ad503a5253b63c069ca2398bafb7f36c22e47a7dea11e0ee22cd`.
Its package source is revision `fef8e64`; the later independent-reference
quadrature correction changes excluded validation tools only. All 580 source
package files passed the distribution audit.

- Linux R 4.6.1 and the available R-devel snapshot (2026-10-10 r90659)
  each completed `R CMD check --as-cran --timings`, including the PDF manual,
  with `Status: OK`: 46,103 passing assertions, zero failures or warnings and
  the same 11 documented skips. Incoming CRAN checks were disabled.
- R 4.1.0 passed isolated installation, 19 core compatibility checks,
  optional-dependency absence, all three numerical campaigns, the public
  checker, bootstrap arithmetic and installed export-failure checks. This is
  not the full optional-package test suite on minimum R.
- The first minimum-R inverse run exposed a reference-integration failure.
  Log-denominator and sine-inner coordinates fixed it without changing
  package calculations, tolerances, tail bounds or the sample-size scan.
  Independent kernel comparisons and fresh minimum/current-R runs passed;
  the original failed runs remain retained.
- Every advertised effect helper layout passed server checks: 20 layouts,
  374 assertions, explicit apply behavior, stale-input refusal and EN/DE
  rebinding. Export transaction tests cover rendering/rename interruption,
  existing destinations and symbolic links. The installed probe also checks
  actual JSON, TXT and PDF replacement.
- Installed browser checks passed history/method/language transitions,
  independent sessions, keyboard operation, stale-result refusal and actual
  downloads. The complete EN/DE desktop/mobile plot matrix was inspected;
  final-package checks repeated representative panels with source-equivalence
  evidence. The 109-page manual and 8-/36-page report PDFs had no text outside
  their physical page bounds; report fonts were embedded.
- The local site passed its build and audit of 119 HTML pages, their local
  references and anchors, plus browser search and mobile navigation.
  Hosted deployment remains a separate gate.
- Historical mathematical replay passed G5's 406 comparisons, G6's 4,575
  comparisons across 604 cases, and G12's 989 calculator cases. G8 retained
  22,062 passing case rows and two independently reproduced sensitivity
  refusals. G6 used two workers, fresh checkpoints and bounded enumeration;
  its measured maximum resident set was about 257 MiB. G12 reused the
  retained high-precision reference table with exact case-ID matching.

### Remaining acceptance and next work

1. Obtain successful hosted Windows, macOS, release/devel and minimum-R jobs
   for the published revision, and verify the deployed Pages HTTP endpoint.
   Local success and queued jobs do not satisfy these gates.
2. Adjudicate the retained native-reference disagreements and the exact
   attained-size compromise/no-solution contract. Preserve typed refusals
   while these contracts are unresolved.
3. Extend scientific coverage beyond representative mode samples: 825 of
   1,511 declared marginal option values remain unsampled. The finite option
   candidate grid also includes input-dependent or infeasible combinations;
   it is not a list of certified designs. Bootstrap undercoverage and
   tetrachoric confidence-interval limits remain explicit.
4. Complete scientific content/source-equivalence judgments and the 19
   lesson-gap dispositions before registering additional exact lessons.
   Further planner bridges and the 31 API effect drawers without UI editors
   are separate product extensions requiring exact mappings and tests.
5. Keep the development version until the required numerical, product and
   platform gates are accepted. A stable release, tag or CRAN submission is
   not implied by this implementation or by merging it to `main`.
