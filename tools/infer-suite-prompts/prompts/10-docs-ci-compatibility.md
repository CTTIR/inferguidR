# P10 Legacy reports, documentation and CI

Role: coordinator with documentation/CI owner. Prerequisites: G07 and G09; consume completed classifier/report contracts.

## Objective

Make the three packages installable, documented and coherent as a suite. Finish legacy reporting compatibility and provide executable validation workflows.

## Legacy reports and schemas

Read the current envelope, report_snapshot, report_document, JSON dispatch, resource and export validators before changing them. Preserve historical readers, schema versions, numerical bounds, row-level status, raw-data omission and existing safe export behavior. Retain the native guid exporter during transition if moving it immediately would create a cycle or lose compatibility.

New report components use the shared calc view. If guid export functions forward to report, that dependency is explicitly optional and no report-to-guid edge exists. Test both installed and absent report cases. Existing exported APIs must have the tested treatment recorded in compatibility.csv; a clean package structure does not justify breaking them.

Verify old JSON read/round-trip identity separately from new calculation provenance. Test source-data mutation detection where promised. Preserve the distinction between an analysis hash, a snapshot identity and presentation settings. Do not execute R code or markup while importing a data-only saved analysis.

## Documentation

Update each package README, DESCRIPTION, NAMESPACE, man pages, NEWS, CITATION as appropriate, installation guidance and examples. Preserve attribution and valid licensing. Documentation must use implemented names and actual supported arguments. No functions described only by a conceptual stub.

Provide three concise runnable workflows: ordinary native-R analysis to table/text/ggplot; matched-method planning and observed analysis using Boschloo; and a classifier example preserving evaluation context. Include migration examples for existing guid users and fresh-session examples using each package independently.

Document dependency tiers, optional backend errors, base-R tests, exact method variants, missingness/pairing policy, adjustment families and ggplot coordinate/facet scope. Recheck upstream links/statuses used in optional integrations. Do not promise CRAN availability for development/private or archived packages.

Ship documentation needed to build/check without a runtime roxygen2/devtools dependency. Use existing development tooling only when it adds value and is declared. Keep implementation internals out of the user workflow unless needed to choose a method or understand a result.

## Build and CI lanes

Create base-R scripts that build provider tarballs, install candidates into isolated libraries and check consumers against those exact candidates. Exercise installed public APIs; sourcing R/*.R is not a release check. Avoid hardcoded paths to the author's library.

Required lanes:

1. **Minimal calculation/report core:** only mandatory declared dependencies, no tidyverse/testing/table/plot packages accidentally visible. Verify plain output and optional-feature error paths.
2. **Legacy inferguidR:** declared current prerequisites, historical fixtures, existing regression coverage and new base-R checks. Baseline defects remain separately listed.
3. **ggplot integration:** supported ggplot2 installed; P09 actually executes.
4. **DescrTab2 import/render and claimed specialist/table backends:** explicit optional package versions/sources and real adapters exercised. The requested DescrTab2 import/render check is required even though runtime installation is optional; use a deliberate pinned source if still archived. Every preserved or newly advertised production backend is promoted to a required verification lane. Optional independent reference comparators can remain separate if not a required existing gate. Unavailable requested support is BLOCKED, never silently green.
5. **Platform/minimum-R:** preserve the existing minimum R unless a documented requirement compels a change. Run the available Linux/Windows/macOS checks and minimum-R profile. Report infrastructure limitations precisely.

A named minimal lane may intentionally use _R_CHECK_FORCE_SUGGESTS_=false; document the environment and skipped checks. It never substitutes for the full supported-feature lanes. No continue-on-error for a required check. Do not require every historical optional comparator to be installed just to test base arithmetic, but preserve its separate evidence lane where claimed.

Control and record R versions, .libPaths(), package versions, locale, seeds and material check environment variables. Do not assume --vanilla alone hides all site libraries. Fail core installation/load errors. Logs preserve exit status and warnings.

## Cross repository authentication

Use exact provider candidate commits or verified source artifacts. A workflow's GITHUB_TOKEN does not automatically read another private repository. Use an existing authorized App/token or approved artifact handoff; never hardcode secrets, create broad credentials or make a repository public to fix authentication. If access is unavailable, complete local suite validation and mark the hosted lane BLOCKED.

Make CI behavior real: do not claim success merely because YAML exists. Hosted execution may legitimately remain NOT_RUN before P12 creates/pushes new repositories. G11 evaluates required local candidates separately; after its PASS, P12 may create/push the remotes to enable hosted checks. A hosted-only access blocker does not create a circular prerequisite. Pin external actions/source inputs according to repository policy; avoid fetching an unreviewed script and executing it as a shortcut.

## Gate G10

Examples and help match actual APIs, old report handling is tested, dependency/build instructions are truthful, and executable CI/local scripts cover the required lanes. Candidate source tarballs contain all necessary resources. CI infrastructure still blocked is explicitly carried to G11/G12; it is not a package pass.

