# P01 Current source discovery and baseline

Role: source auditor. Read-only analysis and isolated baseline builds; no statistical implementation changes.

## Objective

Create the concrete ownership and capability inventory for the checkout actually supplied. The architecture was informed by e91938f77eaa991b69989b121bfd32c124c95db0; use that as a comparison reference, never as an assumed or forced HEAD.

## Actions

1. Read applicable AGENTS.md/repository instructions. Record owner/remotes, branch, HEAD, package version, default branch, tracked modifications and relevant untracked user work. Never print credential-bearing remote URLs. Identify R versions and build tools available.
2. Inspect README, DESCRIPTION, NAMESPACE, validation/roadmap documents, schemas, exported objects, source files, resources, templates and test inventory. Read actual implementations for the critical paths below.
3. Build the original package from a clean detached checkout/worktree at the recorded HEAD into an isolated library. If the user's dirty changes are intentionally part of the task, preserve them separately and document the chosen baseline; never silently discard them.
4. Run existing available baseline checks. Save commands, exit status, warnings, package versions, platform and failures. R or dependency absence is a recorded blocker. Continue inventories even if baseline execution is blocked.
5. Capture small deterministic synthetic fixtures from exported APIs and representative historical JSON documents. Use fresh processes for baseline and later candidate versions. Retain expected values and provenance, not just screenshots or pretty-printed output.

## Critical inspection

- core-method/core-registry/core-engine/core-api and methods-*/num-*/effect-* are largely planning machinery. Verify actual function ownership.
- The current provider discovery scans the guid namespace; identify every lookup or system.file() package reference that a split would break.
- advisor_r_code() exports teaching scripts, often with synthetic execution or plots at the bottom. Inventory their computational functions and external packages; do not source the scripts as production code.
- Identify current observed-test and classifier engines, immutable classes, conditions, status fields, result hashes, envelope/report schema versions and raw-data omission rules.
- Identify test methods versus planning variants. In particular inspect Welch versus pooled t, rank-test approximation policies, exact/asymptotic McNemar, and Boschloo central/minlike/directional variants.
- Verify current reporting exports and the observed-test restrictions; distinguish existing reports from the new general reporting interface.

## Required artifacts

Populate the evidence templates with:

- `ownership.csv`: every exported API and significant internal numerical kernel, final owner, source location, resources, tests and compatibility route.
- `capabilities.csv`: all currently advised/implemented procedures, current execution status, backend, result shape, proposed calc adapter and planning relation. Separate advisor-only, executable-template, implemented and externally checked states. Freeze a mandatory_adapter field: existing callable analyses and every template with an executable analysis implementation are mandatory within their present procedure/settings scope. Merely cited/recommended methods without executable implementation may remain documented unsupported.
- `dependencies.csv`: base/recommended/contributed classification, runtime/test/build roles, proposed disposition and justification.
- `compatibility.csv`: old calls/classes/conditions/schemas/docs and acceptance fixtures.
- `baseline/`: exact source SHA, environment, build logs and baseline outcomes.
- `ARCHITECTURE.md`: small public API proposal, explicit ownership boundaries and first extraction slice.

## Gate G01

The coordinator can account for every current public capability and every analysis template, with no assumption that filenames imply ownership. Baseline execution is either evidenced or explicitly blocked. Existing failures are distinguishable from later regressions. The inventory defines full task scope: do not reduce it later by relabelling previously supported methods unsupported.

Handoff includes actual source changes since the observed reference and any material effect on the plan. No claim of numerical validation is made from source inspection alone.

