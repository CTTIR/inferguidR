# P02 Local repositories and execution state

Role: coordinator/integrator. Prerequisite: G01.

## Objective

Prepare isolated working locations and local package repositories while preserving the original checkout. Make the later remote creation concrete and idempotent.

## Actions

1. Resolve the configured repository paths to absolute locations. Check existence, ownership, symlinks, remote identity and dirty state. If a target directory/repository already exists, inspect it. Reuse only if it is clearly the requested project and doing so preserves user work; an unrelated collision is a blocker for that target.
2. Create an inferguidR feature branch in an isolated worktree or clean checkout at the agreed baseline. Use a unique run suffix if needed. Record the original checkout and working checkout separately. Do not change the user's branch behind their back.
3. Initialize separate local infercalcR and infereportR Git repositories or use verified existing ones. For a true resume, verify the recorded run ID, branch and expected SHAs. For a matching nonempty project with newer/unrelated work, use a new feature branch and reviewable PR; do not silently write its default branch. New package content must be independent of the source .git directory. Preserve applicable copyright and licenses. Use real configured commit identity; do not invent one or disable signing globally. Local commits are expected throughout implementation; candidates must be clean and committed before P11.
4. Prepare valid DESCRIPTION, NAMESPACE, R/, man/, tests/, README, NEWS and licensing for each new package. Package names and case must match the configuration. Keep developmental version claims honest. No mandatory tidyverse/testing/document backend dependencies.
5. Create `tools/infer-suite-migration/` in the coordinating guid worktree. Keep reviewed architecture/contract documents and ownership/capability/dependency/compatibility inventories there. Create the mutable run directory `tools/infer-suite-migration/runs/<actual-run-id>/` and put STATUS.dcf, phase-status.csv, test-results.csv, candidate-manifest.csv, logs, handoffs and inventory snapshots there. Ignore that run-directory subtree in Git and exclude it from source tarballs. This is the external-to-candidate run ledger referred to later; its final SHAs can change without changing a checked commit. Add narrowly scoped ignore rules for temporary libraries, credentials and build intermediates. Never record secrets. Publish a concise evidence summary through the final handoff/PR without altering checked package content.
6. Establish integration rules: coordinator owns shared metadata, contract files and remote operations; workers use assigned worktrees and branches. Each has an exact starting SHA and file scope. Never run parallel writers against one Git index.
7. If GitHub access is available, inspect the verified owner and whether the requested new remote names already exist. This phase does not need to create/push remotes. Distinguish inaccessible from absent; do not create duplicate targets after an ambiguous response.

## Repository operation scope

Phase P12 performs the configured remote creation and pushes after local candidates are reviewable. New repositories default to private. Existing repository visibility stays unchanged. There is no authorization here to delete/archive/rename the original repository, merge a PR, create a public release, publish secrets or force-push.

If remote credentials are unavailable, prepare all local repositories and continue. If filesystem permissions prevent sibling creation, use an allowed workspace location and record the actual path. Do not weaken sandbox or access controls.

## Gate G02

There are three explicit local working paths, collision-safe branches, package skeletons, an initialized run ledger and a writer-ownership plan. Original work is preserved. Remote state is verified or marked pending; no remote success is invented.

