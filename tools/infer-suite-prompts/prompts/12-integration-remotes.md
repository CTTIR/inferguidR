# P12 Final integration and new repositories

Role: coordinator only. Prerequisite: G11 LOCAL PASS for candidate creation/push actions. Hosted status may be NOT_RUN or BLOCKED solely because new repositories/access are not yet available; create/push the verified candidates to enable those checks. A failed or blocked required LOCAL gate prevents candidate pushes. Local preparation and read-only access diagnosis can always continue.

## Objective

Create the requested new repositories, push exact verified content, and leave inferguidR changes reviewable. Preserve all unrelated remote state.

## Before remote mutation

1. Re-read stack.dcf and the user's current instructions. Verify exact owner/repository names, desired visibility and scope. Default new repositories are private; existing visibility never changes under this task.
2. Inspect authentication/permissions without exposing secrets. Check each exact target, distinguishing absent, inaccessible and an existing project. Reuse only a clearly matching authorized project. An unrelated collision stops mutation of that target.
3. Verify the clean commits/source tree IDs already frozen and checked in P11. Do not create a new unvalidated implementation commit here. The provider and both consumers must agree on API versions. Record exact candidate commits and tarball digests in the ignored run directory from P02; avoid a self-referential commit whose own SHA is embedded in itself. Any subsequent code/packaged-content change requires the relevant local gates to run again.
4. Re-read source and target remote HEADs before pushing. If the source default branch moved, integrate changes on the feature branch without rewriting unrelated history, resolve conflicts and rerun affected gates. Never overwrite remote work or force-push to hide divergence.

## Create and push

Use the available official GitHub connector or authenticated CLI. With gh, prefer structured explicit commands such as `gh repo create CTTIR/infercalcR --private --source <verified-local-path> --remote origin`; adapt only after verifying the exact installed CLI and existing remote state. Avoid duplicate origin creation or --push until the intended local branch and candidate SHA are confirmed.

Create infercalcR and infereportR once, with explicit visibility. Initial main pushes to newly created repositories are within scope; existing inferguidR changes go to the configured feature branch only. For a true resumed new repository, verify the recorded run/branch/SHA before continuing. For another matching nonempty existing calc/report project, use a feature branch and reviewable PR; no unapproved default-branch push. Push calc before consumers. Keep consumer dependency configuration reproducible through exact development commits/source artifacts where needed. Do not pretend Remotes metadata alone supplies a package to plain R CMD INSTALL or CRAN.

Use ordinary non-force pushes. If any operation times out or returns ambiguous completion, inspect server state before retrying. Verify repository owner, visibility, URLs, default branches and remote SHAs after success. Do not claim a repository exists from a local git init alone.

Open the inferguidR pull request against its verified default branch. The title and body should explain the existing coupling, the three new ownership boundaries, retained compatibility, base-R dependency policy, concrete validation evidence, cross-repository references and actual remaining blockers. Use a structured body or a temporary body file with literal newlines. Do not put credentials, raw patient data or ephemeral signed URLs in the PR. Do not send email/Slack or unrelated messages.

## Hosted checks and access blockers

Observe actual CI runs where available. Private sibling repositories may require an already-authorized cross-repository credential or source-artifact route. Do not change visibility, invent a token, broaden permissions or mark YAML as passing. If access is blocked, preserve the exact local/remote candidates already completed and identify the one missing permission/configuration with its impact. Respect actual automatic approval review rejections and report their stated reason if they leave work blocked.

If creating/pushing a target is not permitted, finish the local implementation, commits/patches, source tarballs, instructions and evidence. Report remote creation/push as blocked and give the concrete remaining action; do not ask the user to approve work that has not yet been prepared.

## Final handoff

Provide a concise table with each repository URL, branch, package version, verified SHA, pushed SHA and check status. Include the source PR URL, architecture/dependency summary, retained non-base exceptions, base-R test coverage, actual ggplot/backend checks, legacy API/document compatibility and scoped scientific limitations.

Link durable run evidence and local source artifacts where the environment permits. List unresolved blockers distinctly from known preexisting limitations. Explain installation order and authenticated access requirements for private repositories. No source PR merge, public release/tag, CRAN submission, original-repository retirement or visibility change is performed under this stack.

## Gate G12

Repository creation, pushes and PR status are confirmed from actual provider responses. The remote commits match verified candidates. If a remote-only action remains blocked, the final state names it explicitly and supplies complete local deliverables; it does not claim end-to-end completion.

