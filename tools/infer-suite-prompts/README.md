# inferguidR three package migration prompt stack

Version: 2026-10-31.1

This stack instructs a coding agent, running from an existing inferguidR checkout, to perform the agreed split into infercalcR, inferguidR and infereportR. It includes source discovery, implementation, base-R validation, creation of two new GitHub repositories, ordinary pushes and a reviewable inferguidR pull request. The stack itself is the deliverable; its creation does not mean the package migration or remote repository creation has run.

## Start

1. Extract the bundle into the inferguidR repository root. It creates `tools/infer-suite-prompts/` only; it contains no replacement package source files.
2. Inspect `tools/infer-suite-prompts/stack.dcf`. The defaults target `CTTIR/infercalcR` and `CTTIR/infereportR`, create new repositories as **private**, and leave the existing inferguidR default branch and visibility unchanged. Set another visibility only when the user has explicitly chosen it.
3. Give the coding agent this instruction:

```text
Read tools/infer-suite-prompts/START_HERE.md and execute the complete orchestrated migration it defines. Treat the shared rules, phase prompts and stack.dcf as one task. Reinspect the current inferguidR checkout before acting, implement the packages, run the required base-R and integration gates, and create/push the two new repositories using the configured scope. Continue through the dependency graph rather than stopping after a plan. Preserve unrelated work and record any actual access or environment blocker.
```

If only the combined Markdown document is available, the agent may materialize its named file sections into `tools/infer-suite-prompts/`, verify them, then run START_HERE. Configuration is data: never evaluate DCF values as shell or R code.

## Resulting responsibilities

| Package | Canonical owner of |
|---|---|
| infercalcR | Method identity, observed analyses, numerical inference, reusable statistical kernels, diagnostic computations, common result/accessor contract |
| inferguidR | Advisor policy, study assumptions, power models, sample-size/design searches, planning comparisons, guidance/workbench |
| infereportR | Descriptive and inferential reporting, value formatting, manuscript text, ggplot annotations, output adapters |

Mandatory suite dependencies: inferguidR imports infercalcR; infereportR imports infercalcR. The computational package has no dependency on either consumer. Reporting works without the guidance package. The initial target is no mandatory contributed dependency for infercalcR's new computational core, and only infercalcR plus R-distributed base packages for plain reporting. Compatibility-sensitive existing dependencies are handled explicitly, not removed by slogan.

## What base R means here

The user constraint applies both to statistical implementations and to software tests. Prefer the standard `stats` implementations and ordinary lists/data.frames. Write new tests as `.R` scripts using `stopifnot`, `identical`, `isTRUE(all.equal(...))` and standard condition assertions. R-distributed base packages such as stats, utils, tools, graphics, grDevices, grid and methods are allowed. Recommended packages and contributed packages are separately declared. Calling standard R does not imply that the R runtime contains no compiled code.

The requested ggplot feature needs ggplot2 when used; it remains optional for package loading and plain reporting. Existing canonical JSON/SHA-256 behavior must not be replaced by an ad hoc implementation to reduce a dependency count.

## Execution map

| Phase | Work | Depends on |
|---|---|---|
| P01 | Inspect current source, capture baseline, inventory capabilities | none |
| P02 | Prepare isolated repositories, branches and execution state | P01 |
| P03 | Implement and freeze shared contracts and base-R harness | P02 |
| P04 | Extract Boschloo and other existing numerical kernels | P03 |
| P05 | Promote observed analyses through standard-R adapters | P03 |
| P06 | Extract classifier calculations | P03 |
| P07 | Integrate planning, advisor and compatibility wrappers | P04, P05, P06 |
| P08 | Implement report tables, descriptions and text | P03 |
| P09 | Implement and verify ggplot annotation layers | P05, P08 |
| P10 | Complete legacy report handling, documentation and CI | P07, P09 |
| P11 | Independent numerical, compatibility and integration gates | P10 |
| P12 | Integrate exact candidates, create repos, push and hand off | P11 |

`phases.csv` is the machine-readable DAG. Phase numbers are labels, not a requirement to serialize all work. After P03, P04/P05/P06/P08 can be scheduled independently subject to worker capacity and file ownership. One coordinator plus at most three workers is the default. Use sequential role passes if delegation is unavailable.

## Completion

Completion requires working packages, the inventoried functionality assigned and implemented, base-R core checks, actual ggplot and pinned DescrTab2 import/render checks, legacy compatibility evidence, checked source tarballs, and the authorized repository actions. Optional installation does not make validation optional for a claimed backend. A scaffold, a stub registry, a skipped feature check, or generated CI YAML is not completion. G11 is the required local candidate gate; hosted checks are tracked separately so creating the new repositories can enable them. A blocked remote action does not prevent completion of local implementation and evidence; it does prevent claiming remote completion.

The final handoff identifies repositories, branches, exact SHAs, versions, required and optional dependencies, checks actually run, baseline failures, unsupported capabilities and the remaining concrete blockers. Do not claim CRAN submission, a public release, a merge or complete scientific validation.

## Files

- `START_HERE.md`: controlling orchestrator prompt.
- `stack.dcf`: execution defaults readable by `base::read.dcf()`.
- `phases.csv`: phase prerequisites and ownership areas.
- `prompts/00-shared-rules.md`: mandatory cross-cutting rules.
- `prompts/01-*.md` through `12-*.md`: implementation stages.
- `prompts/90-independent-review.md`: reviewer mandate.
- `prompts/99-resume.md`: resuming a long or interrupted run.
- `templates/`: state, evidence and migration-inventory templates.
- `references/source-notes.md`: the inspected baseline and primary sources.

