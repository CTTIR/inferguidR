# inferguidR <img src="man/figures/logo.svg" align="right" height="139" alt="inferguidR hex logo" />

<!-- badges: start -->
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![License: MIT](https://img.shields.io/badge/license-MIT-yellow.svg)](LICENSE.md)
<!-- badges: end -->

**Guided inference, power analysis and diagnostic accuracy in R.**

inferguidR connects explicit statistical questions with validated analysis
specifications, typed results and reproducible exports. It includes a method
advisor, power and sample-size calculations, classification and cutoff tools,
and an optional Shiny workbench.

**Development software.** Implemented methods differ in their validation status.
A successful calculation is not proof of agreement with every external reference
or suitability for a particular study. The package is not a CRAN release.

## Installation

R 4.1.0 or later is required. Install the development version from GitHub:

```r
install.packages("remotes")
remotes::install_github("CTTIR/inferguidR", build_vignettes = TRUE)
```

Building the tutorials needs `knitr`, `rmarkdown` and Pandoc. Set
`build_vignettes = FALSE` if those are unavailable; the tutorial sources remain
readable in this repository. Core calculations do not require a reference
application, browser, database or external service.

For the optional workbench, install `shiny` and `ggplot2`. Some advisor-generated
analysis scripts use additional packages named in their output; inspect the
script and its requirements before running it.

## Quick start

This hypothetical two-group design assumes independent normal observations,
equal variances and a prespecified standardised difference. Its inputs illustrate
the interface; they are not a recommended study design.

```r
library(inferguidR)
spec <- analysis_spec(
  "t_means_two_independent", "a_priori",
  inputs = list(alpha = 0.05, power = 0.80, d = 0.5),
  options = list(tails = "two")
)
result <- solve_design(spec)
result$status
result$solved
result$achieved
result$validation
stopifnot(result$status$code == "ok")
```

Always inspect `status`, `warnings`, `assumptions` and `validation` before using
a number. For this example, `n_total` is the total across both groups; the
`group_sizes` field records the allocation. Unsupported, infeasible or unresolved
calculations carry explicit statuses rather than a usable design estimate.

## Save and reproduce

```r
json <- ig_to_json(spec)
restored <- ig_from_json(json)
stopifnot(identical(ig_hash(spec), ig_hash(restored)))
file <- tempfile(fileext = ".json")
writeLines(json, file)
packageVersion("inferguidR")
```

Keep the specification, result, package revision and `sessionInfo()` together.
Hashes identify content; they do not establish scientific correctness.

## Function guide

| Task | Entry points |
|---|---|
| Inspect available methods | `list_methods()`, `ig_method()` |
| Specify and calculate a design | `analysis_spec()`, `solve_design()`, `calculate_power()` |
| Compare design assumptions | `sweep_design()`, `plot_design_sweep()` |
| Evaluate diagnostic or classifier performance | `classification_spec()`, `evaluate_cutoff()`, `classification_metrics()` |
| Explore operating points | `cutoff_sweep()`, `select_cutoff()` |
| Work through an analysis question | `advisor_session()`, `advisor_question()`, `advisor_answer()` |
| Save specifications and results | `ig_to_json()`, `ig_from_json()`, `ig_hash()` |
| Explore interactively | `launch_app()` |

## Tutorials and workbench

- [Guided analysis and reproducible power calculations](https://github.com/CTTIR/inferguidR/blob/main/vignettes/getting-started.Rmd).
- [Classification metrics and explicit operating points](https://github.com/CTTIR/inferguidR/blob/main/vignettes/classification.Rmd).

After installing with vignettes, use `vignette("getting-started", package =
"inferguidR")` or `vignette("classification", package = "inferguidR")`.

Launch the workbench in an interactive R session:

```r
# install.packages(c("shiny", "ggplot2"))
inferguidR::launch_app()
```

## Troubleshooting and scope

An unknown method, an invalid probability or an inapplicable option is an input
error: inspect `?analysis_spec` and the selected method definition. A non-`ok`
result is a calculation outcome to inspect, not a value to coerce into a number.
Diagnostic examples are synthetic; selecting and evaluating a cutoff on the
same observations does not establish out-of-sample performance.

## Authors, citation and licence

Raban Heller. Use `citation("inferguidR")` for installed-version metadata and
record the exact repository revision for a development analysis. The package is
licensed under [MIT](LICENSE.md). Report reproducible problems with synthetic
inputs through [GitHub issues](https://github.com/CTTIR/inferguidR/issues).

## Use of LLM tools

Where applicable, LLM tools were used for narrowly defined, non-authorial tasks:
copyediting, prose refinement, Markdown/LaTeX formatting, boilerplate scaffolding
(such as CI configurations and build scripts), and code refactoring. Tools included [Chat AI](https://kisski.gwdg.de/leistungen/2-02-llm-service/), the KISSKI
(GWDG) service, and self-hosted **Mistral Small (24B, Apache-2.0)** served through
a local **vLLM architecture** on an Ubuntu Linux workstation with Intel
processors and NVIDIA GPUs. The self-hosted model ran locally, with no data sent
to third parties. Maintainers remain responsible for review, validation and
the final content.
