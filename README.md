# inferguidR <picture class="cttir-readme-logo"><source class="cttir-logo-dark-source" media="(prefers-color-scheme: dark)" srcset="man/figures/logo-with-name-dark.png" /><img src="man/figures/logo-with-name.png" align="right" alt="inferguidR logo" width="120" height="139" /></picture>

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

See the [validation guide](VALIDATION.md) for evidence levels, reproducible
public checks and known limitations, and the [improvement roadmap](ROADMAP.md)
for the numerical, product and release acceptance gates.

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
- [Understanding validation status](https://github.com/CTTIR/inferguidR/blob/main/vignettes/validation.Rmd).

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

## Acknowledgements and sources

inferguidR is part of an inference suite with infercalcR and infereportR. The
suite builds on the work of others. We gratefully acknowledge that work here and
cite it, so that every method can be traced to its source.

### DescrTab2, Institute of Medical Biometry, Heidelberg University

The suite is deeply inspired by **DescrTab2** [1], the R package for
publication-quality descriptive statistics tables developed at the Institute of
Medical Biometry (IMBI), Heidelberg University, by Jan Meis, Lukas Baumann,
Maximilian Pilz and Lukas Sauer, with contributions from Lorenz Uhlmann and
Csilla van Lunteren. Its documented test-choice tree, its validation statement
and its practice of computing each test with an established implementation (R's
stats package, exact2x2, DescTools and nlme) are the model for how the suite
chooses, computes and reports statistical tests. The inferguidR method advisor
maps each of its recommendations to the DescrTab2 test-choice tree
(`inst/advisor/descrtab2-tree.json`). DescrTab2 was archived on CRAN on 18
October 2025; its source code remains available on GitHub. We thank its authors
and IMBI Heidelberg for making their work openly available.

### Methodenberatung, University of Zurich

The methods-consulting pages of the University of Zurich [2] inspired the
guidance in this suite: their decision tree and step-by-step decision assistant
for choosing a statistical method, and their worked examples and explanations of
each method. The inferguidR method advisor maps its questions to that decision
tree (`inst/advisor/uzh-tree.json`) and records where and why it departs from
it, for example where the pages choose between tests with a preliminary test.

### G\*Power

The power and sample-size procedures of inferguidR follow the test inventory of
G\*Power 3.1 [3, 4], and the distribution calculator of infercalcR follows the
G\*Power 3.1.9.7 calculator. Their results are compared with values recorded
from G\*Power 3.1.9.7. G\*Power is neither distributed with nor called by these
packages.

### Software

inferguidR is written in R [5] and imports infercalcR [6], digest [7] and
jsonlite [8]. The optional workbench uses shiny [9] and ggplot2 [10]. Tests
compare Boschloo results with Exact [11] and classifier metrics with pROC [12]
and yardstick [13]. Boschloo's test follows Boschloo (1970) [14]. Use
`citation("inferguidR")` for this package and `citation()` or
`citation("<package>")` for R and the packages you used.

### References

1. Meis J, Baumann L, Pilz M, Sauer L (2022). *DescrTab2: Publication Quality Descriptive Statistics Tables*. R package version 2.1.16. <https://doi.org/10.32614/CRAN.package.DescrTab2>. Source code, including the development version 2.1.26: <https://github.com/imbi-heidelberg/DescrTab2>.
2. Universität Zürich, Methodenberatung. *Datenanalyse mit SPSS*: decision tree (Entscheidbaum), decision assistant (Entscheidassistent) and method pages with examples and explanations. <https://www.methodenberatung.uzh.ch/de/datenanalyse_spss.html>. Accessed 3 October 2026.
3. Faul F, Erdfelder E, Lang A-G, Buchner A (2007). G\*Power 3: A flexible statistical power analysis program for the social, behavioral, and biomedical sciences. *Behavior Research Methods*, 39(2), 175-191. <https://doi.org/10.3758/BF03193146>.
4. Faul F, Erdfelder E, Buchner A, Lang A-G (2009). Statistical power analyses using G\*Power 3.1: Tests for correlation and regression analyses. *Behavior Research Methods*, 41(4), 1149-1160. <https://doi.org/10.3758/BRM.41.4.1149>.
5. R Core Team. *R: A Language and Environment for Statistical Computing*. R Foundation for Statistical Computing, Vienna, Austria. <https://www.R-project.org/>.
6. Heller R. *infercalcR: Statistical Inference and Diagnostic Calculations*. R package. <https://github.com/CTTIR/infercalcR>.
7. Eddelbuettel D. *digest: Create Compact Hash Digests of R Objects*. R package. <https://doi.org/10.32614/CRAN.package.digest>.
8. Ooms J (2014). The jsonlite package: A practical and consistent mapping between JSON data and R objects. arXiv:1403.2805. <https://arxiv.org/abs/1403.2805>.
9. Chang W, Cheng J, Allaire JJ, Sievert C, Schloerke B, Aden-Buie G, Xie Y, Allen J, McPherson J, Dipert A, Borges B. *shiny: Web Application Framework for R*. R package. <https://doi.org/10.32614/CRAN.package.shiny>.
10. Wickham H (2016). *ggplot2: Elegant Graphics for Data Analysis*. Springer, New York. <https://doi.org/10.1007/978-3-319-24277-4>.
11. Calhoun P. *Exact: Unconditional Exact Test*. R package. <https://doi.org/10.32614/CRAN.package.Exact>.
12. Robin X, Turck N, Hainard A, Tiberti N, Lisacek F, Sanchez J-C, Müller M (2011). pROC: an open-source package for R and S+ to analyze and compare ROC curves. *BMC Bioinformatics*, 12, 77. <https://doi.org/10.1186/1471-2105-12-77>.
13. Kuhn M, Vaughan D, Hvitfeldt E. *yardstick: Tidy Characterizations of Model Performance*. R package. <https://doi.org/10.32614/CRAN.package.yardstick>.
14. Boschloo RD (1970). Raised conditional level of significance for the 2 x 2-table when testing the equality of two probabilities. *Statistica Neerlandica*, 24(1), 1-9. <https://doi.org/10.1111/j.1467-9574.1970.tb00104.x>.

## Authors, citation and licence

Raban Heller. Use `citation("inferguidR")` for installed-version metadata and
record the exact repository revision for a development analysis. The package is
licensed under [MIT](LICENSE.md). Report reproducible problems with synthetic
inputs through [GitHub issues](https://github.com/CTTIR/inferguidR/issues).

## Use of LLM tools

Where applicable, LLM tools were used for narrowly defined, non-authorial tasks:
copyediting, prose refinement, Markdown/LaTeX formatting, boilerplate scaffolding
(such as CI configurations and build scripts), and code refactoring. Tools included
[Chat AI](https://kisski.gwdg.de/leistungen/2-02-llm-service/), the KISSKI
(GWDG) service, and self-hosted **Mistral Small (24B, Apache-2.0)** served through
a local **vLLM architecture** on an Ubuntu Linux workstation with Intel
processors and NVIDIA GPUs. Self-hosted use also includes models that were
further trained (fine-tuned) or weight-edited in-house for specific tasks.
All self-hosted models ran locally, with no data sent to third parties.

All scientific claims, methodological choices, analyses, interpretations and
conclusions remain the authors' own. Where applicable, LLM-generated material
must be reviewed and revised before inclusion, and references must be verified
against their DOI, arXiv ID or ISBN.

Local software setup: vLLM with Python/PyTorch, plus R and RStudio.

## Branding

Editable R artwork and regeneration instructions are in [tools/branding](https://github.com/CTTIR/inferguidR/blob/main/tools/branding/README.md).
