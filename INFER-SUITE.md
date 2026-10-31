# Inference suite development guide

`infercalcR` calculates observed analyses, descriptions and diagnostic metrics.
`inferguidR` guides method choice and plans power, sample size and study designs.
`infereportR` presents finalized results. Both consumers import infercalcR;
calculation has no dependency on either consumer, and reporting works without
guidance. This development migration does not establish universal scientific
validation, a public release or CRAN availability.

## Installation

Use the exact revisions and source archives identified by the migration handoff.
With sibling checkouts, install from their source archives in provider-first
order, into a library you control:

```sh
R CMD build infercalcR
R CMD INSTALL --library=/path/to/library infercalcR_0.0.0.9000.tar.gz
R CMD build infereportR
R CMD INSTALL --library=/path/to/library infereportR_0.0.0.9000.tar.gz
R CMD build inferguidR
R CMD INSTALL --library=/path/to/library inferguidR_0.0.0.9000.tar.gz
```

Set `R_LIBS_USER` to that library for builds and sessions. Plain R CMD INSTALL
does not fetch dependencies or honor GitHub Remotes automatically. GitHub
installation requires authenticated access to the private provider/report
repositories. Keep credentials outside source files and command history.

Core calculations and plain reporting do not require ggplot2 or a testing
framework. Legacy identity-bearing APIs require digest/jsonlite; their absence
is an explicit feature error. Existing canonical JSON and SHA-256 are preserved,
including cross-platform number conversion. Specialist native adapters declare
their optional backend and retain the original procedure/settings. Installing
an additional backend does not select a different method.

## Ordinary observed analysis and presentation

```r
spec <- infercalcR::method_spec("t_means_two_independent_welch")
result <- infercalcR::run_inference(spec,
  x = c(5.2, 6.1, 4.9, 7.3, 5.8, 6.5, 7.1, 4.7),
  y = c(3.8, 5.9, 4.2, 5.1, 6.7, 4.8, 3.9, 5.3, 6.0))
infercalcR::result_data(result)
infereportR::report_table(result)
infereportR::report_text(result)
```

The effect direction is x minus y. Prespecify sidedness and confidence level;
no observed normality test chooses the method. Results retain raw numerical
values, population/exclusion metadata and the exact estimand. Rendering does
not repeat the test. `describe_data()` and `report_descriptive()` provide
continuous/categorical tables with explicit denominators and missing counts;
they do not automatically add inferential p-values.

With ggplot2 installed, add explicit original-unit annotations without raw
analysis data in the result:

```r
plot_data <- data.frame(group = rep(c("x", "y"), c(8, 9)),
  value = c(5.2, 6.1, 4.9, 7.3, 5.8, 6.5, 7.1, 4.7,
            3.8, 5.9, 4.2, 5.1, 6.7, 4.8, 3.9, 5.3, 6.0))
positions <- data.frame(row_id = infercalcR::result_data(result)$row_id[1],
  xmin = 1, xmax = 2, y_position = 8)
ggplot2::ggplot(plot_data, ggplot2::aes(group, value)) +
  ggplot2::geom_boxplot() +
  infereportR::stat_inference(result, positions, rows = positions$row_id)
```

Facet-specific results require explicit facet keys. Coordinates control display
only; transformed scales and repeated plot builds never recompute inference.
For multiplicity, adjust the complete declared family in calc before filtering
rows for presentation. Requesting unavailable adjusted p-values is an error.

## Matched Boschloo analysis and planning

```r
control <- infercalcR::boschloo_control()
planned <- inferguidR::boschloo_power(15, 15, .2, .6,
  alpha = .05, variant = "minlike", control = control)
observed <- infercalcR::boschloo_test(3, 15, 9, 15,
  alpha = .05, variant = "minlike", groups = c("control", "treatment"),
  event = "response", control = control)
infereportR::report_table(infercalcR::as_inference(observed))
```

The same rejection engine and numerical controls serve both paths. Central,
minlike and prespecified directional variants are distinct. A numerical upper
bound for the p-value is not a sampling confidence interval. An unresolved
decision stays unresolved. Discrete sample-size searches retain their stated
allocation sequence and finite stability window.

`planning_methods()` records explicit relations between observed procedures
and available planning engines. Pooled-t planning is a comparator for Welch;
rank-test planning uses declared approximations. The existing advisor bridge
remains specifically scoped to Boschloo. An observed estimate never silently
becomes a prospective assumption.

## Classification with evaluation context

```r
d <- data.frame(id = paste0("unit", 1:6),
  truth = c("yes", "no", "yes", "no", "yes", "no"),
  score = c(.9, .2, .7, .4, .6, .3))
spec <- infercalcR::classification_spec(d, truth = "truth",
  positive_label = "yes", scores = "score", id = "id", cutoff = .5,
  score_scale = "probability", evaluation_split = "held_out_test")
result <- infercalcR::evaluate_cutoff(spec)
infereportR::report_table(infercalcR::as_inference(result))
```

This synthetic example records a held-out evaluation declaration; the label
does not prove independence. A cutoff selected on tuning data retains that
origin. Aggregate counts cannot supply a score distribution or ROC curve.
Undefined metrics retain their row status and denominator. Portable reporting
omits patient rows; the original source-data fingerprint continues to protect
reattachment and mutation checks.

## Existing users and saved analyses

Existing inferguidR exports remain available through retained guidance code or
tested compatibility routes. Continue reading saved JSON with
`inferguidR::ig_from_json()`. Historical content and identity are not rewritten.
New results may record the new producing package, so numerical preservation is
checked separately from historical hash preservation.

`infercalcR::as_inference()` converts supported finalized objects to the common
view. It neither selects a new test nor fits a replacement model. For legacy
planning JSON in a fresh session, read with inferguidR first so its accessor
methods are registered; retain a common result for subsequent reporting
without guidance installed. Trusted RDS fixtures are not an untrusted-input
format.

## Verification

From an inferguidR checkout, run:

```sh
Rscript tools/ci/check-infer-suite.R /path/to/infercalcR /path/to/inferguidR \
  /path/to/infereportR /path/to/check-output
```

The runner builds archives and checks installed packages, records environment
and logs, and verifies a genuinely isolated minimal core. Full integration
requires the declared ggplot2, pinned DescrTab2 and specialist backends.
Inspect check warnings and the lane evidence; an optional-installation guard
does not demonstrate that a feature was tested. Existing testthat regressions
are retained until mapped replacement checks have executed. New tests use
ordinary R scripts.

Hosted consumer checks require a credential authorized to read the private
provider repository. A workflow's default token does not automatically have
cross-repository access. Local checks and hosted checks have separate states;
workflow configuration alone is not evidence of a passing run.
