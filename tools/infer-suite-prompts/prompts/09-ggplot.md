# P09 ggplot inference annotation layer

Role: plotting adapter owner. Prerequisites: G05 and G08.

## Objective

Implement `infereportR::stat_inference()` as an easy add-on to ggplot, similar in use to stat_compare_means, consuming finalized results and the common formatter.

## Implementation boundary

Prefer a small list of ggplot2::geom_segment()/geom_text() layers with explicit annotation data and inherit.aes = FALSE. All input normalization and coordinate preparation can use base R. ggplot2 is optional for package loading and plain reporting, but required when this function is called. Do not add ggpubr, rstatix, ggsignif or a testing framework unless a concrete feature justifies it. An explicitly selected ggpubr::stat_pvalue_manual backend can be added later; it is not required for this first layer.

Do not compute a test in a custom ggplot Stat. Scale transformations and out-of-range handling occur before layer statistics; display choices must not alter the scientific calculation. No fit, CI computation, multiplicity adjustment or cutoff optimization at plot construction/build/print time. Do not reevaluate stored calls.

## Required interface behavior

- Accept the shared result object and choose rows by stable IDs.
- Pairwise comparison brackets and text-only omnibus annotations.
- Label presets `p`, `estimate_ci`, `estimate_ci_p`, and optional `stars`, plus a documented custom template route that formats approved fields.
- Explicit raw/adjusted/reported p-value choice. An adjusted request with no adjusted value errors; no fallback.
- Original-data-unit y_position and explicit named group endpoints or numeric xmin/xmax. Document ordering of contrast labels separately from visual endpoint order.
- Preserve supplied facet key columns. Require an explicit facet mapping contract for facet-specific results; missing keys must not silently broadcast a subgroup result across panels. Explain that a layer cannot validate arbitrary later-added facets automatically without a plot-level helper.
- Support the documented horizontal path, e.g. coord_flip, and explicit dodge endpoints. Avoid claiming fully automatic positioning for arbitrary nonlinear coordinates/dodging/free scales.
- A bounded single-panel placement helper is optional. It must never extract raw data to recalculate inference. Do not replace a user's theme, scale or coordinates to force a label into view.
- A clear policy for unavailable/non-estimable results: retain them in report data; warn/error or show an unavailable text when requested. Never replace missing p-values with ns or invent a pairwise estimate from an omnibus test.
- Nonsignificant results are visible by default. Display filtering cannot change stored hypothesis families or adjusted values.

## Required base-R checks with ggplot2 installed

Write ordinary R test scripts using actual public ggplot2 functions. Inspect ggplot_build() output where suitable. No vdiffr/testthat requirement.

1. Annotate a valid result without retained raw analysis data.
2. Rebuild/print a plot repeatedly; stored results and labels remain unchanged.
3. Add scale_y_log10() to appropriate positive-coordinate examples; inferential values remain unchanged.
4. An explicit original y_position of 100 maps to log10 position 2 once, not a double transformation.
5. Distinct facet-key results land in the intended panels. Check duplicates/missing requested keys and no unwanted recycling.
6. Inherited fill/colour/group mappings do not split or change annotation data.
7. Reorder factor levels and raw plotting rows; named comparisons and scientific contrast sign remain correct.
8. Explicit numeric dodge endpoints connect the intended groups.
9. The documented horizontal plot builds with correct text/brackets.
10. Selecting a subset of a larger adjusted family preserves the original adjusted p-values and family metadata.
11. Missing results, ambiguous comparisons, incompatible coordinates and absent optional ggplot2 trigger the documented behavior.
12. The table, text and annotation formatter use the same selected numerical values and precision policy.

Produce a small set of representative plot files for visual inspection: two groups, several comparisons, facets, horizontal, log scale and unavailable inference. Verify legibility as well as structure. Do not treat visual snapshots as proof of statistical correctness.

## Gate G09

The actual ggplot-enabled integration lane ran successfully on installed candidate packages. Minimal core tests also prove package load and plain reporting without ggplot2. A headless session is sufficient: use a noninteractive grDevices::pdf() or another available device. Only an actual missing R/ggplot2/usable-device requirement makes this gate BLOCKED; lack of an interactive display alone does not. Document precisely which layouts are supported and the explicit-coordinate escape hatch.

