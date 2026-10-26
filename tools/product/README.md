# Product evidence from an installed package

These tools use public synthetic inputs and an installed package. They do not
load a development namespace, read the private qualification corpus, or assert
release readiness. Supply a new output directory for each run. To pin the
candidate installation, put its fresh library directory first in `R_LIBS_USER`
and set `INFERGUIDR_EXPECTED_LIBRARY` to that directory. Keep dependency libraries
later in `R_LIBS_USER` or the usual site-library paths. The expected-library
guard refuses an older inferguidR from any fallback dependency library.
`INFERGUIDR_SOURCE_SHA256` optionally records
the source archive hash supplied by the installation process.

```sh
Rscript tools/product/check-bootstrap.R evidence/bootstrap
Rscript tools/product/check-bootstrap.R evidence/bootstrap-negative --negative-control
Rscript tools/product/inventory-advisor.R evidence/advisor
```

The negative control must exit with status 1 and retain a failed summary. It
perturbs one lower interval endpoint by 0.01. A positive exit requires every
recorded arithmetic comparison to pass. Coverage diagnostics are descriptive;
they are never turned into a blanket coverage pass.

The bootstrap checker independently computes sensitivity, false-positive
rate, positive predictive value, accuracy, Brier score and pair-count ROC AUC
from sampled frequency weights. It implements type-7 endpoint interpolation
directly. Three designs use class-stratified, unstratified and whole-cluster
resampling, with 250 replicates and a fixed seed. The checker shares R's
documented multinomial random-number generator with the package; it does not
reuse package draw helpers, metric arithmetic or quantile code. Agreement
checks this finite fixture and seed, not independence of the RNG implementation
or broad statistical validity. Caller RNG restoration is also checked.

`coverage-diagnostic.csv` enumerates binomial sampling distributions for
5, 10 and 20 event units and true sensitivities 0.01, 0.1, 0.5, 0.9 and 0.99.
For each possible observed count it computes the limiting percentile interval
from the binomial empirical bootstrap distribution, then sums the true
sampling probabilities of intervals covering the truth. This is an exact
infinite-replicate diagnostic, not coverage of a finite replicate run or of
all classifier metrics. Boundary and small-sample undercoverage is expected.
The bootstrap remains `development`; these checks do not certify nominal
coverage. Stratification conditions on class totals and cannot estimate
population-prevalence uncertainty. Fixed-score/fixed-cutoff intervals also
omit training and cutoff-selection uncertainty.

The advisor inventory records all 69 variant identities, estimands, actual
template hashes and implementation declarations, route-specific assumptions,
allowed overrides, rejected routes, and the recorded relationship to both
bundled source trees. Source extensions are explicit. Historical call sketches
remain distinguishable from executable templates. The one rejected-only
variance route has no template. An inventory is a review aid: parsing,
source-link availability and executable teaching examples do not establish
general scientific validity, source equivalence or G*Power parity.

The [26 October 2026 review tables](reviews/2026-10-26/README.md) retain the
substantive review of all 69 variants, bibliographic identity checks for all
30 rule references, and scoped candidates for the 19 missing lesson mappings.
They distinguish actual executable templates, historical source judgments and
external references that have not been registered as lessons.

Portable classifier reports omit full score-derived ROC/PR/density tables.
`include_probability_summaries = TRUE` explicitly includes calibration-bin
means and counts plus aggregates at requested net-benefit thresholds, as
report schema 0.4. There are no fabricated curves or intervals. Score ranges
and individual rows are omitted, but small aggregate groups can disclose
sensitive information; sharing requires review of the included content.
Default reports keep their existing privacy boundary and schema versions.

Boschloo power reports also use schema 0.4 for two frozen possible-table
lattices (H0 and H1), up to the 5000-table display limit. Readback checks full
ordered support, conditional binomial probabilities, rejection masses and
the result's recorded rejection-region hash. It does not rerun a Boschloo
region solver. Above the display limit the report retains an explicit absent
figure and its reason. The lattice contains hypothetical tables, not observed
individual records, and does not opt into classifier probability summaries.
General XY power sweeps already retain every requested input and frozen
result, including failed-row reasons and gaps; point figures never interpolate
or create uncertainty bands.
