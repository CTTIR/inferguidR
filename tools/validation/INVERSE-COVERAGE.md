# Independent inverse and minimality checks

`check-inverses.R` adds 144 synthetic cases and 598 scalar comparisons to the
public numerical campaigns. It uses an installed inferguidR package, ordinary
imports, base R and jsonlite. Optional digest records checker hashes. It uses
only the standalone MIT-licensed algorithms in `reference/`; no private files,
captures, development namespace, network or Python are needed.

```sh
R_LIBS=/path/to/candidate-library \
INFERGUIDR_EXPECTED_LIBRARY=/path/to/candidate-library \
INFERGUIDR_SOURCE_SHA256=archive-sha256 \
Rscript --vanilla tools/validation/check-inverses.R /path/to/new-inverse-results
```

The output directory must be new or empty. The library guard fails before
numerical checks when the loaded installation differs from the specified
candidate. The archive hash is caller-supplied provenance, not independently
verified by the R process. Keep the archive and installation receipt with it.

`results.csv` records every comparison and its tolerance. `cases.csv` records
inputs, options, status, solved values, diagnostics, package warnings, errors
and elapsed time. `summary.json` includes package identity, R/platform details,
checker hashes, the expected case count and failure count. Complete registry
inventories are emitted using the same code as the other campaigns. A normal
run exits zero only if all 144 scheduled cases pass. Local single-process,
single-thread runtime was about one minute; other platforms vary.

Append `--negative-control` with a different output directory to perturb the
first observed power by .01. This must produce exactly one failed power
comparison and exit 1. A wrong expected-library path must also fail.

## Independent calculations and design grids

The oracle receives resolved input values, including a reported solve target.
It does not use the engine's power, degrees of freedom, noncentrality, critical
values or solver diagnostics to derive the expected result. It independently
recomputes the forward probability at the solved inputs. Standard central
base-R quantiles and probability primitives are shared in some branches;
alternative tails use independent integrals, formulas or enumeration. This
is algorithmic cross-checking, not independence from every numerical library.

The F reference in `reference/inverse-f.R` integrates the denominator on its
log scale, retaining the exponential Jacobian. Its inner noncentral
chi-square convolution uses a sine coordinate and the corresponding cosine
Jacobian, with the complete integration interval split around the shifted
normal density's center. These coordinates regularize endpoint behavior
across R versions. Inner and outer relative quadrature targets remain
`1e-13` and `1e-11`, with zero absolute target; the two omitted denominator
tails bound truncation error by `2e-14`. No smaller sample size is skipped.

| Families | Inverse evidence |
|---|---|
| t means, slopes/intercepts, point-biserial correlation, single coefficient, ANOVA/ANCOVA | A priori size with independently mapped parameters and noncentral t/F integrals |
| Hotelling, fixed-model R-squared, RM ANOVA and MANOVA | A priori, sensitivity, criterion and compromise; Pillai/O'Brien-Shieh for the selected MANOVA inverse cases |
| Generic t and F | Sensitivity, criterion and compromise |
| Goodness of fit, variances | A priori; variance sensitivity, criterion and compromise, including lower and upper alternatives where supported |
| Exact correlation and random-model R-squared | A priori through separate density integrals under zero null correlation/R-squared |
| Rank methods | A priori for all three methods under normal, Laplace and logistic parent distributions using Lehmann moments |
| Independent/dependent and tetrachoric correlations | Missing inverse modes, with both exact-asymptotic and Bonett-Price tetrachoric power branches |
| Normal proportions | A priori, sensitivity and compromise, with and without the arcsine transform and without continuity correction |
| Binomial/sign, McNemar, Fisher | All four inverse modes through complete count/table enumeration; both McNemar computation branches |
| Boschloo | Hand-enumerated minimal N=4 cases for all four variants; the earlier N=2 design has no rejection at the selected alpha |
| Unconditional nonzero offset | Exact criterion with nominal and attained alpha checked separately |

A-priori cases scan **every earlier admissible integer design** through the
reported size. This does not assume global monotonicity. Equal two-group
allocations use even totals; one-way/RM balanced designs use group-count
multiples; special-effects and within-between MANOVA use their declared
integer grid. The tests explicitly cap returned N at 400. No unbounded
parameter sweep or large table lattice is constructed.

Continuous sensitivity and criterion solutions must satisfy an independent
forward residual; compromise checks beta against nominal alpha for q=1.
The default comparison tolerance is `1e-8 + 1e-6 * abs(reference)`, with
stricter row-specific tolerances for selected finite enumerations and exact
integer-size comparisons. These are ordinary-domain checks, not extreme-tail
relative-accuracy qualifications.

## Discrete boundaries

Binomial, McNemar and Fisher rejection regions are step functions of nominal
alpha. Criterion returns the smallest nominal alpha whose power reaches the
target. Compromise returns the first alpha with beta <= q*alpha if an exact
equality is unattainable. Neither rule is an equality with attained size.

An inverse solver commonly returns an alpha precisely at a jump. Independently
summed null probabilities can differ in their last binary digit; evaluating
both implementations at that one floating-point number can then select
different regions. The checker therefore locates the first attaining alpha
independently by 60 bisection steps over enumerated probabilities, compares
the reported alpha within `2e-12 + 1e-10 * abs(reference)`, and evaluates the
canonical attaining region eight machine epsilons above that reference
boundary. It checks reported power and attained size against that region, and
checks that a relative `1e-7` step to the left fails the target or ratio.
This qualifies the threshold to the stated tolerance, not bit-identical CDF
rounding at a discontinuity. Original unsuccessful development receipts are
retained separately; tolerances are fixed in this driver.

## Coverage and remaining limits

Together, `check-numerics.R`, `check-coverage.R` and this driver contain 728
cases and 2655 comparisons. Their finite sample covers all 243 declared
method/mode pairs: **241 have independent numerical evidence and two have
only intentional exact-compromise refusal evidence**. All 50 methods have
independent numerical representatives. This fills the 117 method/mode gaps
listed for the earlier two campaigns; it does not fill every finite option
combination or continuous parameter domain.

The following remain separate requirements:

- Exact unconditional compromise needs an implemented attained-size solver
  and an explicit unattainable-ratio policy. Its existing guard remains; this
  campaign does not implement or qualify that mode.
- Unequal allocations, very large designs, multiple-root selection, extreme
  tails, all option interactions and all effect-direction/boundary cases need
  their own scoped evidence. Selected full integer scans do not certify every
  design search.
- Historical native-reference disagreements and the preserved G6 memory
  failure remain separate records. No historical result is silently promoted
  because the source hash is unchanged.
- Tetrachoric observed-frequency approximation confidence intervals remain
  blocked by competing estimator conventions. The two tetrachoric **power**
  branches tested here do not resolve that confidence-interval contract.
- Bootstrap coverage, complete utility-domain checks, cross-platform results
  and external-application parity have their own qualification gates.

See [EXTENDED-COVERAGE.md](EXTENDED-COVERAGE.md) and
[NUMERICAL-CONTRACTS.md](NUMERICAL-CONTRACTS.md) for reference provenance and
other campaign boundaries. Method/mode sample coverage is not full numerical
certification or release approval.
