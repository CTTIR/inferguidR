# Numerical contracts and remaining gates

`check-numerics.R` runs a fixed **313-case synthetic campaign** against an
installed package. It extends the compact `check-public.R` smoke campaign.
Every scheduled case is scored, including non-OK outcomes and exceptions.
Rows are written after each case, so an interrupted run retains completed
evidence. A complete run currently has 990 scalar comparisons. These counts
are coverage checks, not a claim that all numerical domains are certified.
The complementary [extended campaign](EXTENDED-COVERAGE.md) adds 271 cases
covering the remaining families, random-model and rank inverse branches,
multivariate parameter mappings and risk/odds-ratio offsets. The table below
describes the scope of this 313-case driver; the extension has its own scope
and inventories.

```sh
R_LIBS=/path/to/candidate-library \
INFERGUIDR_EXPECTED_LIBRARY=/path/to/candidate-library \
Rscript --vanilla tools/validation/check-numerics.R /path/to/new-numerical-results
```

Requirements are the installed package and its ordinary imports. The script
uses base R and jsonlite, with digest for optional source hashes. It does not
use Python, pkgload, captured results, network access, or private files. Each
output directory must be new or empty. The driver uses one process and no
workers. Its enumeration cases have N=60; its reference Poisson enumerator
also has a hard two-million-row cap and a mean-700 recursion limit. These
bounds prevent the historical unbounded exemplary-data allocation pattern.

`INFERGUIDR_SOURCE_SHA256` and `INFERGUIDR_SOURCE_COMMIT` are recorded as
caller-supplied identity, not verified installation provenance. Keep build
and install receipts alongside the outputs. The expected-library guard
prevents accidentally testing another installed candidate.

## Exact unconditional compromise

For these methods, nominal alpha determines a normal or t critical value.
The attained alpha is the probability of rejection at the specified null
point. It is neither the nominal alpha nor a nuisance-parameter supremum.
The intended compromise equation is

`beta = error_ratio * alpha_attained`.

At fixed integer group sizes, moving a deterministic critical threshold
produces finitely many rejection regions. Both sides of this equation are
step functions. Equality therefore need not exist. Solving the generic
nominal-alpha equation, interpolating between attained sizes, or returning
the nearest ratio would change this contract. Randomised rejection at a
boundary would define a different test and is not offered.

The campaign enumerates every threshold region for the synthetic design
n1=4, n2=6, p1=.7, p2=.3, pooled statistic, one-sided rejection, q=1. None
satisfies the equality. `exact-compromise-states.csv` exposes every attained
size, beta, and residual. Both unconditional method variants must continue
to return `reference_pending` and no alpha estimate in exact compromise.
This verifies an intentional refusal; it does not implement the missing
attained-size solver. A future implementation needs a reviewed policy for
unattainable ratios, an explicit outcome for no solution, and independent
tests over all statistic and offset branches. Large-sample compromise uses
its declared continuous approximation and is a separate branch.

## Discrete critical-boundary correction

For n1=8 and n2=12, the six tables (2,0), (0,4), (6,4), (2,8), (8,8),
(6,12) have pooled squared statistic exactly 10/3. Computing proportions
and standard errors separately previously produced four floating-point
values for this one equality class. An inverse alpha search could include
only part of the tied class.

The no-offset pooled, corrected pooled and Mantel-Haenszel statistics now
use signed square roots of reduced rational squared statistics formed
directly from counts. Integer common factors are cancelled before multiplying
by N or N-1. This prevents another rounding split observed near N=10000
when an unreduced numerator product exceeded 2^53.
The criterion regression determines membership by integer cross
multiplication against 10/3, independently of the engine statistic. The
test suite additionally checks all rational equality classes for each of
these three statistics over n1,n2=1,...,8, and their zero-count corners. It
checks the six-table pooled/Mantel equality class at all 500 integer scales
up to N=10000, and complementary tables at unequal N=10000 designs.
This fixes the demonstrated boundary defect without changing nominal-alpha
semantics, offset formulas, or reference-pending statuses. Manually extended
exact designs whose required count products exceed the exactly representable
integer range return `numerical_failure`; callers can explicitly select the
large-sample approximation. Default exact designs through N=10000 stay
within the verified integer range. This is not a proof of exact arithmetic
for arbitrary large floating-point counts.

## Campaign scope and tolerances

| Work item | Checked supported contract | Intentional boundary or remaining work |
|---|---|---|
| N1 | Full finite rejection-state example; refusal by both exact compromise methods | Attained-size compromise solver remains unavailable; equality can be unattainable |
| N2 | Seven no-offset statistics, both tails, unequal designs and effects near 0/1; five nonzero-difference offset statistics; exhaustive a-priori scans; sensitivity directions and all four effect indexes; criterion boundary replay | Risk/odds-ratio offset score branches, wider inverse domains and native-reference ambiguities require broader qualification |
| N3 | GLM six covariate distributions, corrected and uncorrected procedures; small-grid Wald/LR enumeration; four inverse modes for binary/normal GLM; rank ARE and manual-moment Lehmann; two-proportion CC/arcsine cases; fixed-model regression and generic t/F | Synthetic coverage does not requalify every historical input, very large exemplary grids, random-model regression, or every rank parent/moment inverse branch |
| N4 | RM covariance/contrast quadratic forms; independent F probabilities for RM, MANOVA and Hotelling; Hotelling/regression df and noncentralities; generic t/F kernels | Evaluating an independent F tail at reported df/ncp alone does not validate every approximation's parameter mapping; full multivariate/utility domains remain broader than this campaign |
| N5 | Tetrachoric orthant quadrature, Pearson asymptotic variance, Bonett-Price numerical-derivative variance, and power in both tails | Exact cell integration is not finite-sample exact power; small-cell/high-correlation robustness and approximation-mode observed-frequency CI remain separate gates; bootstrap has its own contract/tests |
| N6 | All 50 methods, all 243 registered modes, all 1511 marginal option values and 25660 finite Cartesian option candidates enumerated | Registration and a sampled passing case do not establish coverage of every input or interaction; continuous domains and input-dependent feasibility are explicit |

Exact table probability comparisons use absolute tolerance 2e-12. Tied
integer membership, sample sizes and interface statuses are exact. GLM,
rank, t/F and noncentral distributions use absolute 1e-8 plus relative 1e-6.
Two-proportion and tetrachoric checks use absolute 1e-10 plus relative 1e-8.
Inverse target tolerances are stated per result row. Tolerances are fixed
before comparison. Observed discrepancies are not used to widen them.

The independent oracles use a different enumeration, integration, matrix
calculation, or root-search route from production. They are not wholly
independent of base-R central distributions, standard mathematical formulas,
or floating-point arithmetic. In particular the unconditional table reference
and production share the count-based scalar statistic equation; rejection
enumeration and the rational equality-class regression provide the distinct
checking routes. The F convolution checks distinguish probability kernels
from the separately checked parameter mappings.

## Machine-readable evidence

- `results.csv`: every scalar result, oracle, fixed tolerance and pass flag.
- `cases.csv`: every case, full synthetic inputs/options, status, warning data,
  exception text, elapsed time and pass flag.
- `summary.json`: runtime and package identity, reference hashes, coverage
  counts, failures and overall pass flag.
- `method-mode-inventory.csv`: every registered method/mode, including
  uncovered modes.
- `method-mode-option-inventory.csv`: every marginal declared option value;
  evidence resolves registry defaults and checks known applicability predicates.
- `method-mode-option-combinations.csv`: every finite Cartesian declaration
  candidate. Inactive options are separated in `active_options_json`.
  These are candidate declarations, not a claim of input feasibility. Unknown
  or input-dependent applicability cannot establish checked coverage.
- `registry-declarations.json`: complete installed declarations, retaining
  continuous domains and applicability conditions.
- `exact-compromise-states.csv`: finite rejection-state nonattainment proof.

Inventory evidence statuses are `uncovered`, `sample_failed`,
`sample_checked_not_certified`, and `guard_checked_no_numeric_solution`.
They describe this campaign only. Evidence from another run or an older
candidate is not silently imported into the matrix. A case matching only
one option value establishes a marginal sample, not every combination.

Append `--negative-control` with a different output directory to perturb the
first observed power by +.01 without changing package source. The run must
finish with exit status 1, one failed power comparison, and
`negative_control_injected: true`. Ordinary successful runs exit 0. Setup
failures can occur before a summary is available; per-case failures retain
results and allow later scheduled cases to run.

## Reference provenance

All reference code is MIT-licensed package-project code, provided here as
standalone functions rather than importing test helpers. `unconditional.R`
and `correlations.R` adapt the project's independent count/orthant/moment
test references. `multivariate.R` uses the helper-reference convolution;
`continuous.R` documents its existing provenance. `glm.R` adapts the pure
Gauss-Legendre, explicit-information and bounded enumeration functions from
the independent G6 checker. `proportions.R` adapts the independent FTU
bisection reference. No captured input/output corpus accompanies these files.
