# Extended independent numerical coverage

`check-coverage.R` complements `check-numerics.R` with 271 synthetic cases
and 1067 scalar comparisons. It uses an installed candidate, its ordinary
imports, base R and jsonlite. Optional digest supplies checker hashes. It
needs no captured data, Python, private file, development namespace or network.

```sh
R_LIBS=/path/to/candidate-library \
INFERGUIDR_EXPECTED_LIBRARY=/path/to/candidate-library \
Rscript --vanilla tools/validation/check-coverage.R /path/to/new-coverage-results
```

The output directory must be new or empty. The expected-library guard fails
before numerical checks if another installation is loaded. The driver streams
`results.csv` and `cases.csv`, retaining package warnings and errors. Its
`summary.json` records runtime, package path/version, source hashes, fixed
coverage count and failures. `INFERGUIDR_SOURCE_SHA256` is recorded as
caller-supplied identity; retain the source archive and installation receipt
separately. It also writes the complete registry inventories described in
[NUMERICAL-CONTRACTS.md](NUMERICAL-CONTRACTS.md). Runtime was 52 seconds with
one process and one BLAS thread on the local Linux check; other systems vary.

An ordinary run exits zero only when all scheduled comparisons pass. Append
`--negative-control` with a different output directory to add .01 to the first
observed power. The expected outcome is exit 1, one failed power comparison,
and `negative_control_injected: true`. Setup failures can precede creation
of the summary. All numerical tolerances are fixed in each recorded row.

## Scope and independent calculations

| Branches | Reference route |
|---|---|
| t means, point-biserial correlation, slope/intercept comparisons, single regression coefficient, ANOVA/ANCOVA; forward, sensitivity, criterion and compromise | Input-derived df/noncentrality and normal/chi-square density integration |
| RM ANOVA effect conventions and all four MANOVA approximations under both algorithms | Separate parameter formulas and F convolution |
| Variance ratios, generic chi-square and goodness of fit | Integrated log densities, inverted quantiles and noncentral chi-square convolution |
| Binomial/sign balancing rules, McNemar and Fisher | Direct count/table enumeration and independent tail balancing |
| Boschloo elementary designs | Hand-enumerated n1=n2=2 rejection regions; the separate public checker covers larger exact conditional orderings |
| Exact correlation and random-model regression, exact/approximation branches and inverse modes | Hypergeometric correlation/R-squared densities and Lee F approximation, with independent integration |
| Nonzero risk/odds-ratio score offsets, ordinary and small boundary designs, sensitivity and a-priori | Restricted-null score roots, complete table probabilities and full design scans |
| Three rank methods, normal/logistic/Laplace parents and inverse modes; manual moments | Parent-distribution quadrature, ARE or Lehmann moment approximation |
| Independent/dependent correlations, generic normal and seven effect utilities | Explicit covariance/normal-tail and algebraic formulas |

Manual-moment rank sensitivity is intentionally refused with `invalid_input`
for all three methods. Those three cases verify refusal, not numerical output.
Their criterion and compromise cases are independently checked. The 271
cases also include seven utility cases outside the 50-method power registry.

Combined with the 313-case numerical campaign, the current local sample is
584 cases and 2057 comparisons. Every one of the 50 registered methods has
at least one independent numerical representative. Of 243 method/mode
declarations, 124 have numerical samples, two have only the intentional
exact-compromise guard, and 117 remain uncovered by these two campaigns.
The inventory exposes 1511 marginal option values and 25660 finite option
candidates; registration is not evidence that every combination is feasible
or tested. Historical evidence is separate and is never silently imported.

## Remaining qualification boundaries

- Exact unconditional compromise still needs a reviewed attainable-size and
  no-solution policy. The current guard is required by its discrete contract.
- Offset-score and inverse samples cover selected domains; they do not prove
  every boundary, larger design, root multiplicity or global monotonicity.
- GLM large exemplary-data grids remain outside this bounded campaign.
  Large-domain and native-reference disagreements retain their own evidence.
- Random-model and rank inverse samples do not certify every parameter
  domain, moment-admissibility condition or historical display convention.
- Multivariate formula checks cover the scheduled mappings; tetrachoric
  small-cell robustness, approximation confidence intervals and bootstrap
  coverage require their separate contracts and evidence.
- All-family representative coverage is not full numerical certification,
  external-application parity or cross-platform acceptance.

## Reference provenance

The standalone reference functions are MIT-licensed project code. `regression.R`
contains the independent G5 quadrature/density functions; `rank.R` the G10
parent and moment functions; `manova.R` the G4 alternative parameter mappings;
`discrete.R` the G8 count/tail calculations; `offset.R` the G9 scalar score-root
and table functions; and `variance.R` the G11 density/quadrature test reference.
Only pure algorithms are included: no original drivers, captures, recorded
results or private paths. Other reference files are described in the
[numerical contracts](NUMERICAL-CONTRACTS.md#reference-provenance).

These checks share standard formulas, floating-point arithmetic and some
central base-R probability primitives with the package. Distinct integrations,
root searches and exhaustive enumerations supply algorithmic cross-checks;
they are not independent implementations of every underlying library.
