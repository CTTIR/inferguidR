# P05 Observed analyses using standard R

Role: observed-method adapter owner. Prerequisite: G03. Can be divided into disjoint method-family workers after the coordinator freezes method IDs and result schemas.

## Objective

Turn the currently used observed analyses into callable infercalcR APIs. Prefer standard R implementation and base data handling. Cover the inventory, not just one demonstration; do not advertise an advisor catalogue as a completed execution engine.

## Plan the work from capabilities.csv

Group ownership into disjoint families such as location/rank, categorical/exact, regression/ANOVA/repeated measurements, and association/variance/structure. Enumerate every assigned method, current source/template, supported settings, backend, expected outputs, and evidence. Compare against the mandatory_adapter list frozen in P01/P03. Every executable analysis template in that list needs a production adapter and executed evidence; not previously exported is not an exemption. Methods merely mentioned in advice without executable code can remain documented unsupported. Mandatory work that cannot execute is BLOCKED, not removed from the scope. Do not relabel a previously executable capability to avoid work.

Use `stats` where its procedure matches the approved method. Candidate calls include t.test, wilcox.test, binom.test, prop.test, chisq.test, fisher.test, mcnemar.test, cor.test, lm, glm, oneway.test, aov, kruskal.test, friedman.test, var.test, fligner.test, prcomp and kmeans, but actual support and defaults follow the source inventory and declared semantics. An available R function is not by itself proof that it matches a teaching template.

Keep specialist methods using an explicit optional backend when standard R does not provide the intended method reliably. Do not replace mixed models with fixed models, GEE with GLM, Brunner-Munzel with Wilcoxon, or a profile interval with Wald to meet a dependency target. Existing small R-native algorithms may be extracted after review; avoid inventing a new unvalidated advanced implementation during packaging.

## Method adapter requirements

- Explicit method ID/variant, settings, reference/event coding and contrast sign.
- Exact versus approximate versus simulated behavior; tie and continuity-correction policy. Never allow a silent backend fallback when exact computation fails or is unavailable.
- Null value, sidedness and confidence level. Do not infer a direction from data.
- Defined missingness behavior and counts. Exclusions are recorded; inference and displayed denominators must identify their populations.
- Pairing by explicit subject identifiers for long data, with duplicates/incomplete pairs handled by an explicit policy. Vector pairing requires documented aligned inputs.
- Clear failure/non-estimability behavior for empty/singleton/constant data, singular models, separation, convergence, zero variance and sparse cells.
- Estimates/intervals that target the same stated estimand. Preserve multi-parameter and omnibus results separately from coefficient tests.
- Deterministic backend choice. Installing another optional package must not silently change covariance estimates, degrees of freedom or confidence intervals.
- Explicit stochastic seeds and recorded computation settings for Monte Carlo/permutation/bootstrap procedures; preserve global RNG state when the function promises that behavior.

For `as_inference.htest`, map method-specific fields deliberately. Two group means are not themselves the estimate of their difference, and some htest objects have no effect estimate or interval. Confidence-level attributes, alternative and named null values matter. Require metadata for unrecoverable sample sizes, group labels or analysis populations.

For lm/glm, distinguish extracting an existing object from calculating requested numerical summaries once. All resulting estimates/CIs/contrasts are finalized before reporting. Define interval/covariance policies explicitly; do not let a generic confint dispatch unexpectedly load or require a specialist backend. For logistic GLM, distinguish odds ratios from risk ratios and preserve link/response scale. Avoid automatic variable selection or interpretive causal claims.

Multiplicity belongs in calc: declare full families, retain raw and adjusted values, do not re-adjust imported values, and keep interval multiplicity separate. Defaults must be documented and scientifically appropriate to the actual request, not selected for shorter code.

## Advisor consolidation

Extract callable computation from teaching templates; never source the whole template or run example code on namespace load. Provide P07 with regenerated example bodies that call calc's public APIs while retaining design explanations in guid. Preserve current supported settings and limits until a separately justified change is tested.

## Base-R verification

Test adapters against explicit matched `stats` calls, fixed known examples and independent invariants. Cover group order, subject order, NA/Inf, boundary probabilities, all-equal inputs, event orientation, method options and unsupported settings. Verify no unwanted package attachment, working-directory change or plot side effects.

Standard-R parity is evidence of correct wrapping, not independent revalidation of stats. Preserve distinct levels of numerical evidence. Keep old tests until their scientific checks are mapped and run.

## Gate G05

The completed adapter set matches the frozen mandatory list, and each mandatory adapter has installed-package tests and documentation. Ordinary methods run without contributed calculation/tidying/testing packages; optional exceptions are justified and feature-local, with executed evidence for their claimed support. The exact public signatures and updated capability rows are handed to P07/P09. Unimplemented or unverified mandatory entries prevent G05 PASS.

