# P07 Guidance, planning and compatibility integration

Role: inferguidR integration owner. Prerequisites: G04, G05, G06.

## Objective

Make inferguidR consume infercalcR while retaining its guidance and study-planning identity and its existing public behavior.

## Planning ownership

Retain method-specific power models, prospective assumptions, planning specifications, allocation rules, inverse/sample-size search, sensitivity/criterion/compromise modes, sweeps and design comparisons in guid. Shared distributions/rejection rules may be supplied through public calc capabilities. Preserve current numerical precision and discrete-search contracts.

Replace namespace-local discovery that now crosses package boundaries with explicit providers. Existing scans of guid for .ig_define_methods_* and lookups of .ig_drawer_* require review. Do not simply point get0()/getFromNamespace() into calc. Methods remaining wholly in guid can keep a sound local mechanism until a targeted improvement is needed.

## Method-to-planning mapping

Maintain explicit versioned mappings from observed method specifications to planning engines. Label each as same-procedure under stated assumptions, approximation, simulation, comparator, or unsupported. List required prospective inputs. Keep calculator status separate from evidence status.

Preserve concrete distinctions from the current catalogue: Welch versus pooled t, Welch ANOVA versus pooled F, exact versus asymptotic McNemar, and finite-sample rank inference versus large-sample planning approximations. Boschloo variants and direction/numerical controls must agree with calc. A common method name is not sufficient.

If simulation is supported, require a declared data-generating model/design, repetitions/seed, analysis decision rule and Monte Carlo uncertainty. Failures/non-estimable repetitions have an explicit policy; do not silently remove them from the denominator. Do not add a whole new simulation platform merely to hide missing analytic mappings.

## Advisor and examples

Keep questions, decision policies, rationale, sources, learning content and Shiny guidance in guid. Recommendations refer to calc method IDs and capability metadata. They remain recommendations, not claims that data were analyzed. Generated scripts call supported public calc APIs and show explicit missingness, pairing and method settings. They do not source bundled analyser scripts at runtime.

## Compatibility

Implement tested guid forwarding functions for moved exports. Preserve accepted argument names, return shapes/classes/conditions and documented failure behavior where required by compatibility.csv. Build legacy spec/result wrappers around the single calc engine when necessary. Do not retain a second divergent copy of the numeric implementation.

Register shared accessor methods for power/design/classifier objects. Test fresh-process behavior. Old saved JSON readers remain available; do not change historical schema meaning or hashes. Newly computed provenance may change and must be labelled/versioned deliberately.

Keep existing plot/export entry points usable. Retain legacy native reporting temporarily or add an optional guid-to-report forwarding path only when report does not depend on guid. Do not make installing ordinary calculation depend on the plotting/reporting stack. Existing optional UI test tools may remain until equivalents are justified.

## Verification

Run baseline-versus-candidate checks for representative planning modes, allocations, statuses, unsupported branches, advisor mappings, generated examples, compatibility wrappers and serialized legacy documents. Install exact core candidate first, then guid in an isolated library. Execute current public examples. Record any baseline failure separately.

## Gate G07

Guid owns guidance/planning and uses calc through public APIs. Existing supported calls have tested compatibility. Numerical calculations, discrete search behavior, method identity, legacy document identity and privacy rules remain coherent. New code does not introduce a dependency cycle or infer same-method planning where none was established.

