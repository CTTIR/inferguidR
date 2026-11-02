# P03 Shared contracts and base R verification harness

Role: computational contract owner. Prerequisite: G02. Coordinator integrates shared metadata.

## Objective

Implement and freeze the smallest public interface needed by all three packages. A registry of advertised names without real execution is insufficient.

## Public surface

Use the following names unless a documented current-source conflict requires a coordinator-approved change:

- `method_spec()`: validated identity/options of an observed statistical procedure.
- `run_inference()`: dispatches a supported method on supplied observations or supported aggregate inputs.
- `as_inference()`: explicit adapters for existing supported R results; conversion does not silently fit another analysis.
- `describe_data()`: a small standard-R descriptive provider for numeric/categorical data, explicit variable/group selection, missing counts, denominators and labels, returning kind=descriptive. It performs no hypothesis tests and does not treat ordinal category codes as interval measurements by default. The calc contract owner implements and tests it in this phase so P08's minimal report_descriptive() has a concrete provider.
- `result_data()` and `result_metadata()`: public extraction generics shared by computation, planning and reporting.
- A small capability listing and a constructor/validator for the common result contract.

Use ordinary S3 lists/data.frames. Distinguish observed inference, planning, classification and descriptions through a kind and typed components. Do not force sample-size outputs into artificial p-value/CI columns. The core must have useful print methods without the reporting package.

## Result and method requirements

Define exact schemas for stable result/row IDs, method and variant IDs, estimands, contrasts, nulls/alternatives, event/reference orientation, numerical results, intervals, p-value selection, sample/exclusion counts, status/reasons, warnings and provenance. Retain method-specific payloads where a flat table would lose meaning. Lists for multi-df statistics or multi-parameter hypotheses are preferable to mislabeled scalar fields.

Record confidence level, interval method/scale and pointwise/simultaneous status. Preserve raw and adjusted p-values separately and the declared family. Mark absent quantities explicitly; allow documented unbounded intervals. Distinguish numerical bounds from sampling intervals. Keep status and validation evidence separate.

Design generic dispatch so guid can provide a public view of power results and report can consume that view without importing guid. Prefer new planning objects carrying the shared data-view contract at construction, or a documented registered method when guid is loaded. Specify what happens to legacy deserialized objects in a fresh session; do not rely on an accidentally registered method from a development session.

The method catalogue describes execution semantics. The guid planning catalogue describes available planning engines and the relation to those semantics. Do not put guid function names/callbacks into calc or invent a cyclic plugin architecture.

## Dependency and serialization boundary

Implement this new core with standard R. No required broom, tibble, dplyr, rstatix or testthat. Explicitly audit transitive dependencies. Hashing is not a prerequisite for a plain numeric result unless the approved contract requires it. Do not label a missing or non-SHA identifier as SHA-256.

Keep the existing guid serializers/schema dispatch in their current owner until a separately verified migration is ready. They refer to guid-specific power/advisor/report classes. If reusable primitives are extracted, expose only a narrow supported interface and preserve historical bytes/hash behavior. Any indispensable digest/jsonlite exception must be explicit in dependencies.csv, not concealed by an undeclared requireNamespace call.

## Base R tests

Create top-level tests/*.R against the installed package. A small helper subdirectory is acceptable. The package's mandatory tests must run through R CMD check without testthat. Use fixed synthetic fixtures and explicit expected failure classes. This illustrates the intended style; adapt to the actual frozen API and execute it:

```r
stopifnot(inherits(result, "inference_result"))
stopifnot(is.data.frame(result_data(result)))
stopifnot(isTRUE(all.equal(
  result_data(result)$p.value,
  reference$p.value,
  tolerance = 1e-12
)))
tools::assertError(
  method_spec("a_method_that_is_not_registered"),
  classes = "infercalcR_invalid_input"
)
```

Do not use the illustrative tolerance universally. Give each numerical family an appropriate error policy. Add input validation, unknown fields/options, malformed result rows, legitimate NA/Inf cases, status retention, serialization of trusted fixtures and generic-dispatch tests. Keep source objects unmodified.

## Deliverables and gate G03

Deliver code including describe_data(), man pages, base-R tests, `CONTRACT.md`, a versioned contract fixture set, and the exact implemented signatures. Build/install/check the core skeleton and descriptive provider. The coordinator freezes the contract revision, mandatory method list, tests and dependency allowance. Every downstream agent reads that revision. Report conceptual names as proposed until actual implementation and tests exist.

