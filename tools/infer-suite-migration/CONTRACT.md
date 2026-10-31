# Shared inference result contract 1.0.0

`infercalcR` is the only owner of observed numerical engines and this contract. `inferguidR` imports it for guidance/planning; `infereportR` imports it for finalized-result presentation. Calculation never depends on either consumer. Legacy serializers and native report exporters remain in guidance during the transition.

## Implemented public interfaces

- `method_spec(method_id, alternative = "two.sided", conf.level = .95, missing = "error", options = list())`
- `run_inference(spec, ...)`, with method-specific data/design arguments.
- `new_inference_result(data, metadata, kind = "inference", payload = list())`
- `validate_inference(x)`, `as_inference(x, ...)`, `result_data(x, ...)`, `result_metadata(x, ...)`
- `method_catalogue()` describes executable adapters, backend requirements and option names.
- `describe_data(data, variables = setdiff(names(data), group), group = NULL, labels = NULL, units = NULL, types = NULL)`
- `adjust_inference(x, method = "holm", family_id, rows = NULL)` adjusts a complete declared family once.

## Data and metadata

Ordinary S3 lists contain `data` (data.frame), `metadata` (named list), and `payload` (list of finalized method-specific quantities). All result kinds inherit `inference_result`; kind is inference, planning, classification or descriptive. Analytical precision is preserved. Renderers must not refit, test, compute confidence intervals, adjust multiplicity or evaluate stored calls.

Mandatory row columns are row_id, variable_id, component, estimand, contrast, estimate, statistic, df, conf.low, conf.high, conf.level, interval.method, interval.scale, interval.scope, p.value, p.adjusted, p.reported, p.method, family_id, family_size, null.value, alternative, n, n.excluded, status, reason. Extra facet keys and method-specific columns are allowed. Multi-parameter statistics/df remain in payload if flattening would be misleading. All numerical columns use numeric NA for unavailable values. Documented infinite confidence endpoints are valid. Numerical optimization bounds must use separate columns or payload, never conf.low/conf.high.

Metadata requires method_id, variant_id, result_id, contract_version and kind. Constructors supply variant_id and result_id if absent. IDs are stable labels within a result, not global identifiers or cryptographic hashes. Cross-result joins require explicit distinct result_id values and row_id keys. Metadata carries population, exclusions, settings, orientation, warnings, assumptions, provenance and validation evidence when relevant. Source-specific status/fields stay in payload. Common status is one of ok, undefined, not_applicable, unsupported, invalid_input, numerical_failure, infeasible, reference_pending or unresolved. Original unresolved/partial computation is never promoted to ok.

Raw and adjusted p-values are distinct. Missing adjusted p-values may not fall back to raw values when adjustment is requested. Adjustment records complete family membership before display filtering and does not alter intervals. An omnibus test cannot supply pairwise contrasts. Estimates and intervals must name the same estimand and scale. Unknown sample size remains unknown, never inferred from Welch df.

## Execution and adapter rules

Observed method IDs follow the frozen capabilities.csv inventory (69 advisor identities and three callable helper capabilities; 71 mandatory entries). Source templates are development-time evidence only. Production code must use ordinary declared functions and explicit dispatch, never source/evaluate teaching scripts at runtime. Every callable template procedure and helper retains its existing scope and limits. Full finalized payloads are retained when necessary.

Do not forward alternative/conf.level into procedures that do not have those concepts. Mark inapplicable alternatives NA. `missing = "omit"` does not authorize a silent MAR, MCAR, pairwise or complete-case policy: method-specific policies must be explicitly supported and recorded. Long-form pairs are matched by subject ID; wide inputs require unique IDs. Backend choice is deterministic, and no backend is silently substituted.

Core calculations, descriptions and contract use only base/standard R. Legacy canonical JSON, SHA-256 fingerprints and correctly rounded decimal parsing may require optional digest/jsonlite on their explicit feature paths; preserve their algorithms and error if absent. Guidance retains mandatory digest/jsonlite for historical document compatibility. Record every exception. No runtime installs.

## Legacy and planning access

Guidance registers result_data/result_metadata/as_inference methods for its own planning objects. Consumers operate through the generic, with no reverse dependency. In a fresh process a legacy planning JSON document must first be read with inferguidR::ig_from_json(), which loads its registered methods. For a guid-free process, convert it once with as_inference() and retain the common plain result. Trusted RDS is not a safe decoder for untrusted input. Historical JSON bytes/identity stay unchanged; new package provenance is evaluated separately from numerical parity.

## Initial verification

Installed source-tarball tests contract.R and describe.R passed under R 4.6.1 on macOS arm64. R CMD check of the P03 core completed Status: OK. This verifies the contract/descriptive slice, not the later method inventory or suite completion. Evidence: coordinating run contracts/ directory. Downstream additions require their own installed-candidate tests.
