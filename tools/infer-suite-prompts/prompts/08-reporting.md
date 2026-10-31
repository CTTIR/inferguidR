# P08 Tables, descriptive reporting and manuscript text

Role: infereportR core owner. Prerequisite: G03; integration against completed P05/P06/P07 results is required before G10/G11.

## Objective

Implement a practical reporting package that uses the shared result contract and remains useful with only infercalcR and standard R installed.

## Required core

- `report_table()`: ordinary data.frame by default, with raw analytical columns kept distinct from formatted display columns and metadata accessible through a documented interface.
- `report_text()`: deterministic character output or inline values from finalized numerical results.
- `report_descriptive()`: continuous/categorical summaries, labels/units, group denominators and missing counts; optional supplied inference columns.
- Shared estimate/interval/statistic/p-value formatting helpers used by tables, text and the ggplot layer.
- Clear access to method, population, assumptions, status, missing quantities and provenance.

Use base R data manipulation and formatting. A small structured reporting container is allowed only if necessary to preserve metadata; do not introduce a competing computational result schema. Every result type uses the approved accessor contract.

For the minimal raw-data descriptive route, call calc's describe_data() implemented in P03; verify that prerequisite rather than leaving its ownership implicit. An explicitly selected optional backend may supply another route, but cannot replace the required standard-R route. Do not duplicate an inferential engine. Baseline tables have no automatic p-value column by default. The descriptive population may differ from an analysis subset only when explicitly labelled. Joining supplied inference uses stable variable/result IDs and checks missing/duplicate mappings; never match only on a pretty label.

## Statistical boundary

Rendering finalized results must never fit models, select tests, recompute intervals, run p.adjust, resample or optimize cutoffs. An explicit analysis convenience entry point must delegate to calc before rendering and expose the reusable result. Do not let rendering backend defaults change the test.

Preserve estimand, contrast orientation, method/covariance/interval definition, confidence level, raw/adjusted p-value choice and family scope. Omnibus results remain omnibus. Undefined estimates and p-values remain explicit. Support valid one-sided/unbounded intervals. Use the stored Boschloo bound/decision policy; never call its numerical p_upper a sampling CI or infer a significance star contrary to an unresolved decision.

Text templates describe the supplied estimate and uncertainty neutrally. They do not invent clinical/causal meaning, effect-size adjectives, equivalence from a nonsignificant test, or statistical validation. Scientific interpretation remains separately authored. Locale may affect display only, not stored numeric values or method identity.

## Optional table and document integrations

Support the requested DescrTab2 workflow through a narrow explicit adapter: import a precomputed native result and/or render a native object while preserving what its methods actually mean. Do not let descr() automatically retest an analysis during presentation. `print_p = FALSE` only hides presentation and is not proof that testing was disabled. Use exported, version-verified APIs; do not base package interoperability on third-party `:::` internals. If safe precomputed injection is not publicly supported, support importing/rendering an already computed object and report that boundary clearly.

DescrTab2 was archived on CRAN on 2025-10-18 at the design review. Recheck at implementation time. Pin and document the optional source/version; do not assume an ordinary install.packages call can fetch it. Keep unavailable archived packages out of mandatory Imports; determine CRAN packaging implications explicitly before a future submission.

gtsummary/gt/flextable/officer may be optional formatting/export backends when they add requested value. They are not the default mandatory engine under the user's updated dependency policy. Use precomputed-statistic hooks rather than add_p defaults. Ship plain-data/plain-text output even when none is installed. Do not support many alternative backends merely for breadth.

## Verification

Write base-R tests for full-precision numeric preservation, source-object immutability, shared display values, tiny/zero/NA p-value semantics, raw-versus-adjusted selection, missing adjustment errors, filtered families, different confidence levels, omnibus results, directional/ratio/rank estimands, status handling, label joins, denominators and locale settings.

Verify result-only reporting succeeds without raw data. Plain tables/text work without ggplot2, gtsummary, DescrTab2 or testthat. Missing explicitly requested optional features give actionable errors. Integration claims need a separate installed-backend lane; no result can be marked verified solely because code is guarded by requireNamespace.

## Gate G08

Core tables/text/descriptions are implemented and verified in an isolated minimal library. The requested DescrTab2 import/render adapter has an executed check using a pinned installed version for its declared scope; runtime installation remains optional, but this support claim needs actual validation. If the pinned backend cannot be made available, record the feature/gate as BLOCKED and continue independent work without claiming the complete reporting scope. Other optional adapters are documented and validated whenever advertised. Numerical inference remains centralized, shared formatter contracts are ready for P09, and data/metadata remain coherent across output types.

