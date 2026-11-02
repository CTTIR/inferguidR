# P04 Existing exact engines and shared kernels

Role: numerical extraction owner. Prerequisite: G03. Own assigned calc numerical files; propose guid changes for P07.

## Objective

Extract the existing observed-test machinery and reusable statistical kernels without changing their scientific behavior. Begin with the complete Boschloo slice, then finish the other shared-kernel entries assigned by ownership.csv. This stage is extraction and public-interface construction, not numerical simplification.

## Boschloo first

Read current boschloo-test, ordering, nuisance, region, control, design and method code and their tests. Identify numerical helpers versus planning loops versus rendering. Move observed execution, method variants, ordering/tie rules, nuisance optimization, numerical controls and rejection-region semantics into calc. Keep planning orchestration in guid.

Provide public core capabilities for the region/quantities that planning needs. Expose a constrained supported API; do not make guid call private calc functions. Preserve caching as an implementation detail. If an old cache requires a contributed hash package, consider a correct base-R cache key or no cache only when numerical semantics and practical performance are preserved; do not replace cryptographic document identity.

Preserve central, minlike, less and greater definitions; group/event direction; exact equality classes and tie policy; numerical tolerances; attained-size information; result diagnostics; p_value and its rigorous numerical bound; and reject/do_not_reject/unresolved decisions. p_upper is not conf.high. No invented effect interval.

The guid power integration must use the same rejection rule. Keep the first qualifying allocation and stability over a stated finite window distinct. Do not replace the discrete sample-size search by an assumed monotone root search. Do not claim unbounded stability from a finite search.

## Other shared mathematics

Follow the symbol-level map. Move reusable distributions, observed decision rules and numerical primitives only when needed across packages. Keep method-specific power adapters, allocation rules and inverse planning in guid. Preserve precision fixes and platform-sensitive numeric parsing/rounding. Do not delete unusual-looking numerical code before understanding its recorded failure case.

Separate source ownership from compatibility. There must be one active numerical implementation; guid compatibility wrappers can reconstruct old classes/specifications/provenance from the core output. Provide the exact wrapper requirements to P07.

## Verification

Use base-R tests for fixed known cases, recorded baseline cases, small-table enumeration, group/event reversal, variants, boundary alpha, ties, numeric bounds and unresolved decisions. Compare baseline/candidate in separate processes and libraries. Add independent small-support enumeration when feasible and check rejection membership agrees with the observed decision rule under the same numerical settings.

Separate old-document identity checks from new-result parity. A new producing package changes provenance legitimately; the mathematical results and documented decisions must remain within the declared numerical contract. Historical documents keep their original identity.

Profile representative small and moderate cases to detect a material regression from dependency removal; use base system.time/Rprof where useful. No broad benchmark campaign beyond resolving concrete extraction risk.

## Gate G04

All assigned kernels have one canonical owner, public consumer interfaces, meaningful base-R parity/semantic checks, documented tolerances and an explicit compatibility handoff. Boschloo analysis and planning inputs can be tied to the same variant/rejection rule. Any unresolved numerical discrepancy blocks this gate and dependent integration.

