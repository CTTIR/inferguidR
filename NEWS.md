# inferguidR 0.0.0.9000

Development version; numerical and release qualification remain incomplete.

* Exact unconditional pooled, continuity-corrected pooled and Mantel-Haenszel
  calculations now preserve rational equality classes at discrete critical
  boundaries. Some criterion results can change where the previous calculation
  split mathematically tied tables through floating-point rounding. Exact
  compromise remains guarded pending an attained-size solver contract.
* `advisor_power_spec()` creates explicit same-test planning specifications for
  mapped Boschloo recommendations. Variant, direction, probabilities,
  allocation policy and numerical controls are retained. Rules version
  `2026-10-26.2` requires explicit replay of older saved sessions; existing
  recommendation documents remain readable.
* Classifier probability diagnostics now expose calibration and net-benefit
  panels, explicit decision probabilities and partial-AUC settings. Changed or
  invalid evaluation settings make old results ineligible for report export.
* The planner retains method-specific input drafts when dynamic controls are
  rebuilt, including invalid and explicitly empty entries with their input
  locale. Methods without effect helpers no longer raise an indexing error
  during a method transition. Responsive plot labels and legends remain visible
  at narrow browser widths.
* `create_report_snapshot(include_probability_summaries = TRUE)` can retain
  explicitly selected aggregate probability diagnostics in report schema 0.4.
  Default exports continue to omit these summaries and raw scored rows.
  Boschloo power reports now retain both conditional count lattices, including
  their finite support, rejection membership and probability masses. Import
  and export use the frozen grids without rerunning the region solver.
* Report exports stage and validate a complete sibling file before replacing
  the destination with a checked rename. Failed rendering or commit leaves an
  existing report intact; this follows filesystem rename semantics rather than
  promising crash durability or isolation from concurrent writers.
* Canonical JSON string handling avoids repeated scans while preserving the
  serialized bytes and hashes. Report normalization explicitly preserves NULL
  values across supported R versions.
* Public numerical campaigns, method/mode/option inventories, installed release
  preflights and source/resource audits provide reproducible scoped evidence.
  A bilingual statistical glossary explains interpretation and uncertainty.
* Platform workflows use a verified shared source archive, fresh candidate
  libraries and a separate minimum-R dependency profile. Passing these checks
  does not establish complete external-reference parity or release approval.
