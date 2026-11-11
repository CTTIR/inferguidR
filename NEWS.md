# inferguidR 0.1.0

* `analysis_spec()` gains `engine = "reference"`, which computes a procedure
  with an established R package that implements the same definition instead
  of inferguidR's own forward model; `engine = "base"` stays the default and
  base specifications keep their structure and hash. 35 procedures (t, F, z,
  exact, chi-squared and Wilcoxon families) have reference engines (stats,
  pwr, pwrss, WebPower, MESS, powerMediation, Exact), benchmarked against the
  base engines on the originals' help-page examples and G*Power manual
  examples (`inst/validation/engine-benchmarks.rds`, 102 cases agreeing within
  1e-8). A missing, older or inapplicable original gives an "unsupported"
  result, never a fallback. `list_methods()` lists each procedure's engines,
  and `inst/validation/power-provenance.csv` the supported combinations.
* The README and every vignette now acknowledge and cite DescrTab2
  (Institute of Medical Biometry, Heidelberg University), which deeply
  inspired the suite, the Methodenberatung pages of the University of
  Zurich and G*Power; the README also cites R and the packages used.
* Prepared for a first CRAN submission.
* Observed calculations and the shared numerical engines (Boschloo,
  classifier, distribution and calculator code) now live in infercalcR, which
  inferguidR imports; finalized-result tables, text and plot annotations live
  in infereportR. inferguidR keeps guidance, prospective planning and the
  historical JSON and report formats. All previously exported functions remain
  available; status helpers and common result accessors are re-exported from
  infercalcR.
* The maintainer address is now raban.heller@uni-ulm.de, and `URL` also
  lists the documentation site.
* DESCRIPTION declares `Language: en-GB`; `inst/WORDLIST` records the
  technical terms and the German glossary entries used in the
  documentation.
* Nine exported functions without examples gain them, among them
  `advisor_power_spec()`, `advisor_r_code()`, `ig_method()`,
  `create_report_snapshot()` and `export_report()`.
* The GitHub Actions check matrix adds R oldrel-1 on Linux.
* The 23 slowest test files, which together took about ten minutes,
  skip on CRAN. They still run with `NOT_CRAN=true`, which the Linux
  R-release CI job now sets.

## Earlier development changes

* Add public common planning-result accessors and explicit observed-method
  relations to planning engines.
* Standardize the named v5 package logo across README, documentation pages
  and available application headers, with matching light and dark artwork.


Development version; numerical and release qualification remain incomplete.

* Numbers typed into the calculator and the workbench, and numbers read back
  from canonical JSON, are converted with correct rounding on every platform.
  R's own conversion missed the nearest double for many inputs with 16 or more
  significant digits or a large exponent where long double has no extra
  precision (macOS on Apple silicon), which changed calculator results and
  made retained report hashes fail there.
* The unconstrained theoretical Youden optimum reports the analytic stationary
  point instead of a numerical estimate next to it.
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
