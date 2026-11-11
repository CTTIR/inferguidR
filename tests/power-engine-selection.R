library(inferguidR)
# Engine selection for power procedures: base is the default and keeps every
# historical identity; a reference engine is recorded and never substituted.
fixture <- function(name) file.path("power-engine", name)
spec_args <- list("t_means_one_sample", "post_hoc", inputs = list(alpha = .05, d = .5, n_total = 20),
  options = list(tails = "two"))
base <- do.call(analysis_spec, spec_args)

# Identities recorded on main b900523, before engines existed.
stopifnot(identical(base$spec_hash, "467119d8f8eb4143689192445dafda5d7a633dcf53fe2195d3b3e84c63854b58"),
  is.null(base$engine),
  identical(calculate_power(base)$result_hash, "e9297e4fc51fd4d664ebaf6e1b52456b10d4853802bf1fcf95153f27144630a3"))
saved <- ig_from_json(readLines(fixture("spec-t-one-sample-b900523.json")))
stopifnot(identical(saved$spec_hash, base$spec_hash),
  identical(calculate_power(saved)$result_hash, "e9297e4fc51fd4d664ebaf6e1b52456b10d4853802bf1fcf95153f27144630a3"),
  identical(as.character(ig_to_json(base)), paste(readLines(fixture("spec-t-one-sample-b900523.json")), collapse = "\n")))

# A reference specification is a different, recorded analysis.
reference <- do.call(analysis_spec, c(spec_args, list(engine = "reference")))
stopifnot(identical(reference$engine, "reference"), !identical(reference$spec_hash, base$spec_hash))
tools::assertError(do.call(analysis_spec, c(spec_args, list(engine = "fast"))))

# Procedures without an R equivalent have no reference engine.
manova <- analysis_spec("f_manova_global", "post_hoc", inputs = list(alpha = .05, f2 = .1, n_total = 60,
  n_groups = 3, n_responses = 2), engine = "reference")
manova_result <- calculate_power(manova)
stopifnot(identical(manova_result$status$code, "unsupported"),
  grepl("no reference engine", manova_result$status$reason, fixed = TRUE))

# An original that is missing or older than required is reported, never replaced.
too_old <- list(t_means_one_sample = list(source = "test", branches = list(list(
  when = function(options, inputs) TRUE, package = "stats", functions = "power.t.test", min_version = "999.0",
  power = function(p) stop("must not run"), generic_search = TRUE, call = "stats::power.t.test"))))
resolved <- inferguidR:::.ig_engine_resolve(reference, registry = too_old)
stopifnot(inherits(resolved, "ig_status"), identical(resolved$code, "unsupported"),
  grepl("999.0", resolved$reason, fixed = TRUE))

# An unbalanced a priori search brackets designs the two-sample original cannot
# evaluate; the reference engine reports this instead of mixing engines.
if (requireNamespace("pwr", quietly = TRUE)) {
  unbalanced <- solve_design(analysis_spec("t_means_two_independent", "a_priori",
    inputs = list(d = 0.5, alpha = 0.05, power = 0.95, allocation_ratio = 2), options = list(tails = "one"),
    engine = "reference"))
  stopifnot(identical(unbalanced$status$code, "unsupported"))
  balanced <- solve_design(analysis_spec("t_means_two_independent", "a_priori",
    inputs = list(d = 0.5, alpha = 0.05, power = 0.95, allocation_ratio = 1), options = list(tails = "one"),
    engine = "reference"))
  stopifnot(identical(balanced$status$code, "ok"), identical(balanced$provenance$engine$engine, "reference"))
}
# A one-tailed single-coefficient test has no packaged original.
if (requireNamespace("WebPower", quietly = TRUE)) {
  one_tailed <- calculate_power(analysis_spec("t_mrc_single_coefficient", "post_hoc",
    inputs = list(f2 = 0.1, alpha = 0.05, n_total = 100, n_predictors = 3), options = list(tails = "one"),
    engine = "reference"))
  stopifnot(identical(one_tailed$status$code, "unsupported"))
}

listing <- list_methods()
stopifnot(all(c("engines", "reference") %in% names(listing)),
  all(listing$engines %in% c("base", "base+reference")))
cat("PASS power engine selection, historical identities and unsupported originals\n")
