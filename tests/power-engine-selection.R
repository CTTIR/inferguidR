library(inferguidR)
# Engine selection for power procedures: base is the default and keeps every
# historical identity; a reference engine is recorded and never substituted.
fixture <- function(name) file.path("power-engine", name)
spec_args <- list("t_means_one_sample", "post_hoc", inputs = list(alpha = .05, d = .5, n_total = 20),
  options = list(tails = "two"))
base <- do.call(analysis_spec, spec_args)

# Identities recorded on main b900523, before engines existed. The spec hash
# covers inputs only and is the same on every platform; the result hash covers
# computed doubles, whose last bits differ between platforms, so the saved
# specification must reproduce this platform's base result, and the power the
# value recorded on Linux within rounding.
base_result <- calculate_power(base)
stopifnot(identical(base$spec_hash, "467119d8f8eb4143689192445dafda5d7a633dcf53fe2195d3b3e84c63854b58"),
  is.null(base$engine), abs(base_result$achieved$power - 0.56450441843908394) <= 1e-12)
saved <- ig_from_json(readLines(fixture("spec-t-one-sample-b900523.json")))
stopifnot(identical(saved$spec_hash, base$spec_hash),
  identical(calculate_power(saved)$result_hash, base_result$result_hash),
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

# Option combinations without a matching original are unsupported, not replaced.
if (requireNamespace("WebPower", quietly = TRUE)) {
  vc <- calculate_power(analysis_spec("z_logistic_regression", "post_hoc", inputs = list(odds_ratio = 1.5, p1 = 0.5,
    alpha = 0.05, n_total = 300, r2_other_x = 0, x_mu = 0, x_sigma = 1), options = list(tails = "two",
    computation = "large_sample", test_procedure = "demidenko", variance_correction = "on", x_distribution = "normal"),
    engine = "reference"))
  stopifnot(identical(vc$status$code, "unsupported"))
}
corrected <- calculate_power(analysis_spec("z_proportions_two_independent", "post_hoc", inputs = list(p1 = 0.5,
  p2 = 0.75, alpha = 0.05, n1 = 50, n2 = 50), options = list(tails = "two", continuity_correction = "on"),
  engine = "reference"))
stopifnot(identical(corrected$status$code, "unsupported"))

# Final-review fixes. The recorded result hash is checked where it was captured.
if (identical(Sys.info()[["sysname"]], "Linux") && identical(R.version$arch, "x86_64"))
  stopifnot(identical(base_result$result_hash, "e9297e4fc51fd4d664ebaf6e1b52456b10d4853802bf1fcf95153f27144630a3"))
ref_power <- function(id, inputs, options, mode = "post_hoc") {
  spec <- analysis_spec(id, mode, inputs = inputs, options = options, engine = "reference")
  if (identical(mode, "a_priori")) solve_design(spec) else calculate_power(spec)
}
if (all(vapply(c("pwrss", "MESS"), requireNamespace, NA, quietly = TRUE))) {
  # Binomial balancing rules other than half each side have no matching original.
  minor <- ref_power("exact_proportion_one_sample_binomial", list(g = -0.1, p0 = 0.2, alpha = 0.05, n_total = 23),
    list(tails = "two", alpha_balancing = "minor_tail_half_rest_major"))
  stopifnot(identical(minor$status$code, "unsupported"))
  # The discrete reference engines do not report an attained alpha.
  binomial <- ref_power("exact_proportion_one_sample_binomial", list(g = 0.25, p0 = 0.5, alpha = 0.05, n_total = 50),
    list(tails = "two"))
  stopifnot(identical(binomial$status$code, "ok"), is.na(binomial$achieved$alpha_attained))
  # McNemar with an odds ratio below 1 (the two-tailed test is symmetric in OR and 1 / OR).
  mcnemar_args <- list(list(or = 0.5, pd = 0.375, alpha = 0.05, n_total = 120), list(tails = "two"))
  base_mc <- calculate_power(analysis_spec("exact_proportions_mcnemar", "post_hoc", inputs = mcnemar_args[[1]],
    options = mcnemar_args[[2]]))$achieved$power
  ref_mc <- do.call(ref_power, c(list("exact_proportions_mcnemar"), mcnemar_args))
  stopifnot(identical(ref_mc$status$code, "ok"), abs(ref_mc$achieved$power - base_mc) <= 1e-8 * base_mc)
}
if (requireNamespace("pwrss", quietly = TRUE)) {
  # A smallest design below the original's domain is not reported as minimal.
  floor_case <- ref_power("t_correlation_point_biserial", list(rho = 0.9, alpha = 0.05, power = 0.8),
    list(tails = "one"), mode = "a_priori")
  stopifnot(identical(floor_case$status$code, "unsupported"))
  # The reference engines solve post hoc and a priori designs only.
  sensitivity <- calculate_power(analysis_spec("t_correlation_point_biserial", "sensitivity",
    inputs = list(alpha = 0.05, power = 0.8, n_total = 50), options = list(tails = "two"), engine = "reference"))
  stopifnot(identical(sensitivity$status$code, "unsupported"))
}
# The reproducible call keeps the engine.
call_spec <- analysis_spec("t_means_one_sample", "post_hoc", inputs = list(alpha = .05, d = .5, n_total = 20),
  options = list(tails = "two"), engine = "reference")
call_result <- calculate_power(call_spec)
stopifnot(grepl("engine = \"reference\"", call_result$r_call, fixed = TRUE),
  identical(eval(parse(text = call_result$r_call))$spec_hash, call_spec$spec_hash),
  any(grepl("reference", format(call_spec), fixed = TRUE)))

listing <- list_methods()
stopifnot(all(c("engines", "reference") %in% names(listing)),
  all(listing$engines %in% c("base", "base+reference")))
cat("PASS power engine selection, historical identities and unsupported originals\n")
