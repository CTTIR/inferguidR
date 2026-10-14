# These tests describe the foundation behaviour for methods without an
# engine; they run against the bare registry (see th_registry_only()).
th_registry_only()

two_groups <- function(...) {
  analysis_spec("t_means_two_independent", "a_priori",
                inputs = list(alpha = 0.05, power = 0.8, d = 0.5),
                options = list(tails = "two"), ...)
}

test_that("a registry-only specification is built and labelled as such", {
  s <- two_groups()
  expect_s3_class(s, "power_spec")
  expect_identical(s$target, "n_total")
  expect_false(s$joint_target)
  expect_identical(s$validation_level, "registry_metadata")
  expect_identical(s$unvalidated_inputs, "d")
  expect_identical(names(s$inputs), c("alpha", "power", "d"))
  expect_match(s$spec_hash, "^[0-9a-f]{64}$")
  expect_output(print(s), "registry metadata only")
})

test_that("bad, missing and unoffered modes give invalid_input", {
  expect_identical(th_condition_code(
    analysis_spec("t_means_two_independent", "simulation", inputs = list(alpha = 0.05),
                  options = list(tails = "two"))), "invalid_input")
  expect_identical(th_condition_code(
    analysis_spec("t_means_two_independent", inputs = list(alpha = 0.05))), "invalid_input")
  expect_identical(th_condition_code(
    analysis_spec("f_generic", "a_priori", inputs = list(alpha = 0.05, power = 0.8))), "invalid_input")
  expect_identical(th_condition_code(analysis_spec("t_welch", "a_priori")), "invalid_input")
  expect_identical(th_condition_code(analysis_spec(c("a", "b"), "a_priori")), "invalid_input")
})

test_that("a missing or wrong solve target gives invalid_input", {
  # The inventory lists alternative sensitivity targets (d or p1) for this test.
  args <- list(method = "t_wilcoxon_mann_whitney", mode = "sensitivity",
               inputs = list(alpha = 0.05, power = 0.8, n_total = 40),
               options = th_minimal_options("t_wilcoxon_mann_whitney", "sensitivity"))
  cond <- tryCatch(do.call(analysis_spec, args), inferguidR_invalid_input = function(e) e)
  expect_match(unlist(cond$status$details$problems), "missing solve target", all = FALSE)
  args$target <- "p1"
  expect_identical(do.call(analysis_spec, args)$target, "p1")
  args$target <- "median"
  expect_identical(th_condition_code(do.call(analysis_spec, args)), "invalid_input")
  expect_identical(th_condition_code(two_groups(target = "alpha")), "invalid_input")
  expect_identical(th_condition_code(
    analysis_spec("t_means_one_sample", "compromise", inputs = list(n_total = 20, error_ratio = 1),
                  options = list(tails = "two"), target = "alpha")), "invalid_input")
  comp <- analysis_spec("t_means_one_sample", "compromise",
                        inputs = list(n_total = 20, error_ratio = 1), options = list(tails = "two"))
  expect_identical(comp$target, c("alpha", "power"))
  expect_true(comp$joint_target)
})

test_that("supplying the solve target as an input is rejected", {
  cond <- tryCatch(
    analysis_spec("t_means_one_sample", "post_hoc",
                  inputs = list(alpha = 0.05, n_total = 20, power = 0.8),
                  options = list(tails = "two")),
    inferguidR_invalid_input = function(e) e)
  expect_match(unlist(cond$status$details$problems), "solve target", all = FALSE)
  expect_identical(th_condition_code(
    analysis_spec("t_means_one_sample", "a_priori",
                  inputs = list(alpha = 0.05, power = 0.8, n_total = 20),
                  options = list(tails = "two"))), "invalid_input")
})

test_that("inconsistent duplicate inputs are rejected, consistent ones collapse", {
  cond <- tryCatch(two_groups(alpha = 0.01), inferguidR_invalid_input = function(e) e)
  problems <- unlist(cond$status$details$problems)
  expect_match(problems, "inconsistent duplicate input 'alpha'", all = FALSE)
  expect_false(any(grepl("missing input 'alpha'", problems)))
  expect_identical(two_groups(alpha = 0.05)$spec_hash, two_groups()$spec_hash)
  # Conflicts involving missing values are reported, not crashed on.
  expect_identical(th_condition_code(two_groups(alpha = NA_real_)), "invalid_input")
  s1 <- analysis_spec("t_means_one_sample", "post_hoc", inputs = list(alpha = 0.05, n_total = 20L),
                      n_total = 20, options = list(tails = "two"))
  expect_identical(s1$inputs$n_total, 20L)
  expect_identical(th_condition_code(
    analysis_spec("t_means_one_sample", "post_hoc",
                  inputs = list(alpha = 0.05, n_total = 20, alpha = 0.1),
                  options = list(tails = "two"))), "invalid_input")
  expect_identical(th_condition_code(
    analysis_spec("t_means_one_sample", "post_hoc", inputs = list(alpha = 0.05, n_total = 20),
                  options = list(tails = "two", tails = "one"))), "invalid_input")
})

test_that("domain violations give invalid_input with every problem listed", {
  cond <- tryCatch(
    analysis_spec("t_means_one_sample", "sensitivity",
                  inputs = list(alpha = 1.5, power = 0, n_total = 10.5),
                  options = list(tails = "two")),
    inferguidR_invalid_input = function(e) e)
  problems <- unlist(cond$status$details$problems)
  expect_length(problems, 3L)
  expect_match(problems, "'alpha' must lie in", all = FALSE)
  expect_match(problems, "'power' must lie in", all = FALSE)
  expect_match(problems, "whole number", all = FALSE)
  bad <- list(list(alpha = NA_real_), list(alpha = "0.05"), list(alpha = Inf),
              list(alpha = c(0.05, 0.01)), list(alpha = list(0.05)))
  for (b in bad) {
    expect_identical(th_condition_code(
      analysis_spec("t_means_one_sample", "post_hoc", inputs = c(b, list(n_total = 20)),
                    options = list(tails = "two"))), "invalid_input", info = deparse(b))
  }
  expect_identical(th_condition_code(
    analysis_spec("t_means_two_independent", "a_priori",
                  inputs = list(alpha = 0.05, power = 0.8, d = 0.5, error_ratio = 2),
                  options = list(tails = "two"))), "invalid_input")
  expect_identical(th_condition_code(
    analysis_spec("t_means_one_sample", "post_hoc", inputs = list(alpha = 0.05, n_total = 20, eff = NaN),
                  options = list(tails = "two"))), "invalid_input")
})

test_that("options are validated against codes, defaults and applicability", {
  base <- list(method = "exact_proportion_one_sample_binomial", mode = "post_hoc",
               inputs = list(alpha = 0.05, n_total = 30, p0 = 0.5, p1 = 0.7))
  # Certified inventory: Tail(s) defaults to "One" in the live application,
  # so the default is applied and written out.
  dflt <- do.call(analysis_spec, base)
  expect_identical(dflt$options, list(tails = "one"))
  expect_identical(dflt$defaults_applied$options, "tails")
  # An applicable option without a documented default must still be supplied
  # (the extension's variant has none).
  missing <- tryCatch(
    analysis_spec("boschloo_proportions_two_independent", "post_hoc", inputs = th_minimal_inputs("post_hoc")),
    inferguidR_invalid_input = function(e) e)
  expect_match(unlist(missing$status$details$problems), "option 'variant' must be supplied", all = FALSE)
  two <- do.call(analysis_spec, c(base, list(options = list(tails = "two"))))
  expect_identical(two$options, list(tails = "two", alpha_balancing = "half_each_side"))
  expect_identical(two$defaults_applied$options, "alpha_balancing")
  one <- do.call(analysis_spec, c(base, list(options = list(tails = "one"))))
  expect_identical(one$options, list(tails = "one"))
  inapplicable <- tryCatch(
    do.call(analysis_spec, c(base, list(options = list(tails = "one", alpha_balancing = "half_each_side")))),
    inferguidR_invalid_input = function(e) e)
  expect_match(unlist(inapplicable$status$details$problems), "does not apply", all = FALSE)
  expect_identical(th_condition_code(do.call(analysis_spec, c(base, list(options = list(tails = "both"))))),
                   "invalid_input")
  expect_identical(th_condition_code(do.call(analysis_spec, c(base, list(options = list(tails = "two", colour = "red"))))),
                   "invalid_input")
  explicit <- do.call(analysis_spec, c(base, list(options = list(tails = "two", alpha_balancing = "half_each_side"))))
  expect_identical(explicit$spec_hash, two$spec_hash)
  expect_identical(explicit$defaults_applied$options, character())
})

test_that("presentation settings and numeric options are not accepted as options", {
  opts <- th_minimal_options("t_wilcoxon_mann_whitney", "post_hoc")
  # 'output moments p1, p2' only shows or hides outputs (live capture).
  expect_identical(th_condition_code(
    analysis_spec("t_wilcoxon_mann_whitney", "post_hoc", inputs = list(alpha = 0.05, n_total = 40),
                  options = c(opts, list(output_moments = "on")))), "invalid_input")
  # The plot choice changes critical values under Lehmann; it is inactive
  # under A.R.E. (D43), even though the option remains result-affecting.
  expect_identical(th_condition_code(
    analysis_spec("t_wilcoxon_mann_whitney", "post_hoc", inputs = list(alpha = 0.05, n_total = 40),
                  options = list(method = "are", distribution_plot = "normal_rank_sum"))), "invalid_input")
  rs <- analysis_spec("t_wilcoxon_mann_whitney", "post_hoc", inputs = list(alpha = 0.05, n_total = 40),
                      options = list(method = "lehmann", distribution_plot = "normal_rank_sum"))
  expect_identical(rs$options$distribution_plot, "normal_rank_sum")
  expect_identical(th_condition_code(
    analysis_spec("t_means_two_independent", "a_priori",
                  inputs = list(alpha = 0.05, power = 0.8, d = 0.5),
                  options = list(tails = "two", allocation_ratio = 1))), "invalid_input")
  expect_identical(th_condition_code(
    analysis_spec("t_means_two_independent", "a_priori",
                  inputs = list(alpha = 0.05, power = 0.8, d = 0.5, tails = "two"))), "invalid_input")
  ok <- analysis_spec("t_means_two_independent", "a_priori",
                      inputs = list(alpha = 0.05, power = 0.8, d = 0.5, allocation_ratio = 2),
                      options = list(tails = "two"))
  expect_identical(ok$inputs$allocation_ratio, 2)
  expect_identical(th_condition_code(
    analysis_spec("t_means_two_independent", "a_priori",
                  inputs = list(alpha = 0.05, power = 0.8, d = 0.5, allocation_ratio = 0),
                  options = list(tails = "two"))), "invalid_input")
})

test_that("every registry method and offered mode accepts a minimal registry-only spec", {
  # The 49 inventory procedures and the extensions (th_extension_ids).
  expect_setequal(.ig_list_method_ids(), c(th_inventory_ids, th_extension_ids))
  for (id in .ig_list_method_ids()) {
    for (mode in .ig_get_metadata(id)$modes) {
      s <- th_minimal_spec(id, mode)
      expect_s3_class(s, "power_spec")
      expect_identical(s$validation_level, "registry_metadata")
    }
  }
})

test_that("specifications validated against a method definition write defaults out", {
  th_local_definition()
  s <- analysis_spec("t_means_two_independent", "post_hoc",
                     inputs = list(alpha = 0.05, n_total = 40, d = 0.5),
                     options = list(tails = "two"))
  expect_identical(s$validation_level, "method_definition")
  expect_identical(s$inputs$allocation_ratio, 1)
  expect_identical(s$inputs$n_total, 40L)
  expect_identical(s$options, list(tails = "two", computation = "approximate", balancing = "half"))
  expect_setequal(s$defaults_applied$inputs, "allocation_ratio")
  expect_setequal(s$defaults_applied$options, c("computation", "balancing"))
  expect_identical(s$unvalidated_inputs, character())
  expect_identical(names(s$inputs), c("d", "alpha", "n_total", "allocation_ratio"))
  exact <- update(s, options = list(computation = "exact"))
  expect_identical(exact$inputs$exact_threshold, 1000L)
  expect_identical(th_condition_code(
    update(s, inputs = list(exact_threshold = 10L))), "invalid_input")
})

test_that("definition-based validation rejects unknown, missing, aliased and inconsistent inputs", {
  th_local_definition()
  base <- list(method = "t_means_two_independent", mode = "post_hoc", options = list(tails = "two"))
  run <- function(inputs) do.call(analysis_spec, c(base, list(inputs = inputs)))
  expect_identical(th_condition_code(run(list(alpha = 0.05, n_total = 40, d = 0.5, sd = 1))),
                   "invalid_input")
  missing <- tryCatch(run(list(alpha = 0.05, n_total = 40)), inferguidR_invalid_input = function(e) e)
  expect_match(unlist(missing$status$details$problems), "missing input 'd'", all = FALSE)
  aliased <- run(list(alpha = 0.05, n_total = 40, effect_size_d = 0.5, ratio = 2))
  expect_identical(aliased$inputs$d, 0.5)
  expect_identical(aliased$inputs$allocation_ratio, 2)
  expect_identical(aliased$spec_hash, run(list(alpha = 0.05, n_total = 40, d = 0.5, allocation_ratio = 2))$spec_hash)
  expect_identical(th_condition_code(run(list(alpha = 0.05, n_total = 40, d = 0.5,
                                              allocation_ratio = 1, ratio = 2))), "invalid_input")
  odd <- tryCatch(run(list(alpha = 0.05, n_total = 41, d = 0.5)), inferguidR_invalid_input = function(e) e)
  expect_match(unlist(odd$status$details$problems), "must be even", all = FALSE)
  expect_identical(th_condition_code(
    analysis_spec("t_means_two_independent", "a_priori",
                  inputs = list(alpha = 0.05, power = 0.8, d = 0.5, error_ratio = 1),
                  options = list(tails = "two"))), "invalid_input")
})

test_that("the same question always has the same hash", {
  a <- analysis_spec("t_means_two_independent", "a_priori",
                     inputs = list(d = 0.5, power = 0.8, alpha = 0.05), options = c(tails = "two"))
  b <- analysis_spec("t_means_two_independent", "a_priori", alpha = 0.05, power = 0.8, d = 0.5,
                     options = list(tails = "two"))
  expect_identical(a, b)
  expect_identical(ig_hash(a), two_groups()$spec_hash)
  expect_identical(ig_from_json(ig_to_json(a))$spec_hash, a$spec_hash)
})

test_that("every result-affecting change produces a new hash", {
  h <- two_groups()$spec_hash
  variants <- list(
    alpha_ulp = analysis_spec("t_means_two_independent", "a_priori",
                              inputs = list(alpha = 0.05 * (1 + .Machine$double.eps), power = 0.8, d = 0.5),
                              options = list(tails = "two")),
    tails = analysis_spec("t_means_two_independent", "a_priori",
                          inputs = list(alpha = 0.05, power = 0.8, d = 0.5), options = list(tails = "one")),
    effect = analysis_spec("t_means_two_independent", "a_priori",
                           inputs = list(alpha = 0.05, power = 0.8, d = 0.51), options = list(tails = "two")),
    extra_input = analysis_spec("t_means_two_independent", "a_priori",
                                inputs = list(alpha = 0.05, power = 0.8, d = 0.5, allocation_ratio = 2),
                                options = list(tails = "two")),
    mode = analysis_spec("t_means_two_independent", "sensitivity",
                         inputs = list(alpha = 0.05, power = 0.8, n_total = 128),
                         options = list(tails = "two")),
    method = analysis_spec("t_means_matched_pairs", "a_priori",
                           inputs = list(alpha = 0.05, power = 0.8, d = 0.5), options = list(tails = "two"))
  )
  hashes <- vapply(variants, `[[`, "", "spec_hash")
  expect_false(any(hashes == h))
  expect_false(anyDuplicated(hashes) > 0)
  wmw <- function(target) {
    analysis_spec("t_wilcoxon_mann_whitney", "sensitivity",
                  inputs = list(alpha = 0.05, power = 0.8, n_total = 40),
                  options = th_minimal_options("t_wilcoxon_mann_whitney", "sensitivity"),
                  target = target)
  }
  expect_false(identical(wmw("d")$spec_hash, wmw("p1")$spec_hash))
})

test_that("specifications are immutable", {
  s <- two_groups()
  expect_identical(th_condition_code(s$inputs <- list()), "invalid_input")
  expect_identical(th_condition_code(s[["mode"]] <- "post_hoc"), "invalid_input")
  expect_identical(th_condition_code(s["mode"] <- "post_hoc"), "invalid_input")
  expect_identical(s, two_groups())
})

test_that("update() derives new validated specifications", {
  s <- two_groups()
  s2 <- update(s, inputs = list(alpha = 0.01))
  expect_identical(s2$inputs$alpha, 0.01)
  expect_false(identical(s2$spec_hash, s$spec_hash))
  expect_identical(s$inputs$alpha, 0.05)
  ph <- update(s, mode = "post_hoc", n_total = 128)
  expect_identical(ph$target, "power")
  expect_false("power" %in% names(ph$inputs))
  expect_identical(update(ph, n_total = NULL, power = 0.8, mode = "a_priori")$spec_hash, s$spec_hash)
  expect_identical(th_condition_code(update(s, inputs = list(alpha = 2))), "invalid_input")
  expect_identical(th_condition_code(update(s, alpha = 0.01, inputs = list(alpha = 0.02))),
                   "invalid_input")
  expect_identical(update(s, options = list(tails = "one"))$options$tails, "one")
})

test_that("specifications round-trip through JSON exactly and tampering is detected", {
  s <- analysis_spec("exact_proportion_one_sample_binomial", "post_hoc",
                     inputs = list(alpha = 0.05, n_total = 30, p0 = 0.5, p1 = 0.7),
                     options = list(tails = "two"))
  expect_identical(ig_from_json(ig_to_json(s)), s)
  th_local_definition()
  d <- analysis_spec("t_means_two_independent", "post_hoc",
                     inputs = list(alpha = 0.05, n_total = 40, d = 0.5), options = list(tails = "two"))
  expect_identical(ig_from_json(ig_to_json(d)), d)
  lst <- .ig_spec_to_list(d)
  lst$inputs$alpha <- 0.1
  cond <- tryCatch(ig_from_json(ig_to_json(lst)), inferguidR_invalid_input = function(e) e)
  expect_match(cond$status$reason, "spec_hash does not match")
  lst <- .ig_spec_to_list(d)
  lst$inputs$alpha <- 2
  expect_identical(th_condition_code(ig_from_json(ig_to_json(lst))), "invalid_input")
  txt <- sub('"alpha":0.05', '"alpha":0.05,"alpha":0.01', as.character(ig_to_json(d)), fixed = TRUE)
  expect_identical(th_condition_code(ig_from_json(txt)), "invalid_input")
})

test_that("a reproducible R call recreates the same specification", {
  th_local_definition()
  s <- analysis_spec("t_means_two_independent", "post_hoc",
                     inputs = list(alpha = 0.05, n_total = 40, d = 0.1 + 0.2),
                     options = list(tails = "two"))
  call_txt <- .ig_spec_call(s)
  expect_match(call_txt, "d = 0.30000000000000004", fixed = TRUE)
  # The call text is generated by the package from a validated spec (never
  # from imported data); evaluating it here is the point of the test.
  rebuilt <- eval(parse(text = sub("inferguidR::", "", call_txt, fixed = TRUE)))
  expect_identical(rebuilt$spec_hash, s$spec_hash)
})
