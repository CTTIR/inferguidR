# Boschloo API contract: explicit variants, input validation, typed statuses,
# immutable results, labels and provenance, envelopes and JSON round trips.

cond_class <- function(expr) {
  tryCatch({
    force(expr)
    "no condition"
  }, error = function(e) class(e)[1])
}

test_that("the variant is never implied", {
  expect_identical(cond_class(boschloo_test(1, 5, 2, 5)), "inferguidR_invalid_input")
  expect_identical(cond_class(boschloo_region(5, 5, 0.05)), "inferguidR_invalid_input")
  expect_identical(cond_class(boschloo_power(5, 5, 0.1, 0.5)), "inferguidR_invalid_input")
  expect_identical(cond_class(boschloo_sample_size(0.1, 0.5)), "inferguidR_invalid_input")
  expect_identical(cond_class(boschloo_test(1, 5, 2, 5, variant = "two.sided")),
                   "inferguidR_invalid_input")
})

test_that("malformed inputs signal invalid_input with all problems listed", {
  e <- tryCatch(boschloo_test(6, 5, -1, 5, variant = "less", groups = c("a", "a")),
                error = function(e) e)
  expect_s3_class(e, "inferguidR_invalid_input")
  expect_length(e$status$details$problems, 3L)
  expect_identical(cond_class(boschloo_test(1.5, 5, 2, 5, variant = "less")), "inferguidR_invalid_input")
  expect_identical(cond_class(boschloo_test(1, 0, 0, 5, variant = "less")), "inferguidR_invalid_input")
  expect_identical(cond_class(boschloo_test(1, 5, 2, 5, variant = "less", alpha = 1)),
                   "inferguidR_invalid_input")
  expect_identical(cond_class(boschloo_power(5, 5, 1.2, 0.5, variant = "less")), "inferguidR_invalid_input")
  expect_identical(cond_class(boschloo_control(tie_tolerance = -1)), "inferguidR_invalid_input")
  expect_identical(cond_class(boschloo_control(grid_points = 2)), "inferguidR_invalid_input")
  expect_identical(cond_class(boschloo_test(1, 5, 2, 5, variant = "less", control = list())),
                   "inferguidR_invalid_input")
  expect_identical(cond_class(boschloo_region(2500, 2500, 0.05, "minlike")), "inferguidR_unsupported")
})

test_that("results record the variant, orientation, settings and scope", {
  r <- boschloo_test(3, 12, 8, 9, variant = "less", alpha = 0.05,
                     groups = c("placebo", "drug"), event = "response")
  expect_identical(r$status$code, "ok")
  expect_identical(r$variant$variant_id, "boschloo_lower")
  expect_identical(r$variant$alternative, "p1 < p2")
  expect_identical(r$variant$null, "p1 = p2")
  expect_identical(r$variant$validation_status, "independently_checked")
  expect_identical(r$variant$planning_compatibility, "independently_validated_extension")
  expect_match(r$variant$gpower_relation, "gpower_comparator_different_method")
  expect_match(r$data$orientation, "group 1 = placebo")
  expect_identical(r$settings$tie_tolerance, 1e-7)
  expect_match(r$scope, "^observed data")
  expect_identical(r$decision, "reject")
  expect_false(any(grepl("version", names(r$provenance))))
  expect_output(print(r), "Boschloo test, one-sided")
  reg <- boschloo_region(8, 8, 0.05, "central")
  expect_match(reg$scope, "hypothetical future study")
  expect_output(print(reg), "attained size")
  expect_output(print(boschloo_power(8, 8, 0.1, 0.7, variant = "central")), "power")
  expect_output(print(boschloo_control()), "tie tolerance")
})

test_that("Boschloo results are immutable and reproducible from their R call", {
  # The evaluated strings are the reproducible calls that the package itself
  # generated a line earlier (no imported text is evaluated); evaluating them
  # is the point of the check.
  r <- boschloo_test(3, 12, 8, 9, variant = "central", alpha = 0.05, groups = c("A", "B"), event = "cure")
  expect_identical(cond_class(r$p_value <- 0), "inferguidR_invalid_input")
  expect_identical(eval(parse(text = r$r_call)), r)
  pw <- boschloo_power(10, 14, 0.2, 0.75, alpha = 0.05, variant = "minlike",
                       control = boschloo_control(grid_points = 501L))
  expect_match(pw$r_call, "grid_points = 501L")
  expect_identical(eval(parse(text = pw$r_call)), pw)
  expect_identical(boschloo_region(7, 9, 0.05, "less")$region_hash,
                   boschloo_region(7, 9, 0.05, "less", control = boschloo_control(cache = FALSE))$region_hash)
})

test_that("observed results travel in an analysis envelope and survive JSON", {
  r <- boschloo_test(3, 15, 9, 15, variant = "minlike", alpha = 0.05, groups = c("c", "t"), event = "y")
  env <- analysis_envelope(r$spec, r)
  expect_identical(env$analysis_kind, "observed_test")
  back <- ig_from_json(ig_to_json(env))
  expect_identical(back$analysis_hash, env$analysis_hash)
  expect_identical(back$result$p_value, r$p_value)
  r2 <- boschloo_test(0, 10, 0, 10, variant = "central")
  env2 <- analysis_envelope(r2$spec, r2)
  expect_identical(ig_from_json(ig_to_json(env2))$analysis_hash, env2$analysis_hash)
  # The result belongs to its own specification only.
  expect_identical(cond_class(analysis_envelope(r2$spec, r)), "inferguidR_invalid_input")
})

test_that("failed numerical certification yields NA with reasons, never a number", {
  bad <- boschloo_control(grid_points = 3L, max_levels = 1L)
  r <- boschloo_test(10, 35, 24, 61, variant = "minlike", alpha = 0.05, control = bad)
  expect_identical(r$status$code, "numerical_failure")
  expect_true(is.na(r$p_value) && is.na(r$p_upper) && is.na(r$fisher_p_value))
  expect_null(r$decision)
  expect_true(all(c("p_value", "p_upper", "fisher_p_value") %in% names(r$undefined)))
  expect_true(is.na(r$sides$minlike$ordering_statistic))
  reg <- boschloo_region(30, 30, 0.05, "minlike", control = bad)
  expect_identical(reg$status$code, "numerical_failure")
  expect_true(is.na(reg$size$value) && is.na(reg$n_rejected))
  expect_identical(nrow(reg$rejected), 0L)
  pw <- boschloo_power(30, 30, 0.2, 0.5, variant = "minlike", control = bad)
  expect_true(is.na(pw$power) && is.na(pw$beta))
  expect_true(nzchar(pw$undefined$power))
})

test_that("variant metadata matches the advisor identifiers and vocabularies", {
  v <- boschloo_variants()
  expect_identical(v$variant, c("minlike", "central", "less", "greater"))
  expect_identical(v$variant_id, c("boschloo_minlike", "boschloo_central", "boschloo_lower", "boschloo_upper"))
  expect_identical(v$advisor_variant_id, c("boschloo_two_sided_minlike", "boschloo_two_sided_central",
                                           "boschloo_one_sided", "boschloo_one_sided"))
  expect_identical(v$alternative, c("p1 != p2", "p1 != p2", "p1 < p2", "p1 > p2"))
  expect_true(all(v$planning_compatibility %in% c("exact_gpower_match", "gpower_comparator_different_method",
                                                  "independently_validated_extension", "not_yet_validated",
                                                  "not_applicable")))
  expect_false(any(grepl("gpower_validated|exact_gpower_match", unlist(v))))
  expect_match(v$exact_call[1], "tsmethod = \"square\"")
  expect_match(attr(v, "gpower_note"), "no Boschloo procedure")
})

test_that("baseline-plus-effect input is converted with feasibility checks", {
  pr <- boschloo_probabilities(0.2, difference = 0.15)
  expect_equal(pr$p2, 0.35)
  expect_equal(boschloo_probabilities(0.2, ratio = 2)$p2, 0.4)
  expect_equal(boschloo_probabilities(0.2, odds_ratio = 4)$p2, 0.5)
  expect_equal(boschloo_probabilities(0.2, odds_ratio = 4)$odds_ratio, 4)
  z <- boschloo_probabilities(0, p2 = 0.3)
  expect_true(is.na(z$ratio) && is.na(z$odds_ratio))
  expect_true(all(c("ratio", "odds_ratio") %in% names(z$undefined)))
  expect_identical(cond_class(boschloo_probabilities(0.8, difference = 0.3)), "inferguidR_invalid_input")
  expect_identical(cond_class(boschloo_probabilities(0.2, difference = 0.1, ratio = 2)),
                   "inferguidR_invalid_input")
  expect_identical(cond_class(boschloo_probabilities(0, ratio = 2)), "inferguidR_invalid_input")
})

test_that("the tie tolerance must lie inside the ambiguity band", {
  # Outside the band, ties would be decided without an ambiguity flag: a
  # tolerance of 0 splits the exact tie of (1, 6) and (2, 7) at n = (5, 7)
  # (both have the lower Fisher p-value 1/22 but differ by rounding) and gave
  # p = 0.0174 with status ok; a tolerance above the band merges distinct values.
  expect_identical(cond_class(boschloo_control(tie_tolerance = 0)), "inferguidR_invalid_input")
  expect_identical(cond_class(boschloo_control(tie_tolerance = 1e-13)), "inferguidR_invalid_input")
  expect_identical(cond_class(boschloo_control(tie_tolerance = 1e-5)), "inferguidR_invalid_input")
  expect_identical(cond_class(boschloo_control(ambiguity_band = c(1e-9, 1e-6))), "inferguidR_invalid_input")
  expect_identical(cond_class(boschloo_control(ambiguity_band = c(0, 1e-11))), "inferguidR_invalid_input")
  e <- tryCatch(boschloo_control(tie_tolerance = 1e-5), error = function(e) e)
  expect_match(paste(unlist(e$status$details$problems), collapse = " "), "within `ambiguity_band`")
  # The edges themselves are allowed. With tolerance 0 and the band starting
  # at 0 every rounding gap is flagged, and the exact fallback (N <= 56)
  # restores the tie.
  expect_s3_class(boschloo_control(tie_tolerance = 1e-12), "boschloo_control")
  expect_s3_class(boschloo_control(tie_tolerance = 1e-4, ambiguity_band = c(1e-12, 1e-4)), "boschloo_control")
  zero <- boschloo_control(tie_tolerance = 0, ambiguity_band = c(0, 1e-6))
  r <- boschloo_test(1, 5, 6, 7, variant = "less", control = zero)
  expect_identical(r$status$code, "ok")
  expect_lt(abs(r$p_value - 0.02317922), 5e-9)
  expect_identical(r$sides$lower$ordering_arithmetic, "exact_integer")
  # The generic method has the fixed band [1e-12, 1e-6].
  expect_identical(cond_class(analysis_spec("boschloo_proportions_two_independent", "post_hoc",
                                            inputs = list(p1 = 0.2, p2 = 0.7, alpha = 0.05, n1 = 10, n2 = 10,
                                                          tie_tolerance = 1e-5),
                                            options = list(variant = "less"))),
                   "inferguidR_invalid_input")
})
