# Boschloo plot data: every plotted number equals the calculated one (same
# region, same p-values, same power), and observed and planning quantities
# are labelled.

test_that("the lattice shows exactly the rejection region and the per-table p-values", {
  pd <- boschloo_plotdata_lattice(6, 9, alpha = 0.05, variant = "minlike", p1 = 0.2, p2 = 0.7,
                                  observed = c(1, 7))
  d <- pd$data
  expect_identical(nrow(d), 70L)
  reg <- boschloo_region(6, 9, 0.05, "minlike")
  expect_identical(sort(paste(d$x1, d$x2)[d$in_region]), sort(paste(reg$rejected$x1, reg$rejected$x2)))
  expect_identical(d$in_region, d$boschloo_p_upper <= 0.05)
  for (i in c(1, 17, 40, 70)) {
    expect_equal(d$boschloo_p[i], boschloo_test(d$x1[i], 6, d$x2[i], 9, variant = "minlike")$p_value,
                 tolerance = 1e-12)
  }
  expect_equal(sum(d$prob_alternative), 1, tolerance = 1e-12)
  expect_equal(sum(d$prob_alternative[d$in_region]), pd$reference$power, tolerance = 1e-12)
  expect_equal(pd$reference$power, boschloo_power(6, 9, 0.2, 0.7, variant = "minlike")$power,
               tolerance = 1e-15)
  expect_identical(sum(d$observed), 1L)
  expect_match(pd$scope, "observed one")
  # Marking an observed table changes nothing else.
  pd0 <- boschloo_plotdata_lattice(6, 9, alpha = 0.05, variant = "minlike", p1 = 0.2, p2 = 0.7)
  expect_identical(pd0$data[, names(d) != "observed"], d[, names(d) != "observed"])
  expect_output(print(pd), "boschloo_lattice")
})

test_that("the nuisance envelope peaks at the p-value and at the attained size", {
  for (v in c("minlike", "central", "greater")) {
    r <- boschloo_test(9, 15, 3, 15, variant = v)
    env <- boschloo_plotdata_nuisance(r)
    peaks <- as.vector(tapply(env$data$probability, env$data$side, max))
    if (v == "central") {
      expect_equal(min(1, 2 * min(peaks)), r$p_value, tolerance = 1e-9)
    } else {
      expect_equal(unname(peaks), r$p_value, tolerance = 1e-9)
    }
    expect_true(all(peaks <= r$sides[[1]]$p_upper + 1e-12 | v == "central"))
    expect_match(env$scope, "^observed data")
  }
  reg <- boschloo_region(12, 10, 0.05, "central")
  env <- boschloo_plotdata_nuisance(reg)
  expect_equal(max(env$data$probability), reg$size$value, tolerance = 1e-12)
  expect_lte(max(env$data$probability), 0.05)
  expect_match(env$scope, "hypothetical future study")
  expect_identical(tryCatch(boschloo_plotdata_nuisance(list()), error = function(e) class(e)[1]),
                   "inferguidR_invalid_input")
})

test_that("power curves equal the power function of the same regions", {
  pd <- boschloo_plotdata_power_n(0.25, 0.75, alpha = 0.05, variant = "less", n_range = c(10, 24),
                                  power = 0.8, fisher = TRUE)
  d <- pd$data
  expect_identical(d$n_total, seq(10L, 24L, by = 2L))
  for (k in seq_len(nrow(d))) {
    expect_identical(d$power[k], boschloo_power(d$n1[k], d$n2[k], 0.25, 0.75, alpha = 0.05,
                                                variant = "less")$power)
  }
  expect_true(all(d$power >= d$fisher_power - 1e-12))     # same convention: Boschloo dominates
  pe <- boschloo_plotdata_power_effect(14, 11, p1 = 0.3, p2 = c(0, 0.3, 0.6, 0.9, 1), variant = "central")
  expect_identical(pe$data$power[3], boschloo_power(14, 11, 0.3, 0.6, variant = "central")$power)
  expect_lte(pe$data$power[2], pe$reference$attained_size + 1e-12)
})

test_that("the Fisher comparison classifies every table and confirms dominance", {
  cmp <- boschloo_compare_fisher(12, 12, alpha = 0.05, variant = "central", p1 = 0.2, p2 = 0.7)
  expect_identical(sum(cmp$counts), 169L)
  expect_identical(cmp$counts[["fisher_only"]], 0L)
  expect_true(cmp$dominance_holds)
  expect_gt(cmp$counts[["boschloo_only"]], 0L)
  expect_gte(cmp$power$boschloo, cmp$power$fisher)
  expect_lte(cmp$size$fisher$upper, 0.05)
  # The minlike Fisher ordering is the two-sided p-value of stats::fisher.test.
  cm <- boschloo_compare_fisher(7, 9, alpha = 0.05, variant = "minlike")
  tb <- cm$tables
  for (i in c(3, 20, 41, 64)) {
    m <- matrix(c(tb$x1[i], 7 - tb$x1[i], tb$x2[i], 9 - tb$x2[i]), 2, byrow = TRUE)
    expect_equal(tb$fisher_p[i], stats::fisher.test(m)$p.value, tolerance = 1e-12)
  }
  other <- boschloo_compare_fisher(12, 12, alpha = 0.05, variant = "central", fisher_convention = "minlike")
  expect_true(is.na(other$dominance_holds))
  expect_match(other$fisher_test$convention, "different convention")
  expect_output(print(cmp), "Boschloo only")
})

test_that("the nuisance envelope uses the settings recorded in its input", {
  ctl <- boschloo_control(grid_points = 11L)
  r <- boschloo_test(2, 9, 7, 10, variant = "minlike", control = ctl)
  env <- boschloo_plotdata_nuisance(r)
  expect_equal(max(env$data$probability), r$p_value, tolerance = 1e-9)
  expect_identical(boschloo_plotdata_nuisance(r, control = ctl)$data, env$data)
  e <- tryCatch(boschloo_plotdata_nuisance(r, control = boschloo_control()), error = function(e) e)
  expect_s3_class(e, "inferguidR_invalid_input")
  # After a JSON round trip of the envelope the recorded settings still apply.
  back <- ig_from_json(ig_to_json(analysis_envelope(r$spec, r)))$result
  expect_identical(.bl_control_key(.bl_control_from_record(back$settings)), .bl_control_key(ctl))
})

test_that("uncertified p-values and regions give no curve", {
  bad <- boschloo_control(grid_points = 3L, max_levels = 1L)
  r <- boschloo_test(10, 35, 24, 61, variant = "minlike", control = bad)
  env <- boschloo_plotdata_nuisance(r)
  expect_identical(env$status$code, "numerical_failure")
  expect_identical(nrow(env$data), 0L)
  reg <- boschloo_region(30, 30, 0.05, "minlike", control = bad)
  env2 <- boschloo_plotdata_nuisance(reg)
  expect_identical(env2$status$code, "numerical_failure")
  expect_identical(nrow(env2$data), 0L)
  expect_null(env2$reference$attained_size)
})

test_that("the generic set supremum reproduces region and Fisher sizes", {
  reg <- boschloo_region(9, 13, 0.05, "central")
  member <- logical(10 * 14)
  member[reg$rejected$x1 + 1 + reg$rejected$x2 * 10] <- TRUE
  s <- .bl_set_sup(9, 13, member)
  expect_equal(s$value, reg$size$value, tolerance = 1e-12)
  expect_lte(s$upper, 0.05)
  cmp <- boschloo_compare_fisher(9, 13, alpha = 0.05, variant = "central")
  fs <- .bl_set_sup(9, 13, cmp$tables$fisher_reject)
  expect_equal(fs$value, cmp$size$fisher$value, tolerance = 1e-12)
  expect_identical(.bl_set_sup(9, 13, logical(140))$value, 0)
  expect_identical(.bl_set_sup(9, 13, rep(TRUE, 140))$value, 1)
  expect_error(.bl_set_sup(9, 13, TRUE), "one element per table")
})

test_that("power versus baseline and versus allocation use the same regions", {
  base <- c(0.05, 0.2, 0.4, 0.6)
  pe <- boschloo_plotdata_power_effect(14, 16, p1 = base, p2 = base + 0.3, variant = "minlike")
  expect_identical(pe$data$p1, base)
  expect_equal(pe$data$difference, rep(0.3, 4), tolerance = 1e-12)
  for (k in seq_along(base)) {
    expect_identical(pe$data$power[k], boschloo_power(14, 16, base[k], base[k] + 0.3, variant = "minlike")$power)
  }
  expect_identical(tryCatch(boschloo_plotdata_power_effect(14, 16, p1 = c(0.1, 0.2), p2 = c(0.3, 0.4, 0.5),
                                                           variant = "less"), error = function(e) class(e)[1]),
                   "inferguidR_invalid_input")
  al <- boschloo_plotdata_allocation(0.2, 0.65, n_total = 24, variant = "central", n1 = c(6, 12, 18))
  expect_identical(al$data$n1, c(6L, 12L, 18L))
  expect_identical(al$data$n2, c(18L, 12L, 6L))
  for (k in 1:3) {
    ref <- boschloo_power(al$data$n1[k], al$data$n2[k], 0.2, 0.65, variant = "central")
    expect_identical(al$data$power[k], ref$power)
    expect_identical(al$data$size[k], ref$size$value)
  }
  expect_true(all(al$data$size <= 0.05))
  expect_match(al$scope, "fixed total")
  full <- boschloo_plotdata_allocation(0.3, 0.8, n_total = 9, variant = "less")
  expect_identical(full$data$n1, 1:8)
  expect_identical(tryCatch(boschloo_plotdata_allocation(0.3, 0.8, n_total = 9, variant = "less", n1 = 9),
                            error = function(e) class(e)[1]), "inferguidR_invalid_input")
})

test_that("results state the design they are valid for", {
  r <- boschloo_test(2, 10, 7, 11, variant = "central")
  expect_match(r$variant$design, "two independent binomial samples")
  expect_match(r$variant$design, "McNemar")
  expect_match(boschloo_region(5, 6, 0.05, "less")$variant$design, "not for paired")
})
