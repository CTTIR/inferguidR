# Boschloo planning: allocation sequences, the non-monotone power staircase,
# minimality of the first design, the stability window and typed statuses.

test_that("the first design is the minimum of a full scan of the allocation sequence", {
  cases <- list(list(0.2, 0.7, "central", "n2_ceiling", 1),
                list(0.2, 0.7, "central", "total_nearest", 1),
                list(0.6, 0.1, "greater", "n2_ceiling", 2),
                list(0.1, 0.55, "minlike", "total_nearest", 0.5))
  for (cs in cases) {
    ss <- boschloo_sample_size(cs[[1]], cs[[2]], power = 0.8, alpha = 0.05, variant = cs[[3]],
                               allocation = cs[[4]], allocation_ratio = cs[[5]], stability_window = 4L)
    expect_identical(ss$status$code, "ok")
    m_first <- ss$first$index
    # Full scan of every earlier design with the public power function.
    for (m in seq_len(m_first - 1L)) {
      d <- .bl_design_at(m, cs[[5]], cs[[4]])
      pw <- boschloo_power(d$n1, d$n2, cs[[1]], cs[[2]], alpha = 0.05, variant = cs[[3]])$power
      expect_lt(pw, 0.8)
    }
    pw_first <- boschloo_power(ss$first$n1, ss$first$n2, cs[[1]], cs[[2]], alpha = 0.05,
                               variant = cs[[3]])
    expect_identical(pw_first$power, ss$first$power)
    expect_gte(ss$first$power, 0.8)
    expect_identical(ss$first$size, pw_first$size$value)
    expect_true(ss$minimality$first_is_minimum)
  }
})

test_that("power is saw-toothed in N and the stability window is honoured", {
  ss <- boschloo_sample_size(0.2, 0.7, power = 0.8, alpha = 0.05, variant = "central",
                             allocation = "total_nearest", stability_window = 4L)
  expect_identical(c(ss$first$n1, ss$first$n2), c(15L, 15L))
  expect_identical(c(ss$stable$n1, ss$stable$n2), c(16L, 16L))
  st <- ss$staircase
  expect_false(st$meets_target[st$n_total == 31])           # drops below the target after N = 30
  expect_true(all(st$meets_target[st$n_total >= 32 & st$n_total <= 36]))
  expect_identical(ss$minimality$designs_below_target_after_first, 1L)
  expect_match(ss$notes, "not monotone")
  # The staircase is exactly the power of each listed design.
  for (k in seq_len(nrow(st))) {
    expect_identical(st$power[k], boschloo_power(st$n1[k], st$n2[k], 0.2, 0.7, alpha = 0.05,
                                                 variant = "central")$power)
  }
  pd <- boschloo_plotdata_power_n(0.2, 0.7, alpha = 0.05, variant = "central", n_range = c(20, 40),
                                  allocation = "total_nearest", power = 0.8)
  expect_true(any(diff(pd$data$power) < 0))
})

test_that("allocation policies produce the declared group sizes", {
  d <- .bl_design_at(1:6, 2, "n2_ceiling")
  expect_identical(d$n2, 2L * (1:6))
  d <- .bl_design_at(1:4, 1.5, "n2_ceiling")
  expect_identical(d$n2, c(2L, 3L, 5L, 6L))
  d <- .bl_design_at(1:5, 1, "total_nearest")              # N = 2..6
  expect_identical(d$n1, c(1L, 2L, 2L, 3L, 3L))
  expect_identical(d$n1 + d$n2, 2:6)
  d <- .bl_design_at(1:3, 100, "total_nearest")
  expect_true(all(d$n1 >= 1L & d$n2 >= 1L))
  # Both sequences never decrease a group size (needed for the bound below).
  for (al in c("n2_ceiling", "total_nearest")) {
    for (r in c(0.3, 1, 2.5)) {
      d <- .bl_design_at(1:80, r, al)
      expect_true(all(diff(d$n1) >= 0 & diff(d$n2) >= 0))
    }
  }
})

test_that("the Neyman-Pearson bound is never below the Boschloo power", {
  for (nn in list(c(8, 8), c(12, 20), c(25, 25))) {
    for (pp in list(c(0.2, 0.7), c(0.5, 0.1))) {
      b <- .bl_np_bound(nn[1], nn[2], pp[1], pp[2], 0.05)
      for (v in c("minlike", "central")) {
        expect_lte(boschloo_power(nn[1], nn[2], pp[1], pp[2], alpha = 0.05, variant = v)$power, b + 1e-12)
      }
    }
  }
})

test_that("unanswerable requests return typed statuses without numbers", {
  s1 <- boschloo_sample_size(0.4, 0.2, power = 0.8, variant = "less")
  expect_identical(s1$status$code, "invalid_input")
  expect_null(s1$first)
  s2 <- boschloo_sample_size(0.3, 0.3, power = 0.8, variant = "minlike")
  expect_identical(s2$status$code, "infeasible")
  s3 <- boschloo_sample_size(0.3, 0.4, power = 0.95, variant = "central", n_max = 40L)
  expect_identical(s3$status$code, "infeasible")
  expect_match(s3$status$reason, "Neyman-Pearson")
  expect_null(s3$first)
  s4 <- boschloo_sample_size(0.1, 0.6, power = 0.04, alpha = 0.05, variant = "central")
  expect_identical(s4$status$code, "invalid_input")
  expect_identical(tryCatch(boschloo_sample_size(0.1, 0.6, variant = "central", n_max = 1),
                            error = function(e) class(e)[1]), "inferguidR_invalid_input")
  expect_output(print(s3), "infeasible")
})

test_that("baseline-plus-effect planning uses the converted probabilities", {
  pr <- boschloo_probabilities(0.15, odds_ratio = 9)
  ss <- boschloo_sample_size(pr$p1, pr$p2, power = 0.8, alpha = 0.05, variant = "central",
                             stability_window = 2L)
  expect_identical(ss$status$code, "ok")
  expect_identical(ss$inputs$p2, pr$p2)
  expect_match(ss$scope, "planning")
  expect_output(print(ss), "first design reaching the target")
})
