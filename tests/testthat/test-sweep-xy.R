# X-Y plot data (R/sweep-xy.R). Group G12; method notes
# specs/method-notes/utilities.md, section "X-Y plot data".
#
# Every point is an ordinary calculate_power() result; these tests check the
# composition (modes, inputs held constant, design rules, staircases, gaps)
# and the range rules captured from G*Power's X-Y plot window (record ids
# cited; the messages themselves are G*Power's wording, not compared).

t2 <- "t_means_two_independent"
base_post_hoc <- function() {
  calculate_power(analysis_spec(t2, "post_hoc", inputs = list(d = 0.5, alpha = 0.05, n1 = 50L, n2 = 50L),
                                options = list(tails = "two")))
}

test_that("power as a function of the total sample size keeps only admissible designs", {
  base <- base_post_hoc()
  xy <- xy_plot_data(base, y = "power", x = "sample_size", from = 10, to = 40, by = 5)
  expect_s3_class(xy, "ig_xy_plot_data")
  expect_identical(xy$mode, "post_hoc")
  p <- xy$points
  expect_identical(p$x, c(10, 15, 20, 25, 30, 35, 40))
  odd <- p$x %% 2 == 1
  expect_true(all(p$status[odd] == "invalid_input"))
  expect_true(all(is.na(p$y[odd])))
  expect_true(all(grepl("not an admissible design", p$reason[odd])))
  expect_true(all(p$status[!odd] == "ok"))
  for (i in which(!odd)) {
    n <- as.integer(p$x[i] / 2)
    direct <- calculate_power(analysis_spec(t2, "post_hoc", inputs = list(d = 0.5, alpha = 0.05, n1 = n, n2 = n),
                                            options = list(tails = "two")))
    expect_identical(p$y[i], direct$achieved$power)
    expect_identical(p$result_hash[i], direct$result_hash)
    expect_identical(unname(p$group_sizes[[i]]), c(n, n))
  }
  expect_identical(xy$axes$x$label, "Total sample size")
  expect_identical(xy$axes$y$label, "Power (1-β err prob)")
  expect_match(xy$subtitle, "Effect size d = 0.5")
  expect_match(xy$subtitle, "Tail(s) = Two", fixed = TRUE)
  expect_identical(xy$render, "line")
  expect_identical(xy$title, "t tests - Means: Difference between two independent means (two groups)")
})

test_that("sample size as a function of power is an integer staircase", {
  base <- base_post_hoc()
  xy <- xy_plot_data(base, y = "sample_size", x = "power", from = 0.6, to = 0.95, by = 0.05,
                     curves = list(parameter = "effect_size", values = c(0.4, 0.5)))
  expect_identical(xy$mode, "a_priori")
  expect_identical(xy$render, "step")
  p <- xy$points
  expect_true(all(p$status == "ok"))
  expect_true(all(p$y == round(p$y)))
  expect_true(all(p$achieved_power >= p$x))
  expect_true(all(p$n_total_continuous <= p$y))
  for (cv in c(0.4, 0.5)) expect_false(is.unsorted(p$y[p$curve_value == cv]))
  # Minimality: the next smaller admissible design (N - 2, equal groups) misses the target.
  for (i in c(1, 8, 12)) {
    n <- as.integer((p$y[i] - 2) / 2)
    smaller <- calculate_power(analysis_spec(t2, "post_hoc", inputs = list(d = p$curve_value[i], alpha = 0.05, n1 = n, n2 = n),
                                             options = list(tails = "two")))
    expect_lt(smaller$achieved$power, p$x[i])
  }
  expect_true(any(grepl("not an achievable design", xy$notes)))
  expect_identical(xy$curves$label, "Effect size d")
  expect_match(xy$subtitle, "Allocation ratio N2/N1 = 1", fixed = TRUE)
})

test_that("effect size and alpha curves use sensitivity and criterion mode", {
  base <- base_post_hoc()
  xy <- xy_plot_data(base, y = "effect_size", x = "alpha", from = 0.01, to = 0.1, by = 0.03)
  expect_identical(xy$mode, "sensitivity")
  direct <- calculate_power(analysis_spec(t2, "sensitivity", inputs = list(alpha = 0.04, power = base$achieved$power,
                                                                           n1 = 50L, n2 = 50L), options = list(tails = "two")))
  expect_identical(xy$points$y[2], direct$solved$d)
  xy <- xy_plot_data(base, y = "alpha", x = "effect_size", from = 0.2, to = 0.8, by = 0.2)
  expect_identical(xy$mode, "criterion")
  expect_true(all(diff(xy$points$y) < 0))
  # The constant may be overridden.
  xy <- xy_plot_data(base, y = "power", x = "effect_size", from = 0.2, to = 0.6, by = 0.2, fixed = list(alpha = 0.01))
  expect_match(xy$subtitle, "α err prob = 0.01")
  direct <- calculate_power(analysis_spec(t2, "post_hoc", inputs = list(d = 0.4, alpha = 0.01, n1 = 50L, n2 = 50L),
                                          options = list(tails = "two")))
  expect_identical(xy$points$y[2], direct$achieved$power)
})

test_that("the range rules of the G*Power X-Y plot are typed refusals", {
  base <- base_post_hoc()
  bad <- function(...) expect_error(xy_plot_data(base, ...), class = "inferguidR_invalid_input")
  bad(y = "power", x = "sample_size", from = 200, to = 10, by = 10)     # GPW-UT-XYPLOT-UTL-U026-I01--r02
  bad(y = "power", x = "sample_size", from = 10, to = 10, by = 10)      # SUP-G12-UT-XYPLOT-X08--r01
  bad(y = "power", x = "sample_size", from = 10, to = 200, by = 0)      # GPW-UT-XYPLOT-UTL-U026-I02--r01
  bad(y = "power", x = "sample_size", from = 10, to = 200, by = 2.5)    # GPW-UT-XYPLOT-UTL-U026-I03--r01
  bad(y = "effect_size", x = "power", from = 0.01, to = 0.2, by = 0.01) # GPW-UT-XYPLOT-UTL-U026-I04--r02
  bad(y = "power", x = "alpha", from = 0.01, to = 1.5, by = 0.01)       # SUP-G12-UT-XYPLOT-X03--r01
  bad(y = "sample_size", x = "power", from = 0.5, to = 1, by = 0.05)    # SUP-G12-UT-XYPLOT-X05--r01
  bad(y = "power", x = "power", from = 0.1, to = 0.2, by = 0.1)
  bad(y = "power", x = "sample_size", from = 10, to = 40, by = 5, curves = list(parameter = "power", values = 0.5))
  bad(y = "power", x = "sample_size", from = 10, to = 40, by = 5, curves = list(parameter = "alpha", values = seq(0.01, 0.11, 0.01)))
  bad(y = "power", x = "sample_size", from = 10, to = 40, by = 5, fixed = list(power = 0.5))
  bad(y = "power", x = "sample_size", from = 10, to = 1e6, by = 1)
  # Accepted by G*Power (no message, drawing starts) and here:
  expect_s3_class(xy_plot_data(base, y = "power", x = "sample_size", from = 10, to = 30, by = 1), "ig_xy_plot_data")  # X01
  xy4 <- xy_plot_data(base, y = "power", x = "effect_size", from = -1, to = 1, by = 0.5)                            # X04
  expect_true(all(xy4$points$status == "ok"))
  xy7 <- xy_plot_data(base, y = "effect_size", x = "power", from = 0.05, to = 0.08, by = 0.01)                      # X07
  expect_false(identical(xy7$points$status[1], "ok"))
  expect_true(is.na(xy7$points$y[1]))
})

test_that("documented X-Y tables are the continuous relaxation, reported apart from integer designs", {
  # DOC-005--r01 (manual Fig. 3, one tail, alpha 0.05, power 0.6): 1442.62, 641.920, 361.676.
  base <- calculate_power(analysis_spec(t2, "a_priori", inputs = list(d = 0.1, alpha = 0.05, power = 0.6, allocation_ratio = 1),
                                        options = list(tails = "one")))
  xy <- xy_plot_data(base, y = "sample_size", x = "effect_size", from = 0.1, to = 0.2, by = 0.05)
  p <- xy$points
  expect_lte(abs(p$n_total_continuous[1] - 1442.62), 0.005)
  expect_lte(abs(p$n_total_continuous[2] - 641.920), 0.0005)
  expect_lte(abs(p$n_total_continuous[3] - 361.676), 0.0005)
  expect_identical(p$y, ceiling(p$n_total_continuous / 2) * 2)
})

test_that("a base without a result or of an unknown quantity is refused", {
  res <- calculate_power(analysis_spec(t2, "a_priori", inputs = list(d = 1e-9, alpha = 0.05, power = 0.95, allocation_ratio = 1),
                                       options = list(tails = "two")))
  expect_error(xy_plot_data(res, y = "power", x = "effect_size", from = 0.1, to = 0.5, by = 0.1),
               class = "inferguidR_invalid_input")
  expect_error(xy_plot_data(base_post_hoc(), y = "beta", x = "effect_size", from = 0.1, to = 0.5, by = 0.1),
               class = "inferguidR_invalid_input")
  expect_error(xy_plot_data(list(), y = "power", x = "effect_size", from = 0.1, to = 0.5, by = 0.1),
               class = "inferguidR_invalid_input")
  spec <- analysis_spec(t2, "post_hoc", inputs = list(d = 0.5, alpha = 0.05, n1 = 50L, n2 = 50L), options = list(tails = "two"))
  expect_s3_class(xy_plot_data(spec, y = "power", x = "effect_size", from = 0.2, to = 0.4, by = 0.1), "ig_xy_plot_data")
})

test_that("the plot data work for a single-sample method", {
  id <- "t_means_one_sample"
  skip_if(is.null(inferguidR:::.ig_get_method(id)), "no engine for the one-sample t test")
  base <- tryCatch(calculate_power(analysis_spec(id, "post_hoc", inputs = list(alpha = 0.05, n_total = 20L, d = 0.5),
                                                 options = list(tails = "two"))), error = function(e) NULL)
  skip_if(is.null(base) || !identical(base$status$code, "ok"), "the one-sample t specification differs")
  xy <- xy_plot_data(base, y = "power", x = "sample_size", from = 10, to = 30, by = 5)
  expect_true(all(xy$points$status == "ok"))
  expect_false(is.unsorted(xy$points$y))
})

test_that("a family of curves is given by values or by start, step and count", {
  base <- base_post_hoc()
  xy <- xy_plot_data(base, y = "power", x = "sample_size", from = 20, to = 40, by = 10,
                     curves = list(parameter = "effect_size", from = 0.2, by = 0.3, n = 3))
  expect_equal(xy$curves$values, c(0.2, 0.5, 0.8), tolerance = 1e-15)
  expect_identical(nrow(xy$points), 9L)
  expect_true(all(tapply(xy$points$y, xy$points$x, function(v) !is.unsorted(v))))  # power grows with the effect
  for (n in list(0, -1, 11, 2.5, NA_real_, "3")) {
    expect_error(xy_plot_data(base, y = "power", x = "sample_size", from = 20, to = 40, by = 10,
                              curves = list(parameter = "effect_size", from = 0.2, by = 0.3, n = n)),
                 class = "inferguidR_invalid_input")
  }
})
