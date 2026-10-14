test_that("conditional cells use attained engine probabilities and preserve identity", {
  s <- analysis_spec("exact_proportion_one_sample_binomial", "post_hoc",
    inputs = list(g = .2, p0 = .5, n_total = 20, alpha = .05),
    options = list(tails = "two"))
  r <- calculate_power(s)
  cells <- power_plot_data(r, "decision_table")
  expect_identical(cells$status, "ok")
  expect_identical(cells$spec_hash, r$spec_hash)
  expect_equal(cells$data$probability, c(r$achieved$alpha_attained,
    1 - r$achieved$alpha_attained, r$achieved$power, r$achieved$beta))
  expect_equal(as.numeric(tapply(cells$data$probability, cells$data$hypothesis, sum)), c(1, 1))
  mass <- power_plot_data(r)
  expect_identical(mass$ordinate, "probability_mass")
  expect_equal(as.numeric(tapply(mass$data$y, mass$data$hypothesis, sum)), c(1, 1))
  selected <- mass$data[mass$data$reject, ]
  expect_equal(as.numeric(tapply(selected$y, selected$hypothesis, sum)),
    c(r$achieved$alpha_attained, r$achieved$power), tolerance = 1e-12)
  mosaic <- power_plot_data(r, "mosaic")$data
  expect_equal(mosaic$ymax - mosaic$ymin, cells$data$probability)
})

test_that("normal display preserves the declared law and one-sided direction", {
  s <- analysis_spec("z_generic", "post_hoc",
    inputs = list(mu = -2, sigma = 1.5, alpha = .05), options = list(tails = "one"))
  r <- calculate_power(s)
  p <- power_plot_data(r)
  expect_identical(p$status, "ok")
  h1 <- p$data[p$data$hypothesis == "H1", ]
  expect_equal(h1$y, dnorm(h1$x, -2, 1.5))
  expect_identical(h1$reject, h1$x <= qnorm(.05))
  expect_equal(pnorm(qnorm(.05), -2, 1.5), r$achieved$power)
  expect_error(power_plot_data(r, points = 2), "51 to 2001")
  expect_error(power_plot_data(list()), "power_result")
})

test_that("continuous family adapters retain scaling and noncentrality", {
  laws <- list(list(family = "t", df = 29, ncp = 2),
               list(family = "F", df1 = 3, df2 = 40, ncp = 5),
               list(family = "chisq", df = 4, ncp = 3),
               list(family = "scaled_F", df1 = 10, df2 = 20, scale = 2),
               list(family = "scaled_chisq", df = 12, scale = 1.5))
  expected <- c(dt(2, 29, ncp = 2), df(2, 3, 40, ncp = 5),
                dchisq(2, 4, ncp = 3), df(1, 10, 20) / 2,
                dchisq(2 / 1.5, 12) / 1.5)
  expect_equal(vapply(laws, .ig_plot_law, numeric(1), x = 2, operation = "d"), expected)
  for (law in laws) {
    q <- .ig_plot_law(law, c(.1, .5, .9), "q")
    expect_true(all(is.finite(q)) && all(diff(q) > 0))
  }
})

test_that("two-sided t display includes both critical boundaries", {
  s <- analysis_spec("t_means_one_sample", "post_hoc",
    inputs = list(d = .5, alpha = .05, n_total = 30), options = list(tails = "two"))
  r <- calculate_power(s)
  p <- power_plot_data(r)
  expect_identical(p$status, "ok")
  expect_true(all(r$quantities$critical %in% p$data$x))
  expect_equal(p$data$y[p$data$hypothesis == "H1"],
    dt(p$data$x[p$data$hypothesis == "H1"], 29, ncp = sqrt(30) * .5))
})

test_that("unsupported and inconsistent descriptors never acquire a substitute curve", {
  s <- analysis_spec("z_generic", "post_hoc",
    inputs = list(mu = 2, sigma = 1, alpha = .05), options = list(tails = "one"))
  custom <- function(family, power = .7) .ig_power_result(s, ig_status("ok"), solved = list(power = power),
    achieved = list(alpha_attained = .05, power = power, beta = 1 - power),
    quantities = list(critical = qnorm(.95)),
    method_outputs = list(statistic_distribution = list(
      h0 = list(family = family, mean = 0, sd = 1),
      h1 = list(family = family, mean = 2, sd = 1), rejection_side = "upper")))
  unsupported <- power_plot_data(custom("special_joint_law"))
  expect_identical(unsupported$status, "unavailable")
  expect_equal(nrow(unsupported$data), 0L)
  inconsistent <- power_plot_data(custom("normal"))
  expect_identical(inconsistent$status, "unavailable")
  expect_match(inconsistent$reason, "do not reproduce")
  expect_identical(power_plot_data(custom("normal"), "decision_table")$status, "ok")
  declined <- .ig_power_result(s, ig_status("unsupported", "Test domain unavailable"))
  expect_identical(power_plot_data(declined)$status, "unavailable")
})

test_that("static panels retain the same numeric source in both languages", {
  skip_if_not_installed("ggplot2")
  s <- analysis_spec("t_means_one_sample", "post_hoc",
    inputs = list(d = .5, alpha = .05, n_total = 30), options = list(tails = "two"))
  r <- calculate_power(s)
  for (panel in c("distribution", "decision_table", "mosaic")) {
    en <- plot_power(r, panel, "en")
    de <- plot_power(r, panel, "de")
    expect_identical(attr(en, "inferguidR_plot_data"), attr(de, "inferguidR_plot_data"))
    expect_identical(attr(en, "inferguidR_plot_data")$result_hash, r$result_hash)
    expect_no_error(ggplot2::ggplot_build(en))
    expect_no_error(ggplot2::ggplot_build(de))
  }
})
