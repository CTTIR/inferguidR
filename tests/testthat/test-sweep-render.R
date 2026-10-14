test_that("sweep gaps and stored point identities survive rendering", {
  skip_if_not_installed("ggplot2")
  s <- analysis_spec("t_means_two_independent", "post_hoc",
    inputs = list(d = .5, n1 = 10, n2 = 10, alpha = .05), options = list(tails = "two"))
  z <- xy_plot_data(s, "power", "sample_size", 20, 24, 1)
  seg <- .ig_xy_render_segments(z)
  expect_identical(seg$status, c("ok", "invalid_input", "ok", "invalid_input", "ok"))
  expect_equal(length(unique(na.omit(seg$segment))), 3)
  p <- plot_design_sweep(z)
  expect_equal(length(ggplot2::ggplot_build(p)$data), 1)
  expect_equal(nrow(ggplot2::ggplot_build(p)$data[[1]]), 3)
  expect_identical(attr(p, "inferguidR_plot_data"), z)
  changed <- z
  changed$points$y[1] <- .99
  expect_error(plot_design_sweep(changed), "differ from the stored point result")
  for (y in c("power", "sample_size")) {
    z <- if (y == "power") xy_plot_data(s, y, "sample_size", 20, 24, 2) else
      xy_plot_data(s, y, "power", .70, .74, .02)
    p <- plot_design_sweep(z)
    expect_no_error(ggplot2::ggplot_build(p))
    expect_identical(attr(p, "inferguidR_plot_data"), z)
    expect_equal(length(unique(.ig_xy_render_segments(z)$segment)), 1)
  }
})
