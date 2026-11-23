test_that("power panels print numbers without padding", {
  expect_identical(unname(.ig_power_error_labels(0.05, 0.24604, "en")),
    c("Attained alpha = 0.05", "beta = 0.24604"))
  expect_identical(unname(.ig_power_error_labels(0.05, 0.5, "de")),
    c("Erreichtes alpha = 0,05", "beta = 0,5"))
  skip_if_not_installed("ggplot2")
  r <- calculate_power(analysis_spec("t_means_one_sample", "post_hoc",
    inputs = list(d = .5, alpha = .05, n_total = 30), options = list(tails = "two")))
  labels <- ggplot2::ggplot_build(plot_power(r, "decision_table"))$data[[2]]$label
  expect_identical(labels, trimws(labels))
  expect_true("0.05" %in% labels)
  p <- plot_power(r, "distribution")
  expect_identical(p$theme$legend.box, "vertical")
  rb <- calculate_power(analysis_spec("boschloo_proportions_two_independent", "post_hoc",
    inputs = list(n1 = 4, n2 = 5, p1 = .8, p2 = .2, alpha = .05), options = list(variant = "minlike")))
  expect_match(plot_power(rb, "lattice")$labels$caption, "(common p = 0.5)", fixed = TRUE)
})

test_that("design sweep titles, subtitles and captions fit the default device width", {
  skip_if_not_installed("ggplot2")
  s <- analysis_spec("t_means_two_independent", "post_hoc",
    inputs = list(d = .5, n1 = 20, n2 = 20, alpha = .05), options = list(tails = "two"))
  for (lang in c("en", "de")) {
    p <- plot_design_sweep(xy_plot_data(s, "power", "sample_size", 20, 60, 10), lang)
    expect_true(all(nchar(strsplit(p$labels$title, "\n", fixed = TRUE)[[1]]) <= 55))
    expect_true(all(nchar(strsplit(p$labels$subtitle, "\n", fixed = TRUE)[[1]]) <= 65))
    expect_true(all(nchar(strsplit(p$labels$caption, "\n", fixed = TRUE)[[1]]) <= 85))
    expect_false(grepl("0 points|0 Punkte", p$labels$caption))
  }
  gaps <- plot_design_sweep(xy_plot_data(s, "power", "sample_size", 20, 24, 1), "de")
  expect_match(gaps$labels$caption, "über Lücken", fixed = TRUE)
})

test_that("classification panels use German umlauts and unpadded node values", {
  skip_if_not_installed("ggplot2")
  r <- evaluate_cutoff(classification_examples("scored_synthetic")$spec)
  roc <- plot_classification(r, "roc", "model_b", lang = "de")
  expect_match(roc$labels$caption, "gewählter Schwellenwert", fixed = TRUE)
  expect_false(grepl("Praevalenz|gewaehlt|geaendert|fuer ", roc$labels$caption))
  pr <- plot_classification(r, "pr", "model_b", lang = "de")
  expect_identical(pr$labels$x, "Sensitivität (Recall)")
  network <- plot_classification(r, "metric_network", "model_b")
  values <- vapply(strsplit(ggplot2::ggplot_build(network)$data[[2]]$label, "\n", fixed = TRUE),
    `[[`, "", 2L)
  expect_identical(values, trimws(values))
})
