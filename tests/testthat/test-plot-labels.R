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
  expect_match(gaps$labels$caption, "\u00fcber L\u00fccken", fixed = TRUE)
})

test_that("classification panels use German umlauts and unpadded node values", {
  skip_if_not_installed("ggplot2")
  r <- evaluate_cutoff(classification_examples("scored_synthetic")$spec)
  roc <- plot_classification(r, "roc", "model_b", lang = "de")
  expect_match(roc$labels$caption, "gew\u00e4hlter Schwellenwert", fixed = TRUE)
  expect_false(grepl("Praevalenz|gewaehlt|geaendert|fuer ", roc$labels$caption))
  pr <- plot_classification(r, "pr", "model_b", lang = "de")
  expect_identical(pr$labels$x, "Sensitivit\u00e4t (Recall)")
  network <- plot_classification(r, "metric_network", "model_b")
  values <- vapply(strsplit(ggplot2::ggplot_build(network)$data[[2]]$label, "\n", fixed = TRUE),
    `[[`, "", 2L)
  expect_identical(values, trimws(values))
})

test_that("German classification subtitles translate the reference basis", {
  expect_identical(.ig_cls_basis_subtitle("observed", FALSE), "Basis: observed")
  expect_identical(.ig_cls_basis_subtitle("observed", TRUE), "Bezugsbasis: beobachtet")
  expect_identical(.ig_cls_basis_subtitle("scenario", TRUE), "Bezugsbasis: Szenario")
  skip_if_not_installed("ggplot2")
  r <- evaluate_cutoff(classification_examples("scored_synthetic")$spec)
  expect_identical(plot_classification(r, "mosaic", "model_b", lang = "de")$labels$subtitle,
    "Bezugsbasis: beobachtet")
})

test_that("titles, captions and bottom legends stay inside the default figure width", {
  skip_if_not_installed("ggplot2")
  width <- 7
  # Horizontal extent in inches of a layout element, measured on a 7-inch
  # wide device with the standard PostScript font metrics.
  extent <- function(p, name) {
    grDevices::pdf(NULL, width = width, height = width / 1.618)
    on.exit(grDevices::dev.off())
    gt <- suppressWarnings(ggplot2::ggplotGrob(p))
    i <- which(gt$layout$name == name)
    if (!length(i)) return(0)
    g <- gt$grobs[[i]]
    w <- if (inherits(g, "gtable")) sum(g$widths) else grid::grobWidth(g)
    left <- 0
    if (name != "guide-box-bottom" && gt$layout$l[i] > 1L) {
      left <- grid::convertWidth(sum(gt$widths[seq_len(gt$layout$l[i] - 1L)]), "in",
                                 valueOnly = TRUE)
    }
    left + grid::convertWidth(w, "in", valueOnly = TRUE)
  }
  r <- calculate_power(analysis_spec("t_means_one_sample", "post_hoc",
    inputs = list(d = .5, alpha = .05, n_total = 30), options = list(tails = "two")))
  cr <- evaluate_cutoff(classification_examples("scored_synthetic")$spec)
  s <- analysis_spec("t_means_two_independent", "post_hoc",
    inputs = list(d = .5, n1 = 20, n2 = 20, alpha = .05), options = list(tails = "two"))
  for (lang in c("en", "de")) {
    plots <- list(power_decision_table = plot_power(r, "decision_table", lang = lang),
      power_mosaic = plot_power(r, "mosaic", lang = lang),
      power_distribution = plot_power(r, "distribution", lang = lang),
      design_sweep = plot_design_sweep(xy_plot_data(s, "power", "sample_size", 20, 24, 1), lang))
    for (panel in c("roc", "pr", "mosaic", "metric_network", "metric_pair", "prevalence_curve",
                    "model_comparison", "score_distribution"))
      plots[[panel]] <- plot_classification(cr, panel, "model_b", lang = lang)
    for (name in names(plots)) {
      for (text in c("title", "subtitle", "caption"))
        expect_lte(extent(plots[[name]], text), 0.9 * width, label = paste(lang, name, text))
      expect_lte(extent(plots[[name]], "guide-box-bottom"), 0.85 * width,
                 label = paste(lang, name, "legend"))
    }
  }
})
