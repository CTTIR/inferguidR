probability_panel_spec <- function() {
  classification_spec(data.frame(y = c(0, 1, 0, 1, 1, 0), p = c(.1, .2, .4, .6, .8, .9)),
    truth = "y", positive_label = 1, scores = "p", score_scale = "probability", cutoff = .7)
}

test_that("probability panels use frozen tables and independently checked arithmetic", {
  spec <- probability_panel_spec()
  r <- evaluate_cutoff(spec, threshold_probability = c(.2, .5, .8), calibration_bins = 3)
  cal <- classification_plot_data(r, "calibration")
  expect_equal(cal$data$points$mean_predicted, c(.15, .5, .85), tolerance = 1e-14)
  expect_equal(cal$data$points$observed_rate, rep(.5, 3))
  expect_equal(cal$data$points$n_units, rep(2, 3))
  brier <- cal$data$summaries$value[cal$data$summaries$metric_id == "brier_score"]
  expect_equal(brier, sum((c(.1, .2, .4, .6, .8, .9) - c(0, 1, 0, 1, 1, 0))^2) / 6)
  nb <- classification_plot_data(r, "net_benefit")
  expect_equal(nb$data$points$tp, c(3, 2, 1))
  expect_equal(nb$data$points$fp, c(2, 1, 1))
  expect_equal(nb$data$points$net_benefit, c(3/6 - 2/6/4, 2/6 - 1/6, 1/6 - 4/6))
  expect_equal(nb$data$points$net_benefit_alert_all, .5 - .5 * c(.2/.8, 1, 4))
  expect_equal(nb$data$points$net_benefit_alert_none, rep(0, 3))
  restored <- ig_from_json(ig_to_json(analysis_envelope(spec, r)))$result
  expect_identical(classification_plot_data(restored, "net_benefit"), nb)
  expect_identical(classification_plot_data(restored, "calibration"), cal)
  for (lang in c("en", "de")) for (panel in c("calibration", "net_benefit")) {
    p <- plot_classification(r, panel, lang = lang)
    expect_s3_class(p, "ggplot")
    expect_identical(attr(p, "inferguidR_plot_data"), classification_plot_data(r, panel, lang = lang))
  }
})

test_that("probability diagnostics do not invent scenario or unsupported results", {
  r <- evaluate_cutoff(update(probability_panel_spec(), prevalence = .2), threshold_probability = .5)
  for (panel in c("calibration", "net_benefit")) {
    expect_error(classification_plot_data(r, panel, basis = "scenario"), class = "inferguidR_unsupported")
    for (kind in c("counts_fixed", "binormal_equal_variance"))
      expect_error(classification_plot_data(evaluate_cutoff(classification_examples(kind)$spec), panel),
        class = "inferguidR_unsupported")
  }
  no_threshold <- evaluate_cutoff(probability_panel_spec())
  expect_error(classification_plot_data(no_threshold, "net_benefit"), class = "inferguidR_unsupported")
  d <- data.frame(y = c(0, 1, 1, 0), p = c(.1, .2, .8, .9))
  raw <- classification_spec(d, truth = "y", positive_label = 1, scores = "p", cutoff = .5)
  reversed <- classification_spec(d, truth = "y", positive_label = 1, scores = "p", cutoff = .5,
    score_scale = "probability", direction = "lower_is_positive")
  for (spec in list(raw, reversed)) for (panel in c("calibration", "net_benefit"))
    expect_error(classification_plot_data(evaluate_cutoff(spec, threshold_probability = .5), panel),
      class = "inferguidR_unsupported")
})

test_that("probability controls parse locale and reject malformed thresholds", {
  expect_equal(.ig_classifier_probability_thresholds("0,2; 0,5", "de"), c(.2, .5))
  expect_equal(.ig_classifier_probability_thresholds("0.5; 0.2; 0.5", "en"), c(.2, .5))
  expect_null(.ig_classifier_probability_thresholds("", "en"))
  for (bad in c(".2;", ";.2", ".2;;.4", "0", "1", "NaN", "garbage"))
    expect_error(.ig_classifier_probability_thresholds(bad, "en"), class = "inferguidR_invalid_input")
  expect_error(.ig_classifier_evaluation_settings(list(calibration_bins = "1"), "en"), class = "inferguidR_invalid_input")
})

test_that("evaluation settings invalidate a prior classifier result and report eligibility", {
  shiny::testServer(.ig_classifier_workspace_server, args = list(language = shiny::reactive("en")), {
    session$setInputs(example = "scored_synthetic", selection = "manual", cutoff = "0.3", scenario = FALSE,
      threshold_probability = "0.2; 0.5", calibration_bins = "3")
    session$setInputs(evaluate = 1)
    old <- result()
    expect_true(evaluation_current())
    session$setInputs(`report-freeze` = 1)
    frozen <- report_module$snapshot()
    expect_s3_class(frozen, "report_document")
    session$setInputs(threshold_probability = "0.2; 0.6")
    expect_identical(candidate()$spec_hash, old$spec_hash)
    expect_false(evaluation_current())
    session$setInputs(`report-freeze` = 2)
    expect_identical(report_module$snapshot()$snapshot_hash, frozen$snapshot_hash)
    session$setInputs(evaluate = 2)
    expect_true(evaluation_current())
    expect_equal(result()$settings$threshold_probability, c(.2, .6))
    session$setInputs(calibration_bins = "2")
    expect_false(evaluation_current())
    session$setInputs(evaluate = 3)
    expect_true(evaluation_current())
    session$setInputs(threshold_probability = "0;")
    before_failure <- result()
    session$setInputs(evaluate = 4)
    expect_identical(result(), before_failure)
    expect_false(evaluation_current())
    expect_match(output$message, "previous evaluation")
  })
})
