test_that("large finite calibration-bin requests retain finite observed groups", {
  d <- data.frame(y=c(0,1,0,1,1,0),p=c(.1,.2,.4,.6,.8,.9))
  s <- classification_spec(d,truth="y",positive_label=1,scores="p",score_scale="probability",cutoff=.5)
  for (bins in c(1e10,1e308)) {
    r <- evaluate_cutoff(s,calibration_bins=bins)
    pd <- classification_plot_data(r,"calibration")
    expect_equal(pd$data$points$group,seq_len(6))
    expect_equal(pd$data$points$n_units,rep(1,6))
    expect_equal(pd$data$points$mean_predicted,d$p)
    expect_equal(pd$data$points$observed_rate,d$y)
    expect_true(all(vapply(pd$data$points,function(x) all(is.finite(x)),TRUE)))
    doc <- create_report_snapshot(analysis_envelope(s,r),include_probability_summaries=TRUE)
    expect_identical(ig_from_json(ig_to_json(doc)),doc)
  }
  d$p[2] <- d$p[1]
  s <- classification_spec(d,truth="y",positive_label=1,scores="p",score_scale="probability",cutoff=.5)
  pd <- classification_plot_data(evaluate_cutoff(s,calibration_bins=1e10),"calibration")
  expect_equal(pd$data$points$n_units,c(2,1,1,1,1))
  expect_equal(pd$data$points$mean_predicted,unique(d$p))
})
