test_that("binomial confidence limits agree with independent extreme-tail references", {
  # 70-digit normal-tail inversion, Wilson score roots and binomial tail sums.
  ref <- read.csv(test_path("fixtures", "classifier-confidence-boundary.csv"), colClasses = "character")
  ref[] <- lapply(ref, .ig_parse_decimal)
  for (i in seq_len(nrow(ref))) {
    r <- ref[i, ]
    w <- .ig_cls_wilson(r$x, r$n, r$level)
    cp <- .ig_cls_clopper_pearson(r$x, r$n, r$level)
    expect_lte(max(abs(c(w$lower, w$upper) - c(r$wilson_lower, r$wilson_upper))), 2e-14)
    expect_lte(max(abs(c(cp$lower, cp$upper) - c(r$cp_lower, r$cp_upper))), 2e-14)
  }
})

test_that("the largest valid confidence level survives paired reports and JSON", {
  d <- classification_examples("scored_synthetic")$data
  levels <- c(.95, 1 - 1e-12, 1 - .Machine$double.eps / 2)
  # Independently inverted standard-normal upper tails at exact double inputs.
  critical <- c(1.9599639845400538, 7.130509892879272, 8.292361075813595)
  for (i in seq_along(levels)) {
    s <- classification_spec(d, truth = "outcome", positive_label = "event",
      negative_label = "no_event", scores = c("model_a", "model_b"), id = "id",
      score_scale = "probability", cutoff = .5,
      uncertainty = classification_uncertainty(level = levels[i]))
    r <- evaluate_cutoff(s)
    ci <- classification_table(r, "intervals")
    ok <- ci$status == "ok"
    expect_true(all(is.finite(ci$lower[ok]) & is.finite(ci$upper[ok])))
    pair <- classification_table(r, "comparisons")
    expect_identical(pair$inference_status, "ok")
    expect_equal(pair$lower, pair$difference - critical[i] * pair$se, tolerance = 2e-14)
    expect_equal(pair$upper, pair$difference + critical[i] * pair$se, tolerance = 2e-14)
    x <- create_report_snapshot(analysis_envelope(s, r), report_kind = "comparison")
    expect_identical(ig_to_json(ig_from_json(ig_to_json(x))), ig_to_json(x))
  }
})

test_that("individual DeLong limits remain finite before probability clipping", {
  level <- 1 - .Machine$double.eps / 2
  # A small variance keeps both reference limits in the interior.
  interval <- .ig_cls_delong_interval(list(status = "ok", var = 1e-6, auc = .5), level, FALSE)
  expect_equal(c(interval$lower, interval$upper),
    c(.491707638924186405, .508292361075813595), tolerance = 2e-14)
  expect_true(is.na(interval$note))
})

test_that("paired report confidence levels must equal the frozen specification", {
  for (levels in list(c(.95, .95 + 5e-11), c(1 - 1e-12, 1 - 2e-12))) {
    s <- update(classification_examples("scored_synthetic")$spec,
      uncertainty = classification_uncertainty(level = levels[1]))
    y <- unclass(create_report_snapshot(analysis_envelope(s, evaluate_cutoff(s))))
    z <- y$results[[1]]$comparisons
    z$level <- levels[2]
    critical <- stats::qnorm((1 - z$level) / 2, lower.tail = FALSE)
    z$lower <- z$difference - critical * z$se
    z$upper <- z$difference + critical * z$se
    y$results[[1]]$comparisons <- z
    y$comparisons[[1]] <- z
    y$visualisations <- list()
    y$snapshot_hash <- NULL
    y$snapshot_hash <- .ig_hash_value(y)
    expect_error(ig_from_json(.ig_canonical_json(y)), class = "inferguidR_invalid_input")
  }
})
