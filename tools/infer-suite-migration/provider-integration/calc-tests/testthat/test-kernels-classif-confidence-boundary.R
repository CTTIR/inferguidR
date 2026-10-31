# Numerical regression blocks retained from the historical guidance suite.
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
test_that("individual DeLong limits remain finite before probability clipping", {
  level <- 1 - .Machine$double.eps / 2
  # A small variance keeps both reference limits in the interior.
  interval <- .ig_cls_delong_interval(list(status = "ok", var = 1e-6, auc = .5), level, FALSE)
  expect_equal(c(interval$lower, interval$upper),
    c(.491707638924186405, .508292361075813595), tolerance = 2e-14)
  expect_true(is.na(interval$note))
})
