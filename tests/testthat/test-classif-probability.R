# Probability views: Brier score, calibration and net benefit, checked against
# direct arithmetic and an own Newton fit of the logistic recalibration
# models (the engine uses stats::glm.fit).

pb_data <- function(seed = 3, n = 150, weighted = FALSE) {
  set.seed(seed)
  lp <- rnorm(n, -0.4, 1.1)
  y <- rbinom(n, 1, plogis(0.25 + 0.7 * lp))
  p <- pmin(pmax(round(plogis(lp), 2), 0.01), 0.99)
  data.frame(y = y, p = p, w = if (weighted) sample(1:3, n, replace = TRUE) else 1)
}

pb_newton <- function(p, y, w) {
  x <- log(p / (1 - p))
  mu <- function(e) 1 / (1 + exp(-e))
  th <- c(0, 1)
  for (i in 1:50) {
    m <- mu(th[1] + th[2] * x); v <- w * m * (1 - m)
    g <- c(sum(w * (y - m)), sum(w * (y - m) * x))
    H <- -matrix(c(sum(v), sum(v * x), sum(v * x), sum(v * x^2)), 2)
    th <- th - solve(H, g)
  }
  a <- 0
  for (i in 1:50) {
    m <- mu(a + x)
    a <- a + sum(w * (y - m)) / sum(w * m * (1 - m))
  }
  c(intercept = a, slope = th[2])
}

test_that("Brier score and calibration equal direct computation and an own maximum-likelihood fit", {
  for (wt in c(FALSE, TRUE)) {
    d <- pb_data(weighted = wt)
    spec <- classification_spec(d, truth = "y", positive_label = 1, scores = "p", cutoff = 0.5,
                                score_scale = "probability", weights = if (wt) "w" else NULL)
    cs <- calibration_summary(spec)
    s <- cs$summaries
    v <- function(id) s$value[s$metric_id == id]
    expect_equal(v("brier_score"), sum(d$w * (d$p - d$y)^2) / sum(d$w), tolerance = 1e-14)
    expect_equal(v("calibration_in_the_large"), (sum(d$w * d$p) - sum(d$w * d$y)) / sum(d$w), tolerance = 1e-14)
    ref <- pb_newton(d$p, d$y, d$w)
    expect_equal(v("calibration_intercept"), ref[["intercept"]], tolerance = 1e-9)
    expect_equal(v("calibration_slope"), ref[["slope"]], tolerance = 1e-9)
    # the grouped table covers every unit once; equal probabilities stay in one group
    tb <- cs$table
    expect_equal(sum(tb$n_units), sum(d$w))
    expect_equal(sum(tb$events), sum(d$w * d$y))
    expect_true(all(tb$max_predicted[-nrow(tb)] < tb$min_predicted[-1]))
  }
  # weights act as replication
  d <- pb_data(weighted = TRUE)
  long <- d[rep(seq_len(nrow(d)), d$w), c("y", "p")]
  a <- calibration_summary(classification_spec(d, truth = "y", positive_label = 1, scores = "p", cutoff = 0.5,
                                               score_scale = "probability", weights = "w"))$summaries
  b <- calibration_summary(classification_spec(long, truth = "y", positive_label = 1, scores = "p", cutoff = 0.5,
                                               score_scale = "probability"))$summaries
  expect_equal(a$value, b$value, tolerance = 1e-9)
})

test_that("net benefit uses explicit threshold probabilities and shows alert-all and alert-none references", {
  d <- pb_data()
  spec <- classification_spec(d, truth = "y", positive_label = 1, scores = "p", cutoff = 0.5,
                              score_scale = "probability")
  thr <- c(0.1, 0.25, 0.5)
  nb <- net_benefit(spec, thr)
  n <- nrow(d); P <- sum(d$y)
  for (i in seq_along(thr)) {
    t <- thr[i]
    a <- d$p >= t
    expect_equal(nb$net_benefit[i], sum(a & d$y == 1) / n - sum(a & d$y == 0) / n * t / (1 - t), tolerance = 1e-14)
    expect_equal(nb$net_benefit_alert_all[i], P / n - (n - P) / n * t / (1 - t), tolerance = 1e-14)
    expect_equal(nb$harm_weight[i], t / (1 - t))
  }
  expect_true(all(nb$net_benefit_alert_none == 0))
  expect_identical(nb$rule[1], "alert if p >= 0.1")
  # the comparator ">" is honoured
  spec_gt <- classification_spec(d, truth = "y", positive_label = 1, scores = "p", cutoff = 0.5,
                                 score_scale = "probability", comparator = ">")
  nb_gt <- net_benefit(spec_gt, 0.25)
  a <- d$p > 0.25
  expect_equal(nb_gt$net_benefit, sum(a & d$y == 1) / n - sum(a & d$y == 0) / n / 3, tolerance = 1e-14)
  expect_identical(nb_gt$rule, "alert if p > 0.25")
  # inside evaluate_cutoff(): stored per model, view available, calibration caveat listed
  r <- evaluate_cutoff(spec, threshold_probability = thr)
  expect_equal(classification_table(r, "net_benefit")$net_benefit, nb$net_benefit)
  expect_true(r$models[[1]]$views$net_benefit$available)
  expect_true("net_benefit_calibration" %in% r$limitations$code)
  r0 <- evaluate_cutoff(spec)
  expect_identical(r0$models[[1]]$views$net_benefit$reason_code, "threshold_probability_missing")
  expect_identical(th_condition_code(net_benefit(spec, c(0, 0.2))), "invalid_input")
})

test_that("the score cutoff and the threshold probability are different settings", {
  d <- pb_data()
  spec <- classification_spec(d, truth = "y", positive_label = 1, scores = "p", cutoff = 0.7,
                              score_scale = "probability")
  r <- evaluate_cutoff(spec, threshold_probability = 0.2)
  nb <- classification_table(r, "net_benefit")
  # net benefit is evaluated at p >= 0.2, not at the operating cutoff 0.7
  expect_equal(nb$tp, sum(d$p >= 0.2 & d$y == 1))
  expect_false(nb$tp == r$models[[1]]$operating_point$tp)
})

test_that("calibration and net benefit are unavailable for raw scores and theory", {
  raw <- classification_spec(data.frame(y = c(1, 0, 1, 0), s = c(4, 1, 3, 2)), truth = "y", positive_label = 1,
                             scores = "s", cutoff = 2)
  expect_identical(th_condition_code(calibration_summary(raw)), "unsupported")
  th <- classification_examples("binormal_equal_variance")$spec
  r <- evaluate_cutoff(th)
  expect_identical(r$models[[1]]$views$calibration$reason_code, "not_probability_scale")
  expect_identical(th_condition_code(calibration_summary(th)), "unsupported")
})
