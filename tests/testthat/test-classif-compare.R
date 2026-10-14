# Comparison with pROC and yardstick (Suggests only; decisions D13 and D17).
# Conventions are aligned explicitly before comparing (engine-decisions-aux
# section 3): pROC gets `levels = c(<non-event>, <event>)` and an explicit
# `direction` (never "auto"), and its "youden" is J + 1; yardstick gets an
# explicit `event_level`. These packages are comparators, not independent
# references.

cp_data <- function(seed = 20261003, n0 = 60, n1 = 40) {
  set.seed(seed)
  y <- c(rep(0, n0), rep(1, n1))
  s <- round(c(rnorm(n0), rnorm(n1, 1, 1.3)), 1)
  data.frame(y = y, s = s, s2 = round(s + rnorm(n0 + n1, 0, 0.7), 1), p = plogis(s - 0.3))
}

test_that("ROC AUC, partial AUC (both foci, McClish) and DeLong agree with pROC", {
  skip_if_not_installed("pROC")
  for (seed in c(20261003, 7, 11)) {
    d <- cp_data(seed)
    spec <- classification_spec(d, truth = "y", positive_label = 1, scores = "s", cutoff = 0.5)
    r <- pROC::roc(d$y, d$s, levels = c(0, 1), direction = "<", quiet = TRUE)
    a <- classification_auc(spec, partial_auc = list(focus = "fpr", range = c(0, 0.2)), level = 0.95)
    expect_equal(a$value[a$metric_id == "roc_auc"], as.numeric(pROC::auc(r)), tolerance = 1e-12)
    expect_equal(a$value[a$variant %in% "raw"],
                 as.numeric(pROC::auc(r, partial.auc = c(1, 0.8), partial.auc.focus = "specificity")),
                 tolerance = 1e-12)
    expect_equal(a$value[a$variant %in% "mcclish"],
                 as.numeric(pROC::auc(r, partial.auc = c(1, 0.8), partial.auc.focus = "specificity",
                                      partial.auc.correct = TRUE)), tolerance = 1e-12)
    a2 <- classification_auc(spec, partial_auc = list(focus = "tpr", range = c(0.8, 1)))
    expect_equal(a2$value[a2$variant %in% "raw"],
                 as.numeric(pROC::auc(r, partial.auc = c(1, 0.8), partial.auc.focus = "sensitivity")),
                 tolerance = 1e-12)
    ci <- as.numeric(pROC::ci.auc(r, method = "delong", conf.level = 0.95))
    expect_equal(c(a$ci_lower[1], a$ci_upper[1]), ci[c(1, 3)], tolerance = 1e-10)
  }
  # a lower-is-positive score: pROC direction ">" means positive if predictor <= threshold
  d <- cp_data(5)
  spec <- classification_spec(d, truth = "y", positive_label = 1, scores = "s", cutoff = 0.5,
                              direction = "lower_is_positive")
  r <- pROC::roc(d$y, d$s, levels = c(0, 1), direction = ">", quiet = TRUE)
  expect_equal(classification_auc(spec)$value[1], as.numeric(pROC::auc(r)), tolerance = 1e-12)
})

test_that("the paired DeLong comparison agrees with pROC::roc.test", {
  skip_if_not_installed("pROC")
  d <- cp_data()
  spec <- classification_spec(d, truth = "y", positive_label = 1, scores = c("s", "s2"), cutoff = 0.5)
  cmp <- compare_classifiers(spec)
  r1 <- pROC::roc(d$y, d$s, levels = c(0, 1), direction = "<", quiet = TRUE)
  r2 <- pROC::roc(d$y, d$s2, levels = c(0, 1), direction = "<", quiet = TRUE)
  rt <- pROC::roc.test(r1, r2, method = "delong", paired = TRUE)
  expect_equal(cmp$z, unname(rt$statistic), tolerance = 1e-10)
  expect_equal(cmp$p_value, rt$p.value, tolerance = 1e-10)
  expect_equal(c(cmp$auc_a, cmp$auc_b), unname(rt$estimate), tolerance = 1e-12)
})

test_that("counts at a cutoff and Youden's J agree with pROC::coords after mapping conventions", {
  skip_if_not_installed("pROC")
  d <- cp_data()
  r <- pROC::roc(d$y, d$s, levels = c(0, 1), direction = "<", quiet = TRUE)
  for (cc in c(-0.5, 0.5, 1.2)) {
    res <- evaluate_cutoff(classification_spec(d, truth = "y", positive_label = 1, scores = "s", cutoff = cc))
    op <- res$models[[1]]$operating_point
    co <- pROC::coords(r, x = cc, input = "threshold", ret = c("tp", "fn", "fp", "tn", "youden"),
                       transpose = FALSE)
    expect_equal(c(op$tp, op$fn, op$fp, op$tn), unname(unlist(co[1, c("tp", "fn", "fp", "tn")])))
    m <- classification_table(res, "metrics")
    expect_equal(m$value[m$metric_id == "youden_j"], co$youden - 1, tolerance = 1e-12)
  }
})

test_that("summaries and class metrics agree with yardstick when the event level is explicit", {
  skip_if_not_installed("yardstick")
  d <- cp_data()
  f <- factor(ifelse(d$y == 1, "event", "none"), levels = c("event", "none"))
  spec <- classification_spec(d, truth = "y", positive_label = 1, scores = "p", cutoff = plogis(0.2),
                              score_scale = "probability")
  res <- evaluate_cutoff(spec)
  sm <- classification_table(res, "summaries")
  v <- function(id) sm$value[sm$metric_id == id]
  expect_equal(v("roc_auc"), yardstick::roc_auc_vec(f, d$p, event_level = "first"), tolerance = 1e-12)
  expect_equal(v("average_precision"), yardstick::average_precision_vec(f, d$p, event_level = "first"),
               tolerance = 1e-12)
  expect_equal(v("pr_area_trapezoid"), yardstick::pr_auc_vec(f, d$p, event_level = "first"), tolerance = 1e-12)
  expect_equal(v("brier_score"), yardstick::brier_class_vec(f, d$p, event_level = "first"), tolerance = 1e-12)
  # the wrong event level would silently give 1 - AUC; the engine never guesses it
  expect_equal(1 - v("roc_auc"), yardstick::roc_auc_vec(f, d$p, event_level = "second"), tolerance = 1e-12)
  est <- factor(ifelse(d$p >= plogis(0.2), "event", "none"), levels = c("event", "none"))
  m <- classification_table(res, "metrics")
  mv <- function(id) m$value[m$metric_id == id]
  ys <- c(tpr = yardstick::sens_vec(f, est, event_level = "first"),
          tnr = yardstick::spec_vec(f, est, event_level = "first"),
          ppv = yardstick::ppv_vec(f, est, event_level = "first"),
          npv = yardstick::npv_vec(f, est, event_level = "first"),
          f1 = yardstick::f_meas_vec(f, est, event_level = "first"),
          f2 = yardstick::f_meas_vec(f, est, beta = 2, event_level = "first"),
          f05 = yardstick::f_meas_vec(f, est, beta = 0.5, event_level = "first"),
          youden_j = yardstick::j_index_vec(f, est, event_level = "first"),
          balanced_accuracy = yardstick::bal_accuracy_vec(f, est, event_level = "first"),
          mcc = yardstick::mcc_vec(f, est, event_level = "first"),
          cohen_kappa = yardstick::kap_vec(f, est),
          markedness = yardstick::markedness_vec(f, est, event_level = "first"),
          accuracy = yardstick::accuracy_vec(f, est),
          alert_rate = yardstick::detection_prevalence_vec(f, est, event_level = "first"),
          fpr = yardstick::fall_out_vec(f, est, event_level = "first"),
          fnr = yardstick::miss_rate_vec(f, est, event_level = "first"))
  for (id in names(ys)) expect_equal(mv(id), unname(ys[[id]]), tolerance = 1e-12, info = id)
  # Bayes PPV at an assumed prevalence
  pv <- predictive_values_at_prevalence(mv("tpr"), mv("fpr"), 0.05)
  expect_equal(pv$ppv, yardstick::ppv_vec(f, est, prevalence = 0.05, event_level = "first"), tolerance = 1e-12)
})

test_that("frequency weights agree with yardstick case weights", {
  skip_if_not_installed("yardstick")
  d <- cp_data(3)
  d$w <- rep(c(1, 2, 3, 0), length.out = nrow(d))
  f <- factor(ifelse(d$y == 1, "event", "none"), levels = c("event", "none"))
  cw <- d$w   # numeric case weights (equal to hardhat::frequency_weights() for these metrics)
  spec <- classification_spec(d, truth = "y", positive_label = 1, scores = "p", cutoff = 0.5,
                              score_scale = "probability", weights = "w")
  sm <- classification_table(evaluate_cutoff(spec), "summaries")
  expect_equal(sm$value[sm$metric_id == "roc_auc"],
               yardstick::roc_auc_vec(f, d$p, event_level = "first", case_weights = cw), tolerance = 1e-12)
  expect_equal(sm$value[sm$metric_id == "average_precision"],
               yardstick::average_precision_vec(f, d$p, event_level = "first", case_weights = cw), tolerance = 1e-12)
  # yardstick 1.4.0 transforms Brier case weights by a softmax (brier_ind), so its
  # weighted Brier is not frequency-weight replication; compare with replicated rows
  idx <- rep(seq_len(nrow(d)), d$w)
  expect_equal(sm$value[sm$metric_id == "brier_score"],
               yardstick::brier_class_vec(f[idx], d$p[idx], event_level = "first"), tolerance = 1e-12)
})

test_that("where the comparators return NaN or NA, the engine returns a typed status", {
  skip_if_not_installed("yardstick")
  f <- factor(c("event", "none", "event", "none"), levels = c("event", "none"))
  est <- factor(rep("none", 4), levels = c("event", "none"))
  ppv_y <- suppressWarnings(yardstick::ppv_vec(f, est, event_level = "first"))
  expect_true(is.na(ppv_y))
  m <- classification_metrics(c(tp = 0, fn = 2, fp = 0, tn = 2))
  expect_identical(m$status[m$metric_id == "ppv"], "undefined")
  expect_identical(m$reason_code[m$metric_id == "ppv"], "no_alerts")
})
