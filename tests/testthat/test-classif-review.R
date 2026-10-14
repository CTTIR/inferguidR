# Regression tests for defects found in the phase 10A review
# (specs/reviews/phase10A-classifier.md). Every case below used to end in an
# untyped internal error or in a status that was not true.

rv_no_events <- function(...) {
  classification_spec(data.frame(y = c(0, 0, 0, 0), s = c(0.3, 0.1, 0.5, 0.9)), truth = "y",
                      positive_label = 1, scores = "s", prevalence = 0.05, ...)
}

test_that("a scenario basis without events gives typed results, never an internal error", {
  spec <- rv_no_events(cutoff = 0.4)
  res <- evaluate_cutoff(spec)
  expect_identical(res$status$code, "ok")
  expect_identical(res$models[[1]]$scenario$status, "undefined")
  # sweeps and panels on the scenario basis: every metric of the undefined counts is typed
  sw <- cutoff_sweep(res, metrics = c("tpr", "ppv", "f1"), basis = "scenario")
  expect_true(all(is.na(sw$tp)))
  expect_true(all(sw$ppv_status == "undefined"))
  for (p in c("metric_vs_cutoff", "metric_pair", "sweep_table")) {
    expect_no_error(classification_plot_data(res, p, basis = "scenario"))
  }
  expect_identical(th_condition_code(classification_plot_data(res, "confusion_table", basis = "scenario")),
                   "infeasible")
  expect_error(classification_plot_data(res, "confusion_table", basis = "scenario"), "scenario undefined")
  v <- inferguidR:::.ig_cls_metrics_vec(c(NA, 1), c(NA, 1), c(NA, 1), c(NA, 1))
  expect_identical(v$ppv$status, c("undefined", "ok"))
  expect_identical(v$ppv$reason, c("counts_undefined", NA))
  # a selection on that basis is infeasible for the stated reason
  for (obj in c("max_youden", "max_f1")) {
    r <- select_cutoff(spec, obj, prevalence_basis = "scenario")
    expect_identical(r$status$code, "infeasible")
    sel <- classification_table(r, "selection")
    expect_identical(sel$reason_code, "no_events")
    expect_null(r$models[[1]]$operating_point)
  }
  # a manual cutoff with constraints on that basis: the bounds are unmet, not an error
  m <- evaluate_cutoff(update(spec, rule = cutoff_rule(cutoff = 0.4, constraints = list(ppv_min = 0.2),
                                                       prevalence_basis = "scenario")))
  expect_identical(m$status$code, "ok")
  expect_false(m$models[[1]]$selection$constraints_met)
  expect_identical(classification_table(m, "bounds")$best_status, "undefined")
  # the result survives the JSON round trip
  r <- select_cutoff(spec, "max_f1", prevalence_basis = "scenario")
  back <- ig_from_json(ig_to_json(analysis_envelope(classification_evaluated_spec(r), r)))
  expect_identical(back$result$status$code, "infeasible")
})

test_that("probability views refuse a lower-is-positive direction everywhere, as the probability functions do", {
  d <- data.frame(y = c(1, 1, 1, 0, 0, 0, 0, 0), p = c(0.1, 0.2, 0.6, 0.3, 0.7, 0.8, 0.8, 0.9))
  spec <- classification_spec(d, truth = "y", positive_label = 1, scores = "p", cutoff = 0.5,
                              score_scale = "probability", direction = "lower_is_positive")
  r <- evaluate_cutoff(spec, threshold_probability = c(0.2, 0.5))
  sm <- classification_table(r, "summaries")
  prob <- sm[sm$metric_id %in% c("brier_score", "calibration_in_the_large", "calibration_intercept",
                                 "calibration_slope"), ]
  expect_equal(nrow(prob), 4)
  expect_true(all(prob$status == "not_applicable"))
  expect_true(all(prob$reason_code == "probability_direction_reversed"))
  expect_true(all(is.na(prob$value)))
  expect_false(r$models[[1]]$views$calibration$available)
  expect_identical(r$models[[1]]$views$net_benefit$reason_code, "probability_direction_reversed")
  expect_null(r$models[[1]]$probability)
  expect_identical(nrow(classification_table(r, "net_benefit")), 0L)
  # the ROC AUC in the declared direction is still reported
  expect_equal(sm$value[sm$metric_id == "roc_auc"], 14 / 15)
  expect_identical(th_condition_code(net_benefit(spec, 0.2)), "unsupported")
  expect_identical(th_condition_code(calibration_summary(spec)), "unsupported")
  # a requested bootstrap Brier score is not computed for such a model
  boot <- update(spec, uncertainty = classification_uncertainty(bootstrap_replicates = 100, seed = 1,
                                                                bootstrap_metrics = c("tpr", "brier_score")))
  iv <- classification_table(evaluate_cutoff(boot), "intervals")
  expect_false("brier_score" %in% iv$metric_id)
})

test_that("totals that would overflow the metric formulas are refused with a typed condition", {
  expect_identical(th_condition_code(classification_metrics(c(tp = 3e300, fn = 1e300, fp = 2e300, tn = 7e300),
                                                            count_kind = "expected")), "invalid_input")
  expect_identical(th_condition_code(classification_spec(theory = list(auc = 0.8, sigma = 1, n_total = 1e300),
                                                         cutoff = 0, prevalence = 0.1)), "invalid_input")
  expect_identical(th_condition_code(prevalence_scenario(c(tp = 1, fn = 1, fp = 1, tn = 1), 0.1, n_total = 1e200)),
                   "invalid_input")
  d <- data.frame(y = c(1, 0), s = c(2, 1), w = c(1e300, 1))
  expect_identical(th_condition_code(classification_spec(d, truth = "y", positive_label = 1, scores = "s",
                                                         cutoff = 1, weights = "w")), "invalid_input")
  # large but admissible totals stay exact
  m <- classification_metrics(c(tp = 3e15, fn = 1e15, fp = 2e15, tn = 7e15))
  expect_equal(m$value[m$metric_id == "dor"], 10.5)
})

test_that("far tails of a theoretical model never give an internal error or a number outside the double range", {
  # small sigma and AUC near 0.5: the selection scan reaches subnormal expected counts
  spec <- classification_spec(theory = list(auc = 0.501, sigma = 0.2, n_total = 1000), cutoff = 0.5,
                              prevalence = 1e-4)
  for (obj in c("max_youden", "max_f1", "max_mcc")) {
    r <- select_cutoff(spec, obj, include_curves = FALSE)
    expect_true(r$status$code %in% c("ok", "infeasible"), info = obj)
  }
  r <- select_cutoff(spec, "constraint", constraints = list(alert_rate_max = 0.5), optimise = "alerts_per_tp")
  expect_true(r$status$code %in% c("ok", "infeasible"))
  wide <- classification_spec(theory = list(auc = 0.6, sigma = 6, n_total = 1000), cutoff = 0, prevalence = 0.1)
  r <- select_cutoff(wide, "constraint", constraints = list(fpr_max = 1e-300))
  expect_identical(r$status$code, "ok")
  expect_equal(r$models[[1]]$operating_point$cutoff, qnorm(1e-300, lower.tail = FALSE), tolerance = 1e-9)
  # a ratio that overflows is "not_computed", with its reason; exact zeros keep their typed status
  v <- inferguidR:::.ig_cls_metrics_vec(tp = 1e-310, fn = 0.1, fp = 1, tn = 999)
  expect_identical(v$fp_per_tp$status, "not_computed")
  expect_identical(v$fp_per_tp$reason, "outside_double_range")
  expect_true(is.na(v$fp_per_tp$value))
  expect_identical(inferguidR:::.ig_cls_metrics_vec(0, 0.1, 1, 999)$fp_per_tp$status, "positive_infinity")
  # uniformly tiny expected counts: products underflow, but no status is invented
  tiny <- classification_metrics(c(tp = 1e-200, fn = 1e-200, fp = 1e-200, tn = 1e-200), count_kind = "expected")
  expect_equal(tiny$value[tiny$metric_id == "f1"], 0.5)
  for (id in c("dor", "cohen_kappa", "mcc", "lr_pos")) {
    expect_identical(tiny$status[tiny$metric_id == id], "not_computed", info = id)
    expect_identical(tiny$reason_code[tiny$metric_id == id], "outside_double_range", info = id)
  }
  # TPR = FPR = 0.5 here, so ppv_gap_max_prevalence is genuinely undefined
  expect_identical(tiny$reason_code[tiny$metric_id == "ppv_gap_max_prevalence"], "equal_tpr_fpr")
  expect_true(all(tiny$status[tiny$metric_id != "ppv_gap_max_prevalence"] %in% c("ok", "not_computed")))
  # ordinary counts never hit the range rule
  m <- classification_metrics(c(tp = 0, fn = 3, fp = 0, tn = 5))
  expect_false(any(m$status == "not_computed"))
  expect_identical(m$status[m$metric_id == "dor"], "undefined")
})

test_that("theory optima near the scan edge are judged by looking beyond it", {
  # References: mpmath at 40 digits (crosscheck-style), not the engine.
  rv_th <- function(auc, sigma, pi) {
    classification_spec(theory = list(auc = auc, sigma = sigma, n_total = 1000), cutoff = 0.3, prevalence = pi)
  }
  # sigma = 1000: the grid step is 2.5, and the Youden optimum c = 7.0621526672045 (root of
  # phi(c) = phi((c - d) / sigma) / sigma) lies within two steps of the scan edge
  r <- select_cutoff(rv_th(1e-9, 1000, 0.5), "max_youden")
  expect_identical(r$status$code, "ok")
  expect_equal(r$models[[1]]$operating_point$cutoff, 7.06215266720447, tolerance = 1e-9)
  expect_equal(r$models[[1]]$selection$objective_value, 9.56593835098879e-10, tolerance = 1e-9)
  r <- select_cutoff(rv_th(1e-9, 1000, 0.5), "constraint", constraints = list(fpr_max = 1e-6))
  expect_equal(r$models[[1]]$operating_point$cutoff, 4.7534243088229, tolerance = 1e-10)
  # a plateau that reaches alert all: the limiting rule attains the optimum
  for (s in c(0.001, 0.05)) {
    r <- select_cutoff(rv_th(1e-9, s, 1 - 1e-9), "max_f1")
    expect_identical(r$status$code, "ok")
    expect_identical(r$models[[1]]$operating_point$cutoff_kind, "alert_all")
  }
  # an objective that keeps improving is still refused (PPV -> 1 as the cutoff grows, sigma > 1)
  r <- select_cutoff(rv_th(0.8, 1.5, 0.1), "constraint", optimise = "ppv", constraints = list(fpr_max = 0.5))
  expect_identical(r$models[[1]]$selection$reason_code, "supremum_not_attained")
  # extreme parameters integrate over the standardised event score
  sm <- classification_table(evaluate_cutoff(rv_th(1 - 1e-9, 1000, 0.5),
                                             partial_auc = list(focus = "tpr", range = c(0.5, 1))), "summaries")
  expect_true(all(sm$status[sm$metric_id %in% c("partial_auc", "average_precision")] == "ok"))
  expect_equal(sm$value[sm$metric_id == "average_precision"], 0.999999999500485, tolerance = 1e-11)
  expect_equal(sm$value[sm$metric_id == "partial_auc" & sm$variant == "raw"], 0.499999998985358, tolerance = 1e-9)
})

test_that("selection on a theoretical model needs a prevalence strictly between 0 and 1", {
  for (pv in c(0, 1)) {
    th <- classification_spec(theory = list(auc = 0.8, sigma = 1), cutoff = 0, prevalence = pv)
    expect_identical(evaluate_cutoff(th)$status$code, "ok")
    expect_identical(th_condition_code(select_cutoff(th, "max_f1")), "invalid_input")
    expect_identical(th_condition_code(select_cutoff(th, "cost_utility", costs = c(fn = 1, fp = 1))),
                     "invalid_input")
  }
  th <- classification_spec(theory = list(auc = 0.8, sigma = 1), cutoff = 0, prevalence = 0.3)
  expect_identical(th_condition_code(update(th, prevalence = 1, rule = cutoff_rule(selection = "max_youden",
                                                                                    selection_data = "predeclared_rule"))),
                   "invalid_input")
})
