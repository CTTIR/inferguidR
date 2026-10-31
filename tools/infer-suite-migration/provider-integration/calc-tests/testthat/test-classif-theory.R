# Binormal theory model: S0 ~ N(0, 1), S1 ~ N(d, sigma^2), d = qnorm(AUC) sqrt(1 + sigma^2).
# Reference values below come from numerical integration of the densities
# (stats::integrate over dnorm), not from the engine's analytic tails.

th_spec <- function(auc = 0.8, sigma = 1, cutoff = 1, prevalence = 0.1, n_total = 1000) {
  classification_spec(theory = list(auc = auc, sigma = sigma, n_total = n_total),
                      cutoff = cutoff, prevalence = prevalence)
}

upper_tail_integral <- function(c, mean, sd) {
  stats::integrate(function(z) dnorm(z, mean, sd), c, Inf, rel.tol = 1e-12)$value
}

test_that("the parameterisation reproduces the AUC by double integration", {
  for (auc in c(0.6, 0.8, 0.95)) for (s in c(0.5, 1, 2)) {
    m <- binormal_model(auc, s)
    inner <- function(z) vapply(z, function(zz) upper_tail_integral(zz, m$d, s), 0)
    num <- stats::integrate(function(z) dnorm(z) * inner(z), -Inf, Inf, rel.tol = 1e-10)$value
    expect_equal(num, auc, tolerance = 1e-8)
    expect_equal(m$auc_check, auc, tolerance = 1e-14)
  }
})

test_that("rates at real cutoffs are analytic tails, including far tails", {
  for (s in c(0.5, 1.5)) {
    m <- binormal_model(0.8, s)
    for (cc in c(-2, 0, 1.234, 4, 6)) {
      r <- binormal_rates(cc, 0.8, s)
      expect_equal(r$tpr, upper_tail_integral(cc, m$d, s), tolerance = 1e-9)
      expect_equal(r$fpr, upper_tail_integral(cc, 0, 1), tolerance = 1e-9)
    }
  }
  # relative accuracy deep in the tail (the source's erf approximation lost 3.6e-3 at z = 6)
  expect_equal(binormal_rates(6, 0.8, 1)$fpr, 9.865876450376981e-10, tolerance = 1e-12)
})

test_that("expected counts scale with the model N and rates do not", {
  r1 <- evaluate_cutoff(th_spec(n_total = 1000))
  r2 <- evaluate_cutoff(th_spec(n_total = 50000))
  m1 <- classification_table(r1, "metrics")
  m2 <- classification_table(r2, "metrics")
  rates <- c("tpr", "fpr", "ppv", "npv", "mcc", "lr_pos")
  expect_equal(m1$value[m1$metric_id %in% rates], m2$value[m2$metric_id %in% rates], tolerance = 1e-12)
  expect_equal(m2$value[m2$metric_id == "tp"], 50 * m1$value[m1$metric_id == "tp"], tolerance = 1e-12)
  expect_identical(unique(m1$estimate_kind), "expected_theoretical")
  expect_true(all(m1$note[m1$unit == "count"] == "expected count, not observed units"))
})

test_that("changing AUC or sigma keeps the numerical cutoff", {
  s <- th_spec(auc = 0.8, sigma = 1, cutoff = 1.25)
  s2 <- update(s, theory = list(sigma = 1.5))
  s3 <- update(s2, theory = list(auc = 0.9))
  expect_identical(s2$operating_rule$score_cutoff, 1.25)
  expect_identical(s3$operating_rule$score_cutoff, 1.25)
  r3 <- evaluate_cutoff(s3)
  expect_equal(r3$models[[1]]$operating_point$cutoff, 1.25)
  expect_equal(r3$models[[1]]$operating_point$fp, 900 * pnorm(1.25, lower.tail = FALSE), tolerance = 1e-12)
})

test_that("theoretical partial AUC and average precision equal numerical integration", {
  spec <- th_spec(auc = 0.8, sigma = 1.5, prevalence = 0.1)
  m <- binormal_model(0.8, 1.5)
  res <- evaluate_cutoff(spec, partial_auc = list(focus = "fpr", range = c(0, 0.2)))
  sm <- classification_table(res, "summaries")
  # partial AUC: integral of TPR over FPR in [0, 0.2], written over FPR directly
  tpr_at_fpr <- function(f) pnorm((qnorm(f, lower.tail = FALSE) - m$d) / 1.5, lower.tail = FALSE)
  ref <- stats::integrate(tpr_at_fpr, 0, 0.2, rel.tol = 1e-11)$value
  expect_equal(sm$value[sm$metric_id == "partial_auc" & sm$variant == "raw"], ref, tolerance = 1e-8)
  # average precision: Riemann-Stieltjes sum of PPV over recall on a fine cutoff grid
  cc <- seq(-12, 14, length.out = 400001)
  tpr <- pnorm((cc - m$d) / 1.5, lower.tail = FALSE)
  fpr <- pnorm(cc, lower.tail = FALSE)
  ppv <- 0.1 * tpr / (0.1 * tpr + 0.9 * fpr)
  dr <- -diff(tpr)
  ap_ref <- sum(dr * (head(ppv, -1) + tail(ppv, -1)) / 2)
  expect_equal(sm$value[sm$metric_id == "average_precision"], ap_ref, tolerance = 1e-7)
  expect_identical(sm$status[sm$metric_id == "pr_area_trapezoid"], "not_applicable")
  expect_identical(sm$value[sm$metric_id == "roc_auc"], 0.8)
  r2 <- evaluate_cutoff(spec, partial_auc = list(focus = "tpr", range = c(0.8, 1)))
  s2 <- classification_table(r2, "summaries")
  fpr_at_tpr <- function(t) pnorm(m$d + 1.5 * qnorm(t, lower.tail = FALSE), lower.tail = FALSE)
  ref2 <- stats::integrate(function(t) 1 - fpr_at_tpr(t), 0.8, 1, rel.tol = 1e-11)$value
  expect_equal(s2$value[s2$metric_id == "partial_auc" & s2$variant == "raw"], ref2, tolerance = 1e-8)
})

test_that("theoretical results carry no sampling interval", {
  res <- evaluate_cutoff(th_spec())
  iv <- classification_table(res, "intervals")
  expect_true(nrow(iv) > 0)
  expect_true(all(iv$status == "not_applicable"))
  expect_true(all(iv$reason_code == "theoretical_no_interval"))
  expect_true(all(is.na(iv$lower)))
  expect_identical(res$evaluation$label, "theoretical_model")
  expect_false(grepl("observed", res$evaluation$statement))
})

test_that("the Youden optimum is analytic (d/2 for equal variances) and matches a fine grid", {
  for (s in c(1, 0.6, 1.5, 2)) {
    spec <- th_spec(auc = 0.8, sigma = s)
    m <- binormal_model(0.8, s)
    res <- select_cutoff(spec, "max_youden")
    cut <- res$models[[1]]$operating_point$cutoff
    if (s == 1) expect_equal(cut, m$d / 2, tolerance = 1e-12)
    g <- seq(-6, 8, by = 1e-4)
    J <- pnorm((g - m$d) / s, lower.tail = FALSE) - pnorm(g, lower.tail = FALSE)
    expect_lt(abs(cut - g[which.max(J)]), 2e-4)
    expect_gte(res$models[[1]]$selection$objective_value, max(J) - 1e-12)
  }
  # the source grid gave 0.9773 at sigma 1.5 and AUC 0.8; the optimum is 0.9713
  res <- select_cutoff(th_spec(auc = 0.8, sigma = 1.5), "max_youden")
  expect_equal(res$models[[1]]$operating_point$cutoff, 0.9712919, tolerance = 1e-6)
})

test_that("numerical theory optima (F1, MCC, cost, constraints) match a brute-force grid", {
  spec <- th_spec(auc = 0.85, sigma = 1.3, prevalence = 0.2)
  m <- binormal_model(0.85, 1.3)
  g <- seq(-6, 9, by = 2e-4)
  tpr <- pnorm((g - m$d) / 1.3, lower.tail = FALSE)
  fpr <- pnorm(g, lower.tail = FALSE)
  tp <- 200 * tpr; fp <- 800 * fpr; fn <- 200 - tp; tn <- 800 - fp
  f1 <- 2 * tp / (2 * tp + fp + fn)
  mcc <- (tp * tn - fp * fn) / sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
  cost <- (5 * fn + fp) / 1000
  check <- function(res, target) {
    cut <- res$models[[1]]$operating_point$cutoff
    expect_lt(abs(cut - target), 5e-4)
  }
  check(select_cutoff(spec, "max_f1"), g[which.max(f1)])
  check(select_cutoff(spec, "max_mcc"), g[which.max(mcc)])
  check(select_cutoff(spec, "cost_utility", costs = c(fn = 5, fp = 1)), g[which.min(cost)])
  # highest TPR with FPR <= 0.05: the boundary qnorm(0.95)
  res <- select_cutoff(spec, "constraint", constraints = list(fpr_max = 0.05))
  expect_equal(res$models[[1]]$operating_point$cutoff, qnorm(0.95), tolerance = 1e-9)
  # capacity: alert rate at most 10 %
  res <- select_cutoff(spec, "capacity", capacity = list(max_alert_rate = 0.1))
  op <- res$models[[1]]$operating_point
  expect_equal(op$n_alerts / op$n_total, 0.1, tolerance = 1e-9)
})

test_that("theory selection reports limits, endpoints and unattainable suprema honestly", {
  # Costs that make alerting everyone optimal: with sigma = 2 the likelihood
  # ratio f1/f0 never falls below 0.40 > fp/fn = 0.01, so the expected cost
  # increases with the cutoff and the limiting rule "alert all" is optimal.
  # (With sigma = 1 the likelihood ratio tends to 0 in the left tail, so an
  # interior cutoff near -5.8 is very slightly cheaper than alerting all.)
  spec <- th_spec(auc = 0.7, sigma = 2, prevalence = 0.5)
  res <- select_cutoff(spec, "cost_utility", costs = c(fn = 100, fp = 1))
  op <- res$models[[1]]$operating_point
  expect_identical(op$cutoff_kind, "alert_all")
  expect_true(is.na(op$cutoff))
  expect_equal(op$tp, 500)
  # maximising PPV with sigma > 1: PPV rises towards 1 as the cutoff grows, never attained
  res <- select_cutoff(th_spec(auc = 0.8, sigma = 1.5), "constraint", optimise = "ppv",
                       constraints = list(fpr_max = 0.5))
  expect_identical(res$status$code, "infeasible")
  expect_identical(res$models[[1]]$selection$reason_code, "supremum_not_attained")
  expect_null(res$models[[1]]$operating_point)
  # infeasible constraints
  res <- select_cutoff(th_spec(auc = 0.7), "constraint", constraints = list(tpr_min = 0.95, fpr_max = 0.05))
  expect_identical(res$status$code, "infeasible")
  b <- classification_table(res, "bounds")
  expect_equal(nrow(b), 2)
  expect_true(all(b$n_satisfied > 0))
})

test_that("theory specifications refuse data-only settings", {
  expect_identical(th_condition_code(classification_spec(theory = list(auc = 0.8, sigma = 1), cutoff = 0)), "invalid_input")
  expect_identical(th_condition_code(classification_spec(theory = list(auc = 0.8, sigma = 1), cutoff = 0,
                                                         prevalence = 0.1, evaluation_split = "same_data")), "invalid_input")
  expect_identical(th_condition_code(classification_spec(theory = list(auc = 0.8, sigma = 1), prevalence = 0.1,
                                                         rule = cutoff_rule(cutoff = 0, selection_data = "same_data"))),
                   "invalid_input")
  expect_identical(th_condition_code(classification_spec(theory = list(auc = 1, sigma = 1), cutoff = 0, prevalence = 0.1)),
                   "invalid_input")
  expect_identical(th_condition_code(classification_spec(theory = list(auc = 0.8, sigma = 1), cutoff = 0, prevalence = 0.1,
                                                         direction = "lower_is_positive")), "invalid_input")
  s <- th_spec()
  expect_identical(s$operating_rule$selection_data, "predeclared_rule")
  expect_true("selection_data_theoretical" %in% s$operating_rule$defaults_applied)
  expect_identical(s$operating_rule$cutoff_scale, "z_standardised_to_non_events")
  expect_identical(th_condition_code(cutoff_sweep(classification_spec(counts = c(tp = 1, fn = 1, fp = 1, tn = 1)))),
                   "unsupported")
})

test_that("the theoretical curve keeps the exact operating point and both limiting rules", {
  res <- evaluate_cutoff(th_spec(cutoff = 1.2345))
  sw <- cutoff_sweep(res)
  expect_true(any(sw$selected))
  expect_equal(sw$cutoff[sw$selected], 1.2345)
  expect_true(all(c("alert_none", "alert_all") %in% sw$cutoff_kind))
  expect_identical(sw$count_kind[1], "expected")
  expect_equal(sw$tp[sw$cutoff_kind == "alert_all"], 100)
  expect_equal(sw$fp[sw$cutoff_kind == "alert_none"], 0)
})
