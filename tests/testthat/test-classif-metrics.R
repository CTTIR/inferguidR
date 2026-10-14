# Metric dictionary, identities, typed statuses and prevalence relations.

cm_value <- function(m, id) m$value[m$metric_id == id]
cm_status <- function(m, id) m$status[m$metric_id == id]

test_that("the dictionary documents every computed metric completely and bilingually", {
  doc <- inferguidR:::.ig_cls_dictionary()
  ids <- vapply(doc$metrics, `[[`, "", "id")
  expect_false(anyDuplicated(ids) > 0)
  expect_true(all(grepl("^[a-z][a-z0-9_]*$", ids)))
  need <- c("id", "scope", "group", "label", "formula", "denominator", "unit", "domain",
            "prevalence_dependence", "preferred_direction", "undefined", "source_feature_ids")
  for (e in doc$metrics) {
    expect_true(all(need %in% names(e)), info = e$id)
    expect_true(all(c("short", "name") %in% names(e$label$en)), info = e$id)
    expect_true(all(c("short", "name") %in% names(e$label$de)), info = e$id)
    expect_true(e$group %in% names(doc$groups), info = e$id)
    expect_true(e$preferred_direction %in% c("higher", "lower", "none"), info = e$id)
    for (u in e$undefined) {
      expect_true(u$status %in% names(doc$statuses), info = e$id)
      expect_true(u$reason_code %in% names(doc$reason_codes), info = paste(e$id, u$reason_code))
    }
    if (!is.null(e$denominator)) expect_true(all(c("expression", "en", "de") %in% names(e$denominator)))
  }
  for (rc in doc$reason_codes) expect_true(all(c("en", "de") %in% names(rc)))
  # every operating-point metric is computed by the engine, and nothing else
  computed <- names(inferguidR:::.ig_cls_metrics_vec(1, 1, 1, 1))
  expect_setequal(computed, ids[vapply(doc$metrics, function(e) e$scope == "operating_point", TRUE)])
  # the required families of CLASSIFICATION_CUTOFF are present
  expect_true(all(c("tp", "fn", "fp", "tn", "n_total", "n_events", "n_non_events", "n_alerts",
                    "n_non_alerts", "prevalence", "alert_rate", "tpr", "fnr", "tnr", "fpr", "ppv",
                    "fdr", "npv", "for_rate", "accuracy", "balanced_accuracy", "youden_j",
                    "markedness", "f1", "f05", "f2", "fowlkes_mallows", "jaccard", "mcc",
                    "cohen_kappa", "lr_pos", "lr_neg", "dor", "alerts_per_tp", "fp_per_tp",
                    "fp_per_1000_non_events", "missed_per_1000_total", "roc_auc", "partial_auc",
                    "average_precision", "pr_area_trapezoid", "brier_score", "calibration_slope",
                    "net_benefit", "ppv_gap_max_prevalence") %in% ids))
  # the misleading source label is gone; the renamed quantity has no preferred direction
  gap <- doc$metrics[[match("ppv_gap_max_prevalence", ids)]]
  expect_false(grepl("fastest", gap$label$en$name))
  expect_identical(gap$preferred_direction, "none")
  # denominators appear in the labels of the per-1,000 metrics (audit L-12)
  d <- metric_dictionary()
  expect_match(d$short[d$id == "fp_per_1000_non_events"], "non-events")
  expect_match(d$short[d$id == "missed_per_1000_total"], "units")
  dde <- metric_dictionary("de")
  expect_match(dde$name[dde$id == "ppv"], "Vorhersagewert")
  expect_true(all(c("TP + FP (alerts)") %in% d$denominator[d$id == "ppv"]))
})

test_that("prevalence-dependence flags include the cells and prevalence (audit C-02)", {
  d <- metric_dictionary()
  dep <- d$id[d$prevalence_dependent]
  expect_true(all(c("tp", "fn", "fp", "tn", "prevalence", "ppv", "npv", "f1", "mcc", "accuracy",
                    "alert_rate", "missed_per_1000_total") %in% dep))
  expect_false(any(c("tpr", "fpr", "tnr", "fnr", "youden_j", "lr_pos", "lr_neg", "dor",
                     "fp_per_1000_non_events", "balanced_accuracy", "roc_auc") %in% dep))
})

test_that("complement and sum identities hold", {
  for (tb in list(c(40, 10, 20, 130), c(3, 7, 11, 2), c(25.5, 4.5, 12.25, 57.75))) {
    m <- classification_metrics(c(tp = tb[1], fn = tb[2], fp = tb[3], tn = tb[4]), count_kind = "expected")
    v <- function(id) cm_value(m, id)
    expect_equal(v("tpr") + v("fnr"), 1)
    expect_equal(v("tnr") + v("fpr"), 1)
    expect_equal(v("ppv") + v("fdr"), 1)
    expect_equal(v("npv") + v("for_rate"), 1)
    expect_equal(v("n_events") + v("n_non_events"), v("n_total"))
    expect_equal(v("n_alerts") + v("n_non_alerts"), v("n_total"))
    expect_equal(v("youden_j"), v("tpr") - v("fpr"))
    expect_equal(v("balanced_accuracy"), (1 + v("youden_j")) / 2)
    expect_equal(v("markedness"), v("ppv") + v("npv") - 1)
    expect_equal(v("mcc")^2, v("youden_j") * v("markedness"))
    expect_equal(sign(v("mcc")), sign(v("youden_j")))
    expect_equal(v("dor"), v("lr_pos") / v("lr_neg"))
    expect_equal(v("alerts_per_tp"), 1 / v("ppv"))
    expect_equal(v("fp_per_tp"), v("fdr") / v("ppv"))
    expect_equal(v("fp_per_1000_non_events"), 1000 * v("fpr"))
    expect_equal(v("missed_per_1000_total"), 1000 * v("prevalence") * v("fnr"))
    expect_equal(v("alert_rate"), v("prevalence") * v("tpr") + (1 - v("prevalence")) * v("fpr"))
    expect_equal(v("accuracy"), v("prevalence") * v("tpr") + (1 - v("prevalence")) * v("tnr"))
    expect_equal(v("tp_share"), v("prevalence") * v("tpr"))
    expect_equal(v("fowlkes_mallows"), sqrt(v("ppv") * v("tpr")))
    # F-beta count form equals the PPV/TPR form when both are defined
    for (b in c(0.5, 1, 2)) {
      id <- c("0.5" = "f05", "1" = "f1", "2" = "f2")[[as.character(b)]]
      expect_equal(v(id), (1 + b^2) * v("ppv") * v("tpr") / (b^2 * v("ppv") + v("tpr")))
    }
    # kappa label form (accuracy - pe) / (1 - pe) with a = alert rate
    a <- v("alert_rate"); pi <- v("prevalence")
    pe <- pi * a + (1 - pi) * (1 - a)
    expect_equal(v("cohen_kappa"), (v("accuracy") - pe) / (1 - pe))
    # Bayes transform of PPV and NPV
    expect_equal(v("ppv"), v("tpr") * pi / (v("tpr") * pi + v("fpr") * (1 - pi)))
    expect_equal(v("npv"), v("tnr") * (1 - pi) / (v("tnr") * (1 - pi) + v("fnr") * pi))
  }
})

test_that("zero denominators give typed statuses, never silent numbers", {
  z <- classification_metrics(c(tp = 0, fn = 5, fp = 0, tn = 95))
  expect_identical(cm_status(z, "ppv"), "undefined")
  expect_identical(z$reason_code[z$metric_id == "ppv"], "no_alerts")
  expect_true(is.na(cm_value(z, "ppv")))
  expect_identical(cm_status(z, "mcc"), "undefined")
  # count form of F1 is 0 with events but no alerts (declared in the dictionary)
  expect_identical(cm_value(z, "f1"), 0)
  expect_identical(cm_status(z, "lr_pos"), "undefined")
  expect_identical(z$reason_code[z$metric_id == "lr_pos"], "tpr_and_fpr_zero")
  # +infinity is a status with value NA
  inf <- classification_metrics(c(tp = 7, fn = 3, fp = 0, tn = 10))
  expect_identical(cm_status(inf, "lr_pos"), "positive_infinity")
  expect_true(is.na(cm_value(inf, "lr_pos")))
  expect_identical(cm_status(inf, "dor"), "positive_infinity")
  allin <- classification_metrics(c(tp = 5, fn = 0, fp = 95, tn = 0))
  expect_identical(cm_status(allin, "npv"), "undefined")
  expect_identical(allin$reason_code[allin$metric_id == "npv"], "no_non_alerts")
  expect_identical(cm_status(allin, "lr_neg"), "undefined")
  nothing <- classification_metrics(c(tp = 0, fn = 0, fp = 0, tn = 10))
  expect_identical(cm_status(nothing, "f1"), "undefined")
  expect_identical(cm_status(nothing, "cohen_kappa"), "undefined")
  # every non-ok row has a reason in both languages
  de <- classification_metrics(c(tp = 0, fn = 5, fp = 0, tn = 95), lang = "de")
  expect_true(all(!is.na(de$reason[de$status != "ok"])))
  expect_match(de$reason[de$metric_id == "ppv"], "keine Alarme")
  expect_match(de$denominator[de$metric_id == "ppv"], "Alarme")
})

test_that("count validation separates observed, weighted and expected counts", {
  expect_identical(th_condition_code(classification_metrics(c(tp = 1.5, fn = 1, fp = 1, tn = 1))), "invalid_input")
  expect_silent(classification_metrics(c(tp = 1.5, fn = 1, fp = 1, tn = 1), count_kind = "expected"))
  expect_identical(th_condition_code(classification_metrics(c(tp = -1, fn = 1, fp = 1, tn = 1))), "invalid_input")
  expect_identical(th_condition_code(classification_metrics(c(tp = 0, fn = 0, fp = 0, tn = 0))), "invalid_input")
  expect_identical(th_condition_code(classification_metrics(c(a = 1, b = 2))), "invalid_input")
  m <- classification_metrics(c(tp = 1.5, fn = 1, fp = 1, tn = 1), count_kind = "expected")
  expect_identical(unique(m$estimate_kind), "expected")
})

test_that("ppv_gap_max_prevalence maximises PPV - pi and is not the steepest PPV rise", {
  ppv_at <- function(pi, a, b) a * pi / (a * pi + b * (1 - pi))
  pairs <- expand.grid(a = c(0.3, 0.6, 0.9, 1), b = c(0.01, 0.05, 0.2))
  pairs <- pairs[pairs$a > pairs$b, ]
  res <- ppv_gap_max_prevalence(pairs$a, pairs$b)
  expect_true(all(res$status == "ok"))
  for (i in seq_len(nrow(pairs))) {
    a <- pairs$a[i]; b <- pairs$b[i]
    pstar <- res$value[i]
    # stable form equals the source's expression
    expect_equal(pstar, (sqrt(a * b) - b) / (a - b), tolerance = 1e-12)
    # numerical maximiser of PPV(pi) - pi
    opt <- optimize(function(p) ppv_at(p, a, b) - p, c(0, 1), maximum = TRUE, tol = 1e-12)
    expect_equal(pstar, opt$maximum, tolerance = 1e-6)
    # derivative of PPV is 1 there; PPV = 1 - pi there
    h <- 1e-6
    expect_equal((ppv_at(pstar + h, a, b) - ppv_at(pstar - h, a, b)) / (2 * h), 1, tolerance = 1e-6)
    expect_equal(ppv_at(pstar, a, b), 1 - pstar, tolerance = 1e-12)
    # PPV is concave: its slope is largest at pi -> 0 and equals LR+ = a/b there
    slope0 <- (ppv_at(1e-9, a, b) - ppv_at(0, a, b)) / 1e-9
    expect_equal(slope0, a / b, tolerance = 1e-5)
    grid <- seq(0.001, 0.999, length.out = 200)
    slopes <- diff(ppv_at(grid, a, b))
    expect_true(all(diff(slopes) < 0))
  }
  # domain statuses (audit section 4)
  edge <- ppv_gap_max_prevalence(tpr = c(0.5, 0.6, 0.1, 0, 1, NA), fpr = c(0.5, 0, 0.4, 0, 0.04, 0.1))
  expect_identical(edge$status, c("undefined", "undefined", "not_applicable", "undefined", "ok", "undefined"))
  expect_identical(edge$reason_code[1:4], c("equal_tpr_fpr", "fpr_zero_supremum", "tpr_below_fpr", "tpr_and_fpr_zero"))
  expect_equal(edge$value[5], sqrt(0.04) / (1 + sqrt(0.04)))
  expect_true(all(is.na(edge$value[-5])))
  expect_identical(th_condition_code(ppv_gap_max_prevalence(1.2, 0.1)), "invalid_input")
})

test_that("Bayes transforms reproduce the source's worked example and their domains", {
  pv <- predictive_values_at_prevalence(0.2, 0.01, c(0.1, 0.2))
  expect_equal(round(pv$ppv, 3), c(0.690, 0.833))
  need <- prevalence_for_ppv_target(0.2, 0.01, 0.9)
  expect_equal(round(need$prevalence_needed, 3), 0.310)
  expect_equal(need$prevalence_needed, 9 * 0.01 / (0.2 + 9 * 0.01))
  # PPV 90 % is reachable at TPR 0.04, FPR 0.10 for pi >= need (audit L-05)
  need2 <- prevalence_for_ppv_target(0.04, 0.10, 0.9)
  expect_identical(need2$status, "ok")
  expect_equal(predictive_values_at_prevalence(0.04, 0.10, need2$prevalence_needed)$ppv, 0.9)
  edge <- prevalence_for_ppv_target(tpr = c(0, 0.5), fpr = c(0.1, 0), target = 0.9)
  expect_identical(edge$reason_code, c("tpr_zero_unreachable", "fpr_zero_any_positive_prevalence"))
  z <- predictive_values_at_prevalence(0, 0, c(0, 0.5))
  expect_true(all(z$ppv_status == "undefined"))
  expect_identical(th_condition_code(prevalence_for_ppv_target(0.5, 0.1, 1)), "invalid_input")
  # at prevalence 0 PPV is undefined unless there are false positives
  p0 <- predictive_values_at_prevalence(0.5, 0.1, 0)
  expect_identical(p0$ppv, 0)
  p0b <- predictive_values_at_prevalence(0.5, 0, 0)
  expect_identical(p0b$ppv_status, "undefined")
})

test_that("prevalence scenarios keep TPR/FPR and recompute predictive values", {
  sc <- prevalence_scenario(c(tp = 80, fn = 20, fp = 10, tn = 90), prevalence = c(0.02, 0.5), n_total = 1000)
  v <- function(p, id) sc$value[sc$scenario_prevalence == p & sc$metric_id == id]
  expect_equal(v(0.02, "tpr"), 0.8)
  expect_equal(v(0.02, "fpr"), 0.1)
  expect_equal(v(0.02, "ppv"), 0.8 * 0.02 / (0.8 * 0.02 + 0.1 * 0.98))
  expect_equal(v(0.02, "tp"), 16)
  expect_equal(v(0.02, "fp"), 98)
  expect_equal(v(0.5, "ppv"), 0.8 / 0.9)
  expect_identical(unique(sc$estimate_kind), "expected_scenario")
  expect_match(attr(sc, "assumption"), "not a claim about real population shift")
  # no events: the scenario is undefined, never computed with a rate of 0 (audit C-04)
  und <- prevalence_scenario(c(tp = 0, fn = 0, fp = 10, tn = 90), prevalence = 0.1)
  expect_identical(und$scenario_status, "undefined")
  expect_identical(und$reason_code, "no_events")
  expect_identical(th_condition_code(prevalence_scenario(c(tp = 1, fn = 1, fp = 1, tn = 1), 1.5)), "invalid_input")
})
