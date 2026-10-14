# Export consistency: the result record, its JSON round trip, the analysis
# envelope and report snapshot built from it, value formatting and the
# absence of any G*Power validation claim.

ex_results <- function() {
  ex <- classification_examples("scored_synthetic")
  list(
    scored = evaluate_cutoff(update(ex$spec, prevalence = 0.05),
                             partial_auc = list(focus = "fpr", range = c(0, 0.2)),
                             threshold_probability = c(0.1, 0.3)),
    selected = select_cutoff(ex$spec, "max_mcc"),
    infeasible = select_cutoff(ex$spec, "constraint", constraints = list(tpr_min = 1, fpr_max = 0)),
    counts = evaluate_cutoff(classification_examples("counts_fixed")$spec),
    theory = select_cutoff(classification_examples("binormal_unequal_variance")$spec, "max_youden",
                           partial_auc = list(focus = "tpr", range = c(0.8, 1)))
  )
}

test_that("results pair with their evaluated specification and round-trip through JSON unchanged", {
  for (nm in names(rs <- ex_results())) {
    r <- rs[[nm]]
    spec <- classification_evaluated_spec(r)
    expect_identical(spec$spec_hash, r$spec_hash, info = nm)
    env <- analysis_envelope(spec, r)
    js <- ig_to_json(env)
    back <- ig_from_json(js)
    expect_identical(back$analysis_hash, env$analysis_hash, info = nm)
    expect_identical(back$result_hash, env$result_hash, info = nm)
    expect_identical(as.character(ig_to_json(back)), as.character(js), info = nm)
    # the accessors read the rebuilt result exactly as the original
    for (tb in c("metrics", "intervals", "summaries", "selection", "limitations")) {
      expect_equal(classification_table(back$result, tb), classification_table(r, tb), info = paste(nm, tb))
    }
    # a rebuilt result has no attached specification; the envelope holds it
    expect_identical(th_condition_code(classification_evaluated_spec(back$result)), "invalid_input")
    expect_identical(back$spec$spec_hash, spec$spec_hash)
  }
})

test_that("a result cannot be paired with another specification or modified unnoticed", {
  rs <- ex_results()
  ex <- classification_examples("scored_synthetic")
  expect_identical(th_condition_code(analysis_envelope(ex$spec, rs$selected)), "invalid_input")
  env <- analysis_envelope(classification_evaluated_spec(rs$counts), rs$counts)
  js <- sub("\"tp\":40", "\"tp\":41", as.character(ig_to_json(env)), fixed = TRUE)
  expect_identical(th_condition_code(ig_from_json(js)), "invalid_input")
  expect_error(rs$counts$status <- list(code = "ok"))
})

test_that("a rebuilt scored specification evaluates again only with the matching data", {
  ex <- classification_examples("scored_synthetic")
  r <- evaluate_cutoff(ex$spec)
  back <- ig_from_json(ig_to_json(analysis_envelope(ex$spec, r)))
  expect_identical(th_condition_code(evaluate_cutoff(back$spec)), "invalid_input")
  again <- evaluate_cutoff(attach_classification_data(back$spec, ex$data))
  expect_identical(again$models, r$models)
  wrong <- ex$data
  wrong$model_a[1] <- 0.5
  expect_identical(th_condition_code(attach_classification_data(back$spec, wrong)), "invalid_input")
})

test_that("results enter report snapshots; final reports need ok results", {
  rs <- ex_results()
  env <- analysis_envelope(classification_evaluated_spec(rs$scored), rs$scored)
  pd <- classification_plot_data(rs$scored, "roc")
  snap <- inferguidR:::.ig_new_report_snapshot(env, report_kind = "model_cutoff", plot_data = list(pd),
                                               created_at = "2026-10-05T10:00:00.000Z")
  expect_identical(snap$completeness, "complete")
  back <- ig_from_json(ig_to_json(snap))
  expect_identical(back$snapshot_hash, snap$snapshot_hash)
  expect_identical(back$envelopes[[1]]$analysis_hash, env$analysis_hash)
  bad <- analysis_envelope(classification_evaluated_spec(rs$infeasible), rs$infeasible)
  expect_identical(th_condition_code(inferguidR:::.ig_new_report_snapshot(bad, report_kind = "model_cutoff")),
                   "invalid_input")
  inc <- inferguidR:::.ig_new_report_snapshot(bad, report_kind = "incomplete_settings")
  expect_identical(inc$completeness, "incomplete")
})

test_that("statuses are honest: no G*Power claim, no number without an ok status", {
  for (r in ex_results()) {
    v <- r$validation
    expect_identical(v$gpower_reference_status, "not_applicable")
    expect_false(any(grepl("gpower_validated", unlist(v))))
    expect_true(all(unlist(v$components) %in% c("independently_checked", "development")))
    expect_identical(v$components$bootstrap, "development")
    expect_false(grepl("gpower_validated", ig_to_json(r)))
    for (m in r$models) {
      for (tb in list(m$metrics, m$summaries, m$scenario$metrics)) {
        if (is.null(tb)) next
        expect_true(all(is.finite(tb$value[tb$status == "ok"])))
        expect_true(all(is.na(tb$value[tb$status != "ok"])))
        expect_false(anyNA(tb$reason_code[tb$status != "ok"]))
      }
    }
  }
  # the files named in the validation record exist in the source tree
  root <- th_pkg_root()
  if (!is.null(root) && dir.exists(file.path(root, "validation"))) {
    v <- ex_results()$counts$validation
    expect_true(file.exists(file.path(root, v$independent_checker)))
    expect_true(file.exists(file.path(root, v$checker_summary)))
    expect_true(file.exists(file.path(root, v$cross_check_summary)))
  }
})

test_that("results print a readable summary and record a reproducible call", {
  r <- ex_results()$scored
  txt <- format(r)
  expect_match(txt[1], "classification_result")
  expect_true(any(grepl("rule:   alert if score >= 0.3", txt, fixed = TRUE)))
  expect_true(any(grepl("ROC AUC", txt)))
  expect_true(any(grepl("limitation: prevalence_scenario", txt)))
  expect_output(print(r), "spec hash")
  expect_match(r$r_call, "^inferguidR::evaluate_cutoff\\(spec, partial_auc = ")
  expect_match(r$r_call, r$spec_hash, fixed = TRUE)
  expect_identical(th_condition_code(classification_table(list(), "metrics")), "invalid_input")
  de <- classification_table(ex_results()$counts, "intervals", lang = "de")
  expect_true(all(!is.na(de$reason[de$status != "ok"])))
})

test_that("display formatting keeps values, marks expected counts and never leaves a dash unexplained", {
  m <- classification_metrics(c(tp = 0, fn = 5, fp = 0, tn = 95))
  f <- format_metric_values(m)
  expect_identical(f$value, m$value)
  expect_identical(f$display[f$metric_id == "ppv"], "\u2013")
  expect_match(f$display_reason[f$metric_id == "ppv"], "no alerts")
  expect_true(all(!is.na(f$display_reason[f$display == "\u2013"])))
  expect_identical(f$display[f$metric_id == "tn"], "95")
  expect_identical(f$display[f$metric_id == "tnr"], "100.0%")
  inf <- format_metric_values(classification_metrics(c(tp = 7, fn = 3, fp = 0, tn = 10)))
  expect_identical(inf$display[inf$metric_id == "lr_pos"], "\u221e")
  expect_match(inf$display_reason[inf$metric_id == "lr_pos"], "FPR")
  expect_identical(inf$display[inf$metric_id == "lr_neg"], "0.3")
  sc <- prevalence_scenario(c(tp = 8, fn = 2, fp = 9, tn = 81), prevalence = 0.02, n_total = 1000)
  fs <- format_metric_values(sc, lang = "de")
  expect_identical(fs$display[fs$metric_id == "tp"], "\u224816,0")
  expect_identical(fs$display[fs$metric_id == "fp"], "\u224898,0")
  expect_match(fs$display[fs$metric_id == "ppv"], ",")
  expect_match(fs$display[fs$metric_id == "ppv"], "\u00a0%$")
  tiny <- format_metric_values(data.frame(metric_id = "fpr", value = 1e-9, status = "ok"))
  expect_false(tiny$display == "0.000%")
  big <- format_metric_values(classification_metrics(c(tp = 500, fn = 1, fp = 1, tn = 5000)))
  expect_identical(big$display[big$metric_id == "n_total"], "5,502")
  expect_identical(big$display[big$metric_id == "dor"], "2,500,000")
  sm <- format_metric_values(classification_table(ex_results()$theory, "summaries"))
  expect_true(all(nzchar(sm$display)))
  expect_identical(th_condition_code(format_metric_values(data.frame(x = 1))), "invalid_input")
  expect_identical(th_condition_code(format_metric_values(m, digits = -1)), "invalid_input")
})
