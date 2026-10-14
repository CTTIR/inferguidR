test_that("power workspace commits through the R API and keeps draft edits separate", {
  skip_if_not_installed("shiny")
  lang <- shiny::reactiveVal("en")
  shiny::testServer(.ig_power_workspace_server, args = list(language = lang), {
    session$flushReact()
    session$setInputs(mode = "post_hoc", option_t_means_one_sample_tails = "two", field_t_means_one_sample_d = "0.5", field_t_means_one_sample_alpha = "0.05",
                      field_t_means_one_sample_n_total = "30", research = FALSE, calculate = 1)
    expect_null(state())
    session$setInputs(research = TRUE, calculate = 2)
    direct <- calculate_power(analysis_spec("t_means_one_sample", "post_hoc",
      inputs = list(d = .5, alpha = .05, n_total = 30), options = list(tails = "two")))
    expect_identical(state()$committed$result_hash, direct$result_hash)
    session$setInputs(field_t_means_one_sample_d = "0.7")
    expect_true(state()$dirty)
    expect_identical(state()$committed$result_hash, direct$result_hash)
    session$setInputs(calculate = 3)
    expect_false(state()$dirty)
    expect_gt(state()$committed$result$achieved$power, direct$achieved$power)
    before <- state()$committed
    session$setInputs(field_t_means_one_sample_n_total = "not a number", calculate = 4)
    expect_identical(state()$committed, before)
    expect_match(message(), "Enter one number")
    expect_match(output$status, "previous result remains displayed", fixed = TRUE)
    lang("de")
    session$flushReact()
    expect_identical(state()$committed, before)
    expect_match(output$method_selector$html, "f_generic")
    expect_match(output$method_selector$html, "boschloo_proportions_two_independent")
  })
})

test_that("classifier workspace uses synthetic API data without recomputing on edits", {
  skip_if_not_installed("shiny")
  shiny::testServer(.ig_classifier_workspace_server, args = list(language = shiny::reactive("en")), {
    session$flushReact()
    session$setInputs(example = "counts_fixed", cutoff = "0.3", evaluate = 1)
    direct <- evaluate_cutoff(classification_examples("counts_fixed")$spec)
    expect_equal(result()$models[[1]]$operating_point, direct$models[[1]]$operating_point)
    before <- result()
    session$setInputs(cutoff = "0.9")
    expect_identical(result(), before)
    expect_false(grepl("previous evaluation", output$message))
    session$setInputs(example = "scored_synthetic")
    expect_match(output$message, "previous evaluation")
  })
})

test_that("the bundled workbench creates an app without launching a server", {
  skip_if_not_installed("shiny")
  skip_if_not_installed("ggplot2")
  expect_s3_class(workbench_app(), "shiny.appobj")
})

test_that("method switching keeps drafts separate and uses each family's canonical inputs", {
  skip_if_not_installed("shiny")
  shiny::testServer(.ig_power_workspace_server, args = list(language = shiny::reactive("en")), {
    session$flushReact()
    fields <- function(method, values) do.call(session$setInputs,
      stats::setNames(values, vapply(names(values), function(nm) .ig_workbench_field_id(method, nm), "")))
    session$setInputs(method = "f_generic", mode = "post_hoc", research = TRUE)
    fields("f_generic", list(df1 = "3", df2 = "40", lambda = "5", alpha = "0.05"))
    session$setInputs(calculate = 1)
    direct <- calculate_power(analysis_spec("f_generic", "post_hoc",
      inputs = list(df1 = 3, df2 = 40, lambda = 5, alpha = .05)))
    expect_identical(state()$committed$result_hash, direct$result_hash)
    session$setInputs(method = "t_means_two_independent")
    expect_true(state()$dirty)
    expect_identical(state()$committed$result_hash, direct$result_hash)
    fields("t_means_two_independent", list(d = "0.5", alpha = "0.05", n1 = "20", n2 = "25"))
    session$setInputs(option_t_means_two_independent_tails = "two", calculate = 2)
    direct2 <- calculate_power(analysis_spec("t_means_two_independent", "post_hoc",
      inputs = list(d = .5, alpha = .05, n1 = 20, n2 = 25), options = list(tails = "two")))
    expect_identical(state()$committed$result_hash, direct2$result_hash)
    session$setInputs(method = "f_generic", calculate = 3)
    expect_identical(state()$committed$result_hash, direct$result_hash)
    expect_identical(draft()$fields$df1, "3")
  })
})

test_that("all registry methods render editable controls in both locales", {
  skip_if_not_installed("shiny")
  lang <- shiny::reactiveVal("en")
  shiny::testServer(.ig_power_workspace_server, args = list(language = lang), {
    session$flushReact()
    for (locale in c("en", "de")) {
      lang(locale)
      for (id in list_methods()$id) {
        session$setInputs(method = id, mode = "post_hoc")
        expect_true(nzchar(output$fields$html), info = paste(id, locale))
        expect_match(output$method_selector$html, id, fixed = TRUE)
      }
    }
    expect_null(state())
  })
})

test_that("classifier model and panel selection use the committed evaluation", {
  skip_if_not_installed("shiny")
  shiny::testServer(.ig_classifier_workspace_server, args = list(language = shiny::reactive("en")), {
    session$flushReact()
    session$setInputs(example = "scored_synthetic", cutoff = "0.3", evaluate = 1)
    before <- result()
    session$setInputs(model = "model_b", curve_panel = "pr")
    expect_identical(panel(), classification_plot_data(before, "confusion_table", "model_b"))
    expect_identical(curve_data(), classification_plot_data(before, "pr", "model_b"))
    session$setInputs(cutoff = "0.7", curve_panel = "roc")
    expect_identical(result(), before)
    expect_identical(curve_data(), classification_plot_data(before, "roc", "model_b"))
    session$setInputs(example = "counts_fixed", evaluate = 2)
    expect_false(identical(selected_model(), "model_b"))
    expect_error(curve_data(), "Fixed counts do not define a curve")
    expect_identical(panel(), classification_plot_data(result(), "confusion_table"))
  })
})

test_that("classifier scenario applies only on explicit evaluation and basis choice", {
  skip_if_not_installed("shiny")
  shiny::testServer(.ig_classifier_workspace_server, args = list(language = shiny::reactive("en")), {
    session$flushReact()
    session$setInputs(example = "scored_synthetic", cutoff = "0.3", scenario = TRUE, prevalence = "0.05", evaluate = 1)
    session$setInputs(basis = "scenario", curve_panel = "pr")
    before <- result()
    expect_identical(selected_basis(), "scenario")
    expect_equal(curve_data()$data$baseline$precision, .05)
    expect_identical(panel(), classification_plot_data(before, "confusion_table", basis = "scenario"))
    session$setInputs(prevalence = "0.2")
    expect_identical(result(), before)
    expect_equal(curve_data()$data$baseline$precision, .05)
    session$setInputs(example = "binormal_equal_variance", evaluate = 2)
    expect_identical(selected_basis(), "observed")
  })
})


test_that("classifier additional panels expose their authoritative tables", {
  skip_if_not_installed("shiny")
  shiny::testServer(.ig_classifier_workspace_server, args = list(language = shiny::reactive("en")), {
    session$flushReact()
    session$setInputs(example = "scored_synthetic", cutoff = "0.3", scenario = TRUE, prevalence = "0.05", evaluate = 1)
    session$setInputs(basis = "scenario", curve_panel = "model_comparison")
    expect_identical(curve_data()$basis, "observed")
    expect_equal(length(unique(curve_data()$data$curves$model_id)), 2)
    session$setInputs(curve_panel = "prevalence_curve")
    expect_match(output$data_selector$html, "rates")
    session$setInputs(data_table = "rates")
    expect_match(output$curve_data, "tpr")
    session$setInputs(example = "counts_fixed", evaluate = 2, curve_panel = "mosaic")
    expect_identical(curve_data()$basis, "scenario")
    expect_equal(sum(curve_data()$data$rectangles$area_share, na.rm = TRUE), 1)
  })
})

test_that("classifier metric controls and complete tables share the selected evaluation", {
 skip_if_not_installed("shiny")
 shiny::testServer(.ig_classifier_workspace_server,args=list(language=shiny::reactive("en")), {
  session$flushReact()
  session$setInputs(example="scored_synthetic",cutoff="0.3",evaluate=1)
  before <- result()
  session$setInputs(curve_panel="metric_pair",metric_x="lr_pos",metric_y="ppv",log_x=TRUE)
  expect_identical(curve_data(),classification_plot_data(before,"metric_pair",metric_x="lr_pos",metric_y="ppv"))
  expect_match(output$metric_controls$html,"Log X")
  session$setInputs(curve_panel="metric_vs_cutoff",metric_family="summary")
  expect_match(output$metric_controls$html,"Likelihood ratios")
  session$setInputs(curve_panel="current_table")
  expect_equal(nrow(displayed_table()),nrow(classification_plot_data(before,"current_table")$data$metrics))
  expect_match(output$primary_table,"denominator")
  raw <- classification_plot_data(before,"current_table")$data$metrics
  expect_equal(as.numeric(displayed_table()$value),raw$value,tolerance=0)
  session$setInputs(curve_panel="sweep_table")
  expect_equal(nrow(displayed_table()),nrow(cutoff_sweep(before, model_id="model_a")))
  expect_identical(result(),before)
 })
})
