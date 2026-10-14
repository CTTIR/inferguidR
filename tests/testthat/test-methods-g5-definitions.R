# Group G5 (correlation and regression): method definitions, registry
# agreement, statuses and typed refusals. R/methods-correlation-regression.R.

g5_ids <- c("t_correlation_point_biserial", "t_linreg_one_group_slope", "t_linreg_two_groups_slopes",
            "t_linreg_two_groups_intercepts", "t_mrc_single_coefficient", "f_mrc_r2_deviation_from_zero",
            "f_mrc_r2_increase", "exact_correlation_bivariate_normal", "exact_mrc_random_model")

g5_status <- function(method, mode, inputs, options = list()) {
  spec <- tryCatch(analysis_spec(method, mode, inputs = inputs, options = options), inferguidR_error = function(e) e)
  if (inherits(spec, "inferguidR_error")) return(spec$status$code)
  calculate_power(spec)$status$code
}

test_that("the nine G5 procedures are registered in their registry families with all five modes", {
  for (id in g5_ids) {
    def <- .ig_get_method(id)
    meta <- .ig_get_metadata(id)
    expect_false(is.null(def), info = id)
    expect_identical(def$family, meta$family, info = id)
    expect_setequal(def$modes, meta$modes)
    expect_identical(def$status, "independently_checked", info = id)
    expect_true(file.exists(test_path("..", "..", def$checker)) || !dir.exists(test_path("..", "..", "validation")), info = id)
    expect_identical(validation_status(id)$gpower_reference$status, "reference_pending", info = id)
  }
})

test_that("definition options match the certified registry choices and defaults", {
  for (id in g5_ids) {
    def <- .ig_get_method(id)
    meta <- .ig_get_metadata(id)
    choice <- Filter(function(o) identical(o$value_type, "choice"), meta$options)
    expect_setequal(names(def$options), vapply(choice, function(o) o$option_id, ""))
    for (o in choice) {
      d <- def$options[[o$option_id]]
      expect_setequal(d$values, o$codes)
      expect_identical(d$default, o$default, info = paste(id, o$option_id))
    }
  }
  # Continuous registry options are parameters of the definitions.
  expect_identical(.ig_get_method("t_linreg_two_groups_slopes")$params$allocation_ratio$default, 1)
  expect_identical(.ig_get_method("exact_correlation_bivariate_normal")$params$exact_n_threshold$default, 10000)
})

test_that("G*Power refusals are typed refusals without numbers (TD-0010 cases of the captures)", {
  # target power below alpha (e.g. GPW-T-PBIS-APR-B001-I03, GPW-EX-MRCRND-SEN-B002-I02)
  expect_identical(g5_status("t_correlation_point_biserial", "a_priori", list(rho = .3, alpha = .05, power = .04)), "invalid_input")
  expect_identical(g5_status("exact_mrc_random_model", "sensitivity", list(rho2_h0 = 0, alpha = .05, power = .04, n_total = 95, n_predictors = 5)),
                   "invalid_input")
  # zero effect in a priori mode (GPW-T-PBIS-APR-B001-I02, GPW-T-INTER2-APR-B001-I02, GPW-EX-CORRBN-APR-B001-I02)
  expect_identical(g5_status("t_correlation_point_biserial", "a_priori", list(rho = 0, alpha = .05, power = .8)), "invalid_input")
  expect_identical(g5_status("t_linreg_two_groups_intercepts", "a_priori", list(delta_intercept = 0, alpha = .05, power = .8, sd_residual = 1,
                                                                                 sd_x1 = 1, sd_x2 = 1, mean_x1 = 0, mean_x2 = 0)), "invalid_input")
  expect_identical(g5_status("exact_correlation_bivariate_normal", "a_priori", list(rho_h1 = 0, rho_h0 = 0, alpha = .05, power = .8)), "invalid_input")
  expect_identical(g5_status("f_mrc_r2_deviation_from_zero", "a_priori", list(f2 = 0, alpha = .05, power = .8, n_predictors = 5)), "invalid_input")
  # domain refusals (alpha 1.5, N 0, sigma 0, |rho| > 1, negative f2)
  expect_identical(g5_status("t_linreg_one_group_slope", "post_hoc", list(slope_h1 = .3, sd_x = 1, sd_y = 1, alpha = 1.5, n_total = 50)), "invalid_input")
  expect_identical(g5_status("t_mrc_single_coefficient", "compromise", list(f2 = .15, error_ratio = 1, n_total = 0, n_predictors = 5)), "invalid_input")
  expect_identical(g5_status("t_linreg_two_groups_slopes", "post_hoc", list(delta_slope = .5, alpha = .05, n1 = 50, n2 = 50, sd_residual = 0,
                                                                             sd_x1 = 1, sd_x2 = 1)), "invalid_input")
  expect_identical(g5_status("exact_correlation_bivariate_normal", "post_hoc", list(rho_h1 = 1.2, alpha = .05, n_total = 84)), "invalid_input")
  expect_identical(g5_status("t_mrc_single_coefficient", "post_hoc", list(f2 = -.15, alpha = .05, n_total = 100, n_predictors = 5)), "invalid_input")
  # cross-parameter constraints
  expect_identical(g5_status("f_mrc_r2_increase", "post_hoc", list(f2 = .1, alpha = .05, n_total = 100, n_tested = 6, n_predictors = 5)),
                   "invalid_input")
  expect_identical(g5_status("t_mrc_single_coefficient", "post_hoc", list(f2 = .1, alpha = .05, n_total = 6, n_predictors = 5)), "invalid_input")
  expect_identical(g5_status("t_linreg_one_group_slope", "post_hoc", list(slope_h1 = 1.2, sd_x = 1, sd_y = 1, alpha = .05, n_total = 50)),
                   "invalid_input")
})

test_that("options apply only where the inventory shows them", {
  # effect direction: sensitivity only
  expect_identical(g5_status("t_linreg_one_group_slope", "post_hoc", list(slope_h1 = .3, sd_x = 1, sd_y = 1, alpha = .05, n_total = 50),
                             list(effect_direction = "slope_lt_h0")), "invalid_input")
  # the exact-distribution threshold belongs to the exact option
  expect_identical(g5_status("exact_correlation_bivariate_normal", "post_hoc", list(rho_h1 = .3, alpha = .05, n_total = 84, exact_n_threshold = 50),
                             list(computation = "fisher_z_approx")), "invalid_input")
  # allocation ratio: a priori only
  expect_identical(g5_status("t_linreg_two_groups_slopes", "post_hoc", list(delta_slope = .5, alpha = .05, n1 = 50, n2 = 50, sd_residual = 1,
                                                                             sd_x1 = 1, sd_x2 = 1, allocation_ratio = 2)), "invalid_input")
})
