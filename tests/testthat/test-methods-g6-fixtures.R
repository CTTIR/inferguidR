# G*Power reference fixtures for z_logistic_regression and z_poisson_regression
# (group G6): a curated, fast subset of the comparison that
# validation/tools/compare_fixtures.R runs on every fixture of the two methods
# (validation/compare/z_logistic_regression/, validation/compare/z_poisson_regression/).
#
# Every expected string below is the value G*Power 3.1.9.7 displayed in the
# cited fixture record (emulator_capture: G*Power under Wine, D22/D23;
# documented_example: the printed manual, secondary evidence from older
# builds). The specification arguments are those of the method's mapping
# (validation/compare/<method_id>/mapping.R). Match rules as in
# validation/tolerances.md: integers exactly; reals within the closed rounding
# interval of the displayed digits plus 8 eps max(1, |x|). The last test
# checks, in a development checkout, that the embedded strings equal the
# fixture files. Records chosen: ordinary, edge and difficult cases of every
# mode, procedure family and several covariate distributions, among those
# that agree; disagreements are triaged in the comparison summaries, not here.

g6_display_ok <- function(ours, shown) {
  if (!is.finite(ours)) return(FALSE)
  s <- trimws(shown)
  if (grepl("^-?[0-9]+$", s)) return(ours == as.numeric(s))
  x <- as.numeric(s)
  mant <- sub("[eE].*$", "", s)
  ex <- if (grepl("[eE]", s)) as.numeric(sub("^.*[eE]", "", s)) else 0
  dec <- if (grepl(".", mant, fixed = TRUE)) nchar(sub("^[^.]*\\.", "", mant)) else 0
  half <- 0.5 * 10^(ex - dec)
  allowance <- 8 * .Machine$double.eps * max(1, abs(x))
  ours >= x - half - allowance && ours <= x + half + allowance
}

g6_outputs <- function(res) {
  q <- res$quantities
  a <- res$achieved
  s <- res$solved
  mo <- res$method_outputs
  chi <- grepl("^enumeration", mo$procedure %||% "")
  crit <- q$critical
  cz <- if (length(crit) == 1L) crit else if (isTRUE(mo$beta1 < 0) && !identical(mo$procedure, "hsieh")) min(crit) else max(crit)
  out <- list("Total sample size" = s$n_total, "Actual power" = a$power, "Power (1-\u03b2 err prob)" = a$power,
              "\u03b1 err prob" = s$alpha, "\u03b2 err prob" = a$beta,
              "Odds ratio" = s$odds_ratio %||% mo$odds_ratio, "Pr(Y=1|X=1) H1" = s$p2 %||% mo$p2,
              "Exp(\u03b21)" = s$exp_beta1)
  if (chi) {
    out[["Critical \u03c7\u00b2"]] <- crit
    out[["Df"]] <- q$df
    out[["Noncentrality parameter \u03bb"]] <- q$ncp
  } else {
    out[["Critical z"]] <- cz
  }
  out
}

g6_fx_spec <- function(fx) {
  args <- list(method = fx$method, mode = fx$mode, inputs = fx$inputs, options = fx$options)
  if (!is.null(fx$target)) args$target <- fx$target
  do.call(analysis_spec, args)
}

# Sensitivity through the method's own search (the core hook that runs it
# from calculate_power() is change request CR-G6-01).
g6_sensitivity_result <- function(fx) {
  spec <- g6_fx_spec(fx)
  model <- .ig_get_method(fx$method)$power
  p <- c(spec$inputs, spec$options)
  s <- .ig_glm_sensitivity_search(model)(p, spec$target)
  p[[spec$target]] <- s$value
  fw <- model(p)
  list(search = s, solved = stats::setNames(list(s$value), spec$target),
       achieved = list(power = fw$power, beta = fw$beta),
       quantities = list(critical = fw$critical, df = fw$df, ncp = fw$ncp), method_outputs = fw$outputs)
}

g6_fixtures <- list(
  list(id = "GPW-Z-LOGIT-APR-B081-O01--r01", method = "z_logistic_regression", mode = "a_priori", target = NULL,
       inputs = list(odds_ratio = 1.5, p1 = 0.5, alpha = 0.05, power = 0.8, r2_other_x = 0, x_mu = 0, x_sigma = 1),
       options = list(computation = "large_sample", test_procedure = "demidenko", variance_correction = "on", effect_size_input_mode = "odds_ratio", tails = "one", x_distribution = "normal"),
       shown = list("Critical z" = "1.6448536", "Total sample size" = "163", "Actual power" = "0.8012292")),
  list(id = "GPW-Z-LOGIT-APR-B213-O01--r01", method = "z_logistic_regression", mode = "a_priori", target = NULL,
       inputs = list(odds_ratio = 1.5, p1 = 0.5, alpha = 0.05, power = 0.95, r2_other_x = 0, x_mu = 0, x_sigma = 1),
       options = list(computation = "large_sample", test_procedure = "demidenko", variance_correction = "on", effect_size_input_mode = "odds_ratio", tails = "two", x_distribution = "normal"),
       shown = list("Critical z" = "1.9599640", "Total sample size" = "337", "Actual power" = "0.9500770")),
  list(id = "GPW-Z-LOGIT-APR-B213-O02--r01", method = "z_logistic_regression", mode = "a_priori", target = NULL,
       inputs = list(odds_ratio = 1.5, p1 = 0.5, alpha = 0.05, power = 0.95, r2_other_x = 0.1, x_mu = 0, x_sigma = 1),
       options = list(computation = "large_sample", test_procedure = "demidenko", variance_correction = "on", effect_size_input_mode = "odds_ratio", tails = "two", x_distribution = "normal"),
       shown = list("Critical z" = "1.9599640", "Total sample size" = "375", "Actual power" = "0.9503607")),
  list(id = "GPW-Z-LOGIT-APR-B218-O01--r01", method = "z_logistic_regression", mode = "a_priori", target = NULL,
       inputs = list(odds_ratio = 1.5, p1 = 0.5, alpha = 0.05, power = 0.95, r2_other_x = 0, x_mu = 0, x_sigma = 1),
       options = list(computation = "large_sample", test_procedure = "hsieh", effect_size_input_mode = "odds_ratio", tails = "two", x_distribution = "normal"),
       shown = list("Critical z" = "1.9599640", "Total sample size" = "317", "Actual power" = "0.9504862")),
  list(id = "GPW-Z-LOGIT-APR-B016-O01--r01", method = "z_logistic_regression", mode = "a_priori", target = NULL,
       inputs = list(odds_ratio = 1.5, p1 = 0.5, alpha = 0.05, power = 0.8, r2_other_x = 0, x_pi = 0.5),
       options = list(computation = "large_sample", test_procedure = "demidenko", variance_correction = "on", effect_size_input_mode = "odds_ratio", tails = "one", x_distribution = "binomial"),
       shown = list("Critical z" = "1.6448536", "Total sample size" = "613", "Actual power" = "0.8003552")),
  list(id = "GPW-Z-LOGIT-APR-B132-O01--r01", method = "z_logistic_regression", mode = "a_priori", target = NULL,
       inputs = list(odds_ratio = 1.5, p1 = 0.5, alpha = 0.05, power = 0.8, r2_other_x = 0, x_var_h0 = 4, x_var_h1 = 4.17),
       options = list(computation = "large_sample", test_procedure = "demidenko", variance_correction = "on", effect_size_input_mode = "odds_ratio", tails = "one", x_distribution = "manual"),
       shown = list("Critical z" = "1.6448536", "Total sample size" = "155", "Actual power" = "0.8008208")),
  list(id = "GPW-Z-LOGIT-APR-B066-O01--r01", method = "z_logistic_regression", mode = "a_priori", target = NULL,
       inputs = list(enumeration_threshold = 100000L, odds_ratio = 1.5, p1 = 0.5, alpha = 0.05, power = 0.8, r2_other_x = 0, x_mu = 0, x_sigma = 1),
       options = list(computation = "enumeration", test_statistic = "wald", effect_size_input_mode = "odds_ratio", tails = "one", x_distribution = "normal"),
       shown = list("Noncentrality parameter \u03bb" = "6.1831728", "Critical \u03c7\u00b2" = "2.7055435", Df = "1", "Total sample size" = "171", "Actual power" = "0.8000527")),
  list(id = "GPW-Z-LOGIT-APR-B066-D01--r01", method = "z_logistic_regression", mode = "a_priori", target = NULL,
       inputs = list(enumeration_threshold = 100000L, odds_ratio = 1.1, p1 = 0.3, alpha = 0.001, power = 0.99, r2_other_x = 0.5, x_mu = 0, x_sigma = 1),
       options = list(computation = "enumeration", test_statistic = "wald", effect_size_input_mode = "odds_ratio", tails = "one", x_distribution = "normal"),
       shown = list("Noncentrality parameter \u03bb" = "29.3408463", "Critical \u03c7\u00b2" = "9.5495357", Df = "1", "Total sample size" = "30925", "Actual power" = "0.9900037")),
  list(id = "GPW-Z-LOGIT-APR-B001-O01--r01", method = "z_logistic_regression", mode = "a_priori", target = NULL,
       inputs = list(enumeration_threshold = 100000L, odds_ratio = 1.5, p1 = 0.5, alpha = 0.05, power = 0.8, r2_other_x = 0, x_pi = 0.5),
       options = list(computation = "enumeration", test_statistic = "wald", effect_size_input_mode = "odds_ratio", tails = "one", x_distribution = "binomial"),
       shown = list("Noncentrality parameter \u03bb" = "6.1904259", "Critical \u03c7\u00b2" = "2.7055435", Df = "1", "Total sample size" = "615", "Actual power" = "0.8004605")),
  list(id = "GPW-Z-LOGIT-APR-B051-D01--r01", method = "z_logistic_regression", mode = "a_priori", target = NULL,
       inputs = list(enumeration_threshold = 100000L, odds_ratio = 1.1, p1 = 0.3, alpha = 0.001, power = 0.99, r2_other_x = 0.5, x_mu = 0, x_sigma = 0.5),
       options = list(computation = "enumeration", test_statistic = "likelihood_ratio", effect_size_input_mode = "odds_ratio", tails = "one", x_distribution = "lognormal"),
       shown = list("Noncentrality parameter \u03bb" = "29.3398650", "Critical \u03c7\u00b2" = "9.5495357", Df = "1", "Total sample size" = "79397", "Actual power" = "0.9900013")),
  list(id = "GPW-Z-LOGIT-POS-B076-O01--r01", method = "z_logistic_regression", mode = "post_hoc", target = NULL,
       inputs = list(odds_ratio = 1.5, p1 = 0.5, alpha = 0.05, n_total = 200L, r2_other_x = 0, x_mu = 0, x_sigma = 1),
       options = list(computation = "large_sample", test_procedure = "demidenko", variance_correction = "off", effect_size_input_mode = "odds_ratio", tails = "one", x_distribution = "normal"),
       shown = list("Critical z" = "1.6448536", "Power (1-\u03b2 err prob)" = "0.8558525")),
  list(id = "GPW-Z-LOGIT-POS-B061-O01--r01", method = "z_logistic_regression", mode = "post_hoc", target = NULL,
       inputs = list(odds_ratio = 1.5, p1 = 0.5, alpha = 0.05, n_total = 200L, r2_other_x = 0, x_mu = 0, x_sigma = 1),
       options = list(computation = "large_sample", test_procedure = "demidenko", variance_correction = "on", effect_size_input_mode = "odds_ratio", tails = "one", x_distribution = "lognormal"),
       shown = list("Critical z" = "1.6448536", "Power (1-\u03b2 err prob)" = "0.9756963")),
  list(id = "GPW-Z-LOGIT-POS-B041-O01--r01", method = "z_logistic_regression", mode = "post_hoc", target = NULL,
       inputs = list(odds_ratio = 1.5, p1 = 0.5, alpha = 0.05, n_total = 200L, r2_other_x = 0, x_lambda = 1),
       options = list(computation = "large_sample", test_procedure = "demidenko", variance_correction = "on", effect_size_input_mode = "odds_ratio", tails = "one", x_distribution = "exponential"),
       shown = list("Critical z" = "1.6448536", "Power (1-\u03b2 err prob)" = "0.8017698")),
  list(id = "GPW-Z-LOGIT-POS-B106-O01--r01", method = "z_logistic_regression", mode = "post_hoc", target = NULL,
       inputs = list(odds_ratio = 1.5, p1 = 0.5, alpha = 0.05, n_total = 200L, r2_other_x = 0, x_lambda = 1),
       options = list(computation = "large_sample", test_procedure = "demidenko", variance_correction = "on", effect_size_input_mode = "odds_ratio", tails = "one", x_distribution = "poisson"),
       shown = list("Critical z" = "1.6448536", "Power (1-\u03b2 err prob)" = "0.8404797")),
  list(id = "GPW-Z-LOGIT-POS-B126-O01--r01", method = "z_logistic_regression", mode = "post_hoc", target = NULL,
       inputs = list(odds_ratio = 1.5, p1 = 0.5, alpha = 0.05, n_total = 200L, r2_other_x = 0, x_left = 0, x_right = 1),
       options = list(computation = "large_sample", test_procedure = "demidenko", variance_correction = "on", effect_size_input_mode = "odds_ratio", tails = "one", x_distribution = "uniform"),
       shown = list("Critical z" = "1.6448536", "Power (1-\u03b2 err prob)" = "0.2042623")),
  list(id = "GPW-Z-LOGIT-POS-B264-O01--r01", method = "z_logistic_regression", mode = "post_hoc", target = NULL,
       inputs = list(odds_ratio = 1.5, p1 = 0.5, alpha = 0.05, n_total = 200L, r2_other_x = 0, x_var_h0 = 4, x_var_h1 = 4.17),
       options = list(computation = "large_sample", test_procedure = "demidenko", variance_correction = "on", effect_size_input_mode = "odds_ratio", tails = "two", x_distribution = "manual"),
       shown = list("Critical z" = "1.9599640", "Power (1-\u03b2 err prob)" = "0.8067257")),
  list(id = "GPW-Z-LOGIT-POS-B153-O01--r01", method = "z_logistic_regression", mode = "post_hoc", target = NULL,
       inputs = list(odds_ratio = 1.5, p1 = 0.5, alpha = 0.05, n_total = 200L, r2_other_x = 0, x_pi = 0.5),
       options = list(computation = "large_sample", test_procedure = "hsieh", effect_size_input_mode = "odds_ratio", tails = "two", x_distribution = "binomial"),
       shown = list("Critical z" = "1.9599640", "Power (1-\u03b2 err prob)" = "0.2941273")),
  list(id = "GPW-Z-LOGIT-POS-B218-O01--r01", method = "z_logistic_regression", mode = "post_hoc", target = NULL,
       inputs = list(odds_ratio = 1.5, p1 = 0.5, alpha = 0.05, n_total = 200L, r2_other_x = 0, x_mu = 0, x_sigma = 1),
       options = list(computation = "large_sample", test_procedure = "hsieh", effect_size_input_mode = "odds_ratio", tails = "two", x_distribution = "normal"),
       shown = list("Critical z" = "1.9599640", "Power (1-\u03b2 err prob)" = "0.8178250")),
  list(id = "GPW-Z-LOGIT-POS-B066-O01--r01", method = "z_logistic_regression", mode = "post_hoc", target = NULL,
       inputs = list(enumeration_threshold = 100000L, odds_ratio = 1.5, p1 = 0.5, alpha = 0.05, n_total = 200L, r2_other_x = 0, x_mu = 0, x_sigma = 1),
       options = list(computation = "enumeration", test_statistic = "wald", effect_size_input_mode = "odds_ratio", tails = "one", x_distribution = "normal"),
       shown = list("Noncentrality parameter \u03bb" = "7.2446808", "Critical \u03c7\u00b2" = "2.7055435", Df = "1", "Power (1-\u03b2 err prob)" = "0.8523977")),
  list(id = "GPW-Z-LOGIT-POS-B111-D01--r01", method = "z_logistic_regression", mode = "post_hoc", target = NULL,
       inputs = list(enumeration_threshold = 100000L, odds_ratio = 1.1, p1 = 0.3, alpha = 1e-04, n_total = 3000L, r2_other_x = 0.5, x_left = 0, x_right = 1),
       options = list(computation = "enumeration", test_statistic = "wald", effect_size_input_mode = "odds_ratio", tails = "one", x_distribution = "uniform"),
       shown = list("Noncentrality parameter \u03bb" = "0.2427738", "Critical \u03c7\u00b2" = "13.8310836", Df = "1", "Power (1-\u03b2 err prob)" = "0.0006396893")),
  list(id = "GPW-Z-LOGIT-POS-B091-O01--r01", method = "z_logistic_regression", mode = "post_hoc", target = NULL,
       inputs = list(enumeration_threshold = 100000L, odds_ratio = 1.5, p1 = 0.5, alpha = 0.05, n_total = 200L, r2_other_x = 0, x_lambda = 1),
       options = list(computation = "enumeration", test_statistic = "wald", effect_size_input_mode = "odds_ratio", tails = "one", x_distribution = "poisson"),
       shown = list("Noncentrality parameter \u03bb" = "6.5312908", "Critical \u03c7\u00b2" = "2.7055435", Df = "1", "Power (1-\u03b2 err prob)" = "0.8188091")),
  list(id = "GPW-Z-LOGIT-POS-B046-E01--r01", method = "z_logistic_regression", mode = "post_hoc", target = NULL,
       inputs = list(enumeration_threshold = 100000L, odds_ratio = 0.5, p1 = 0.5, alpha = 0.05, n_total = 20L, r2_other_x = 0.9, x_mu = 0, x_sigma = 2),
       options = list(computation = "enumeration", test_statistic = "wald", effect_size_input_mode = "odds_ratio", tails = "one", x_distribution = "lognormal"),
       shown = list("Noncentrality parameter \u03bb" = "0.2656023", "Critical \u03c7\u00b2" = "2.7055435", Df = "1", "Power (1-\u03b2 err prob)" = "0.1447240")),
  list(id = "GPW-Z-LOGIT-CRI-B016-O01--r01", method = "z_logistic_regression", mode = "criterion", target = NULL,
       inputs = list(odds_ratio = 1.5, p1 = 0.5, power = 0.8, n_total = 200L, r2_other_x = 0, x_pi = 0.5),
       options = list(computation = "large_sample", test_procedure = "demidenko", variance_correction = "on", effect_size_input_mode = "odds_ratio", tails = "one", x_distribution = "binomial"),
       shown = list("Critical z" = "0.5809800", "\u03b1 err prob" = "0.2806270")),
  list(id = "GPW-Z-LOGIT-CRI-B061-O01--r01", method = "z_logistic_regression", mode = "criterion", target = NULL,
       inputs = list(odds_ratio = 1.5, p1 = 0.5, power = 0.8, n_total = 200L, r2_other_x = 0, x_mu = 0, x_sigma = 1),
       options = list(computation = "large_sample", test_procedure = "demidenko", variance_correction = "on", effect_size_input_mode = "odds_ratio", tails = "one", x_distribution = "lognormal"),
       shown = list("Critical z" = "2.3885539", "\u03b1 err prob" = "0.008457413")),
  list(id = "GPW-Z-LOGIT-CRI-B066-O01--r01", method = "z_logistic_regression", mode = "criterion", target = NULL,
       inputs = list(enumeration_threshold = 100000L, odds_ratio = 1.5, p1 = 0.5, power = 0.8, n_total = 200L, r2_other_x = 0, x_mu = 0, x_sigma = 1),
       options = list(computation = "enumeration", test_statistic = "wald", effect_size_input_mode = "odds_ratio", tails = "one", x_distribution = "normal"),
       shown = list("Noncentrality parameter \u03bb" = "7.2446808", "Critical \u03c7\u00b2" = "3.4224379", Df = "1", "\u03b1 err prob" = "0.0321580")),
  list(id = "GPW-Z-LOGIT-CRI-B006-O01--r01", method = "z_logistic_regression", mode = "criterion", target = NULL,
       inputs = list(enumeration_threshold = 100000L, odds_ratio = 1.5, p1 = 0.5, power = 0.8, n_total = 200L, r2_other_x = 0, x_pi = 0.5),
       options = list(computation = "enumeration", test_statistic = "likelihood_ratio", effect_size_input_mode = "odds_ratio", tails = "one", x_distribution = "binomial"),
       shown = list("Noncentrality parameter \u03bb" = "2.0237560", "Critical \u03c7\u00b2" = "0.4195521", Df = "1", "\u03b1 err prob" = "0.2585803")),
  list(id = "GPW-Z-LOGIT-COM-B016-O01--r01", method = "z_logistic_regression", mode = "compromise", target = NULL,
       inputs = list(odds_ratio = 1.5, p1 = 0.5, error_ratio = 1, n_total = 200L, r2_other_x = 0, x_pi = 0.5),
       options = list(computation = "large_sample", test_procedure = "demidenko", variance_correction = "on", effect_size_input_mode = "odds_ratio", tails = "one", x_distribution = "binomial"),
       shown = list("Critical z" = "0.7110082", "\u03b1 err prob" = "0.2385396", "\u03b2 err prob" = "0.2385396", "Power (1-\u03b2 err prob)" = "0.7614604")),
  list(id = "GPW-Z-LOGIT-COM-B173-O01--r01", method = "z_logistic_regression", mode = "compromise", target = NULL,
       inputs = list(odds_ratio = 1.5, p1 = 0.5, error_ratio = 1, n_total = 200L, r2_other_x = 0, x_lambda = 1),
       options = list(computation = "large_sample", test_procedure = "demidenko", variance_correction = "on", effect_size_input_mode = "odds_ratio", tails = "two", x_distribution = "exponential"),
       shown = list("Critical z" = "1.4632425", "\u03b1 err prob" = "0.1434010", "\u03b2 err prob" = "0.1434010", "Power (1-\u03b2 err prob)" = "0.8565990")),
  list(id = "GPW-Z-LOGIT-COM-B066-O01--r01", method = "z_logistic_regression", mode = "compromise", target = NULL,
       inputs = list(enumeration_threshold = 100000L, odds_ratio = 1.5, p1 = 0.5, error_ratio = 1, n_total = 200L, r2_other_x = 0, x_mu = 0, x_sigma = 1),
       options = list(computation = "enumeration", test_statistic = "wald", effect_size_input_mode = "odds_ratio", tails = "one", x_distribution = "normal"),
       shown = list("Noncentrality parameter \u03bb" = "7.2446808", "Critical \u03c7\u00b2" = "1.8113956", Df = "1", "\u03b1 err prob" = "0.0891705", "\u03b2 err prob" = "0.0891705", "Power (1-\u03b2 err prob)" = "0.9108295")),
  list(id = "GPW-Z-LOGIT-COM-B138-D01--r01", method = "z_logistic_regression", mode = "compromise", target = NULL,
       inputs = list(enumeration_threshold = 100000L, odds_ratio = 1.1, p1 = 0.3, error_ratio = 0.01, n_total = 3000L, r2_other_x = 0.5, x_pi = 0.5),
       options = list(computation = "enumeration", test_statistic = "likelihood_ratio", effect_size_input_mode = "odds_ratio", tails = "two", x_distribution = "binomial"),
       shown = list("Noncentrality parameter \u03bb" = "0.7286077", "Critical \u03c7\u00b2" = "0.0003163361", Df = "1", "\u03b1 err prob" = "0.9858097", "\u03b2 err prob" = "0.009858097", "Power (1-\u03b2 err prob)" = "0.9901419")),
  list(id = "GPW-Z-LOGIT-APR-B001-I03--r01", method = "z_logistic_regression", mode = "a_priori", target = NULL,
       inputs = list(enumeration_threshold = 100000L, odds_ratio = 1.5, p1 = 0.5, alpha = 0.05, power = 0.04, r2_other_x = 0, x_pi = 0.5),
       options = list(computation = "enumeration", test_statistic = "wald", effect_size_input_mode = "odds_ratio", tails = "one", x_distribution = "binomial"),
       shown = list("Noncentrality parameter \u03bb" = "?", "Critical \u03c7\u00b2" = "?", Df = "?", "Total sample size" = "?", "Actual power" = "?"), message = "Power must be greater than or equal to alpha."),
  list(id = "GPW-Z-LOGIT-POS-B001-I02--r01", method = "z_logistic_regression", mode = "post_hoc", target = NULL,
       inputs = list(enumeration_threshold = 100000L, odds_ratio = 1.5, p1 = 1.2, alpha = 0.05, n_total = 200L, r2_other_x = 0, x_pi = 0.5),
       options = list(computation = "enumeration", test_statistic = "wald", effect_size_input_mode = "odds_ratio", tails = "one", x_distribution = "binomial"),
       shown = list("Noncentrality parameter \u03bb" = "?", "Critical \u03c7\u00b2" = "?", Df = "?", "Power (1-\u03b2 err prob)" = "?"), message = "Proportion 1 is outside the admissible range (1e-06, 1-1e-06)."),
  list(id = "GPW-Z-LOGIT-SEN-B017-O01--r01", method = "z_logistic_regression", mode = "sensitivity", target = "odds_ratio",
       inputs = list(alpha = 0.05, p1 = 0.5, power = 0.8, n_total = 200L, r2_other_x = 0, x_pi = 0.5),
       options = list(computation = "large_sample", test_procedure = "demidenko", variance_correction = "on", effect_size_input_mode = "odds_ratio", tails = "one", effect_direction = "p2_le_p1", x_distribution = "binomial"),
       shown = list("Critical z" = "-1.6448536", "Odds ratio" = "0.4851898", "Actual power" = "0.8000000")),
  list(id = "GPW-Z-LOGIT-SEN-B018-O01--r01", method = "z_logistic_regression", mode = "sensitivity", target = "p2",
       inputs = list(alpha = 0.05, p1 = 0.5, power = 0.8, n_total = 200L, r2_other_x = 0, x_pi = 0.5),
       options = list(computation = "large_sample", test_procedure = "demidenko", variance_correction = "on", effect_size_input_mode = "two_probabilities", tails = "one", effect_direction = "p2_le_p1", x_distribution = "binomial"),
       shown = list("Critical z" = "-1.6448536", "Pr(Y=1|X=1) H1" = "0.3266854", "Actual power" = "0.8000000")),
  list(id = "GPW-Z-LOGIT-SEN-B012-O01--r01", method = "z_logistic_regression", mode = "sensitivity", target = "odds_ratio",
       inputs = list(alpha = 0.05, p1 = 0.5, power = 0.8, n_total = 200L, r2_other_x = 0, x_pi = 0.5),
       options = list(computation = "large_sample", test_procedure = "demidenko", variance_correction = "off", effect_size_input_mode = "odds_ratio", tails = "one", effect_direction = "p2_le_p1", x_distribution = "binomial"),
       shown = list("Critical z" = "-1.6448536", "Odds ratio" = "0.4832796", "Actual power" = "0.8000000")),
  list(id = "GPW-Z-LOGIT-SEN-B044-O01--r01", method = "z_logistic_regression", mode = "sensitivity", target = "odds_ratio",
       inputs = list(alpha = 0.05, p1 = 0.5, power = 0.8, n_total = 200L, r2_other_x = 0, x_lambda = 1),
       options = list(computation = "large_sample", test_procedure = "demidenko", variance_correction = "on", effect_size_input_mode = "odds_ratio", tails = "one", effect_direction = "p2_ge_p1", x_distribution = "exponential"),
       shown = list("Critical z" = "1.6448536", "Odds ratio" = "1.4982545", "Actual power" = "0.8000000")),
  list(id = "GPW-Z-LOGIT-SEN-B062-D01--r01", method = "z_logistic_regression", mode = "sensitivity", target = "odds_ratio",
       inputs = list(alpha = 0.001, p1 = 0.3, power = 0.99, n_total = 100000L, r2_other_x = 0.5, x_mu = 0, x_sigma = 0.5),
       options = list(computation = "large_sample", test_procedure = "demidenko", variance_correction = "on", effect_size_input_mode = "odds_ratio", tails = "one", effect_direction = "p2_le_p1", x_distribution = "lognormal"),
       shown = list("Critical z" = "-3.0902323", "Odds ratio" = "0.9131431", "Actual power" = "0.9900000")),
  list(id = "GPW-Z-LOGIT-SEN-B024-O01--r01", method = "z_logistic_regression", mode = "sensitivity", target = "odds_ratio",
       inputs = list(alpha = 0.05, p1 = 0.5, power = 0.8, n_total = 200L, r2_other_x = 0, x_pi = 0.5),
       options = list(computation = "large_sample", test_procedure = "hsieh", effect_size_input_mode = "odds_ratio", tails = "one", effect_direction = "p2_ge_p1", x_distribution = "binomial"),
       shown = list("Critical z" = "1.6448536", "Odds ratio" = "2.0514562", "Actual power" = "0.8000000")),
  # Unattainable target (power-maximising effect, manual 29.3): the maximum is flat and the
  # effect agrees to about 6 digits only (triage L08), so only the outputs at it are compared.
  list(id = "GPW-Z-LOGIT-SEN-B012-E01--r01", method = "z_logistic_regression", mode = "sensitivity", target = "odds_ratio",
       inputs = list(alpha = 0.05, p1 = 0.05, power = 0.95, n_total = 20L, r2_other_x = 0.9, x_pi = 0.02),
       options = list(computation = "large_sample", test_procedure = "demidenko", variance_correction = "off", effect_size_input_mode = "odds_ratio", tails = "one", effect_direction = "p2_le_p1", x_distribution = "binomial"),
       shown = list("Critical z" = "-1.6448536", "Actual power" = "0.0535486")),
  list(id = "DOC-150--r01", method = "z_logistic_regression", mode = "a_priori", target = NULL,
       inputs = list(odds_ratio = 1.5, p1 = 0.5, alpha = 0.05, power = 0.95, r2_other_x = 0, x_mu = 0, x_sigma = 1),
       options = list(computation = "large_sample", test_procedure = "hsieh", effect_size_input_mode = "odds_ratio", tails = "two", x_distribution = "normal"),
       shown = list("Critical z" = "1.959964", "Total sample size" = "317", "Actual power" = "0.950486")),
  list(id = "DOC-151--r01", method = "z_logistic_regression", mode = "a_priori", target = NULL,
       inputs = list(odds_ratio = 1.5, p1 = 0.5, alpha = 0.05, power = 0.95, r2_other_x = 0, x_mu = 0, x_sigma = 1),
       options = list(computation = "large_sample", test_procedure = "demidenko", variance_correction = "on", effect_size_input_mode = "odds_ratio", tails = "two", x_distribution = "normal"),
       shown = list("Total sample size" = "337")),
  list(id = "DOC-152--r01", method = "z_logistic_regression", mode = "a_priori", target = NULL,
       inputs = list(odds_ratio = 1.5, p1 = 0.5, alpha = 0.05, power = 0.95, r2_other_x = 0, x_mu = 0, x_sigma = 1),
       options = list(computation = "large_sample", test_procedure = "demidenko", variance_correction = "off", effect_size_input_mode = "odds_ratio", tails = "two", x_distribution = "normal"),
       shown = list("Total sample size" = "355")),
  list(id = "DOC-156--r01", method = "z_logistic_regression", mode = "a_priori", target = NULL,
       inputs = list(p2 = 0.1, p1 = 0.05, alpha = 0.05, power = 0.95, r2_other_x = 0, x_pi = 0.5),
       options = list(computation = "large_sample", test_procedure = "demidenko", variance_correction = "off", effect_size_input_mode = "two_probabilities", tails = "two", x_distribution = "binomial"),
       shown = list("Total sample size" = "1498")),
  list(id = "DOC-157--r01", method = "z_logistic_regression", mode = "a_priori", target = NULL,
       inputs = list(enumeration_threshold = 10000L, p2 = 0.1, p1 = 0.05, alpha = 0.05, power = 0.95, r2_other_x = 0, x_pi = 0.5),
       options = list(computation = "enumeration", test_statistic = "wald", effect_size_input_mode = "two_probabilities", tails = "two", x_distribution = "binomial"),
       shown = list("Total sample size" = "1498")),
  list(id = "DOC-158--r01", method = "z_logistic_regression", mode = "a_priori", target = NULL,
       inputs = list(p2 = 0.1, p1 = 0.05, alpha = 0.05, power = 0.95, r2_other_x = 0, x_pi = 0.2),
       options = list(computation = "large_sample", test_procedure = "demidenko", variance_correction = "on", effect_size_input_mode = "two_probabilities", tails = "two", x_distribution = "binomial"),
       shown = list("Total sample size" = "2158")),
  list(id = "DOC-159--r01", method = "z_logistic_regression", mode = "a_priori", target = NULL,
       inputs = list(p2 = 0.1, p1 = 0.05, alpha = 0.05, power = 0.95, r2_other_x = 0, x_pi = 0.8),
       options = list(computation = "large_sample", test_procedure = "demidenko", variance_correction = "on", effect_size_input_mode = "two_probabilities", tails = "two", x_distribution = "binomial"),
       shown = list("Total sample size" = "2368")),
  list(id = "GPW-Z-POIS-APR-B010-O01--r01", method = "z_poisson_regression", mode = "a_priori", target = NULL,
       inputs = list(exp_beta1 = 1.3, alpha = 0.05, power = 0.95, base_rate = 0.85, mean_exposure = 1, r2_other_x = 0, x_pi = 0.5),
       options = list(computation = "large_sample", test_procedure = "demidenko", variance_correction = "on", tails = "one", x_distribution = "binomial"),
       shown = list("Critical z" = "1.6448536", "Total sample size" = "649", "Actual power" = "0.9500244")),
  list(id = "GPW-Z-POIS-APR-B013-O01--r01", method = "z_poisson_regression", mode = "a_priori", target = NULL,
       inputs = list(exp_beta1 = 1.3, alpha = 0.05, power = 0.95, base_rate = 0.85, mean_exposure = 1, r2_other_x = 0, x_pi = 0.5),
       options = list(computation = "large_sample", test_procedure = "signorini", variance_correction = "on", tails = "one", x_distribution = "binomial"),
       shown = list("Critical z" = "1.6448536", "Total sample size" = "697", "Actual power" = "0.9501214")),
  list(id = "GPW-Z-POIS-APR-B055-O01--r01", method = "z_poisson_regression", mode = "a_priori", target = NULL,
       inputs = list(exp_beta1 = 1.3, alpha = 0.05, power = 0.8, base_rate = 0.85, mean_exposure = 1, r2_other_x = 0, x_mu = 0, x_sigma = 1),
       options = list(computation = "large_sample", test_procedure = "demidenko", variance_correction = "on", tails = "one", x_distribution = "normal"),
       shown = list("Critical z" = "1.6448536", "Total sample size" = "103", "Actual power" = "0.8030738")),
  list(id = "GPW-Z-POIS-APR-B092-O01--r01", method = "z_poisson_regression", mode = "a_priori", target = NULL,
       inputs = list(exp_beta1 = 1.3, alpha = 0.05, power = 0.8, base_rate = 0.85, mean_exposure = 1, r2_other_x = 0, x_var_h0 = 4, x_var_h1 = 4.17),
       options = list(computation = "large_sample", test_procedure = "demidenko", variance_correction = "on", tails = "one", x_distribution = "manual"),
       shown = list("Critical z" = "1.6448536", "Total sample size" = "370", "Actual power" = "0.8006361")),
  list(id = "GPW-Z-POIS-APR-B085-O01--r01", method = "z_poisson_regression", mode = "a_priori", target = NULL,
       inputs = list(exp_beta1 = 1.3, alpha = 0.05, power = 0.8, base_rate = 0.85, mean_exposure = 1, r2_other_x = 0, x_left = 0, x_right = 1),
       options = list(computation = "large_sample", test_procedure = "demidenko", variance_correction = "on", tails = "one", x_distribution = "uniform"),
       shown = list("Critical z" = "1.6448536", "Total sample size" = "1112", "Actual power" = "0.8001702")),
  list(id = "GPW-Z-POIS-POS-B148-O01--r01", method = "z_poisson_regression", mode = "post_hoc", target = NULL,
       inputs = list(exp_beta1 = 0.90483741803596, alpha = 0.05, n_total = 200L, base_rate = 1.64872127070013, mean_exposure = 1, r2_other_x = 0, x_mu = 0, x_sigma = 1),
       options = list(computation = "large_sample", test_procedure = "demidenko", variance_correction = "on", tails = "two", x_distribution = "normal"),
       shown = list("Critical z" = "-1.9599640", "Power (1-\u03b2 err prob)" = "0.4445930")),
  list(id = "GPW-Z-POIS-POS-B151-O01--r01", method = "z_poisson_regression", mode = "post_hoc", target = NULL,
       inputs = list(exp_beta1 = 0.90483741803596, alpha = 0.05, n_total = 200L, base_rate = 1.64872127070013, mean_exposure = 1, r2_other_x = 0, x_mu = 0, x_sigma = 1),
       options = list(computation = "large_sample", test_procedure = "signorini", variance_correction = "on", tails = "two", x_distribution = "normal"),
       shown = list("Critical z" = "-1.9599640", "Power (1-\u03b2 err prob)" = "0.4426539")),
  list(id = "GPW-Z-POIS-POS-B100-O01--r01", method = "z_poisson_regression", mode = "post_hoc", target = NULL,
       inputs = list(exp_beta1 = 0.818730753077982, alpha = 0.05, n_total = 200L, base_rate = 1.64872127070013, mean_exposure = 1, r2_other_x = 0, x_pi = 0.2),
       options = list(computation = "large_sample", test_procedure = "demidenko", variance_correction = "off", tails = "two", x_distribution = "binomial"),
       shown = list("Critical z" = "-1.9599640", "Power (1-\u03b2 err prob)" = "0.2677840")),
  list(id = "GPW-Z-POIS-POS-B103-O01--r01", method = "z_poisson_regression", mode = "post_hoc", target = NULL,
       inputs = list(exp_beta1 = 0.818730753077982, alpha = 0.05, n_total = 200L, base_rate = 1.64872127070013, mean_exposure = 1, r2_other_x = 0, x_pi = 0.2),
       options = list(computation = "large_sample", test_procedure = "demidenko", variance_correction = "on", tails = "two", x_distribution = "binomial"),
       shown = list("Critical z" = "-1.9599640", "Power (1-\u03b2 err prob)" = "0.2543882")),
  list(id = "GPW-Z-POIS-POS-B106-O01--r01", method = "z_poisson_regression", mode = "post_hoc", target = NULL,
       inputs = list(exp_beta1 = 0.818730753077982, alpha = 0.05, n_total = 200L, base_rate = 1.64872127070013, mean_exposure = 1, r2_other_x = 0, x_pi = 0.2),
       options = list(computation = "large_sample", test_procedure = "signorini", variance_correction = "on", tails = "two", x_distribution = "binomial"),
       shown = list("Critical z" = "-1.9599640", "Power (1-\u03b2 err prob)" = "0.3208871")),
  list(id = "GPW-Z-POIS-POS-B118-O02--r01", method = "z_poisson_regression", mode = "post_hoc", target = NULL,
       inputs = list(exp_beta1 = 0.49658530379141, alpha = 0.05, n_total = 200L, base_rate = 1.64872127070013, mean_exposure = 1, r2_other_x = 0, x_lambda = 3),
       options = list(computation = "large_sample", test_procedure = "demidenko", variance_correction = "on", tails = "two", x_distribution = "exponential"),
       shown = list("Critical z" = "-1.9599640", "Power (1-\u03b2 err prob)" = "0.9189405")),
  list(id = "GPW-Z-POIS-POS-B121-O02--r01", method = "z_poisson_regression", mode = "post_hoc", target = NULL,
       inputs = list(exp_beta1 = 0.49658530379141, alpha = 0.05, n_total = 200L, base_rate = 1.64872127070013, mean_exposure = 1, r2_other_x = 0, x_lambda = 3),
       options = list(computation = "large_sample", test_procedure = "signorini", variance_correction = "on", tails = "two", x_distribution = "exponential"),
       shown = list("Critical z" = "-1.9599640", "Power (1-\u03b2 err prob)" = "0.9517969")),
  list(id = "GPW-Z-POIS-POS-B185-O01--r01", method = "z_poisson_regression", mode = "post_hoc", target = NULL,
       inputs = list(exp_beta1 = 1.3, alpha = 0.05, n_total = 200L, base_rate = 0.85, mean_exposure = 1, r2_other_x = 0, x_var_h0 = 4, x_var_h1 = 4.17),
       options = list(computation = "large_sample", test_procedure = "demidenko", variance_correction = "on", tails = "two", x_distribution = "manual"),
       shown = list("Critical z" = "1.9599640", "Power (1-\u03b2 err prob)" = "0.4420244")),
  list(id = "GPW-Z-POIS-POS-B186-O01--r01", method = "z_poisson_regression", mode = "post_hoc", target = NULL,
       inputs = list(exp_beta1 = 1.3, alpha = 0.05, n_total = 200L, base_rate = 0.85, mean_exposure = 1, r2_other_x = 0, x_var_h0 = 4, x_var_h1 = 4.17),
       options = list(computation = "large_sample", test_procedure = "signorini", variance_correction = "on", tails = "two", x_distribution = "manual"),
       shown = list("Critical z" = "1.9599640", "Power (1-\u03b2 err prob)" = "0.4592292")),
  list(id = "GPW-Z-POIS-POS-B067-O01--r01", method = "z_poisson_regression", mode = "post_hoc", target = NULL,
       inputs = list(exp_beta1 = 1.3, alpha = 0.05, n_total = 200L, base_rate = 0.85, mean_exposure = 1, r2_other_x = 0, x_lambda = 1),
       options = list(computation = "large_sample", test_procedure = "demidenko", variance_correction = "off", tails = "one", x_distribution = "poisson"),
       shown = list("Critical z" = "1.6448536", "Power (1-\u03b2 err prob)" = "0.9980533")),
  list(id = "GPW-Z-POIS-CRI-B013-O01--r01", method = "z_poisson_regression", mode = "criterion", target = NULL,
       inputs = list(exp_beta1 = 1.3, power = 0.8, n_total = 200L, base_rate = 0.85, mean_exposure = 1, r2_other_x = 0, x_pi = 0.5),
       options = list(computation = "large_sample", test_procedure = "signorini", variance_correction = "on", tails = "one", x_distribution = "binomial"),
       shown = list("Critical z" = "0.9188274", "\u03b1 err prob" = "0.1790929")),
  list(id = "GPW-Z-POIS-CRI-B055-O01--r01", method = "z_poisson_regression", mode = "criterion", target = NULL,
       inputs = list(exp_beta1 = 1.3, power = 0.8, n_total = 200L, base_rate = 0.85, mean_exposure = 1, r2_other_x = 0, x_mu = 0, x_sigma = 1),
       options = list(computation = "large_sample", test_procedure = "demidenko", variance_correction = "on", tails = "one", x_distribution = "normal"),
       shown = list("Critical z" = "2.6385676", "\u03b1 err prob" = "0.004162855")),
  list(id = "GPW-Z-POIS-CRI-B085-O01--r01", method = "z_poisson_regression", mode = "criterion", target = NULL,
       inputs = list(exp_beta1 = 1.3, power = 0.8, n_total = 200L, base_rate = 0.85, mean_exposure = 1, r2_other_x = 0, x_left = 0, x_right = 1),
       options = list(computation = "large_sample", test_procedure = "demidenko", variance_correction = "on", tails = "one", x_distribution = "uniform"),
       shown = list("Critical z" = "0.2139690", "\u03b1 err prob" = "0.4152856")),
  list(id = "GPW-Z-POIS-COM-B010-O01--r01", method = "z_poisson_regression", mode = "compromise", target = NULL,
       inputs = list(exp_beta1 = 1.3, error_ratio = 1, n_total = 200L, base_rate = 0.85, mean_exposure = 1, r2_other_x = 0, x_pi = 0.5),
       options = list(computation = "large_sample", test_procedure = "demidenko", variance_correction = "on", tails = "one", x_distribution = "binomial"),
       shown = list("Critical z" = "0.9131686", "\u03b1 err prob" = "0.1805769", "\u03b2 err prob" = "0.1805769", "Power (1-\u03b2 err prob)" = "0.8194231")),
  list(id = "GPW-Z-POIS-COM-B055-O01--r01", method = "z_poisson_regression", mode = "compromise", target = NULL,
       inputs = list(exp_beta1 = 1.3, error_ratio = 1, n_total = 200L, base_rate = 0.85, mean_exposure = 1, r2_other_x = 0, x_mu = 0, x_sigma = 1),
       options = list(computation = "large_sample", test_procedure = "demidenko", variance_correction = "on", tails = "one", x_distribution = "normal"),
       shown = list("Critical z" = "1.7400944", "\u03b1 err prob" = "0.0409212", "\u03b2 err prob" = "0.0409212", "Power (1-\u03b2 err prob)" = "0.9590788")),
  list(id = "GPW-Z-POIS-COM-B092-O01--r01", method = "z_poisson_regression", mode = "compromise", target = NULL,
       inputs = list(exp_beta1 = 1.3, error_ratio = 1, n_total = 200L, base_rate = 0.85, mean_exposure = 1, r2_other_x = 0, x_var_h0 = 4, x_var_h1 = 4.17),
       options = list(computation = "large_sample", test_procedure = "demidenko", variance_correction = "on", tails = "one", x_distribution = "manual"),
       shown = list("Critical z" = "0.9179461", "\u03b1 err prob" = "0.1793236", "\u03b2 err prob" = "0.1793236", "Power (1-\u03b2 err prob)" = "0.8206764")),
  list(id = "GPW-Z-POIS-POS-B001-I02--r01", method = "z_poisson_regression", mode = "post_hoc", target = NULL,
       inputs = list(enumeration_threshold = 100000L, exp_beta1 = 1.3, alpha = 0.05, n_total = 200L, base_rate = -0.5, mean_exposure = 1, r2_other_x = 0, x_pi = 0.5),
       options = list(computation = "enumeration", test_statistic = "wald", tails = "one", x_distribution = "binomial"),
       shown = list("Noncentrality parameter \u03bb" = "?", "Critical \u03c7\u00b2" = "?", Df = "?", "Power (1-\u03b2 err prob)" = "?"), message = "The base rate exp(b0) must be > 0."),
  list(id = "GPW-Z-POIS-SEN-B053-O01--r01", method = "z_poisson_regression", mode = "sensitivity", target = "exp_beta1",
       inputs = list(alpha = 0.05, power = 0.8, n_total = 200L, base_rate = 0.85, mean_exposure = 1, r2_other_x = 0, x_mu = 0, x_sigma = 1),
       options = list(computation = "large_sample", test_procedure = "demidenko", variance_correction = "off", tails = "one", effect_direction = "beta1_le_0", x_distribution = "normal"),
       shown = list("Critical z" = "-1.6448536", "Exp(\u03b21)" = "0.8277796")),
  list(id = "GPW-Z-POIS-SEN-B054-O01--r01", method = "z_poisson_regression", mode = "sensitivity", target = "exp_beta1",
       inputs = list(alpha = 0.05, power = 0.8, n_total = 200L, base_rate = 0.85, mean_exposure = 1, r2_other_x = 0, x_mu = 0, x_sigma = 1),
       options = list(computation = "large_sample", test_procedure = "demidenko", variance_correction = "off", tails = "one", effect_direction = "beta1_ge_0", x_distribution = "normal"),
       shown = list("Critical z" = "1.6448536", "Exp(\u03b21)" = "1.2080510")),
  list(id = "GPW-Z-POIS-SEN-B011-O01--r01", method = "z_poisson_regression", mode = "sensitivity", target = "exp_beta1",
       inputs = list(alpha = 0.05, power = 0.8, n_total = 200L, base_rate = 0.85, mean_exposure = 1, r2_other_x = 0, x_pi = 0.5),
       options = list(computation = "large_sample", test_procedure = "demidenko", variance_correction = "on", tails = "one", effect_direction = "beta1_le_0", x_distribution = "binomial"),
       shown = list("Critical z" = "-1.6448536", "Exp(\u03b21)" = "0.6532347")),
  list(id = "GPW-Z-POIS-SEN-B086-O01--r01", method = "z_poisson_regression", mode = "sensitivity", target = "exp_beta1",
       inputs = list(alpha = 0.05, power = 0.8, n_total = 200L, base_rate = 0.85, mean_exposure = 1, r2_other_x = 0, x_left = 0, x_right = 1),
       options = list(computation = "large_sample", test_procedure = "demidenko", variance_correction = "on", tails = "one", effect_direction = "beta1_le_0", x_distribution = "uniform"),
       shown = list("Critical z" = "-1.6448536", "Exp(\u03b21)" = "0.4467300")),
  list(id = "DOC-160--r01", method = "z_poisson_regression", mode = "a_priori", target = NULL,
       inputs = list(exp_beta1 = 1.3, alpha = 0.05, power = 0.95, base_rate = 0.85, mean_exposure = 1, r2_other_x = 0, x_pi = 0.5),
       options = list(computation = "large_sample", test_procedure = "signorini", variance_correction = "on", tails = "one", x_distribution = "binomial"),
       shown = list("Critical z" = "1.644854", "Total sample size" = "697", "Actual power" = "0.950121")),
  list(id = "DOC-161--r01", method = "z_poisson_regression", mode = "a_priori", target = NULL,
       inputs = list(exp_beta1 = 1.3, alpha = 0.05, power = 0.95, base_rate = 0.85, mean_exposure = 1, r2_other_x = 0, x_pi = 0.5),
       options = list(computation = "large_sample", test_procedure = "demidenko", variance_correction = "on", tails = "one", x_distribution = "binomial"),
       shown = list("Total sample size" = "649")),
  list(id = "DOC-162--r01", method = "z_poisson_regression", mode = "a_priori", target = NULL,
       inputs = list(enumeration_threshold = 10000L, exp_beta1 = 1.3, alpha = 0.05, power = 0.95, base_rate = 0.85, mean_exposure = 1, r2_other_x = 0, x_pi = 0.5),
       options = list(computation = "enumeration", test_statistic = "likelihood_ratio", tails = "one", x_distribution = "binomial"),
       shown = list("Total sample size" = "649")),
  list(id = "DOC-163--r01", method = "z_poisson_regression", mode = "a_priori", target = NULL,
       inputs = list(exp_beta1 = 1.3, alpha = 0.05, power = 0.95, base_rate = 0.85, mean_exposure = 1, r2_other_x = 0, x_pi = 0.5),
       options = list(computation = "large_sample", test_procedure = "demidenko", variance_correction = "off", tails = "one", x_distribution = "binomial"),
       shown = list("Total sample size" = "655")),
  list(id = "DOC-164--r01", method = "z_poisson_regression", mode = "a_priori", target = NULL,
       inputs = list(enumeration_threshold = 10000L, exp_beta1 = 1.3, alpha = 0.05, power = 0.95, base_rate = 0.85, mean_exposure = 1, r2_other_x = 0, x_pi = 0.5),
       options = list(computation = "enumeration", test_statistic = "wald", tails = "one", x_distribution = "binomial"),
       shown = list("Total sample size" = "655")),
  list(id = "DOC-165--r01", method = "z_poisson_regression", mode = "a_priori", target = NULL,
       inputs = list(exp_beta1 = 1.3, alpha = 0.05, power = 0.95, base_rate = 0.85, mean_exposure = 1, r2_other_x = 0.1, x_pi = 0.5),
       options = list(computation = "large_sample", test_procedure = "signorini", variance_correction = "on", tails = "one", x_distribution = "binomial"),
       shown = list("Total sample size" = "774"))
)

g6_check_fixture <- function(fx) {
  if (!is.null(fx$message)) {
    # G*Power refused with an 'Input error' (TD-0010): the engine refuses too,
    # with a typed status and no number.
    res <- tryCatch(calculate_power(g6_fx_spec(fx)), inferguidR_error = function(e) e)
    code <- if (inherits(res, "inferguidR_error")) res$status$code else res$status$code
    expect_true(code %in% c("invalid_input", "infeasible"), info = paste(fx$id, fx$message, code))
    return(invisible())
  }
  res <- if (identical(fx$mode, "sensitivity")) g6_sensitivity_result(fx) else calculate_power(g6_fx_spec(fx))
  if (!identical(fx$mode, "sensitivity")) expect_identical(res$status$code, "ok", info = fx$id)
  ours <- g6_outputs(res)
  for (lab in names(fx$shown)) {
    v <- ours[[lab]]
    expect_true(!is.null(v) && g6_display_ok(as.double(v), fx$shown[[lab]]),
                info = sprintf("%s %s: G*Power %s, engine %s", fx$id, lab, fx$shown[[lab]],
                               if (is.null(v)) "none" else format(v, digits = 12)))
  }
}

test_that("a priori fixtures are reproduced (large sample, enumeration, documented examples)", {
  for (fx in Filter(function(f) identical(f$mode, "a_priori"), g6_fixtures)) g6_check_fixture(fx)
})

test_that("post hoc fixtures are reproduced", {
  for (fx in Filter(function(f) identical(f$mode, "post_hoc"), g6_fixtures)) g6_check_fixture(fx)
})

test_that("criterion and compromise fixtures are reproduced", {
  for (fx in Filter(function(f) f$mode %in% c("criterion", "compromise"), g6_fixtures)) g6_check_fixture(fx)
})

test_that("sensitivity fixtures are reproduced by the method's sensitivity search", {
  for (fx in Filter(function(f) identical(f$mode, "sensitivity"), g6_fixtures)) g6_check_fixture(fx)
})

test_that("calculate_power() runs the method sensitivity search (CR-G6-01)", {
  for (fx in Filter(function(f) identical(f$mode, "sensitivity"), g6_fixtures)) {
    res <- calculate_power(g6_fx_spec(fx))
    expect_identical(res$status$code, "ok", info = fx$id)
    ours <- g6_outputs(res)
    for (lab in names(fx$shown)) expect_true(g6_display_ok(as.double(ours[[lab]]), fx$shown[[lab]]), info = paste(fx$id, lab))
  }
})

test_that("the embedded G*Power strings equal the fixture files (development checkout)", {
  root <- normalizePath(testthat::test_path("..", ".."), mustWork = FALSE)
  fx_root <- file.path(root, "validation", "fixtures")
  skip_if_not(dir.exists(fx_root), "validation fixtures not available (installed package)")
  for (fx in g6_fixtures) {
    kind <- if (startsWith(fx$id, "DOC-")) "documented_example" else "emulator_capture"
    f <- file.path(fx_root, kind, paste0(fx$id, ".json"))
    expect_true(file.exists(f), info = fx$id)
    rec <- jsonlite::parse_json(paste(readLines(f, warn = FALSE, encoding = "UTF-8"), collapse = "\n"))
    expect_identical(rec$method$method_id, fx$method, info = fx$id)
    expect_identical(rec$method$analysis_type, fx$mode, info = fx$id)
    outs <- Filter(function(o) !startsWith(o$name, "message"), rec$outputs)
    shown <- stats::setNames(lapply(outs, `[[`, "displayed"), vapply(outs, `[[`, "", "label_as_displayed"))
    expect_identical(shown[names(fx$shown)], fx$shown, info = fx$id)
    if (!is.null(fx$message)) {
      msg <- Filter(function(o) startsWith(o$name, "message"), rec$outputs)
      expect_identical(msg[[1]]$displayed, fx$message, info = fx$id)
    }
  }
})
