#!/usr/bin/env Rscript
# Copyright (c) 2026 Raban Heller. MIT licence; see ../../LICENSE.md.
# Bounded inverse and integer-design checks using standalone mathematical references.
args <- commandArgs(TRUE)
if (length(args) < 1L || length(args) > 2L || !nzchar(args[1]) ||
    (length(args) == 2L && args[2] != "--negative-control")) stop("usage: Rscript check-inverses.R NEW_OUTPUT_DIRECTORY [--negative-control]")
negative <- length(args) == 2L; injected <- FALSE
out <- path.expand(args[1])
if (dir.exists(out) && length(list.files(out, all.files = TRUE, no.. = TRUE))) stop("output directory must be empty")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
out <- normalizePath(out, mustWork = TRUE)
stopifnot(requireNamespace("inferguidR", quietly = TRUE), requireNamespace("jsonlite", quietly = TRUE))
pkg <- normalizePath(find.package("inferguidR"), winslash = "/", mustWork = TRUE)
stopifnot(file.exists(file.path(pkg, "Meta", "package.rds")))
expected <- Sys.getenv("INFERGUIDR_EXPECTED_LIBRARY")
if (nzchar(expected) && pkg != normalizePath(file.path(expected, "inferguidR"), winslash = "/", mustWork = TRUE)) stop("wrong installed candidate")
here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))))
ref <- function(file) { e <- new.env(parent = globalenv()); sys.source(file.path(here, "reference", file), e); e }
cont <- ref("continuous.R"); multi <- ref("multivariate.R")
multi$ref_pchisq <- cont$ref_pchisq
inverse_f <- ref("inverse-f.R")
glm <- ref("glm.R"); rank <- ref("rank.R"); manova <- ref("manova.R")
discrete <- ref("discrete.R"); offset <- ref("offset.R"); regression <- ref("regression.R")
variance <- ref("variance.R"); correlations <- ref("correlations.R"); proportions <- ref("proportions.R"); unconditional <- ref("unconditional.R")
proportions$ck_zq <- glm$ck_zq; proportions$ck_Phi <- glm$ck_Phi
started <- Sys.time(); rows <- list(); case_rows <- list(); current <- NULL
json <- function(x) as.character(jsonlite::toJSON(x, auto_unbox = TRUE, null = "null", digits = 17))
append_csv <- function(row, path) utils::write.table(row, path, sep = ",", row.names = FALSE,
  col.names = !file.exists(path), append = file.exists(path), qmethod = "double", na = "")
record <- function(quantity, observed, reference, absolute = 1e-8, relative = 1e-6, oracle = "independent mathematical reference") {
  if (negative && !injected && identical(quantity, "power") && is.numeric(observed) && length(observed) == 1L && is.finite(observed)) {
    observed <- observed + .01; injected <<- TRUE
  }
  valid <- is.numeric(observed) && length(observed) == 1L && is.finite(observed) &&
    is.numeric(reference) && length(reference) == 1L && is.finite(reference)
  error <- if (valid) abs(observed - reference) else NA_real_
  row <- data.frame(case = current$id, method = current$method, mode = current$mode, quantity = quantity,
    observed = if (valid) observed else NA_real_, reference = if (valid) reference else NA_real_,
    absolute_error = error, absolute_tolerance = absolute, relative_tolerance = relative,
    passed = isTRUE(valid && error <= absolute + relative * abs(reference)), oracle = oracle)
  rows[[length(rows) + 1L]] <<- row
  append_csv(row, file.path(out, "results.csv")); invisible(NULL)
}
condition <- function(quantity, value) record(quantity, as.double(isTRUE(value)), 1, 0, 0, "declared contract")
run <- function(id, method, mode, inputs, options = list(), check, operation = NULL, target = NULL, expected_status = "ok") {
  current <<- list(id = id, method = method, mode = mode)
  before <- length(rows); error <- ""; status <- "not_returned"; warnings <- NULL; resolved <- NULL; diagnostics <- NULL; t0 <- proc.time()[["elapsed"]]
  tryCatch({
    r <- if (is.null(operation)) inferguidR::calculate_power(inferguidR::analysis_spec(method, mode, inputs, options, target = target)) else operation()
    status <- r$status$code
    warnings <- r$warnings
    resolved <- r$solved; diagnostics <- r$diagnostics
    condition("expected_status", identical(status, expected_status))
    if (identical(status, "ok") && expected_status == "ok") check(r)
  }, error = function(e) {
    error <<- conditionMessage(e)
    status <<- if (!is.null(e$status$code)) e$status$code else "error"
    condition("expected_refusal", expected_status != "ok" && identical(status, expected_status))
  })
  last <- do.call(rbind, rows[seq.int(before + 1L, length(rows))])
  row <- data.frame(case = id, method = method, mode = mode, inputs_json = json(inputs), options_json = json(options),
    evidence_kind = if (expected_status == "ok") "independent_numeric" else "guard_contract", status = status, passed = all(last$passed),
    elapsed_seconds = proc.time()[["elapsed"]] - t0, warnings_json = json(warnings), solved_json = json(resolved), diagnostics_json = json(diagnostics), error = error)
  case_rows[[length(case_rows) + 1L]] <<- row
  append_csv(row, file.path(out, "cases.csv"))
  if (length(case_rows) %% 10L == 0L) cat(length(case_rows), "cases completed\n")
}

# The reference model receives resolved inputs only, never engine df, ncp,
# power, critical values or solver diagnostics. Central quantiles are standard
# base-R primitives; alternative probabilities use the separate integrals.
t_power <- function(df, ncp, alpha, tails = "two") {
  critical <- cont$ref_qt_upper(if (tails == "two") alpha / 2 else alpha, df)
  if (tails == "two") cont$ref_pt(critical, df, ncp, FALSE) + cont$ref_pt(-critical, df, ncp) else cont$ref_pt(critical, df, abs(ncp), FALSE)
}
f_power <- function(df1, df2, ncp, alpha) {
  if (!is.finite(df1) || !is.finite(df2) || min(df1, df2) <= 0) return(NA_real_)
  inverse_f$ref_inverse_pf(stats::qf(alpha, df1, df2, lower.tail = FALSE), df1, df2, ncp, FALSE)
}
models <- list()
add <- function(id, input, effect, family, parameters, options = list(), minimum = 2, step = 1, modes = "a_priori") {
  models[[id]] <<- list(input = input, effect = effect, family = family, parameters = parameters,
                       options = options, minimum = minimum, step = step, modes = modes)
}
all_inverse <- c("a_priori", "sensitivity", "criterion", "compromise")
add("t_means_one_sample", list(d = .55, n_total = 24), "d", "t", function(p) c(p$n_total - 1, p$d * sqrt(p$n_total)))
add("t_means_matched_pairs", list(dz = .55, n_total = 24), "dz", "t", function(p) c(p$n_total - 1, p$dz * sqrt(p$n_total)))
add("t_means_two_independent", list(d = .65, n1 = 22, n2 = 31), "d", "t", function(p) c(p$n1 + p$n2 - 2, p$d * sqrt(p$n1 * p$n2 / (p$n1 + p$n2))))
add("t_correlation_point_biserial", list(rho = .35, n_total = 35), "rho", "t", function(p) c(p$n_total - 2, p$rho * sqrt(p$n_total / (1 - p$rho^2))))
add("t_linreg_one_group_slope", list(slope_h1 = .35, slope_h0 = .05, sd_x = 1.2, sd_y = 1.5, n_total = 35), "slope_h1", "t",
  function(p) c(p$n_total - 2, sqrt(p$n_total) * p$sd_x * (p$slope_h1 - p$slope_h0) / (p$sd_y * sqrt(1 - (p$slope_h1 * p$sd_x / p$sd_y)^2))))
add("t_linreg_two_groups_slopes", list(delta_slope = .65, sd_residual = 1.2, sd_x1 = .8, sd_x2 = 1.3, n1 = 22, n2 = 31), "delta_slope", "t",
  function(p) c(p$n1 + p$n2 - 4, p$delta_slope / (p$sd_residual * sqrt(1 / (p$n1 * p$sd_x1^2) + 1 / (p$n2 * p$sd_x2^2)))))
add("t_linreg_two_groups_intercepts", list(delta_intercept = .8, sd_residual = 1.2, sd_x1 = .8, sd_x2 = 1.3, mean_x1 = -.3, mean_x2 = .5, n1 = 22, n2 = 31), "delta_intercept", "t",
  function(p) c(p$n1 + p$n2 - 4, p$delta_intercept / (p$sd_residual * sqrt((1 + p$mean_x1^2 / p$sd_x1^2) / p$n1 + (1 + p$mean_x2^2 / p$sd_x2^2) / p$n2))))
add("t_mrc_single_coefficient", list(f2 = .18, n_total = 35, n_predictors = 4), "f2", "t", function(p) c(p$n_total - p$n_predictors - 1, sqrt(p$f2 * p$n_total)))
add("f_anova_oneway", list(f = .35, n_total = 48, n_groups = 3), "f", "F", function(p) c(p$n_groups - 1, p$n_total - p$n_groups, p$f^2 * p$n_total))
add("f_anova_special_main_interactions", list(f = .35, n_total = 48, n_groups = 6, numerator_df = 2), "f", "F", function(p) c(p$numerator_df, p$n_total - p$n_groups, p$f^2 * p$n_total))
add("f_ancova_main_interactions", list(f = .35, n_total = 48, n_groups = 6, numerator_df = 2, n_covariates = 2), "f", "F", function(p) c(p$numerator_df, p$n_total - p$n_groups - p$n_covariates, p$f^2 * p$n_total))

models$t_means_two_independent$minimum <- 4; models$t_means_two_independent$step <- 2
models$t_correlation_point_biserial$minimum <- 3
models$t_linreg_one_group_slope$minimum <- 3
for (id in c("t_linreg_two_groups_slopes", "t_linreg_two_groups_intercepts")) { models[[id]]$minimum <- 6; models[[id]]$step <- 2 }
models$t_mrc_single_coefficient$minimum <- 6
models$f_anova_oneway$minimum <- 6; models$f_anova_oneway$step <- 3
models$f_anova_special_main_interactions$minimum <- 7
models$f_ancova_main_interactions$minimum <- 9
add("f_hotelling_t2_one_group", list(delta = .8, n_total = 30, n_response = 3), "delta", "F",
    function(p) c(p$n_response, p$n_total - p$n_response, p$n_total * p$delta^2), minimum = 4, modes = all_inverse)
add("f_hotelling_t2_two_groups", list(delta = 1.1, n1 = 20, n2 = 20, n_response = 3), "delta", "F",
    function(p) c(p$n_response, p$n1 + p$n2 - p$n_response - 1, p$n1 * p$n2 / (p$n1 + p$n2) * p$delta^2), minimum = 6, step = 2, modes = all_inverse)
add("f_mrc_r2_deviation_from_zero", list(f2 = .3, n_total = 40, n_predictors = 3), "f2", "F",
    function(p) c(p$n_predictors, p$n_total - p$n_predictors - 1, p$n_total * p$f2), minimum = 5, modes = all_inverse)
add("f_mrc_r2_increase", list(f2 = .3, n_total = 40, n_predictors = 4, n_tested = 2), "f2", "F",
    function(p) c(p$n_tested, p$n_total - p$n_predictors - 1, p$n_total * p$f2), minimum = 6, modes = all_inverse)
add("t_generic", list(delta = 2.5, df = 14), "delta", "t", function(p) c(p$df, p$delta), modes = all_inverse[-1])
add("f_generic", list(lambda = 9, df1 = 3, df2 = 24), "lambda", "F", function(p) c(p$df1, p$df2, p$lambda), modes = all_inverse[-1])
for (which in c("between", "within", "within_between")) local({
  effect <- which
  p <- list(f = .4, n_total = 36, k = 3, m = 4, rho = .3)
  if (effect != "between") p$epsilon <- .8
  parameters <- function(p) {
    eps <- if (effect == "between") 1 else p$epsilon
    df1 <- if (effect == "between") p$k - 1 else (p$m - 1) * eps * if (effect == "within") 1 else p$k - 1
    df2 <- (p$n_total - p$k) * if (effect == "between") 1 else (p$m - 1) * eps
    lambda <- p$f^2 * p$n_total * p$m * if (effect == "between") 1 / (1 + (p$m - 1) * p$rho) else eps / (1 - p$rho)
    c(df1, df2, lambda)
  }
  add(paste0("f_rm_anova_", effect), p, "f", "F", parameters, list(effect_size_specification = "gpower30"), minimum = 6, step = 3, modes = all_inverse)
})
for (id in c("f_manova_global", "f_manova_special", "f_rm_manova_between", "f_rm_manova_within", "f_rm_manova_within_between")) local({
  method <- id; repeated <- grepl("f_rm_", method)
  p <- if (repeated) list(f = .4, n_total = 36, k = 3, m = 4) else list(f2 = .2, n_total = 36, n_groups = 3, n_responses = 3)
  if (method == "f_manova_special") p$n_predictors <- 2
  opt <- list(approximation = "pillai_v", algorithm = "obrien_shieh_1999")
  if (repeated) opt$use_mean_correlation <- "off"
  parameters <- function(p) {
    q <- if (method == "f_rm_manova_within") 1 else 2
    dimension <- if (method == "f_rm_manova_between") 1 else 3
    f2 <- if (repeated) p$f^2 else p$f2
    terms <- manova$ref_terms("pillai_v", p$n_total - 3, dimension, q, f2)
    c(terms$df1, terms$df2, manova$ref_lambda("obrien_shieh_1999", p$n_total, terms, f2))
  }
  add(method, p, if (repeated) "f" else "f2", "F", parameters, opt, minimum = if (method %in% c("f_manova_special", "f_rm_manova_within_between")) 5 else 6, step = if (method %in% c("f_manova_special", "f_rm_manova_within_between")) 1 else 3, modes = all_inverse)
})
# Full enumeration of each declared, bounded integer design grid verifies
# global minimality without inferring monotonicity from two points.
for (id in names(models)) for (mode in models[[id]]$modes) {
  model <- models[[id]]; p <- c(model$input, list(alpha = .05)); opt <- model$options
  if (model$family == "t") opt$tails <- "two"
  two <- !is.null(p$n1)
  if (mode == "a_priori") { p$n_total <- p$n1 <- p$n2 <- NULL; p$power <- .6; if (two) p$allocation_ratio <- 1 }
  if (mode == "sensitivity") { p[[model$effect]] <- NULL; p$power <- .6 }
  if (mode == "criterion") { p$alpha <- NULL; p$power <- .6 }
  if (mode == "compromise") { p$alpha <- NULL; p$error_ratio <- 1 }
  run(paste("inverse-continuous", id, mode, sep = "-"), id, mode, p, opt, function(r) {
    resolved <- p
    if (mode == "sensitivity") {
      resolved[[model$effect]] <- r$solved[[model$effect]]
      if (is.null(resolved[[model$effect]])) resolved[[model$effect]] <- r$method_outputs[[model$effect]]
    }
    if (mode %in% c("criterion", "compromise")) resolved$alpha <- r$solved$alpha
    if (mode == "a_priori") {
      resolved$n_total <- r$solved$n_total
      if (two) { resolved$n1 <- unname(r$group_sizes[1]); resolved$n2 <- unname(r$group_sizes[2]) }
    }
    power <- function(v) {
      a <- model$parameters(v)
      if (model$family == "t") t_power(a[1], a[2], v$alpha) else f_power(a[1], a[2], a[3], v$alpha)
    }
    reference <- power(resolved)
    record("power", r$achieved$power, reference)
    if (mode == "a_priori") {
      n <- resolved$n_total; condition("bounded_integer_design", n <= 400 && n == floor(n))
      if (n <= 400) {
        sizes <- seq(model$minimum, n, by = model$step)
        powers <- vapply(sizes, function(k) { q <- resolved; q$n_total <- k; if (two) { q$n1 <- k / 2; q$n2 <- k / 2 }; power(q) }, 0.0)
        condition("all_grid_probabilities_finite", all(is.finite(powers)))
        hit <- which(is.finite(powers) & powers >= .6)[1]
        record("first_admissible_sample_size", n, sizes[hit], 0, 0, "exhaustive independent scan of every earlier admissible integer design")
      }
    } else if (mode == "compromise") record("nominal_error_ratio", 1 - reference, resolved$alpha) else record("reference_target", reference, .6)
  })
}

# Remaining inverse families. A profile supplies its independently derived
# forward probability and exact admissible sample-size sequence.
profile <- function(label, id, input, effect, power, options = list(), modes = all_inverse,
                    minimum = 2, step = 1, discrete_model = FALSE, target_power = .6, sensitivity_options = list()) {
  for (mode in modes) {
    p <- c(input, list(alpha = .05)); two <- !is.null(p$n1)
    opt <- if (mode == "sensitivity") c(options, sensitivity_options) else options
    if (mode == "a_priori") { p$n_total <- p$n1 <- p$n2 <- NULL; p$power <- target_power; if (two) p$allocation_ratio <- 1 }
    if (mode == "sensitivity") { p[[effect]] <- NULL; p$power <- target_power }
    if (mode == "criterion") { p$alpha <- NULL; p$power <- target_power }
    if (mode == "compromise") { p$alpha <- NULL; p$error_ratio <- 1 }
    run(paste(label, id, mode, sep = "-"), id, mode, p, opt, function(r) {
      resolved <- p
      if (mode == "sensitivity") {
        resolved[[effect]] <- r$solved[[effect]]
        if (is.null(resolved[[effect]])) resolved[[effect]] <- r$method_outputs[[effect]]
      }
      if (mode %in% c("criterion", "compromise")) resolved$alpha <- r$solved$alpha
      if (mode == "a_priori") {
        resolved$n_total <- r$solved$n_total
        if (two) { resolved$n1 <- unname(r$group_sizes[1]); resolved$n2 <- unname(r$group_sizes[2]) }
      }
      # Nominal-alpha solutions at a discrete jump sit on floating-point
      # boundaries. Locate the same boundary independently, rather than
      # comparing two differently rounded CDFs at exactly one binary value.
      if (discrete_model && mode %in% c("criterion", "compromise")) {
        probability <- function(a) { q <- resolved; q$alpha <- a; power(q)$power }
        objective <- if (mode == "criterion") function(a) probability(a) - target_power else function(a) a - (1 - probability(a))
        lo <- 1e-12; hi <- 1 - 1e-12
        condition("independent_alpha_root_bracketed", objective(lo) < 0 && objective(hi) >= 0)
        for (iteration in seq_len(60L)) {
          mid <- (lo + hi) / 2
          if (objective(mid) >= 0) hi <- mid else lo <- mid
        }
        record("first_attaining_nominal_alpha", resolved$alpha, hi, 2e-12, 1e-10,
               "independent bisection of fully enumerated rejection-region probabilities")
        resolved$alpha <- min(1 - 1e-12, hi + 8 * .Machine$double.eps * max(hi, 1e-5))
      }
      ans <- power(resolved)
      reference <- if (is.list(ans)) ans$power else ans
      value <- function(q) { a <- power(q); if (is.list(a)) a$power else a }
      record("power", r$achieved$power, reference)
      if (is.list(ans) && !is.null(ans$alpha)) record("alpha_attained", r$achieved$alpha_attained, ans$alpha)
      if (mode == "a_priori") {
        n <- resolved$n_total; condition("bounded_integer_design", n <= 400 && n == floor(n))
        if (n <= 400) {
          sizes <- seq(minimum, n, by = step)
          probabilities <- vapply(sizes, function(k) { q <- resolved; q$n_total <- k; if (two) { q$n1 <- k / 2; q$n2 <- k / 2 }; value(q) }, 0.0)
          condition("all_grid_probabilities_finite", all(is.finite(probabilities)))
          record("first_admissible_sample_size", n, sizes[which(probabilities >= target_power)[1]], 0, 0,
                 "exhaustive independent scan of every earlier admissible integer design")
        }
      } else if (!discrete_model || mode == "sensitivity") {
        if (mode == "compromise") record("nominal_error_ratio", 1 - reference, resolved$alpha) else record("reference_target", reference, target_power)
      } else {
        # For a stepwise region the contract is the first attaining boundary,
        # not equality. A left perturbation resolves a jump or continuous root.
        prior <- resolved; prior$alpha <- resolved$alpha * (1 - 1e-7)
        if (mode == "criterion") {
          condition("target_reached", reference >= target_power - 1e-8)
          condition("left_boundary_below_target", value(prior) < target_power)
        } else {
          condition("nominal_ratio_attained", 1 - reference <= resolved$alpha + 1e-8)
          condition("left_boundary_ratio_unattained", 1 - value(prior) > prior$alpha)
        }
      }
    })
  }
}
profile("chi-square", "chi2_goodness_of_fit_contingency", list(w = .5, n_total = 35, df = 4), "w",
        function(p) cont$ref_pchisq(stats::qchisq(p$alpha, p$df, lower.tail = FALSE), p$df, p$n_total * p$w^2, FALSE), modes = "a_priori", minimum = 2)
for (id in c("chi2_variance_one_sample", "f_variances_two_groups")) for (direction in c("upper", "lower")) local({
  method <- id; pair <- method == "f_variances_two_groups"
  p <- c(list(ratio = if (direction == "upper") 2.5 else .35), if (pair) list(n1 = 18, n2 = 18) else list(n_total = 24))
  opt <- list(tails = "two")
  sensitivity_options <- if (pair) list(effect_direction = if (direction == "upper") "var1_ge_var0" else "var1_le_var0") else list()
  power <- if (pair) function(p) variance$g11r_power_f2(p$ratio, p$n1, p$n2, p$alpha, "two") else function(p) variance$g11r_power_var1(p$ratio, p$n_total, p$alpha, "two")
  profile(paste0("variance-", direction), method, p, "ratio", power, opt, minimum = if (pair) 4 else 2, step = if (pair) 2 else 1,
          modes = if (!pair && direction == "lower") all_inverse[all_inverse != "sensitivity"] else all_inverse, sensitivity_options = sensitivity_options)
})
for (id in c("exact_correlation_bivariate_normal", "exact_mrc_random_model")) local({
  method <- id; corr <- id == "exact_correlation_bivariate_normal"
  p <- if (corr) list(rho_h0 = 0, rho_h1 = .7, n_total = 18) else list(rho2_h0 = 0, rho2_h1 = .5, n_total = 20, n_predictors = 2)
  power <- function(p) {
    n <- p$n_total
    if (corr) {
      critical_t <- stats::qt(p$alpha / 2, n - 2, lower.tail = FALSE)
      critical <- critical_t / sqrt(critical_t^2 + n - 2)
      regression$g5r_pr(-critical, p$rho_h1, n, TRUE) + regression$g5r_pr(critical, p$rho_h1, n, FALSE)
    } else {
      low <- stats::qbeta(p$alpha / 2, 1, (n - 3) / 2)
      high <- stats::qbeta(p$alpha / 2, 1, (n - 3) / 2, lower.tail = FALSE)
      regression$g5r_pr2(low, p$rho2_h1, 2, n, TRUE) + regression$g5r_pr2(high, p$rho2_h1, 2, n, FALSE)
    }
  }
  profile("random-design", method, p, if (corr) "rho_h1" else "rho2_h1", power,
          c(list(tails = "two"), if (corr) list(computation = "exact") else list(distribution = "exact")), modes = "a_priori", minimum = if (corr) 3 else 4)
})
for (method in c("t_wilcoxon_signed_rank_one_sample", "t_wilcoxon_signed_rank_matched_pairs", "t_wilcoxon_mann_whitney")) for (parent in c("normal", "laplace", "logistic")) local({
  id <- method; kind <- if (id == "t_wilcoxon_mann_whitney") "rank_sum" else "signed_rank"
  moments <- rank$ref_moments(kind, parent, .7)
  effect <- if (id == "t_wilcoxon_signed_rank_matched_pairs") "dz" else "d"
  p <- if (kind == "rank_sum") list(n1 = 18, n2 = 18) else list(n_total = 28); p[[effect]] <- .7
  profile(paste0("rank-", parent), id, p, effect,
          function(p) rank$ref_lehmann(kind, if (kind == "rank_sum") c(p$n1, p$n2) else p$n_total, moments, p$alpha, "two", TRUE),
          list(method = "lehmann", parent_distribution = parent, tails = "two", continuity_correction = "on"),
          modes = "a_priori", minimum = 2, step = if (kind == "rank_sum") 2 else 1)
})
profile("independent-correlation", "z_correlations_two_independent", list(q = .6, n1 = 20, n2 = 20), "q",
        function(p) glm$ck_zpower(p$q / sqrt(1 / (p$n1 - 3) + 1 / (p$n2 - 3)), 1, p$alpha, "two")$power,
        list(tails = "two"), modes = "a_priori", minimum = 8, step = 2)
for (common in c(TRUE, FALSE)) local({
  shared <- common
  id <- if (shared) "z_correlations_two_dependent_common_index" else "z_correlations_two_dependent_no_common_index"
  p <- if (shared) list(rho_ab = .35, rho_ac = -.1, rho_bc = .25, n_total = 80) else
    list(rho_ab = .35, rho_cd = -.1, rho_ac = .1, rho_ad = .1, rho_bc = .1, rho_bd = .1, n_total = 80)
  power <- function(p) {
    effect <- if (shared) p$rho_ac else p$rho_cd
    R <- if (shared) matrix(c(1,p$rho_ab,effect,p$rho_ab,1,p$rho_bc,effect,p$rho_bc,1),3) else
      matrix(c(1,p$rho_ab,p$rho_ac,p$rho_ad,p$rho_ab,1,p$rho_bc,p$rho_bd,p$rho_ac,p$rho_bc,1,effect,p$rho_ad,p$rho_bd,effect,1),4)
    H0 <- R; pair <- if (shared) c(1,3) else c(3,4)
    H0[pair[1],pair[2]] <- H0[pair[2],pair[1]] <- p$rho_ab
    c0 <- correlations$g7_cov_z(H0,c(1,2),pair); c1 <- correlations$g7_cov_z(R,c(1,2),pair)
    s0 <- sqrt((2-2*c0)/(p$n_total-3))
    correlations$g7_zpower((atanh(effect)-atanh(p$rho_ab))/s0,sqrt((1-c1)/(1-c0)),p$alpha,"two",-1)
  }
  profile("dependent-correlation", id, p, if (shared) "rho_ac" else "rho_cd", power, list(tails = "two"), minimum = 4)
})
for (branch in c("exact_brown_benedetti", "approx_bonett_price")) local({
  computation <- branch
  power <- function(p) {
    c0 <- correlations$g7_cells(p$rho_h0, p$marginal_x, p$marginal_y)
    c1 <- correlations$g7_cells(p$rho_h1, p$marginal_x, p$marginal_y)
    if (computation == "exact_brown_benedetti") {
      s0 <- correlations$g7_se_pearson(c0,p$rho_h0,p$n_total); s1 <- correlations$g7_se_pearson(c1,p$rho_h1,p$n_total)
      shift <- p$rho_h1-p$rho_h0
    } else {
      b0 <- correlations$ref_bp(c0,p$marginal_x,p$marginal_y,p$n_total); b1 <- correlations$ref_bp(c1,p$marginal_x,p$marginal_y,p$n_total)
      s0 <- b0$se; s1 <- b1$se; shift <- b1$r-b0$r
    }
    correlations$g7_zpower(shift/s0,s1/s0,p$alpha,"two",1)
  }
  profile(paste0("tetrachoric-", branch), "z_tetrachoric_correlation",
          list(rho_h0 = .1, rho_h1 = .55, marginal_x = .4, marginal_y = .6, n_total = 80), "rho_h1", power,
          list(tails = "two", computation_of_correlation = computation), minimum = 2)
})
for (arcsine in c(FALSE, TRUE)) local({
  transform <- arcsine
  profile(paste0("normal-proportion-", transform), "z_proportions_two_independent", list(p1 = .2, p2 = .65, n1 = 20, n2 = 20), "p2",
          function(p) proportions$ref_z2p_power(p$p1,p$p2,p$n1,p$n2,p$alpha,"two",FALSE,transform),
          list(tails = "two", continuity_correction = "off", arcsin_transform = if (transform) "on" else "off"),
          modes = c("a_priori", "sensitivity", "compromise"), minimum = 2, step = 2)
})

# Discrete criterion/compromise follow the declared smallest attaining
# nominal-alpha contract. Reference enumeration checks both sides of a jump.
for (method in c("exact_proportion_one_sample_binomial", "exact_proportion_sign_test", "exact_generic_binomial")) local({
  id <- method; null <- if (id == "exact_proportion_sign_test") .5 else .35
  p <- if (id == "exact_generic_binomial") list(p1 = null, p2 = null + .3) else c(list(g = .3), if (id != "exact_proportion_sign_test") list(p0 = null))
  p$n_total <- 18
  profile("binomial", id, p, if (id == "exact_generic_binomial") "p2" else "g",
          function(p) discrete$ref_binom_test(p$n_total,null,if (id == "exact_generic_binomial") p$p2 else null+p$g,p$alpha,"two","half_each_side"),
          list(tails = "two", alpha_balancing = "half_each_side"), minimum = 1, discrete_model = TRUE)
})
for (computation in c("exact_unconditional", "fast_approx")) local({
  branch <- computation
  power <- function(p) {
    fun <- if (branch == "exact_unconditional") discrete$ref_mcnemar_exact else discrete$ref_mcnemar_fast
    fun(p$n_total,p$or,p$pd,p$alpha,"two","half_each_side")
  }
  profile(paste0("mcnemar-", branch), "exact_proportions_mcnemar", list(or = 4, pd = .5, n_total = 30), "or", power,
          list(tails = "two", alpha_balancing = "half_each_side", computation = branch), minimum = 1, discrete_model = TRUE)
})
profile("fisher", "exact_proportions_fisher", list(p1 = .8, p2 = .2, n1 = 8, n2 = 8), "p1",
        function(p) discrete$ref_fisher_exact(p$n1,p$n2,p$p1,p$p2,p$alpha,"two"),
        list(tails = "two", computation = "exact"), minimum = 2, step = 2, discrete_model = TRUE, sensitivity_options = list(effect_direction = "p1_ge_p2"))
for (variant in c("central", "minlike", "less", "greater")) {
  alpha <- if (variant %in% c("central", "minlike")) .13 else .07
  p1 <- if (variant == "greater") .9 else .1; p2 <- 1-p1
  run(paste0("boschloo-smallest-", variant), "boschloo_proportions_two_independent", "a_priori",
      list(p1 = p1, p2 = p2, alpha = alpha, power = .6, allocation_ratio = 1, n_max = 4), list(variant = variant, allocation = "n2_ceiling"), function(r) {
        expected <- if (variant %in% c("central", "minlike")) .9^4 + .1^4 else .9^4
        record("power", r$achieved$power, expected, 1e-12, 0, "hand-enumerated extreme-table probabilities at n1=n2=2")
        record("first_admissible_sample_size", r$solved$n_total, 4, 0, 0, "n1=n2=1 has no rejection at alpha<1/4; n1=n2=2 reaches target")
      })
}
# Nonzero-offset exact criterion: check the independently enumerated jump and
# keep nominal and attained alpha separate.
run("offset-exact-criterion", "exact_proportions_unconditional_offset", "criterion",
    list(p1_h0 = .4, p1_h1 = .8, p2 = .3, n1 = 8, n2 = 10, power = .6),
    list(effect_size_index = "proportion_p1", computation = "exact", statistic = "z_pooled", tails = "two"), function(r) {
      a <- r$solved$alpha
      now <- unconditional$g9_offset_power("z_pooled",8,10,.4,.8,.3,a,"two")
      before <- unconditional$g9_offset_power("z_pooled",8,10,.4,.8,.3,a*(1-1e-7),"two")
      record("power",r$achieved$power,now[["power"]],1e-11,0)
      record("alpha_attained",r$achieved$alpha_attained,now[["alpha"]],1e-11,0)
      condition("target_reached",now[["power"]]>=.6)
      condition("left_boundary_below_target",before[["power"]]<.6)
    })

cases <- do.call(rbind, case_rows); results <- do.call(rbind, rows)
source(file.path(here, "inventory.R"))
inventory <- write_numerical_inventory(pkg, cases, out)
summary <- list(scope = "bounded synthetic inverse and integer-minimality campaign; not complete domain certification",
  package_version = as.character(utils::packageVersion("inferguidR")), package_path = pkg,
  r_version = R.version.string, platform = R.version$platform,
  source_sha256_caller_supplied = Sys.getenv("INFERGUIDR_SOURCE_SHA256"),
  started_utc = format(started, "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"), finished_utc = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
  cases = nrow(cases), comparisons = nrow(results), planned_cases = 144L, planned_coverage_complete = nrow(cases) == 144L,
  failed_comparisons = sum(!results$passed), negative_control = negative, negative_control_injected = injected,
  checked_methods = length(unique(cases$method[cases$passed & cases$evidence_kind == "independent_numeric"])), inventory = inventory,
  passed = all(results$passed) && nrow(cases) == 144L && (!negative || injected), session_info = capture.output(sessionInfo()))
if (requireNamespace("digest", quietly = TRUE)) {
  files <- c(file.path(here, "check-inverses.R"), file.path(here, "inventory.R"), list.files(file.path(here, "reference"), full.names = TRUE))
  summary$checker_sha256 <- stats::setNames(lapply(files, function(f) digest::digest(file = f, algo = "sha256")), substring(files, nchar(here) + 2L))
}
jsonlite::write_json(summary, file.path(out, "summary.json"), auto_unbox = TRUE, pretty = TRUE, null = "null")
cat(nrow(cases), "cases", nrow(results), "comparisons", sum(!results$passed), "failed comparisons\n")
quit(status = if (summary$passed) 0L else 1L)
