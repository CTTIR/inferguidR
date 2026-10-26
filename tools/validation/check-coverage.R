#!/usr/bin/env Rscript
# Copyright (c) 2026 Raban Heller. MIT licence; see ../../LICENSE.md.
# Extended synthetic campaign. This is separate from the bounded CI smoke run.
args <- commandArgs(TRUE)
if (length(args) < 1L || length(args) > 2L || !nzchar(args[1]) ||
    (length(args) == 2L && args[2] != "--negative-control")) stop("usage: Rscript check-coverage.R NEW_OUTPUT_DIRECTORY [--negative-control]")
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
glm <- ref("glm.R"); rank <- ref("rank.R"); manova <- ref("manova.R")
discrete <- ref("discrete.R"); offset <- ref("offset.R"); regression <- ref("regression.R")
variance <- ref("variance.R"); correlations <- ref("correlations.R")
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
  before <- length(rows); error <- ""; status <- "not_returned"; warnings <- NULL; t0 <- proc.time()[["elapsed"]]
  tryCatch({
    r <- if (is.null(operation)) inferguidR::calculate_power(inferguidR::analysis_spec(method, mode, inputs, options, target = target)) else operation()
    status <- r$status$code
    warnings <- r$warnings
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
    elapsed_seconds = proc.time()[["elapsed"]] - t0, warnings_json = json(warnings), error = error)
  case_rows[[length(case_rows) + 1L]] <<- row
  append_csv(row, file.path(out, "cases.csv"))
  if (length(case_rows) %% 10L == 0L) cat(length(case_rows), "cases completed\n")
}
t_power <- function(df, ncp, alpha, tails = "two") {
  c <- cont$ref_qt_upper(if (tails == "two") alpha / 2 else alpha, df)
  if (tails == "two") cont$ref_pt(c, df, ncp, FALSE) + cont$ref_pt(-c, df, ncp) else cont$ref_pt(c, df, abs(ncp), FALSE)
}
check_t <- function(r, df, ncp, alpha, tails = "two") {
  record("df", unname(r$quantities$df), df, 1e-12, 0)
  record("ncp", r$quantities$ncp, ncp, 1e-10, 1e-10, "design-derived noncentrality")
  record("power", r$achieved$power, t_power(df, ncp, alpha, tails), oracle = "normal/chi-square density integral")
}
check_f <- function(r, df1, df2, lambda, alpha) {
  record("df1", unname(r$quantities$df[1]), df1, 1e-10, 1e-12)
  record("df2", unname(r$quantities$df[2]), df2, 1e-10, 1e-12)
  record("ncp", r$quantities$ncp, lambda, 1e-10, 1e-10, "independently derived parameter mapping")
  critical <- r$quantities$critical
  record("critical_null_tail", multi$ref_pf(critical, df1, df2, 0, FALSE), alpha, oracle = "F convolution through central chi-square denominator")
  record("power", r$achieved$power, multi$ref_pf(critical, df1, df2, lambda, FALSE))
}

# Continuous t/F models. The oracle derives df and ncp from inputs, not from
# reported engine parameters; inverse outputs enter only as the solved input.
models <- list()
add <- function(id, input, effect, family, parameters, options = list()) {
  models[[id]] <<- list(input = input, effect = effect, family = family, parameters = parameters, options = options)
}
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
for (id in names(models)) for (mode in c("post_hoc", "sensitivity", "criterion", "compromise")) {
  model <- models[[id]]; p <- c(model$input, list(alpha = .05)); options <- model$options
  if (model$family == "t") options$tails <- "two"
  if (mode == "sensitivity") { p[[model$effect]] <- NULL; p$power <- .7 }
  if (mode == "criterion") { p$alpha <- NULL; p$power <- .7 }
  if (mode == "compromise") { p$alpha <- NULL; p$error_ratio <- 1 }
  run(paste("continuous", id, mode, sep = "-"), id, mode, p, options, function(r) {
    resolved <- p
    if (mode == "sensitivity") {
      resolved[[model$effect]] <- r$solved[[model$effect]]
      if (is.null(resolved[[model$effect]])) resolved[[model$effect]] <- r$method_outputs[[model$effect]]
    }
    if (mode %in% c("criterion", "compromise")) resolved$alpha <- r$solved$alpha
    pars <- model$parameters(resolved)
    if (model$family == "t") check_t(r, pars[1], pars[2], resolved$alpha) else check_f(r, pars[1], pars[2], pars[3], resolved$alpha)
    if (mode %in% c("sensitivity", "criterion")) record("target", r$achieved$power, .7, 1e-8, 1e-6)
    if (mode == "compromise") record("nominal_error_ratio", r$achieved$beta, resolved$alpha, 1e-8, 1e-6)
  })
}

# RM and MANOVA parameter mappings, all approximation and algorithm choices.
for (which in c("between", "within", "within_between")) {
  id <- paste0("f_rm_anova_", which)
  p <- list(f = .25, n_total = 60, k = 3, m = 4, rho = .3, alpha = .05)
  if (which != "between") p$epsilon <- .8
  for (specification in c("gpower30", "cohen1988", "spss")) {
    input <- p
    if (specification != "gpower30") input$rho <- NULL
    run(paste("rm-anova", which, specification, sep = "-"), id, "post_hoc", input, list(effect_size_specification = specification), function(r) {
      df1 <- if (which == "between") 2 else (if (which == "within") 3 else 6) * .8
      df2 <- if (which == "between") 57 else 57 * 3 * .8
      lambda <- if (specification == "spss") .25^2 * df2 else if (specification == "cohen1988") .25^2 * 60 * (if (which == "between") 1 else .8) else
        .25^2 * 60 * 4 * (if (which == "between") 1 / 1.9 else .8 / .7)
      check_f(r, df1, df2, lambda, .05)
    })
  }
}
for (id in c("f_manova_global", "f_manova_special", "f_rm_manova_between", "f_rm_manova_within", "f_rm_manova_within_between")) {
  repeated <- grepl("f_rm_", id)
  input <- if (repeated) list(f = .25, n_total = 60, k = 3, m = 4, alpha = .05) else
    list(f2 = .08, n_total = 60, n_groups = 3, n_responses = 3, alpha = .05)
  if (id == "f_manova_special") input$n_predictors <- 2
  for (approx in c("pillai_v", "wilks_u", "hotelling_lawley_t1", "hotelling_lawley_t2")) for (algorithm in c("obrien_shieh_1999", "muller_peterson_1984")) {
    options <- list(approximation = approx, algorithm = algorithm)
    if (repeated) options$use_mean_correlation <- "off"
    run(paste("manova", id, approx, algorithm, sep = "-"), id, "post_hoc", input, options, function(r) {
      p <- if (id == "f_rm_manova_between") 1 else 3
      q <- if (id == "f_rm_manova_within") 1 else 2
      f2 <- if (repeated) .25^2 else .08
      terms <- manova$ref_terms(approx, 57, p, q, f2)
      lambda <- manova$ref_lambda(algorithm, 60, terms, f2)
      check_f(r, terms$df1, terms$df2, lambda, .05)
    })
  }
}

# Distribution-based variance models, evaluated by density quadrature.
for (id in c("f_variances_two_groups", "chi2_variance_one_sample")) for (ratio in c(.55, 1, 1.7)) for (tails in c("one", "two")) {
  p <- c(list(ratio = ratio, alpha = .05), if (id == "f_variances_two_groups") list(n1 = 18, n2 = 27) else list(n_total = 23))
  run(paste("variance", id, ratio, tails, sep = "-"), id, "post_hoc", p, list(tails = tails), function(r) {
    ans <- if (id == "f_variances_two_groups") variance$g11r_power_f2(ratio, 18, 27, .05, tails) else variance$g11r_power_var1(ratio, 23, .05, tails)
    record("power", r$achieved$power, ans$power, oracle = "integrated log density and independently inverted quantile")
    record("beta", r$achieved$beta, ans$beta)
  })
}
for (id in c("chi2_goodness_of_fit_contingency", "chi2_generic")) for (mode in c("post_hoc", "sensitivity", "criterion", "compromise")) {
  p <- if (id == "chi2_generic") list(lambda = 7, df = 4, alpha = .05) else list(w = .35, n_total = 60, df = 4, alpha = .05)
  if (mode == "sensitivity") { p$lambda <- p$w <- NULL; p$power <- .7 }
  if (mode == "criterion") { p$alpha <- NULL; p$power <- .7 }
  if (mode == "compromise") { p$alpha <- NULL; p$error_ratio <- 1 }
  run(paste("chi-square", id, mode, sep = "-"), id, mode, p, list(), function(r) {
    alpha <- if (mode %in% c("criterion", "compromise")) r$solved$alpha else .05
    lambda <- if (id == "chi2_generic") { if (mode == "sensitivity") r$solved$lambda else 7 } else
      (if (mode == "sensitivity") r$solved$w else .35)^2 * 60
    record("critical_null_tail", cont$ref_pchisq(r$quantities$critical, 4, 0, FALSE), alpha)
    pw <- cont$ref_pchisq(r$quantities$critical, 4, lambda, FALSE)
    record("power", r$achieved$power, pw)
    if (mode %in% c("criterion", "sensitivity")) record("reference_target", pw, .7)
    if (mode == "compromise") record("reference_error_ratio", 1 - pw, alpha)
  })
}

# Exact discrete models: every table/count is included, not only rejected ones.
for (id in c("exact_proportion_one_sample_binomial", "exact_proportion_sign_test", "exact_generic_binomial")) for (tails in c("one", "two")) for (balancing in if (tails == "one") "half_each_side" else c("half_each_side", "minor_tail_half_rest_major", "half_each_side_then_increase")) {
  p0 <- if (id == "exact_proportion_sign_test") .5 else .35; p1 <- p0 + .25
  p <- if (id == "exact_generic_binomial") list(p1 = p0, p2 = p1) else c(list(g = .25), if (id != "exact_proportion_sign_test") list(p0 = p0))
  p <- c(p, list(n_total = 17, alpha = .05))
  opt <- list(tails = tails)
  if (tails == "two") opt$alpha_balancing <- balancing
  run(paste("binomial", id, tails, balancing, sep = "-"), id, "post_hoc", p, opt, function(r) {
    ans <- discrete$ref_binom_test(17, p0, p1, .05, tails, balancing)
    record("power", r$achieved$power, ans$power, 1e-12, 0, "enumerated binomial probabilities and independent tail region")
    record("alpha_attained", r$achieved$alpha_attained, ans$alpha, 1e-12, 0)
  })
}
for (computation in c("exact_unconditional", "fast_approx")) for (tails in c("one", "two")) {
  opt <- list(computation = computation, tails = tails)
  if (tails == "two") opt$alpha_balancing <- "half_each_side"
  run(paste("mcnemar", computation, tails, sep = "-"), "exact_proportions_mcnemar", "post_hoc",
    list(or = 2.5, pd = .4, n_total = 25, alpha = .05), opt, function(r) {
      ans <- if (computation == "exact_unconditional") discrete$ref_mcnemar_exact(25, 2.5, .4, .05, tails, "half_each_side") else discrete$ref_mcnemar_fast(25, 2.5, .4, .05, tails, "half_each_side")
      record("power", r$achieved$power, ans$power, 1e-12, 0, "binomial discordant-count mixture / conditional enumeration")
    })
}
for (tails in c("one", "two")) {
  run(paste0("fisher-", tails), "exact_proportions_fisher", "post_hoc", list(p1 = .25, p2 = .7, n1 = 6, n2 = 9, alpha = .05),
    list(tails = tails, computation = "exact"), function(r) {
      ans <- discrete$ref_fisher_exact(6, 9, .25, .7, .05, tails)
      record("power", r$achieved$power, ans$power, 1e-12, 0, "full-table hypergeometric ordering and binomial probability sums")
    })
}
for (variant in c("central", "minlike", "less", "greater")) {
  alpha <- if (variant %in% c("central", "minlike")) .13 else .07
  run(paste0("boschloo-elementary-", variant), "boschloo_proportions_two_independent", "post_hoc",
    list(p1 = .25, p2 = .65, n1 = 2, n2 = 2, alpha = alpha), list(variant = variant), function(r) {
      low <- (1 - .25)^2 * .65^2; up <- .25^2 * (1 - .65)^2
      pw <- if (variant == "less") low else if (variant == "greater") up else low + up
      record("power", r$achieved$power, pw, 1e-12, 0, "elementary extreme-table polynomial; null supremum at p=1/2")
    })
}

# Exact correlation and random-design R-squared: hypergeometric density
# integrals, including exact and three-moment approximation options.
for (id in c("exact_correlation_bivariate_normal", "exact_mrc_random_model")) for (tails in c("one", "two")) for (branch in c("exact", "approx")) for (mode in c("post_hoc", "sensitivity", "criterion", "compromise")) {
  corr <- id == "exact_correlation_bivariate_normal"
  p <- if (corr) list(rho_h0 = .1, rho_h1 = .5, n_total = 18, alpha = .05) else list(rho2_h0 = .05, rho2_h1 = .3, n_total = 24, n_predictors = 2, alpha = .05)
  opt <- list(tails = tails)
  if (corr) opt$computation <- if (branch == "exact") "exact" else "fisher_z_approx" else opt$distribution <- if (branch == "exact") "exact" else "three_moment_f"
  if (mode == "sensitivity") { p$rho_h1 <- p$rho2_h1 <- NULL; p$power <- .7 }
  if (mode == "criterion") { p$alpha <- NULL; p$power <- .7 }
  if (mode == "compromise") { p$alpha <- NULL; p$error_ratio <- 1 }
  run(paste("correlation-distribution", id, tails, branch, mode, sep = "-"), id, mode, p, opt, function(r) {
    alpha <- if (mode %in% c("criterion", "compromise")) r$solved$alpha else .05
    effect <- if (mode == "sensitivity") r$method_outputs[[if (corr) "rho_h1" else "rho2_h1"]] else if (corr) .5 else .3
    if (corr && branch == "approx") {
      ans <- glm$ck_zpower((atanh(effect) - atanh(.1)) * sqrt(15), 1, alpha, tails)$power
    } else {
      cdf <- if (corr) function(x, h, lower) regression$g5r_pr(x, h, 18, lower) else if (branch == "exact")
        function(x, h, lower) regression$g5r_pr2(x, h, 2, 24, lower) else function(x, h, lower) regression$g5r_pr2_lee(x, h, 2, 24, lower)
      critical <- r$quantities$critical
      null <- if (corr) .1 else .05
      if (tails == "two") {
        record("lower_null_tail", cdf(min(critical), null, TRUE), alpha / 2)
        record("upper_null_tail", cdf(max(critical), null, FALSE), alpha / 2)
        ans <- cdf(min(critical), effect, TRUE) + cdf(max(critical), effect, FALSE)
      } else {
        record("critical_null_tail", cdf(critical, null, FALSE), alpha)
        ans <- cdf(critical, effect, FALSE)
      }
    }
    record("power", r$achieved$power, ans, oracle = "independent correlation/R-squared density integral or moment-matched F integral")
    if (mode %in% c("sensitivity", "criterion")) record("reference_target", ans, .7)
    if (mode == "compromise") record("reference_error_ratio", 1 - ans, alpha)
  })
}

# Score statistics for nonzero risk/odds-ratio offsets. Restricted nuisance
# estimates come from independent score-equation bisection, not cubic formulas.
for (scale in c("risk_ratio", "odds_ratio")) {
  statistics <- if (scale == "risk_ratio") c("score_gart_nam", "score_farrington_manning", "score_miettinen_nurminen") else c("score_farrington_manning", "score_miettinen_nurminen")
  for (statistic in statistics) for (tails in c("one", "two")) for (design in list(c(3, 5), c(1, 4))) {
    p <- list(p2 = .3, n1 = design[1], n2 = design[2], alpha = .05)
    p[[paste0(scale, "_h0")]] <- 1.3; p[[paste0(scale, "_h1")]] <- 2.1
    convert <- function(value) if (scale == "risk_ratio") value * .3 else value * .3 / (.7 + value * .3)
    null <- list(scale = if (scale == "risk_ratio") "ratio" else "odds_ratio", value = 1.3)
    run(paste("offset-score", scale, statistic, tails, paste(design, collapse = "x"), sep = "-"), "exact_proportions_unconditional_offset", "post_hoc", p,
      list(effect_size_index = scale, statistic = statistic, tails = tails, computation = "exact"), function(r) {
        ans <- offset$ref_exact(statistic, design[1], design[2], convert(2.1), convert(1.3), .3, null, .05, tails, 1)
        record("power", r$achieved$power, ans[["power"]], 1e-11, 0, "restricted-likelihood score bisection and full table enumeration")
        record("alpha_attained", r$achieved$alpha_attained, ans[["alpha"]], 1e-11, 0)
      })
  }
}

# Rank parent/moment inverse options. Parent moments are calculated using
# separate quadrature or closed forms, then the published Lehmann moments.
for (scale in c("risk_ratio", "odds_ratio")) {
  statistics <- if (scale == "risk_ratio") c("score_gart_nam", "score_farrington_manning", "score_miettinen_nurminen") else c("score_farrington_manning", "score_miettinen_nurminen")
  for (statistic in statistics) for (tails in c("one", "two")) for (mode in c("sensitivity", "a_priori")) {
    p <- list(p2 = .2, alpha = .05, power = .6)
    p[[paste0(scale, "_h0")]] <- 1.3
    if (mode == "sensitivity") { p$n1 <- 24; p$n2 <- 28 } else
      p[[paste0(scale, "_h1")]] <- if (scale == "risk_ratio") 4 else 16
    p10 <- if (scale == "risk_ratio") .26 else .26 / 1.06
    null <- list(scale = if (scale == "risk_ratio") "ratio" else "odds_ratio", value = 1.3)
    opt <- list(effect_size_index = scale, statistic = statistic, tails = tails, computation = "exact")
    run(paste("offset-inverse", scale, statistic, tails, mode, sep = "-"), "exact_proportions_unconditional_offset", mode, p, opt, function(r) {
      if (mode == "sensitivity") {
        effect <- r$solved$p1_h1
        ans <- offset$ref_exact(statistic, 24, 28, effect, p10, .2, null, .05, tails, 1)
        record("reference_target", ans[["power"]], .6, 1e-8, 1e-6, "independent score-equation bisection and full-table probability")
        record("power", r$achieved$power, ans[["power"]], 1e-11, 0)
      } else {
        n <- unname(r$group_sizes[1])
        condition("bounded_design", n <= 80)
        if (n <= 80) {
          powers <- vapply(seq_len(n), function(k) offset$ref_exact(statistic, k, k, .8, p10, .2, null, .05, tails, 1)[["power"]], 0.0)
          record("first_total_size", r$solved$n_total, 2 * which(powers >= .6)[1], 0, 0, "exhaustive independent scan of every earlier equal allocation")
          record("power", r$achieved$power, powers[n], 1e-11, 0)
        }
      }
    })
  }
}

# Rank parent/moment inverse options. Parent moments are calculated using
# separate quadrature or closed forms, then the published Lehmann moments.
for (id in c("t_wilcoxon_signed_rank_one_sample", "t_wilcoxon_signed_rank_matched_pairs", "t_wilcoxon_mann_whitney")) for (parent in c("normal", "laplace", "logistic")) for (mode in c("post_hoc", "sensitivity", "criterion", "compromise")) {
  kind <- if (id == "t_wilcoxon_mann_whitney") "rank_sum" else "signed_rank"
  effect_name <- if (id == "t_wilcoxon_signed_rank_matched_pairs") "dz" else "d"
  p <- c(if (kind == "rank_sum") list(n1 = 18, n2 = 25) else list(n_total = 28), list(alpha = .05))
  p[[effect_name]] <- .55
  if (mode == "sensitivity") { p[[effect_name]] <- NULL; p$power <- .7 }
  if (mode == "criterion") { p$alpha <- NULL; p$power <- .7 }
  if (mode == "compromise") { p$alpha <- NULL; p$error_ratio <- 1 }
  options <- list(method = "lehmann", parent_distribution = parent, tails = "two", continuity_correction = "on")
  run(paste("rank-parent", id, parent, mode, sep = "-"), id, mode, p, options, function(r) {
    d <- if (mode == "sensitivity") r$solved[[effect_name]] else .55
    alpha <- if (mode %in% c("criterion", "compromise")) r$solved$alpha else .05
    moments <- rank$ref_moments(kind, parent, abs(d))
    ans <- rank$ref_lehmann(kind, if (kind == "rank_sum") c(18, 25) else 28, moments, alpha, "two", TRUE)
    record("power", r$achieved$power, ans$power, 1e-9, 1e-7, "independent parent quadrature and Lehmann moment formula")
    if (mode %in% c("sensitivity", "criterion")) record("reference_target", ans$power, .7)
    if (mode == "compromise") record("reference_error_ratio", ans$beta, alpha)
  })
}

for (id in c("t_wilcoxon_signed_rank_one_sample", "t_wilcoxon_signed_rank_matched_pairs", "t_wilcoxon_mann_whitney")) for (mode in c("sensitivity", "criterion", "compromise")) {
  kind <- if (id == "t_wilcoxon_mann_whitney") "rank_sum" else "signed_rank"
  moments <- if (kind == "rank_sum") c(.68, .5, .47) else c(.7, .7, .55)
  n <- if (kind == "rank_sum") c(12, 17) else 12
  target_power <- rank$ref_lehmann(kind, n, moments, .05, "two", FALSE)$power
  p <- c(as.list(setNames(moments, c("p1", "p2", "p3"))), if (kind == "rank_sum") list(n1 = n[1], n2 = n[2]) else list(n_total = n), list(alpha = .05))
  if (mode == "sensitivity") { p$p1 <- NULL; p$power <- target_power }
  if (mode == "criterion") { p$alpha <- NULL; p$power <- target_power }
  if (mode == "compromise") { p$alpha <- NULL; p$error_ratio <- 1 }
  options <- list(method = "lehmann", lehmann_input_method = "manual_moments", tails = "two", continuity_correction = "off")
  run(paste("rank-manual", id, mode, sep = "-"), id, mode, p, options, function(r) {
    m <- moments
    if (mode == "sensitivity") m[1] <- r$solved$p1
    alpha <- if (mode %in% c("criterion", "compromise")) r$solved$alpha else .05
    ans <- rank$ref_lehmann(kind, n, m, alpha, "two", FALSE)
    record("power", r$achieved$power, ans$power, 1e-10, 1e-8, "manual-moment Lehmann model evaluated independently")
    if (mode %in% c("sensitivity", "criterion")) record("reference_target", ans$power, target_power, 1e-8, 1e-6)
    if (mode == "compromise") record("reference_error_ratio", ans$beta, alpha, 1e-8, 1e-6)
  }, target = if (mode == "sensitivity") "p1" else NULL,
  expected_status = if (mode == "sensitivity") "invalid_input" else "ok")
}

# Normal correlation tests and generic normal parameters.
for (id in c("z_correlations_two_independent", "z_generic")) for (mode in c("post_hoc", "sensitivity", "criterion", "compromise")) {
  independent <- id == "z_correlations_two_independent"
  p <- if (independent) list(q = .4, n1 = 25, n2 = 35, alpha = .05) else list(mu = 2, sigma = 1.2, alpha = .05)
  effect_name <- if (independent) "q" else "mu"
  if (mode == "sensitivity") { p[[effect_name]] <- NULL; p$power <- .7 }
  if (mode == "criterion") { p$alpha <- NULL; p$power <- .7 }
  if (mode == "compromise") { p$alpha <- NULL; p$error_ratio <- 1 }
  run(paste("normal", id, mode, sep = "-"), id, mode, p, list(tails = "two"), function(r) {
    effect <- if (mode == "sensitivity") r$solved[[effect_name]] else if (independent) .4 else 2
    alpha <- if (mode %in% c("criterion", "compromise")) r$solved$alpha else .05
    ans <- glm$ck_zpower(if (independent) effect / sqrt(1 / 22 + 1 / 32) else effect, if (independent) 1 else 1.2, alpha, "two")
    record("power", r$achieved$power, ans$power, 1e-10, 1e-8, "Fisher-z variance and own normal quadrature")
    if (mode %in% c("sensitivity", "criterion")) record("reference_target", ans$power, .7)
    if (mode == "compromise") record("reference_error_ratio", ans$beta, alpha)
  })
}
for (common in c(TRUE, FALSE)) {
  id <- if (common) "z_correlations_two_dependent_common_index" else "z_correlations_two_dependent_no_common_index"
  p <- if (common) list(rho_ab = .35, rho_ac = .15, rho_bc = .25) else
    list(rho_ab = .35, rho_cd = .15, rho_ac = .1, rho_ad = .1, rho_bc = .1, rho_bd = .1)
  for (tails in c("one", "two")) run(paste("dependent", common, tails, sep = "-"), id, "post_hoc", c(p, list(n_total = 80, alpha = .05)), list(tails = tails), function(r) {
    R <- if (common) matrix(c(1,.35,.15,.35,1,.25,.15,.25,1),3) else
      matrix(c(1,.35,.1,.1,.35,1,.1,.1,.1,.1,1,.15,.1,.1,.15,1),4)
    H0 <- R
    if (common) H0[1,3] <- H0[3,1] <- .35 else H0[3,4] <- H0[4,3] <- .35
    pair <- if (common) c(1,3) else c(3,4)
    c0 <- correlations$g7_cov_z(H0,c(1,2),pair); c1 <- correlations$g7_cov_z(R,c(1,2),pair)
    s0 <- sqrt((2-2*c0)/77)
    pw <- correlations$g7_zpower((atanh(.15)-atanh(.35))/s0,sqrt((1-c1)/(1-c0)),.05,tails,-1)
    record("power", r$achieved$power, pw, 1e-10, 1e-8, "Isserlis/delta-method covariance of sample correlations")
  })
}

# Selected effect utilities: raw-input formulas, independent of the planning
# engines and of the drawer's reported derived quantities.
utilities <- list(
  list("t_means_two_independent", "common_sd", list(mean_1 = 4, mean_2 = 1, sd = 2), 1.5),
  list("t_means_one_sample", "from_means", list(mean_h1 = 4, mean_h0 = 1, sd = 2), 1.5),
  list("f_variances_two_groups", "from_variances", list(var0 = 2, var1 = 3), 1.5),
  list("chi2_variance_one_sample", "from_variances", list(var0 = 3, var1 = 2), 2/3),
  list("t_correlation_point_biserial", "r_squared", list(r_squared = .16), .4),
  list("f_mrc_r2_deviation_from_zero", "r_squared", list(rho2 = .2), .25),
  list("chi2_goodness_of_fit_contingency", "cell_probabilities", list(p0 = c(.2,.3,.5), p1 = c(.3,.4,.3)), sqrt(sum((c(.3,.4,.3)-c(.2,.3,.5))^2/c(.2,.3,.5)))))
for (u in utilities) run(paste("effect", u[[1]], u[[2]], sep = "-"), u[[1]], "effect_utility", u[[3]], list(input_mode = u[[2]]),
  function(r) record("value", r$value, u[[4]], 1e-12, 0, "raw-input algebraic effect definition"),
  operation = function() do.call(inferguidR::effect_size, c(list(method = u[[1]], input_mode = u[[2]]), u[[3]])))

cases <- do.call(rbind, case_rows); results <- do.call(rbind, rows)
source(file.path(here, "inventory.R"))
inventory <- write_numerical_inventory(pkg, cases, out)
summary <- list(scope = "extended synthetic mathematical coverage; no complete domain certification",
  package_version = as.character(utils::packageVersion("inferguidR")), package_path = pkg,
  r_version = R.version.string, platform = R.version$platform,
  source_sha256_caller_supplied = Sys.getenv("INFERGUIDR_SOURCE_SHA256"),
  started_utc = format(started, "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"), finished_utc = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
  cases = nrow(cases), comparisons = nrow(results), planned_cases = 271L, planned_coverage_complete = nrow(cases) == 271L,
  failed_comparisons = sum(!results$passed), negative_control = negative, negative_control_injected = injected,
  checked_methods = length(unique(cases$method[cases$passed & cases$evidence_kind == "independent_numeric"])), inventory = inventory,
  passed = all(results$passed) && nrow(cases) == 271L && (!negative || injected), session_info = capture.output(sessionInfo()))
if (requireNamespace("digest", quietly = TRUE)) {
  files <- c(file.path(here, "check-coverage.R"), file.path(here, "inventory.R"), list.files(file.path(here, "reference"), full.names = TRUE))
  summary$checker_sha256 <- stats::setNames(lapply(files, function(f) digest::digest(file = f, algo = "sha256")), substring(files, nchar(here) + 2L))
}
jsonlite::write_json(summary, file.path(out, "summary.json"), auto_unbox = TRUE, pretty = TRUE, null = "null")
cat(nrow(cases), "cases", nrow(results), "comparisons", sum(!results$passed), "failed comparisons\n")
quit(status = if (summary$passed) 0L else 1L)
