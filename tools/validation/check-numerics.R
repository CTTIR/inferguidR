#!/usr/bin/env Rscript
# Copyright (c) 2026 Raban Heller. MIT licence; see ../../LICENSE.md.

# Bounded, synthetic mathematical checks against an installed package.
# Each case is written immediately; a failing case never suppresses later cases.
args <- commandArgs(TRUE)
if (length(args) < 1L || length(args) > 2L ||
    (length(args) == 2L && args[2L] != "--negative-control")) {
  stop("usage: Rscript check-numerics.R NEW_OUTPUT_DIRECTORY [--negative-control]")
}
out <- path.expand(args[1L])
if (dir.exists(out) && length(list.files(out, all.files = TRUE, no.. = TRUE))) stop("output directory must be empty")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
out <- normalizePath(out, mustWork = TRUE)
stopifnot(requireNamespace("jsonlite", quietly = TRUE), requireNamespace("inferguidR", quietly = TRUE))
pkg <- normalizePath(find.package("inferguidR"), winslash = "/", mustWork = TRUE)
stopifnot(file.exists(file.path(pkg, "Meta", "package.rds")))
expected <- Sys.getenv("INFERGUIDR_EXPECTED_LIBRARY")
if (nzchar(expected) && pkg != normalizePath(file.path(expected, "inferguidR"), winslash = "/", mustWork = TRUE)) {
  stop("loaded package differs from INFERGUIDR_EXPECTED_LIBRARY")
}
filearg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
stopifnot(length(filearg) == 1L)
here <- dirname(normalizePath(sub("^--file=", "", filearg)))
for (f in c("continuous.R", "unconditional.R", "glm.R", "correlations.R", "multivariate.R", "proportions.R")) {
  source(file.path(here, "reference", f))
}
started <- Sys.time()
negative <- length(args) == 2L
injected <- FALSE
results_path <- file.path(out, "results.csv")
cases_path <- file.path(out, "cases.csv")
comparisons <- failures <- cases <- 0L
current <- NULL
case_rows <- list()
json <- function(x) as.character(jsonlite::toJSON(x, auto_unbox = TRUE, null = "null", digits = 17))
append_csv <- function(row, path) utils::write.table(row, path, sep = ",", row.names = FALSE,
  col.names = !file.exists(path), append = file.exists(path), qmethod = "double", na = "")
record <- function(quantity, observed, reference, absolute = 1e-10, relative = 1e-8, oracle = "independent mathematical reference") {
  if (negative && !injected && identical(quantity, "power") && is.numeric(observed) && length(observed) == 1L && is.finite(observed)) {
    observed <- observed + 0.01
    injected <<- TRUE
  }
  valid <- is.numeric(observed) && is.numeric(reference) && length(observed) == 1L &&
    length(reference) == 1L && is.finite(observed) && is.finite(reference)
  error <- if (valid) abs(observed - reference) else NA_real_
  passed <- isTRUE(valid && error <= absolute + relative * abs(reference))
  comparisons <<- comparisons + 1L
  if (!passed) failures <<- failures + 1L
  append_csv(data.frame(case = current$id, method = current$method, mode = current$mode,
    quantity = quantity, observed = if (valid) observed else NA_real_, reference = if (valid) reference else NA_real_,
    absolute_error = error, absolute_tolerance = absolute, relative_tolerance = relative,
    passed = passed, oracle = oracle), results_path)
  invisible(passed)
}
condition <- function(quantity, value, oracle = "declared interface contract") record(quantity, as.numeric(isTRUE(value)), 1, 0, 0, oracle)
solve <- function(method, mode, inputs, options) inferguidR::calculate_power(inferguidR::analysis_spec(method, mode, inputs, options))
run_case <- function(id, method, mode, inputs, options = list(), check, expected_status = "ok", kind = "independent_numeric") {
  current <<- list(id = id, method = method, mode = mode)
  before <- failures
  t0 <- proc.time()[["elapsed"]]
  error <- ""
  status <- "not_returned"
  warnings <- ""
  tryCatch({
    r <- solve(method, mode, inputs, options)
    status <- r$status$code
    warnings <- json(r$warnings)
    condition("expected_status", identical(status, expected_status))
    if (identical(status, expected_status)) check(r)
  }, error = function(e) {
    error <<- conditionMessage(e)
    condition("case_completed_without_error", FALSE)
  })
  cases <<- cases + 1L
  row <- data.frame(case = id, method = method, mode = mode, inputs_json = json(inputs), options_json = json(options),
    evidence_kind = kind, status = status, passed = failures == before, elapsed_seconds = proc.time()[["elapsed"]] - t0,
    warnings_json = warnings, error = error)
  case_rows[[length(case_rows) + 1L]] <<- row
  append_csv(row, cases_path)
  if (cases %% 25L == 0L) cat(sprintf("completed %d cases; %d failed comparisons\n", cases, failures))
  invisible(NULL)
}

# N1/N2: exact table enumeration. This statistic is not an exact-size test.
unc <- "exact_proportions_unconditional"
stats <- c("z_pooled", "z_pooled_cc", "z_unpooled", "z_unpooled_cc", "mantel_haenszel", "likelihood_ratio", "t_test")
for (st in stats) for (tl in c("one", "two")) for (v in list(c(4, 6, .7, .3), c(7, 3, .2, .6), c(3, 5, 1 - 1e-8, .4), c(5, 3, 1e-8, .6))) {
  inp <- list(n1 = v[1], n2 = v[2], p1 = v[3], p2 = v[4], alpha = .05)
  opt <- list(effect_size_index = "proportion_p1", computation = "exact", statistic = st, tails = tl)
  id <- sprintf("unconditional-forward-%s-%s-%d-%g", st, tl, v[1], v[3])
  run_case(id, unc, "post_hoc", inp, opt, function(r) {
    ref <- g9_ref_power(st, v[1], v[2], v[3], v[4], .05, tl)
    record("power", r$achieved$power, ref[["power"]], 2e-12, 0, "full independent binomial table enumeration")
    record("beta", r$achieved$beta, ref[["beta"]], 2e-12, 0)
    record("alpha_attained", r$achieved$alpha_attained, ref[["alpha"]], 2e-12, 0)
  })
}
for (st in c("z_pooled", "z_pooled_cc", "z_unpooled", "z_unpooled_cc", "t_test")) for (tl in c("one", "two")) {
  run_case(paste("unconditional-offset", st, tl, sep = "-"), "exact_proportions_unconditional_offset", "post_hoc",
    list(p1_h0 = .4, p1_h1 = .8, p2 = .3, n1 = 5, n2 = 7, alpha = .05),
    list(effect_size_index = "proportion_p1", computation = "exact", statistic = st, tails = tl), function(r) {
      ref <- g9_offset_power(st, 5, 7, .4, .8, .3, .05, tl)
      record("power", r$achieved$power, ref[["power"]], 2e-12, 0, "full table enumeration with nonzero null difference")
      record("beta", r$achieved$beta, ref[["beta"]], 2e-12, 0)
      record("alpha_attained", r$achieved$alpha_attained, ref[["alpha"]], 2e-12, 0)
    })
}
for (st in c("z_pooled", "z_unpooled_cc", "t_test")) for (ratio in c(.5, 2)) {
  inp <- list(p1 = .85, p2 = .25, alpha = .05, power = .7, allocation_ratio = ratio)
  opt <- list(effect_size_index = "proportion_p1", computation = "exact", statistic = st, tails = "two")
  run_case(paste("unconditional-minimality", st, ratio, sep = "-"), unc, "a_priori", inp, opt, function(r) {
    n1 <- 1:80; n2 <- floor(ratio * n1 + .5)
    keep <- n2 >= 1 & (st != "t_test" | n1 + n2 >= 3)
    n1 <- n1[keep]; n2 <- n2[keep]
    pw <- vapply(seq_along(n1), function(j) g9_ref_power(st, n1[j], n2[j], .85, .25, .05, "two")[["power"]], 0.0)
    first <- which(pw >= .7)[1L]
    condition("reference_scan_bracketed", !is.na(first))
    record("smallest_total_sample_size", r$solved$n_total, n1[first] + n2[first], 0, 0)
    record("power", r$achieved$power, pw[first], 2e-12, 0)
  })
}
for (direction in c("p_h1_ge_p_h0", "p_h1_le_p_h0")) for (index in c("proportion_p1", "difference", "risk_ratio", "odds_ratio")) {
  opt <- list(effect_size_index = index, computation = "exact", statistic = "z_pooled", tails = "one", effect_direction = direction)
  run_case(paste("unconditional-sensitivity", direction, index, sep = "-"), unc, "sensitivity",
    list(p2 = .5, n1 = 30, n2 = 20, alpha = .05, power = .7), opt, function(r) {
      p1 <- r$solved$p1
      ref <- g9_ref_power("z_pooled", 30, 20, p1, .5, .05, "one")
      record("reference_power_at_solution", ref[["power"]], .7, 2e-8, 0)
      delta <- if (direction == "p_h1_ge_p_h0") -1e-6 else 1e-6
      earlier <- g9_ref_power("z_pooled", 30, 20, p1 + delta, .5, .05, "one")[["power"]]
      condition("smaller_effect_below_target", earlier < .7)
    })
}
for (st in c("z_pooled", "z_pooled_cc", "mantel_haenszel", "z_unpooled_cc", "t_test")) {
  opt <- list(effect_size_index = "proportion_p1", computation = "exact", statistic = st, tails = "two")
  run_case(paste0("unconditional-criterion-", st), unc, "criterion",
    list(p1 = .8, p2 = .3, n1 = 8, n2 = 12, power = .7), opt, function(r) {
      ref <- g9_ref_power(st, 8, 12, .8, .3, r$solved$alpha, "two")
      record("power", r$achieved$power, ref[["power"]], 2e-12, 0)
      condition("target_reached_at_discrete_boundary", ref[["power"]] >= .7)
      earlier <- g9_ref_power(st, 8, 12, .8, .3, r$solved$alpha * (1 - 1e-8), "two")[["power"]]
      condition("immediately_lower_alpha_below_target", earlier < .7)
    })
}

# Enumerate every deterministic rejection region available as the critical
# threshold moves. Equality beta=q*attained_alpha need not be attainable.
# These states include a superset of the thresholds permitted by alpha in (0,1).
g <- expand.grid(x1 = 0:4, x2 = 0:6)
z <- g9_ref_stat("z_pooled", g$x1, g$x2, 4, 6)
p0 <- as.vector(outer(g9_pmf(4, .3), g9_pmf(6, .3)))
p1 <- as.vector(outer(g9_pmf(4, .7), g9_pmf(6, .3)))
cuts <- sort(unique(c(-Inf, z, Inf)))
states <- do.call(rbind, lapply(cuts, function(cut) {
  member <- z >= cut
  data.frame(critical = cut, alpha_attained = sum(p0[member]), beta = sum(p1[!member]))
}))
states$error_ratio <- 1
states$residual <- states$beta - states$error_ratio * states$alpha_attained
utils::write.csv(states, file.path(out, "exact-compromise-states.csv"), row.names = FALSE)
for (id in c(unc, "exact_proportions_unconditional_offset")) {
  inp <- if (id == unc) list(p1 = .7, p2 = .3) else list(p1_h0 = .4, p1_h1 = .7, p2 = .3)
  run_case(paste0("compromise-guard-", id), id, "compromise", c(inp, list(n1 = 4, n2 = 6, error_ratio = 1)),
    list(effect_size_index = "proportion_p1", computation = "exact"), function(r) {
      condition("no_invented_alpha", is.null(r$solved$alpha))
      condition("no_deterministic_exact_ratio_in_example", all(abs(states$residual) > 1e-12), "exhaustive finite rejection-state enumeration")
    }, expected_status = "reference_pending", kind = "guard_contract")
}

# N3: bounded GLM quadrature. No exemplary data set is materialised.
pars <- list(binomial = list(pi = .3), normal = list(mu = .5, sigma = 1.2),
  exponential = list(rate = 3), poisson = list(lambda = .7),
  uniform = list(a = -1, b = 2), lognormal = list(mu = 0, sigma = .4))
public_pars <- list(binomial = list(x_pi = .3), normal = list(x_mu = .5, x_sigma = 1.2),
  exponential = list(x_lambda = 3), poisson = list(x_lambda = .7),
  uniform = list(x_left = -1, x_right = 2), lognormal = list(x_mu = 0, x_sigma = .4))
for (family in c("logistic", "poisson")) for (xd in names(pars)) for (vc in c("on", "off")) for (tl in c("one", "two")) {
  procedures <- if (family == "logistic") c("demidenko", "hsieh") else c("demidenko", "signorini")
  for (procedure in procedures) {
    if (family == "logistic" && procedure == "hsieh" && vc == "off") next
    method <- paste0("z_", family, "_regression")
    opt <- list(x_distribution = xd, variance_correction = vc, tails = tl, test_procedure = procedure)
    if (family == "logistic" && procedure == "hsieh") opt$variance_correction <- NULL
    common <- list(n_total = 200, alpha = .05, r2_other_x = .1)
    inp <- c(if (family == "logistic") list(p1 = .3, odds_ratio = .8) else list(base_rate = .85, exp_beta1 = .8), common, public_pars[[xd]])
    cs <- list(p1 = .3, or = .8, base = .85, eb1 = .8, n = 200, alpha = .05, r2 = .1, t = 1,
      xd = xd, par = pars[[xd]], proc = procedure, vc = vc, tails = tl)
    run_case(paste("glm", family, xd, vc, tl, procedure, sep = "-"), method, "post_hoc", inp, opt, function(r) {
      ref <- if (family == "logistic") ck_logistic_power(cs) else ck_poisson_power(cs)
      record("power", r$achieved$power, ref$power, 1e-8, 1e-6, "own Gauss-Legendre information integral / closed-form moments")
      record("beta", r$achieved$beta, ref$beta, 1e-8, 1e-6)
    })
  }
}

for (family in c("logistic", "poisson")) for (xd in names(pars)) for (statistic in c("wald", "likelihood_ratio")) for (tl in c("one", "two")) {
  inp <- c(if (family == "logistic") list(p1 = .3, odds_ratio = .8) else list(base_rate = .85, exp_beta1 = .8),
    list(n_total = 60, alpha = .05, r2_other_x = .1), public_pars[[xd]])
  opt <- list(computation = "enumeration", x_distribution = xd, test_statistic = statistic, tails = tl)
  cs <- list(p1 = .3, or = .8, base = .85, eb1 = .8, n = 60, alpha = .05, r2 = .1, t = 1,
    xd = xd, par = pars[[xd]], stat = statistic, tails = tl)
  run_case(paste("glm-enumeration", family, xd, statistic, tl, sep = "-"), paste0("z_", family, "_regression"),
    "post_hoc", inp, opt, function(r) {
      lambda <- ck_enum_lambda(cs, family)
      ref <- ck_enum_power(lambda, .05, tl)
      record("lambda", r$quantities$ncp, lambda, 1e-8, 1e-6, "bounded exemplary-grid information matrix / Poisson IRLS row fit")
      record("power", r$achieved$power, ref$power, 1e-8, 1e-6)
    })
}
for (family in c("logistic", "poisson")) for (xd in c("binomial", "normal")) for (mode in c("a_priori", "sensitivity", "criterion", "compromise")) {
  inp <- c(if (family == "logistic") list(p1 = .3, odds_ratio = 1.5) else list(base_rate = .85, exp_beta1 = 1.5),
    list(n_total = 150, alpha = .05, r2_other_x = .1, power = .7), public_pars[[xd]])
  if (mode == "a_priori") inp$n_total <- NULL
  if (mode == "sensitivity") { inp$odds_ratio <- NULL; inp$exp_beta1 <- NULL }
  if (mode == "criterion") inp$alpha <- NULL
  if (mode == "compromise") { inp$alpha <- inp$power <- NULL; inp$error_ratio <- 1 }
  opt <- list(x_distribution = xd, test_procedure = "demidenko", variance_correction = "on", tails = "two")
  cs <- list(p1 = .3, or = 1.5, base = .85, eb1 = 1.5, n = 150, alpha = .05, r2 = .1, t = 1,
    xd = xd, par = pars[[xd]], proc = "demidenko", vc = "on", tails = "two")
  run_case(paste("glm-inverse", family, xd, mode, sep = "-"), paste0("z_", family, "_regression"), mode, inp, opt, function(r) {
    model <- if (family == "logistic") ck_logistic_power else ck_poisson_power
    if (mode == "a_priori") cs$n <- r$solved$n_total
    if (mode == "sensitivity") { cs$or <- r$solved$odds_ratio; cs$eb1 <- r$solved$exp_beta1 }
    if (mode %in% c("criterion", "compromise")) cs$alpha <- r$solved$alpha
    ref <- model(cs)
    record("power", r$achieved$power, ref$power, 1e-8, 1e-6)
    if (mode == "a_priori") {
      condition("target_reached", ref$power >= .7)
      cs$n <- cs$n - 1
      condition("predecessor_below_target", model(cs)$power < .7,
        "analytic monotonicity of normal power in sqrt(N), with fixed information variances")
    } else if (mode == "compromise") record("beta_equals_nominal_alpha", ref$beta, cs$alpha, 1e-8, 1e-6) else
      record("reference_target_power", ref$power, .7, 1e-8, 1e-6)
  })
}

# N3 two-proportion normal approximations, including unequal allocation.
for (tl in c("one", "two")) for (cc in c("on", "off")) for (arcsine in c("on", "off")) {
  opt <- list(tails = tl, continuity_correction = cc, arcsin_transform = arcsine)
  for (v in list(c(.2, .55, 40, 120), c(.7, .3, 15, 20), c(.5, .51, 2, 3))) {
    run_case(paste("z-proportions", tl, cc, arcsine, paste(v, collapse = "-"), sep = "-"),
      "z_proportions_two_independent", "post_hoc", list(p1 = v[1], p2 = v[2], n1 = v[3], n2 = v[4], alpha = .05), opt, function(r) {
        ref <- ref_z2p_power(v[1], v[2], v[3], v[4], .05, tl, cc == "on", arcsine == "on")
        record("power", r$achieved$power, ref$power, 1e-10, 1e-8, "FTU relation inverted by bisection and own normal quadrature")
      })
  }
  # The unequal-allocation CC branches have historically had jump targets.
  # A returned solution must meet the reference target; a genuine unattained
  # jump must remain an explicit numerical_failure, never a fabricated root.
  run_case(paste("z-proportions-criterion", tl, cc, arcsine, sep = "-"), "z_proportions_two_independent", "criterion",
    list(p1 = .2, p2 = .55, n1 = 40, n2 = 120, power = .8), opt, function(r) {
      ref <- ref_z2p_power(.2, .55, 40, 120, r$solved$alpha, tl, cc == "on", arcsine == "on")
      record("reference_target_power", ref$power, .8, 1e-8, 0)
    })
}

# N3 rank procedures: independent efficiency integrals and moment formulas.
rankids <- c("t_wilcoxon_signed_rank_one_sample", "t_wilcoxon_signed_rank_matched_pairs", "t_wilcoxon_mann_whitney")
for (method in rankids) for (parent in c("normal", "laplace", "logistic")) for (tl in c("one", "two")) {
  density <- switch(parent, normal = function(x) dnorm(x), laplace = function(x) exp(-sqrt(2) * abs(x)) / sqrt(2),
    logistic = function(x) dlogis(x, scale = sqrt(3) / pi))
  efficiency <- 12 * integrate(function(x) density(x)^2, -Inf, Inf, rel.tol = 1e-12)$value^2
  two <- method == rankids[3]
  inp <- c(list(alpha = .05), if (two) list(d = .5, n1 = 20, n2 = 30) else list(d = .5, n_total = 30))
  if (method == rankids[2]) names(inp)[names(inp) == "d"] <- "dz"
  run_case(paste("rank-are", method, parent, tl, sep = "-"), method, "post_hoc", inp,
    list(method = "are", parent_distribution = parent, tails = tl), function(r) {
      df <- if (two) 50 * efficiency - 2 else 30 * efficiency - 1
      ncp <- .5 * sqrt(if (two) 20 * 30 * efficiency / 50 else 30 * efficiency)
      crit <- ref_qt_upper(if (tl == "two") .025 else .05, df)
      ref <- ref_pt(crit, df, ncp, FALSE) + if (tl == "two") ref_pt(-crit, df, ncp) else 0
      record("power", r$achieved$power, ref, 1e-8, 1e-6, "density-square efficiency integral and normal/chi-square t integral")
    })
}
for (tl in c("one", "two")) for (cc in c("on", "off")) {
  run_case(paste("rank-lehmann", tl, cc, sep = "-"), rankids[3], "post_hoc",
    list(p1 = .62, p2 = .47, p3 = .42, n1 = 5, n2 = 40, alpha = .05),
    list(method = "lehmann", lehmann_input_method = "manual_moments", tails = tl, continuity_correction = cc), function(r) {
      e0 <- 5 * 40 / 2; v0 <- 5 * 40 * 46 / 12
      e1 <- 5 * 40 * .62; v1 <- 5 * 40 * (.62 * .38 + 39 * (.47 - .62^2) + 4 * (.42 - .62^2))
      zc <- ck_zq(if (tl == "two") .025 else .05)
      correction <- if (cc == "on") .5 else 0
      ref <- ck_Q((e0 + zc * sqrt(v0) + correction - e1) / sqrt(v1))
      if (tl == "two") ref <- ref + ck_Phi((e0 - zc * sqrt(v0) - correction - e1) / sqrt(v1))
      record("power", r$achieved$power, ref, 1e-10, 1e-8, "Lehmann moments with independently integrated normal tails")
    })
}

# N3/N4: independent distribution kernels, including regression and RM models.
fcases <- list(
  list("f_rm_anova_within", list(f = .27, alpha = .05, n_total = 60, k = 3, m = 4, rho = .35)),
  list("f_rm_anova_between", list(f = .27, alpha = .05, n_total = 60, k = 3, m = 4, rho = .35)),
  list("f_rm_manova_within", list(f = .27, alpha = .05, n_total = 60, k = 3, m = 4, rho = .35)),
  list("f_manova_global", list(f2 = .08, alpha = .05, n_total = 75, n_groups = 3, n_responses = 2)),
  list("f_hotelling_t2_one_group", list(delta = .5, alpha = .05, n_total = 40, n_response = 3)),
  list("f_hotelling_t2_two_groups", list(delta = .5, alpha = .05, n1 = 30, n2 = 20, n_response = 3)),
  list("f_mrc_r2_deviation_from_zero", list(f2 = .15, alpha = .05, n_total = 80, n_predictors = 3)),
  list("f_mrc_r2_increase", list(f2 = .15, alpha = .05, n_total = 80, n_predictors = 5, n_tested = 2)),
  list("f_generic", list(lambda = 4, df1 = 3, df2 = 40, alpha = .05)))
for (cs in fcases) run_case(paste0("f-integral-", cs[[1]]), cs[[1]], "post_hoc", cs[[2]], list(), function(r) {
  df <- unname(r$quantities$df); crit <- r$quantities$critical; lambda <- r$quantities$ncp
  record("critical_null_tail", ref_pf(crit, df[1], df[2], 0, FALSE), .05, 1e-8, 1e-6,
    "normal-square convolution integrated over central chi-square denominator")
  record("power_from_reported_parameters", r$achieved$power, ref_pf(crit, df[1], df[2], lambda, FALSE), 1e-8, 1e-6)
  if (cs[[1]] == "f_hotelling_t2_one_group") {
    record("ncp", lambda, 40 * .5^2, 1e-12, 0, "Hotelling quadratic-form noncentrality")
    record("df1", df[1], 3, 0, 0); record("df2", df[2], 37, 0, 0)
  }
  if (cs[[1]] == "f_hotelling_t2_two_groups") {
    record("ncp", lambda, (30 * 20 / 50) * .5^2, 1e-12, 0, "Hotelling quadratic-form noncentrality")
    record("df1", df[1], 3, 0, 0); record("df2", df[2], 46, 0, 0)
  }
  if (cs[[1]] %in% c("f_mrc_r2_deviation_from_zero", "f_mrc_r2_increase")) {
    record("ncp", lambda, .15 * 80, 1e-12, 0, "fixed-design regression noncentrality N*f2")
    record("df1", df[1], if (cs[[1]] == "f_mrc_r2_increase") 2 else 3, 0, 0)
    record("df2", df[2], if (cs[[1]] == "f_mrc_r2_increase") 74 else 76, 0, 0)
  }
  if (cs[[1]] %in% c("f_rm_anova_within", "f_rm_manova_within", "f_rm_anova_between")) {
    S <- matrix(.35, 4, 4); diag(S) <- 1
    v <- 1:4 - 2.5; mu <- .27 * v / sqrt(mean(v^2))
    C <- sapply(1:3, function(j) { v <- c(rep(1, j), -j, rep(0, 3 - j)); v / sqrt(sum(v^2)) })
    reference <- if (cs[[1]] == "f_rm_anova_between") 60 * .27^2 * 16 / sum(S) else
      60 * as.numeric(t(t(C) %*% mu) %*% base::solve(t(C) %*% S %*% C, t(C) %*% mu))
    record("ncp_population_quadratic_form", lambda, reference, 1e-10, 1e-10, "explicit covariance matrix and contrast quadratic form")
  }
})
for (cs in list(c(3, .5), c(40, 4), c(120, -3))) {
  run_case(paste0("generic-t-", cs[1]), "t_generic", "post_hoc", list(delta = cs[2], df = cs[1], alpha = .05),
    list(tails = "two"), function(r) {
      crit <- ref_qt_upper(.025, cs[1])
      record("power", r$achieved$power, ref_pt(crit, cs[1], cs[2], FALSE) + ref_pt(-crit, cs[1], cs[2]), 1e-8, 1e-6)
    })
}

# N5: selected tetrachoric cells and Pearson asymptotic SE, not finite-N exactness.
for (v in list(c(.3, 0, .5, .5, 200), c(-.4, .1, .3, .7, 120), c(.6, .2, .2, .4, 80))) for (tl in c("one", "two")) for (computation in c("exact_brown_benedetti", "approx_bonett_price")) {
  run_case(paste("tetrachoric", paste(v, collapse = "-"), tl, computation, sep = "-"), "z_tetrachoric_correlation", "post_hoc",
    list(rho_h1 = v[1], rho_h0 = v[2], marginal_x = v[3], marginal_y = v[4], n_total = v[5], alpha = .05),
    list(tails = tl, computation_of_correlation = computation), function(r) {
      c0 <- g7_cells(v[2], v[3], v[4]); c1 <- g7_cells(v[1], v[3], v[4])
      se0 <- g7_se_pearson(c0, v[2], v[5]); se1 <- g7_se_pearson(c1, v[1], v[5])
      est0 <- v[2]; est1 <- v[1]
      if (computation == "approx_bonett_price") {
        bp0 <- ref_bp(c0, v[3], v[4], v[5]); bp1 <- ref_bp(c1, v[3], v[4], v[5])
        se0 <- bp0$se; se1 <- bp1$se; est0 <- bp0$r; est1 <- bp1$r
      }
      ref <- g7_zpower((est1 - est0) / se0, se1 / se0, .05, tl, sign(est1 - est0))
      record("se_h0", r$method_outputs$se_h0, se0, 1e-10, 1e-8, "orthant quadrature and Pearson closed-form asymptotic variance")
      record("se_h1", r$method_outputs$se_h1, se1, 1e-10, 1e-8)
      record("power", r$achieved$power, ref, 1e-10, 1e-8)
    })
}

# Complete declaration inventory. Evidence attaches only to actual case options,
# with defaults resolved from the installed registry. A passing row is a sample,
# never a claim that every numerical input or option interaction is covered.
source(file.path(here, "inventory.R"))
inventory <- write_numerical_inventory(pkg, do.call(rbind, case_rows), out)
coverage_complete <- cases == 313L
if (!coverage_complete) failures <- failures + 1L
summary <- list(scope = "bounded synthetic numerical campaign; no whole-method certification",
  started_utc = format(started, "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"), finished_utc = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
  package_version = as.character(utils::packageVersion("inferguidR")), package_path = pkg,
  r_version = R.version.string, platform = R.version$platform,
  source_sha256_caller_supplied = Sys.getenv("INFERGUIDR_SOURCE_SHA256"), source_commit_caller_supplied = Sys.getenv("INFERGUIDR_SOURCE_COMMIT"),
  negative_control = negative, negative_control_injected = injected, cases = cases, comparisons = comparisons,
  planned_cases = 313L, planned_coverage_complete = coverage_complete,
  failed_comparisons = failures, inventory = inventory, passed = failures == 0L && (!negative || injected),
  session_info = capture.output(sessionInfo()))
if (requireNamespace("digest", quietly = TRUE)) {
  files <- c(file.path(here, "check-numerics.R"), file.path(here, "inventory.R"), list.files(file.path(here, "reference"), full.names = TRUE))
  summary$checker_sha256 <- stats::setNames(lapply(files, function(f) digest::digest(file = f, algo = "sha256")), substring(files, nchar(here) + 2L))
}
jsonlite::write_json(summary, file.path(out, "summary.json"), auto_unbox = TRUE, pretty = TRUE, null = "null", digits = 17)
cat(sprintf("%d cases, %d comparisons, %d failed comparisons; %d inventoried methods\n", cases, comparisons, failures, inventory$methods))
quit(status = if (summary$passed) 0L else 1L)
