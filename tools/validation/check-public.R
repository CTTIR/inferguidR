#!/usr/bin/env Rscript
# Copyright (c) 2026 Raban Heller. MIT licence; see ../../LICENSE.md.

# This driver deliberately loads an installed package. It never loads the
# development namespace or reads a private fixture/capture directory.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1L || length(args) > 2L || !nzchar(args[1L]) ||
    (length(args) == 2L && !identical(args[2L], "--negative-control"))) {
  stop("usage: Rscript check-public.R OUTPUT_DIRECTORY [--negative-control]", call. = FALSE)
}
out <- path.expand(args[1L])
if (dir.exists(out) && length(list.files(out, all.files = TRUE, no.. = TRUE))) {
  stop("output directory must be new or empty; previous evidence will not be overwritten", call. = FALSE)
}
if (!dir.exists(out) && !dir.create(out, recursive = TRUE)) {
  stop("cannot create output directory", call. = FALSE)
}
out <- normalizePath(out, mustWork = TRUE)
if (!requireNamespace("jsonlite", quietly = TRUE)) stop("jsonlite is required", call. = FALSE)
script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (length(script_arg) != 1L) stop("run this file with Rscript", call. = FALSE)
script <- normalizePath(sub("^--file=", "", script_arg), mustWork = TRUE)
here <- dirname(script)
negative <- length(args) == 2L
started <- Sys.time()
rows <- list()
failure <- NULL
injected <- FALSE
source_sha <- Sys.getenv("INFERGUIDR_SOURCE_SHA256", unset = "")
expected_library <- Sys.getenv("INFERGUIDR_EXPECTED_LIBRARY", unset = "")

record <- function(case, method, quantity, observed, reference,
                   absolute = 0, relative = 0, oracle, detail = "") {
  if (negative && identical(case, "t-forward-01") && identical(quantity, "power")) {
    observed <- observed + 0.01
    injected <<- TRUE
    detail <- paste(detail, "Negative control: observed power perturbed by +0.01.")
  }
  valid <- is.numeric(observed) && length(observed) == 1L && is.finite(observed) &&
    is.numeric(reference) && length(reference) == 1L && is.finite(reference)
  error <- if (valid) abs(observed - reference) else NA_real_
  tolerance <- if (valid) absolute + relative * abs(reference) else NA_real_
  rows[[length(rows) + 1L]] <<- data.frame(
    case = case, method = method, quantity = quantity,
    observed = if (valid) as.double(observed) else NA_real_,
    reference = if (valid) as.double(reference) else NA_real_,
    absolute_error = error, absolute_tolerance = absolute, relative_tolerance = relative,
    passed = isTRUE(valid && error <= tolerance), oracle = oracle, detail = trimws(detail),
    stringsAsFactors = FALSE
  )
  invisible(NULL)
}
record_condition <- function(case, method, quantity, condition, detail = "") {
  record(case, method, quantity, as.double(isTRUE(condition)), 1,
         oracle = "declared interface contract", detail = detail)
}
require_ok <- function(result, case) {
  if (!identical(result$status$code, "ok")) {
    stop(sprintf("%s returned status %s", case, result$status$code), call. = FALSE)
  }
  invisible(result)
}

run_checks <- function() {
  if (!requireNamespace("inferguidR", quietly = TRUE)) {
    stop("install the candidate inferguidR package before running this checker", call. = FALSE)
  }
  if (nzchar(expected_library)) {
    expected_path <- normalizePath(file.path(expected_library, "inferguidR"),
                                   winslash = "/", mustWork = FALSE)
    loaded_path <- normalizePath(find.package("inferguidR"), winslash = "/", mustWork = TRUE)
    if (!identical(loaded_path, expected_path)) {
      stop("loaded inferguidR is not from INFERGUIDR_EXPECTED_LIBRARY", call. = FALSE)
    }
  }
  if (!file.exists(file.path(find.package("inferguidR"), "Meta", "package.rds"))) {
    stop("inferguidR must be installed; a source/development namespace is not a candidate installation", call. = FALSE)
  }
  if (nzchar(source_sha) && !grepl("^[[:xdigit:]]{64}$", source_sha)) {
    stop("INFERGUIDR_SOURCE_SHA256 must be a 64-character SHA-256 or unset", call. = FALSE)
  }
  source(file.path(here, "reference", "continuous.R"), local = TRUE)
  source(file.path(here, "reference", "boschloo.R"), local = TRUE)
  solve <- function(method, mode, inputs, options = list()) {
    inferguidR::calculate_power(inferguidR::analysis_spec(method, mode, inputs, options))
  }

  # Forward t cases include unequal group sizes and both effect directions.
  tm <- "t_means_two_independent"
  tc <- data.frame(d = c(.5, -.7, .2), n1 = c(64, 20, 10), n2 = c(64, 30, 30),
                   alpha = c(.05, .01, .05), tails = c("two", "one", "two"))
  for (i in seq_len(nrow(tc))) {
    cs <- tc[i, ]; id <- sprintf("t-forward-%02d", i)
    r <- solve(tm, "post_hoc", list(d = cs$d, n1 = cs$n1, n2 = cs$n2, alpha = cs$alpha),
               list(tails = cs$tails))
    require_ok(r, id)
    ref <- ref_power_t2(cs$d, cs$n1, cs$n2, cs$alpha, cs$tails)
    oracle <- "normal-CDF integral over chi-square density; reference-integral quantile"
    record(id, tm, "power", r$achieved$power, ref$power, 1e-8, 1e-6, oracle)
    record(id, tm, "df", r$quantities$df, ref$df, oracle = oracle)
    record_condition(id, tm, "critical_count", length(r$quantities$critical) == length(ref$critical))
    for (j in seq_along(ref$critical)) {
      record(id, tm, paste0("critical_", j), r$quantities$critical[j], ref$critical[j],
             1e-8, oracle = oracle)
    }
  }

  # Independent exhaustive scan of every earlier equal-allocation design.
  # Small effects would make this intentionally slow reference unsuitable for CI.
  ap <- data.frame(d = c(1.2, -1), tails = c("two", "one"), alpha = c(.05, .05), power = c(.8, .8))
  for (i in seq_len(nrow(ap))) {
    cs <- ap[i, ]; id <- sprintf("t-minimality-%02d", i)
    r <- solve(tm, "a_priori", list(d = cs$d, alpha = cs$alpha, power = cs$power),
               list(tails = cs$tails))
    require_ok(r, id)
    scan <- vapply(2:60, function(n) ref_power_t2(cs$d, n, n, cs$alpha, cs$tails)$power, 0.0)
    first <- which(scan >= cs$power)[1L]
    if (is.na(first)) stop("reference sample-size scan did not bracket the target", call. = FALSE)
    n <- (2:60)[first]
    oracle <- "exhaustive reference-integral scan of n1=n2=2,...,60"
    record(id, tm, "n_total", r$solved$n_total, 2 * n, oracle = oracle)
    record(id, tm, "power", r$achieved$power, scan[first], 1e-8, 1e-6, oracle)
    record_condition(id, tm, "all_earlier_designs_below_target",
                     all(scan[seq_len(first - 1L)] < cs$power))
  }
  id <- "t-sensitivity-01"
  r <- solve(tm, "sensitivity", list(n1 = 20, n2 = 30, alpha = .05, power = .8), list(tails = "two"))
  require_ok(r, id)
  record(id, tm, "reference_power_at_solved_d", ref_power_t2(r$solved$d, 20, 30, .05, "two")$power,
         .8, 1e-8, oracle = "normal-CDF integral evaluated at the returned effect")

  cm <- "chi2_generic"
  cc <- data.frame(lambda = c(10, .5, 400), df = c(5, 1, 1000))
  for (i in seq_len(nrow(cc))) {
    cs <- cc[i, ]; id <- sprintf("chi2-forward-%02d", i)
    r <- solve(cm, "post_hoc", list(lambda = cs$lambda, df = cs$df, alpha = .05))
    require_ok(r, id)
    critical <- r$quantities$critical
    oracle <- "convolution of squared shifted normal and central chi-square"
    record(id, cm, "critical_null_tail", ref_pchisq(critical, cs$df, 0, FALSE), .05, 1e-8, 1e-6, oracle)
    record(id, cm, "power", r$achieved$power, ref_pchisq(critical, cs$df, cs$lambda, FALSE), 1e-8, 1e-6, oracle)
    record(id, cm, "beta", r$achieved$beta, ref_pchisq(critical, cs$df, cs$lambda, TRUE), 1e-8, 1e-6, oracle)
  }

  bm <- "exact_proportion_one_sample_binomial"
  bc <- data.frame(n = c(12, 20, 35), p0 = c(.3, .5, .6), g = c(.2, -.2, .1))
  for (i in seq_len(nrow(bc))) {
    cs <- bc[i, ]; id <- sprintf("binomial-forward-%02d", i)
    r <- solve(bm, "post_hoc", list(n_total = cs$n, p0 = cs$p0, g = cs$g, alpha = .05), list(tails = "one"))
    require_ok(r, id)
    x <- 0:cs$n
    coefficient <- ck_choose(cs$n)
    null <- coefficient * cs$p0^x * (1 - cs$p0)^(cs$n - x)
    p1 <- cs$p0 + cs$g
    alt <- coefficient * p1^x * (1 - p1)^(cs$n - x)
    tail <- if (cs$g > 0) rev(cumsum(rev(null))) else cumsum(null)
    member <- tail <= .05
    oracle <- "Pascal coefficients, polynomial probabilities and enumerated one-sided null tail"
    record(id, bm, "attained_alpha", r$achieved$alpha_attained, sum(null[member]), 1e-12, oracle = oracle)
    record(id, bm, "power", r$achieved$power, sum(alt[member]), 1e-12, oracle = oracle)
    record(id, bm, "beta", r$achieved$beta, sum(alt[!member]), 1e-12, oracle = oracle)
  }

  # Exact rational ordering and an independent polynomial nuisance maximisation.
  # The reference maximum is numerical, not a rigorous interval certificate.
  python <- Sys.getenv("PYTHON", unset = "")
  if (!nzchar(python)) python <- Sys.which("python3")
  if (!nzchar(python)) python <- Sys.which("python")
  if (!nzchar(python)) stop("Python 3 is required for the exact-rational Boschloo reference", call. = FALSE)
  ranks_dir <- file.path(out, "boschloo-ranks")
  dir.create(ranks_dir)
  py_output <- system2(python, c(shQuote(file.path(here, "reference", "exact_classes.py")),
                                shQuote(ranks_dir), "5:7", "4:6"), stdout = TRUE, stderr = TRUE)
  writeLines(py_output, file.path(out, "python-output.txt"))
  py_status <- attr(py_output, "status")
  if (!is.null(py_status) && py_status != 0L) stop("exact-rational ordering generation failed", call. = FALSE)
  for (nn in list(c(5, 7), c(4, 6))) {
    d <- ck_design(nn[1], nn[2])
    ranks <- ck_read_ranks(ranks_dir, nn[1], nn[2])
    for (variant in c("less", "greater", "minlike", "central")) {
      id <- paste0("boschloo-", nn[1], "-", nn[2], "-", variant)
      oracle <- "exact-rational conditional ordering; polynomial grid plus golden-section supremum"
      ref_p <- ck_pvalues(d, ranks, variant)
      for (k in seq_along(d$x1)) {
        observed <- inferguidR::boschloo_test(d$x1[k], nn[1], d$x2[k], nn[2], variant = variant)
        require_ok(observed, paste0(id, " table ", d$x1[k], "/", d$x2[k]))
        record(id, "boschloo_proportions_two_independent", paste0("p_value_", d$x1[k], "_", d$x2[k]),
               observed$p_value, ref_p[k], 1e-9, oracle = oracle)
        record_condition(id, "boschloo_proportions_two_independent", paste0("upper_bound_", d$x1[k], "_", d$x2[k]),
                         ref_p[k] <= observed$p_upper + 1e-12)
      }
      alpha <- if (variant %in% c("less", "greater")) .025 else .05
      record_condition(id, "boschloo_proportions_two_independent", "no_ambiguous_reference_decision",
                       all(abs(ref_p - alpha) >= 1e-8),
                       "An alpha-adjacent reference value makes this case inconclusive and fails the check.")
      expected_member <- ref_p <= alpha
      region <- inferguidR::boschloo_region(nn[1], nn[2], alpha = alpha, variant = variant)
      require_ok(region, id)
      actual_member <- logical(length(d$x1))
      actual_member[region$rejected$x1 + 1L + region$rejected$x2 * (nn[1] + 1L)] <- TRUE
      record_condition(id, "boschloo_proportions_two_independent", "rejection_region",
                       identical(actual_member, expected_member))
      record(id, "boschloo_proportions_two_independent", "attained_size", region$size$value,
             ck_size(d, expected_member), 1e-9, oracle = oracle)
      for (pp in list(c(.1, .6), c(.7, .2), c(.5, .5), c(0, 1))) {
        pw <- inferguidR::boschloo_power(nn[1], nn[2], pp[1], pp[2], alpha = alpha, variant = variant)
        require_ok(pw, id)
        record(id, "boschloo_proportions_two_independent", paste0("power_", pp[1], "_", pp[2]),
               pw$power, ck_power(d, expected_member, pp[1], pp[2]), 1e-12, oracle = oracle)
      }
    }
  }

  # A reference check must also preserve a deliberately unsupported branch.
  id <- "exact-compromise-refusal"
  r <- solve("exact_proportions_unconditional", "compromise",
             list(p1 = .65, p2 = .4, error_ratio = 1, n1 = 10, n2 = 12),
             list(effect_size_index = "proportion_p1", computation = "exact"))
  record_condition(id, "exact_proportions_unconditional", "reference_pending",
                   identical(r$status$code, "reference_pending"), paste("Observed status:", r$status$code))
  record_condition(id, "exact_proportions_unconditional", "no_alpha_estimate", is.null(r$solved$alpha))
}

tryCatch(run_checks(), error = function(e) failure <<- conditionMessage(e))
results <- if (length(rows)) do.call(rbind, rows) else data.frame()
utils::write.csv(results, file.path(out, "results.csv"), row.names = FALSE)
if (negative && !injected && is.null(failure)) failure <- "negative control was not injected"
if (is.null(failure) && (length(unique(results$case)) != 21L || nrow(results) != 761L)) {
  failure <- "planned coverage incomplete: expected 21 cases and 761 scalar comparisons"
}
passed <- length(rows) > 0L && is.null(failure) && all(results$passed)
available <- requireNamespace("inferguidR", quietly = TRUE)
refs <- c("reference/continuous.R", "reference/boschloo.R", "reference/exact_classes.py")
hashes <- if (requireNamespace("digest", quietly = TRUE)) {
  stats::setNames(lapply(c("check-public.R", refs), function(f) {
    path <- file.path(here, f)
    if (!file.exists(path)) return(NA_character_)
    tryCatch(digest::digest(file = path, algo = "sha256"), error = function(e) NA_character_)
  }),
                  c("check-public.R", refs))
} else list()
summary <- list(
  scope = "selected synthetic independent mathematical checks; not complete method or G*Power certification",
  started_utc = format(started, "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
  finished_utc = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
  elapsed_seconds = as.double(difftime(Sys.time(), started, units = "secs")),
  package_version = if (available) as.character(utils::packageVersion("inferguidR")) else NULL,
  package_path = if (available) find.package("inferguidR") else NULL,
  expected_library = if (nzchar(expected_library)) expected_library else NULL,
  source_archive_sha256 = if (nzchar(source_sha)) source_sha else NULL,
  source_identity_origin = "Caller-supplied archive hash, not independently verified by this script; package path/version recorded separately.",
  runtime = R.version.string, platform = R.version$platform,
  library_paths = .libPaths(), checker_sha256 = hashes,
  case_count = if (nrow(results)) length(unique(results$case)) else 0L,
  comparison_count = nrow(results),
  expected_case_count = 21L, expected_comparison_count = 761L,
  passed = if (nrow(results)) sum(results$passed) else 0L,
  failed = if (nrow(results)) sum(!results$passed) else 0L,
  negative_control = negative, negative_control_injected = injected,
  fatal_error = failure, status = if (passed) "PASS" else "FAIL"
)
jsonlite::write_json(summary, file.path(out, "summary.json"), auto_unbox = TRUE, pretty = TRUE, null = "null", na = "null")
writeLines(capture.output(utils::sessionInfo()), file.path(out, "session-info.txt"))
cat(sprintf("%s: %d cases; %d comparisons passed, %d failed%s\n", summary$status, summary$case_count,
            summary$passed, summary$failed, if (is.null(failure)) "" else paste0("; ", failure)))
quit(save = "no", status = if (passed) 0L else 1L)
