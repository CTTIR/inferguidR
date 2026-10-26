#!/usr/bin/env Rscript

# Installed-package acceptance checks; use --minimal with optional packages absent.
args <- commandArgs(trailingOnly = TRUE)
if (!length(args) || length(args) > 2L ||
    (length(args) == 2L && args[[2L]] != "--minimal")) {
  stop("Usage: Rscript tools/release/check-package.R OUTPUT_DIRECTORY [--minimal]", call. = FALSE)
}
out <- args[[1L]]
minimal <- length(args) == 2L
if (dir.exists(out) && length(list.files(out, all.files = TRUE, no.. = TRUE)))
  stop("OUTPUT_DIRECTORY must be new or empty", call. = FALSE)
dir.create(out, recursive = TRUE, showWarnings = FALSE)
if (!dir.exists(out)) stop("Cannot create output directory", call. = FALSE)
out <- normalizePath(out, winslash = "/", mustWork = TRUE)
started <- format(Sys.time(), tz = "UTC", usetz = TRUE)
checks <- character()
exports <- aliases <- loaded_on_startup <- optional_present <- character()
package_path <- NULL
check <- function(ok, description) {
  if (!isTRUE(ok)) stop(description, call. = FALSE)
  checks <<- c(checks, description)
}
unsupported <- function(expr) {
  tryCatch({ force(expr); FALSE }, error = function(e)
    inherits(e, "inferguidR_unsupported") && identical(e$status$code, "unsupported"))
}
expected <- Sys.getenv("INFERGUIDR_EXPECTED_LIBRARY")
failure <- tryCatch({
  check(nzchar(expected), "INFERGUIDR_EXPECTED_LIBRARY is required for candidate identity")
  package_path <- normalizePath(find.package("inferguidR"), winslash = "/", mustWork = TRUE)
  check(identical(package_path, normalizePath(file.path(expected, "inferguidR"),
    winslash = "/", mustWork = TRUE)), "Loaded package matches the candidate library")

  # Set a documented user option and a seed before loading to detect overrides.
  options(inferguidR.validate_json_schema = FALSE)
  set.seed(9417)
  seed <- .Random.seed
  before_options <- options()
  before_wd <- getwd()
  before_files <- list.files(all.files = TRUE, no.. = TRUE)
  before_namespaces <- loadedNamespaces()
  suppressPackageStartupMessages(library(inferguidR))
  loaded_on_startup <- setdiff(loadedNamespaces(), before_namespaces)
  check(identical(getwd(), before_wd), "Startup preserves the working directory")
  check(identical(list.files(all.files = TRUE, no.. = TRUE), before_files),
        "Startup creates no files in the working directory")
  check(identical(.Random.seed, seed), "Startup preserves the random-number state")
  check(identical(getOption("inferguidR.validate_json_schema"), FALSE),
        "Startup preserves an explicitly configured option")
  changed <- union(names(before_options), names(options()))
  changed <- changed[!vapply(changed, function(n) identical(before_options[[n]], getOption(n)), TRUE)]
  check(!length(changed), "Startup does not change global options")

  exports <- sort(getNamespaceExports("inferguidR"))
  rd <- tools::Rd_db("inferguidR")
  aliases <- unique(unlist(lapply(rd, function(doc) unlist(lapply(doc, function(node)
    if (identical(attr(node, "Rd_tag"), "\\alias")) as.character(node) else NULL)))))
  check(!length(setdiff(exports, aliases)), "Every exported API has an installed help alias")
  check(all(vapply(exports, function(n) exists(n, asNamespace("inferguidR"), inherits = FALSE), TRUE)),
        "Every exported API resolves in the namespace")
  description <- utils::packageDescription("inferguidR")
  check(identical(description$License, "MIT + file LICENSE"), "Installed licence declaration is explicit")
  check(file.exists(file.path(package_path, "LICENSE")), "Installed copyright file is present")

  optional <- c("shiny", "ggplot2", "jsonvalidate", "pROC", "Exact", "exact2x2", "WebPower")
  optional_present <- optional[vapply(optional, requireNamespace, TRUE, quietly = TRUE)]
  check(!any(optional %in% loaded_on_startup), "Startup does not load optional application or validation packages")
  if (minimal) check(!length(optional_present), "Selected optional packages are absent in the minimal profile")

  spec <- analysis_spec("t_means_two_independent", "a_priori",
    inputs = list(alpha = .05, power = .8, d = .5), options = list(tails = "two"))
  result <- solve_design(spec)
  check(identical(result$status$code, "ok") && result$solved$n_total == 128L,
        "Core design works with the active dependency profile")
  check(identical(ig_hash(result), ig_hash(ig_from_json(ig_to_json(result)))),
        "Core result round-trips with the active dependency profile")
  if (!"shiny" %in% optional_present)
    check(unsupported(workbench_app()), "Missing Shiny is reported as typed unsupported")
  if (!"ggplot2" %in% optional_present)
    check(unsupported(plot_power(result)), "Missing ggplot2 is reported as typed unsupported")

  report <- create_report_snapshot(analysis_envelope(spec, result), title = "Dependency profile check")
  for (format in c("json", "txt")) {
    destination <- file.path(out, paste0("report.", format))
    exported <- export_report(report, destination, format)
    check(file.exists(destination) && file.info(destination)$size > 0,
          paste("Offline", format, "report export succeeds"))
    check(identical(exported$sha256, digest::digest(file = destination, algo = "sha256")),
          paste("Exported", format, "report hash matches its bytes"))
  }
  imported <- ig_from_json(paste(readLines(file.path(out, "report.json"), warn = FALSE), collapse = "\n"))
  check(identical(ig_hash(imported), ig_hash(report)), "Exported report round-trips without optional validation packages")

  scored <- classification_spec(data.frame(y = c(0, 1, 0, 1), p = c(.1, .4, .6, .9)),
    truth = "y", positive_label = 1, scores = "p", score_scale = "probability", cutoff = .5)
  evaluated <- evaluate_cutoff(scored, calibration_bins = 2, threshold_probability = c(.25, .75))
  probability_report <- create_report_snapshot(
    analysis_envelope(classification_evaluated_spec(evaluated), evaluated),
    include_probability_summaries = TRUE)
  check(identical(probability_report$schema_version, "inferguidR.report_document/0.4.0"),
        "Opt-in probability report uses its explicit schema version")
  check(isFALSE(probability_report$privacy$raw_rows_included), "Probability report omits raw scored rows")
  check(identical(ig_hash(probability_report), ig_hash(ig_from_json(ig_to_json(probability_report)))),
        "Opt-in probability aggregates round-trip in the active dependency profile")

  advisor <- advisor_replay(c(q.aim = "compare", q.stage = "planning", q.outcome_type = "binary",
    q.unit = "individuals", q.structure = "independent", q.independence_check = "no",
    q.n_factors = "one", q.n_groups = "two", q.covariates = "no",
    q.binary_margins = "group_sizes_fixed", q.direction = "two_sided", q.missingness = "none"))
  mapped <- advisor_power_spec(advisor, list(p1 = .2, p2 = .7, alpha = .13, n1 = 2, n2 = 2),
                              mode = "post_hoc")
  mapped_result <- calculate_power(mapped)
  expected_power <- (1 - .2)^2 * .7^2 + .2^2 * (1 - .7)^2
  check(identical(mapped_result$status$code, "ok") &&
          abs(mapped_result$solved$power - expected_power) < 1e-9,
        "Advisor bridge agrees with the elementary two-extreme-table design")
  lattice_report <- create_report_snapshot(analysis_envelope(mapped, mapped_result))
  lattice_figures <- Filter(function(v) identical(v$plot_type, "power_lattice"),
                            lattice_report$visualisations)
  check(identical(lattice_report$schema_version, "inferguidR.report_document/0.4.0") &&
          length(lattice_figures) == 2L &&
          all(vapply(lattice_figures, function(v) length(v$data) == 9L, TRUE)),
        "Boschloo report retains both nine-cell conditional lattices")
  check(identical(ig_hash(lattice_report), ig_hash(ig_from_json(ig_to_json(lattice_report)))),
        "Frozen lattice report round-trips in the active dependency profile")
  check(isFALSE(lattice_report$privacy$probability_aggregates_included),
        "Lattice-only reports do not opt into scored-data probability aggregates")
  NULL
}, error = function(e) conditionMessage(e))
receipt <- list(scope = "installed API/help, startup and dependency-profile checks; not numerical or release certification",
  started_utc = started, finished_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
  r_version = R.version.string, package_path = package_path, expected_library = expected,
  minimal_profile = minimal,
  exported_api_count = length(exports), exports = as.list(exports),
  undocumented_exports = as.list(setdiff(exports, aliases)),
  namespaces_loaded_on_startup = as.list(loaded_on_startup),
  selected_optional_packages_present = as.list(optional_present),
  passed = length(checks), checks = as.list(checks), success = is.null(failure), failure = failure)
jsonlite::write_json(receipt, file.path(out, "receipt.json"), pretty = TRUE, auto_unbox = TRUE, null = "null")
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
if (!is.null(failure)) stop(failure, call. = FALSE)
cat(sprintf("Installed release preflight: %d checks passed; %d exported APIs documented\n", length(checks), length(exports)))
