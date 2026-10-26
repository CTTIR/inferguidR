#!/usr/bin/env Rscript

# Core compatibility only: this deliberately needs no suggested packages.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L || !nzchar(args[[1L]])) {
  stop("Usage: Rscript tools/ci/minimum-r-smoke.R OUTPUT_DIRECTORY", call. = FALSE)
}
out <- args[[1L]]
if (dir.exists(out) && length(list.files(out, all.files = TRUE, no.. = TRUE))) {
  stop("OUTPUT_DIRECTORY must be new or empty; previous evidence is retained", call. = FALSE)
}
dir.create(out, recursive = TRUE, showWarnings = FALSE)
if (!dir.exists(out)) stop("Cannot create OUTPUT_DIRECTORY", call. = FALSE)

checks <- character()
check <- function(condition, description) {
  if (!isTRUE(condition)) stop(description, call. = FALSE)
  checks <<- c(checks, description)
}
started <- format(Sys.time(), tz = "UTC", usetz = TRUE)
failure <- tryCatch({
  check(getRversion() >= "4.1.0", "R meets the declared minimum")
  expected_library <- Sys.getenv("INFERGUIDR_EXPECTED_LIBRARY")
  if (nzchar(expected_library)) {
    expected <- normalizePath(file.path(expected_library, "inferguidR"), winslash = "/", mustWork = TRUE)
    actual <- normalizePath(find.package("inferguidR"), winslash = "/", mustWork = TRUE)
    check(identical(actual, expected), "The package comes from the expected candidate library")
  }
  library(inferguidR)
  methods <- list_methods()
  check(nrow(methods) == 50L && !anyDuplicated(methods$id),
        "All 50 expected method identifiers load without duplicates")
  check(all(c("t_means_two_independent", "boschloo_proportions_two_independent") %in% methods$id),
        "The core inventory and extension are present")
  check(all(methods$implementation_status == "implemented"),
        "Every current inventory entry has a registered engine and solver")

  design <- analysis_spec("t_means_two_independent", "a_priori",
                          inputs = list(alpha = 0.05, power = 0.8, d = 0.5),
                          options = list(tails = "two"))
  result <- solve_design(design)
  check(identical(result$status$code, "ok"), "The two-group design returns an ok result")
  check(identical(result$solved$n_total, 128L), "The standard design uses 128 total observations")
  reference <- stats::power.t.test(n = 64, delta = 0.5, sd = 1, sig.level = 0.05,
                                   type = "two.sample", alternative = "two.sided", strict = TRUE)$power
  check(abs(result$achieved$power - reference) <= 1e-10,
        "The standard design power agrees with the base-R calculation")
  check(identical(ig_hash(design), ig_hash(ig_from_json(ig_to_json(design)))),
        "Analysis specification JSON preserves the hash")
  check(identical(ig_hash(result), ig_hash(ig_from_json(ig_to_json(result)))),
        "Power result JSON preserves the hash")
  invalid <- tryCatch({
    analysis_spec("t_means_two_independent", "a_priori",
                  inputs = list(alpha = 2, power = 0.8, d = 0.5))
    FALSE
  }, error = function(e) TRUE)
  check(invalid, "An invalid probability is rejected")

  classified <- evaluate_cutoff(classification_spec(counts = c(tp = 8, fn = 2, fp = 1, tn = 9)))
  check(identical(classified$status$code, "ok"), "Count-based classification returns an ok result")
  tab <- classification_table(classified, "metrics")
  expected <- c(tpr = 0.8, tnr = 0.9, ppv = 8 / 9, accuracy = 17 / 20)
  for (id in names(expected)) {
    row <- tab[tab$metric_id == id, , drop = FALSE]
    check(nrow(row) == 1L && row$status == "ok" && abs(row$value - expected[[id]]) <= 1e-14,
          paste("Count-based metric agrees with arithmetic:", id))
  }
  envelope <- analysis_envelope(classification_evaluated_spec(classified), classified)
  check(identical(ig_hash(envelope), ig_hash(ig_from_json(ig_to_json(envelope)))),
        "Classification envelope JSON preserves the hash")

  advisor <- advisor_session(language = "en")
  check(!is.null(advisor_question(advisor)), "The advisor can present its initial question")
  saved <- advisor_save(advisor)
  check(identical(saved, advisor_save(advisor_restore(saved))),
        "Advisor session save and restore preserves canonical JSON")
  NULL
}, error = function(e) conditionMessage(e))

writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
receipt <- list(
  scope = "installed core compatibility; not the full regression suite or numerical certification",
  started_utc = started,
  finished_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
  r_version = R.version.string,
  package_version = if (requireNamespace("inferguidR", quietly = TRUE)) as.character(utils::packageVersion("inferguidR")) else NULL,
  package_path = if (requireNamespace("inferguidR", quietly = TRUE)) find.package("inferguidR") else NULL,
  passed = length(checks),
  checks = as.list(checks),
  success = is.null(failure),
  failure = failure
)
jsonlite::write_json(receipt, file.path(out, "receipt.json"), pretty = TRUE,
                     auto_unbox = TRUE, null = "null")
if (!is.null(failure)) stop(failure, call. = FALSE)
cat(sprintf("Core compatibility: %d checks passed on %s\n", length(checks), R.version.string))
