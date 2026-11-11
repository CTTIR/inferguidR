library(inferguidR)
# Benchmarks of the base power engines against values captured from the
# established R packages (inst/validation/engine-benchmarks.rds). Inputs come
# from the originals' help-page examples and from G*Power manual examples;
# every record names its source. Agreement is a relative difference of at most
# 1e-8; other comparisons need an approved row in
# inst/validation/engine-deviations.csv. When an original is installed, its
# current output is also compared with the stored values.
records <- readRDS(system.file("validation", "engine-benchmarks.rds", package = "inferguidR", mustWork = TRUE))
deviations <- utils::read.csv(system.file("validation", "engine-deviations.csv", package = "inferguidR",
  mustWork = TRUE), colClasses = "character")
stopifnot(identical(names(deviations), c("method_id", "case_id", "field", "kind", "rule", "explanation",
  "approved_by", "approved_on")), length(records) > 0L)
full <- identical(Sys.getenv("INFER_SUITE_FULL"), "true")
registry <- inferguidR:::.ig_engine_info

agree <- function(a, b) {
  if (is.null(a) || is.null(b)) return(is.null(a) && is.null(b))
  if (is.character(a) || is.character(b)) return(identical(as.character(a), as.character(b)))
  a <- as.double(unlist(a)); b <- as.double(unlist(b))
  if (length(a) != length(b)) return(FALSE)
  all(mapply(function(a, b) if (is.na(b)) is.na(a) else if (b == 0) isTRUE(a == 0) else
    isTRUE(abs(a - b) / abs(b) <= 1e-8), a, b))
}
run_engine <- function(record, engine) {
  spec <- analysis_spec(record$method_id, record$mode, inputs = record$inputs, options = record$options,
    engine = engine)
  if (identical(record$mode, "a_priori")) solve_design(spec) else calculate_power(spec)
}
values_of <- function(result, fields) {
  stats::setNames(lapply(fields, function(field) {
    if (identical(field, "status")) result$status$code else if (identical(field, "power")) result$achieved$power else
      if (!is.null(result$solved[[field]])) result$solved[[field]] else result$quantities[[field]]
  }), fields)
}

for (record in records) {
  label <- paste(record$method_id, record$case_id)
  base <- values_of(run_engine(record, "base"), record$fields)
  for (field in record$fields) {
    rule <- deviations[deviations$method_id == record$method_id & deviations$field == field &
      deviations$case_id %in% c(record$case_id, "*"), , drop = FALSE]
    if (nrow(rule)) {
      if (!nzchar(rule$approved_by[[1L]])) stop("deviation awaiting owner approval: ", label, " ", field)
      # Rules are R expressions from the package's own owner-approved CSV, not user input.
      ok <- eval(parse(text = rule$rule[[1L]]), c(list(a = base[[field]], b = record$values[[field]]),
        record$inputs, record$options))
    } else ok <- agree(base[[field]], record$values[[field]])
    if (!isTRUE(ok)) stop("base engine differs from ", record$package, " ", record$version, ": ", label, " ", field)
  }
  entry <- registry[[record$method_id]]
  live <- !is.null(entry) && requireNamespace(record$package, quietly = TRUE)
  if (live) {
    # A changed original is reported, and fails only the full lane.
    live_run <- tryCatch(run_engine(record, "reference"), error = function(e) e)
    if (inherits(live_run, "error") || !identical(live_run$status$code, "ok")) {
      note <- paste0("live reference failed with ", record$package, " ", utils::packageVersion(record$package), ": ",
        label, ": ", if (inherits(live_run, "error")) conditionMessage(live_run) else live_run$status$reason)
      if (full) stop(note) else cat("NOTE", note, "\n")
    } else {
      reference <- values_of(live_run, record$fields)
      same <- vapply(record$fields, function(field) agree(reference[[field]], record$values[[field]]), NA)
      if (!all(same)) {
        note <- paste0("live reference differs from stored ", record$package, " ", record$version, ": ", label, " ",
          paste(record$fields[!same], collapse = ", "))
        if (full) stop(note) else cat("NOTE", note, "\n")
      }
    }
  } else {
    if (full) stop("full lane requires the original for ", label)
    cat("SKIP live reference:", record$package, "not installed for", label, "\n")
  }
  cat("PASS", label, "\n")
}
