#!/usr/bin/env Rscript
# Capture power benchmark values from the original packages (development only).
#
# Usage, from the package root with the development version installed:
#   Rscript tools/engines/capture-power-benchmarks.R <method_id> ... [--replace <method_id>]
#
# Each case runs the reference engine, which computes the procedure with the
# original package, and stores the compared result fields with the package
# version and call. A stored record is never changed silently: if a new capture
# differs, the script stops unless --replace names that method.
library(inferguidR)
source(file.path("tools", "engines", "power-benchmark-cases.R"))
args <- commandArgs(TRUE)
replace <- args[which(args == "--replace") + 1L]
methods <- unique(c(setdiff(args, c("--replace", replace)), replace))
methods <- methods[!is.na(methods)]
if (!length(methods)) stop("name at least one method_id")
unknown <- setdiff(methods, names(power_benchmark_cases))
if (length(unknown)) stop("no benchmark cases for: ", paste(unknown, collapse = ", "))

store <- file.path("inst", "validation", "engine-benchmarks.rds")
records <- if (file.exists(store)) readRDS(store) else list()
key <- function(r) paste(r$method_id, r$case_id)
values_of <- function(result, fields) {
  stats::setNames(lapply(fields, function(field) {
    if (identical(field, "status")) result$status$code else if (identical(field, "power")) result$achieved$power else
      if (!is.null(result$solved[[field]])) result$solved[[field]] else result$quantities[[field]]
  }), fields)
}

for (method_id in methods) {
  for (item in power_benchmark_cases[[method_id]]()) {
    spec <- analysis_spec(method_id, item$mode, inputs = item$inputs, options = item$options, engine = "reference")
    result <- if (identical(item$mode, "a_priori")) solve_design(spec) else calculate_power(spec)
    if (!identical(result$status$code, "ok"))
      stop("reference engine did not compute ", method_id, " ", item$case_id, ": ", result$status$reason)
    engine <- result$provenance$engine
    record <- list(method_id = method_id, case_id = item$case_id, source = item$source, mode = item$mode,
      inputs = item$inputs, options = item$options, package = engine$package,
      version = engine$version, call = engine$fun, fields = item$fields, values = values_of(result, item$fields),
      captured_on = format(Sys.Date()), r_version = R.version.string)
    old <- Filter(function(r) identical(key(r), key(record)), records)
    if (length(old) && !isTRUE(all.equal(old[[1L]]$values, record$values, tolerance = 0)) && !method_id %in% replace)
      stop("captured values differ from the stored record for ", key(record), "; rerun with --replace ", method_id)
    records <- c(Filter(function(r) !identical(key(r), key(record)), records), list(record))
    cat("captured", key(record), "with", engine$package, engine$version, "\n")
  }
}
records <- records[order(vapply(records, key, ""))]
saveRDS(records, store, version = 2)
cat("stored", length(records), "records in", store, "\n")
