# Power results ------------------------------------------------------------------
#
# A power_result holds the canonical input (the full spec), the typed status,
# solved and achieved quantities, group sizes, df/ncp/critical values,
# nominal and attained error rates, diagnostics, provenance, the validation
# state of the method and a reproducible R call. Guard rails enforced by the
# constructor:
#
# * a result whose status is not "ok" contains no finite number in its solved
#   values, achieved values, quantities, group sizes or method outputs;
# * an "ok" result has a finite solved value for every solve target;
# * every NA quantity has a reason in `undefined` (undefined is never 0).
#
# `result_hash` covers everything except itself and `runtime` (R/package
# versions), so the same computation hashes identically on any machine.

.ig_result_numeric_sections <- c("solved", "achieved", "quantities", "method_outputs")

# Internal constructor used by the engine dispatcher and the mode solvers.
.ig_power_result <- function(spec, status, solved = list(), achieved = list(),
                             group_sizes = NULL, quantities = list(), undefined = list(),
                             method_outputs = list(), diagnostics = list(),
                             provenance = list(), validation = NULL,
                             assumptions = character(), warnings = character(),
                             r_call = NULL, runtime = NULL) {
  if (!inherits(spec, "power_spec")) stop("`spec` must be a power_spec", call. = FALSE)
  if (!is_ig_status(status)) stop("`status` must be an ig_status", call. = FALSE)
  sections <- list(solved = solved, achieved = achieved, quantities = quantities,
                   undefined = undefined, method_outputs = method_outputs,
                   diagnostics = diagnostics, provenance = provenance)
  for (nm in names(sections)) {
    s <- sections[[nm]]
    if (is.null(s)) s <- list()
    if (!is.list(s) || (length(s) && !.ig_is_named(s))) {
      stop(sprintf("`%s` must be a named list", nm), call. = FALSE)
    }
    sections[[nm]] <- s
  }
  for (nm in c("solved", "achieved", "quantities")) {
    bad <- names(sections[[nm]])[!vapply(sections[[nm]], function(v) {
      is.atomic(v) && (is.numeric(v) || all(is.na(v)))
    }, logical(1))]
    if (length(bad)) stop(sprintf("`%s` values must be numeric: %s", nm, paste(bad, collapse = ", ")),
                          call. = FALSE)
  }
  if (!all(vapply(sections$undefined, .ig_is_string, logical(1)))) {
    stop("`undefined` must map quantity names to reasons (strings)", call. = FALSE)
  }
  if (!is.null(group_sizes)) {
    if (!is.numeric(group_sizes) || is.null(names(group_sizes)) || anyNA(names(group_sizes)) ||
        !all(nzchar(names(group_sizes)))) {
      stop("`group_sizes` must be a named numeric vector", call. = FALSE)
    }
    fin <- group_sizes[!is.na(group_sizes)]
    if (any(fin < 0 | fin != round(fin))) stop("`group_sizes` must be whole numbers >= 0", call. = FALSE)
    group_sizes <- stats::setNames(as.integer(round(group_sizes)), names(group_sizes))
  }
  if (!is.character(assumptions) || anyNA(assumptions)) stop("`assumptions` must be character", call. = FALSE)
  if (!is.character(warnings) || anyNA(warnings)) stop("`warnings` must be character", call. = FALSE)

  # Truthfulness guard: no number without an ok status.
  if (status$code != "ok") {
    leak <- c(.ig_finite_paths(sections$solved, "solved"),
              .ig_finite_paths(sections$achieved, "achieved"),
              .ig_finite_paths(sections$quantities, "quantities"),
              .ig_finite_paths(sections$method_outputs, "method_outputs"),
              if (!is.null(group_sizes) && any(!is.na(group_sizes))) "group_sizes")
    if (length(leak)) {
      stop(sprintf("a result with status '%s' must not report numbers (found: %s)",
                   status$code, paste(leak, collapse = ", ")), call. = FALSE)
    }
  } else {
    for (tg in spec$target) {
      v <- sections$solved[[tg]]
      if (!is.numeric(v) || length(v) != 1L || !is.finite(v)) {
        stop(sprintf("an ok result needs a finite solved value for target '%s'", tg), call. = FALSE)
      }
    }
  }
  # Undefined values need reasons.
  for (nm in c("solved", "achieved", "quantities")) {
    for (q in names(sections[[nm]])) {
      if (anyNA(sections[[nm]][[q]]) && status$code == "ok" && is.null(sections$undefined[[q]])) {
        stop(sprintf("undefined quantity '%s' needs a reason in `undefined`", q), call. = FALSE)
      }
    }
  }

  # Storage types: solve targets follow their parameter domain, other numbers
  # are doubles (so JSON round trips reproduce the object exactly).
  sections$solved <- .ig_coerce_solved(spec, sections$solved)
  for (nm in c("achieved", "quantities", "method_outputs", "diagnostics")) {
    sections[[nm]] <- .ig_doubles(sections[[nm]])
  }
  # Canonical member order (as in canonical JSON) for exact round trips.
  sections <- lapply(sections, function(s) .ig_named_or_empty(.ig_sort_names(s)))
  if (!is.null(group_sizes)) group_sizes <- group_sizes[order(names(group_sizes), method = "radix")]
  if (is.null(validation)) validation <- .ig_validation_record(spec$method_id)
  validation <- .ig_sort_names(validation)
  if (is.null(r_call)) r_call <- paste0("inferguidR::calculate_power(", .ig_spec_call(spec), ")")
  if (is.null(runtime)) runtime <- .ig_runtime()
  runtime <- .ig_sort_names(runtime)

  out <- list(
    schema_version = .ig_schema_versions[["power_result"]],
    object_type = "power_result",
    status = status,
    method_id = spec$method_id,
    mode = spec$mode,
    target = spec$target,
    spec = spec,
    spec_hash = spec$spec_hash,
    solved = sections$solved,
    achieved = sections$achieved,
    group_sizes = group_sizes,
    quantities = sections$quantities,
    undefined = sections$undefined,
    method_outputs = sections$method_outputs,
    diagnostics = sections$diagnostics,
    provenance = sections$provenance,
    validation = validation,
    assumptions = assumptions,
    warnings = warnings,
    r_call = r_call,
    runtime = runtime,
    result_hash = NA_character_
  )
  out$result_hash <- .ig_result_hash(out)
  structure(out, class = c("power_result", "ig_frozen"))
}

.ig_result_hash <- function(x) {
  lst <- .ig_result_to_list(x)
  lst$result_hash <- NULL
  lst$runtime <- NULL
  .ig_hash_value(lst)
}

.ig_runtime <- function() {
  list(package = "inferguidR",
       package_version = .ig_package_version(),
       r_version = paste(R.version$major, R.version$minor, sep = "."))
}

.ig_package_version <- function() {
  v <- tryCatch(as.character(utils::packageVersion("inferguidR")), error = function(e) NA_character_)
  if (is.na(v)) "unknown" else v
}

# Paths of finite numbers inside a (nested) list.
.ig_finite_paths <- function(x, path) {
  if (is.null(x)) return(character())
  if (is.list(x)) {
    nms <- names(x)
    out <- character()
    for (i in seq_along(x)) {
      sub <- if (is.null(nms)) paste0(path, "[", i, "]") else paste0(path, "$", nms[i])
      out <- c(out, .ig_finite_paths(x[[i]], sub))
    }
    return(out)
  }
  if (is.numeric(x) && any(is.finite(x))) return(path)
  character()
}

# Convert integers to doubles recursively (not inside group sizes).
.ig_doubles <- function(x) {
  if (is.list(x)) {
    nms <- names(x)
    out <- lapply(x, .ig_doubles)
    names(out) <- nms
    return(out)
  }
  if (is.integer(x) && !is.factor(x)) {
    nm <- names(x)
    x <- as.double(x)
    names(x) <- nm
  }
  x
}

.ig_coerce_solved <- function(spec, solved) {
  if (!length(solved)) return(solved)
  def <- .ig_get_method(spec$method_id)
  for (nm in names(solved)) {
    v <- solved[[nm]]
    integer_target <- if (!is.null(def) && !is.null(def$params[[nm]])) {
      identical(def$params[[nm]]$domain$type, "integer")
    } else {
      identical(nm, "n_total")
    }
    solved[[nm]] <- if (integer_target && all(is.na(v) | v == round(v))) as.integer(round(v)) else as.double(v)
  }
  solved
}

# Result for a valid request that has no implementation yet.
.ig_unsupported_result <- function(spec, reason, entry) {
  .ig_power_result(spec, ig_status("unsupported", reason), provenance = list(entry = entry))
}

# Rebuild a result with extra warnings (results are immutable).
.ig_result_add_warnings <- function(res, warnings) {
  .ig_power_result(res$spec, res$status, solved = res$solved, achieved = res$achieved,
                   group_sizes = res$group_sizes, quantities = res$quantities,
                   undefined = res$undefined, method_outputs = res$method_outputs,
                   diagnostics = res$diagnostics, provenance = res$provenance,
                   validation = res$validation, assumptions = res$assumptions,
                   warnings = unique(c(res$warnings, warnings)), r_call = res$r_call,
                   runtime = res$runtime)
}

# Serialisation -----------------------------------------------------------------------

.ig_result_to_list <- function(x) {
  list(
    schema_version = x$schema_version,
    object_type = x$object_type,
    status = .ig_status_to_list(x$status),
    method_id = x$method_id,
    mode = x$mode,
    target = I(x$target),
    spec = .ig_spec_to_list(x$spec),
    spec_hash = x$spec_hash,
    solved = .ig_named_or_empty(x$solved),
    achieved = .ig_named_or_empty(x$achieved),
    group_sizes = if (is.null(x$group_sizes)) NULL else as.list(x$group_sizes),
    quantities = .ig_named_or_empty(x$quantities),
    undefined = .ig_named_or_empty(x$undefined),
    method_outputs = .ig_named_or_empty(x$method_outputs),
    diagnostics = .ig_named_or_empty(x$diagnostics),
    provenance = .ig_named_or_empty(x$provenance),
    validation = x$validation,
    assumptions = I(x$assumptions),
    warnings = I(x$warnings),
    r_call = x$r_call,
    runtime = x$runtime,
    result_hash = x$result_hash
  )
}

# Turn parsed JSON arrays of scalars back into atomic vectors (null -> NA).
.ig_simplify_parsed <- function(x) {
  if (!is.list(x)) return(x)
  if (is.null(names(x)) && length(x) &&
      all(vapply(x, function(e) is.null(e) || (is.atomic(e) && length(e) == 1L), logical(1)))) {
    vals <- lapply(x, function(e) if (is.null(e)) NA else e)
    return(unlist(vals, use.names = FALSE))
  }
  nms <- names(x)
  out <- lapply(x, function(e) if (is.null(e)) NULL else .ig_simplify_parsed(e))
  names(out) <- nms
  # lapply keeps NULL entries; preserve them
  out
}

.ig_result_from_list <- function(x) {
  need <- c("schema_version", "object_type", "status", "spec", "spec_hash", "result_hash")
  miss <- setdiff(need, names(x))
  if (length(miss)) .ig_abort_invalid(sprintf("missing field '%s'", miss), "malformed power_result")
  spec <- .ig_spec_from_list(x$spec)
  if (!identical(spec$spec_hash, x$spec_hash)) {
    .ig_abort("invalid_input", "power_result spec_hash does not match its embedded specification")
  }
  status <- tryCatch(.ig_status_from_list(x$status), error = function(e) {
    .ig_abort("invalid_input", paste("invalid status record:", conditionMessage(e)))
  })
  num_section <- function(s) {
    s <- .ig_simplify_parsed(s %||% list())
    out <- lapply(s, function(v) if (is.null(v) || (is.logical(v) && all(is.na(v)))) NA_real_ else v)
    names(out) <- names(s)
    out
  }
  group_sizes <- NULL
  if (!is.null(x$group_sizes)) {
    gs <- x$group_sizes
    group_sizes <- vapply(gs, function(v) if (is.null(v)) NA_real_ else as.double(v), double(1))
    names(group_sizes) <- names(gs)
  }
  validation <- x$validation
  if (is.list(validation)) {
    for (nm in c("documented_example_ids", "fixture_ids")) {
      validation[[nm]] <- I(.ig_flatten(validation[[nm]], "character") %||% character())
    }
    if (is.null(validation$independent_checker)) validation$independent_checker <- NA_character_
  }
  res <- tryCatch(
    .ig_power_result(
      spec, status,
      solved = num_section(x$solved), achieved = num_section(x$achieved),
      group_sizes = group_sizes, quantities = num_section(x$quantities),
      undefined = x$undefined %||% list(),
      method_outputs = .ig_simplify_parsed(x$method_outputs %||% list()),
      diagnostics = .ig_simplify_parsed(x$diagnostics %||% list()),
      provenance = x$provenance %||% list(),
      validation = validation,
      assumptions = .ig_flatten(x$assumptions, "character") %||% character(),
      warnings = .ig_flatten(x$warnings, "character") %||% character(),
      r_call = x$r_call, runtime = x$runtime),
    error = function(e) .ig_abort("invalid_input", paste("invalid power_result:", conditionMessage(e)))
  )
  if (!identical(res$result_hash, x$result_hash)) {
    .ig_abort("invalid_input", "result_hash does not match the result content (the document was modified)")
  }
  res
}

# Display ----------------------------------------------------------------------------------

#' @export
format.power_result <- function(x, ...) {
  shown <- function(vals) {
    if (!length(vals)) return(NULL)
    paste(vapply(names(vals), function(nm) {
      v <- vals[[nm]]
      txt <- if (is.numeric(v)) paste(ifelse(is.na(v), "NA", format(v, digits = 7)), collapse = ", ") else
        paste(format(v), collapse = ", ")
      paste0(nm, " = ", txt)
    }, ""), collapse = "; ")
  }
  lines <- c(sprintf("<power_result> %s | %s", x$method_id, x$mode),
             sprintf("  status:      %s", format(x$status)))
  if (x$status$code == "ok") {
    lines <- c(lines,
               sprintf("  solved:      %s", shown(x$solved)),
               if (length(x$achieved)) sprintf("  achieved:    %s", shown(x$achieved)),
               if (!is.null(x$group_sizes)) sprintf("  group sizes: %s", shown(as.list(x$group_sizes))),
               if (length(x$quantities)) sprintf("  quantities:  %s", shown(x$quantities)),
               if (length(x$undefined)) sprintf("  undefined:   %s",
                                                paste(names(x$undefined), unlist(x$undefined), sep = ": ", collapse = "; ")))
  } else {
    lines <- c(lines, "  no numerical result")
  }
  c(lines,
    sprintf("  validation:  %s; G*Power reference: %s",
            x$validation$method_validation_status, x$validation$gpower_reference_status),
    if (length(x$warnings)) paste0("  warning: ", x$warnings),
    sprintf("  spec hash:   %s", x$spec_hash))
}

#' @export
print.power_result <- function(x, ...) {
  cat(format(x), sep = "\n")
  invisible(x)
}
