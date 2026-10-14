# Typed statuses and typed conditions -------------------------------------
#
# Every calculation outcome in inferguidR carries one of six status codes.
# A status is a small classed list (code, reason, details). Errors that stop
# a constructor are signalled as conditions whose class names the status
# code, so callers can catch `inferguidR_invalid_input` and read the status
# object from the condition.

#' Status codes used by inferguidR
#'
#' The six typed outcome codes shared by specifications, results and
#' solvers:
#'
#' * `ok`: the requested quantity was computed and passed its checks.
#' * `infeasible`: the request is valid but has no solution (for example a
#'   target power that cannot be reached within the admissible designs).
#' * `unsupported`: the request is valid but the capability is not
#'   implemented (for example a method without an engine yet).
#' * `invalid_input`: the request violates a domain, mode or consistency
#'   rule.
#' * `numerical_failure`: a numerical routine failed or its checks did not
#'   pass; no number is reported.
#' * `reference_pending`: the computation depends on reference evidence that
#'   does not exist yet.
#'
#' @return A character vector of the six codes.
#' @examples
#' ig_status_codes()
#' @family status helpers
#' @export
ig_status_codes <- function() {
  c("ok", "infeasible", "unsupported", "invalid_input",
    "numerical_failure", "reference_pending")
}

#' Create a typed status
#'
#' @param code One of [ig_status_codes()].
#' @param reason A human-readable explanation. Required (non-empty) for every
#'   code except `"ok"`.
#' @param details A named list of structured details (for example the list of
#'   problems found during validation). Values must be serialisable data.
#' @return An object of class `ig_status`.
#' @examples
#' ig_status("unsupported", "engine for t_generic not implemented yet")
#' @family status helpers
#' @export
ig_status <- function(code, reason = "", details = list()) {
  if (!is.character(code) || length(code) != 1L || is.na(code) ||
      !(code %in% ig_status_codes())) {
    stop("`code` must be one of: ", paste(ig_status_codes(), collapse = ", "),
         call. = FALSE)
  }
  if (is.null(reason)) reason <- ""
  if (!is.character(reason) || length(reason) != 1L || is.na(reason)) {
    stop("`reason` must be a single string", call. = FALSE)
  }
  if (code != "ok" && !nzchar(trimws(reason))) {
    stop("a status other than 'ok' needs a non-empty `reason`", call. = FALSE)
  }
  if (is.null(details)) details <- list()
  if (!is.list(details) || (length(details) > 0L && !.ig_is_named(details))) {
    stop("`details` must be a named list", call. = FALSE)
  }
  structure(list(code = code, reason = reason,
                 details = .ig_named_or_empty(.ig_sort_names(details))),
            class = "ig_status")
}

#' Test for a status object
#'
#' @param x Any object.
#' @return `TRUE` if `x` is an `ig_status`.
#' @family status helpers
#' @export
is_ig_status <- function(x) inherits(x, "ig_status")

#' @export
format.ig_status <- function(x, ...) {
  if (identical(x$code, "ok") && !nzchar(x$reason)) return("ok")
  sprintf("%s: %s", x$code, x$reason)
}

#' @export
print.ig_status <- function(x, ...) {
  cat("<ig_status> ", format(x), "\n", sep = "")
  problems <- x$details$problems
  if (length(problems)) {
    cat(paste0("  - ", unlist(problems), collapse = "\n"), "\n", sep = "")
  }
  invisible(x)
}

# Rebuild a status from parsed JSON (a plain list).
.ig_status_from_list <- function(x) {
  if (!is.list(x) || is.null(x$code)) {
    stop("status record must be an object with a `code`", call. = FALSE)
  }
  ig_status(x$code, if (is.null(x$reason)) "" else x$reason,
            if (is.null(x$details)) list() else x$details)
}

.ig_status_to_list <- function(x) {
  list(code = x$code, reason = x$reason,
       details = .ig_named_or_empty(x$details))
}

# Typed conditions ---------------------------------------------------------

# Signal an error whose class names the status code. The status object is
# attached as `condition$status`.
.ig_abort <- function(code, reason, details = list(), call = NULL) {
  status <- ig_status(code, reason, details)
  msg <- reason
  problems <- details$problems
  if (length(problems) > 1L) {
    msg <- paste0(reason, "\n", paste0("* ", unlist(problems), collapse = "\n"))
  }
  cond <- structure(
    class = c(paste0("inferguidR_", code), "inferguidR_error", "error", "condition"),
    list(message = msg, call = call, status = status)
  )
  stop(cond)
}

# Raise invalid_input with a list of problems (character vector).
.ig_abort_invalid <- function(problems, context = "invalid input") {
  problems <- unique(as.character(problems))
  reason <- if (length(problems) == 1L) {
    paste0(context, ": ", problems)
  } else {
    sprintf("%s: %d problems found", context, length(problems))
  }
  .ig_abort("invalid_input", reason, list(problems = as.list(problems)))
}

# Small shared predicates ---------------------------------------------------

.ig_is_named <- function(x) {
  nms <- names(x)
  !is.null(nms) && !anyNA(nms) && all(nzchar(nms))
}

.ig_is_string <- function(x) {
  is.character(x) && length(x) == 1L && !is.na(x) && nzchar(x)
}

.ig_is_flag <- function(x) is.logical(x) && length(x) == 1L && !is.na(x)

# A named empty list serialises as `{}` instead of `[]`.
.ig_empty_object <- function() structure(list(), names = character(0))

.ig_named_or_empty <- function(x) {
  if (is.null(x) || length(x) == 0L) .ig_empty_object() else x
}

# Canonical R-side form of open (free-form) data sections: names sorted in
# C-locale order at every level. Canonical JSON sorts keys the same way, so
# objects rebuilt from JSON are identical to the originals. Empty lists are
# left as they are: a named empty list is a JSON object, an unnamed one an
# array (constructors use .ig_named_or_empty() for object-valued sections).
.ig_sort_names <- function(x) {
  if (!is.list(x) || (is.object(x) && !inherits(x, "AsIs"))) return(x)
  if (length(x) == 0L) return(x)
  nms <- names(x)
  if (!is.null(nms)) x <- x[order(nms, method = "radix")]
  for (i in seq_along(x)) {
    if (!is.null(x[[i]])) x[i] <- list(.ig_sort_names(x[[i]]))
  }
  x
}

`%||%` <- function(x, y) if (is.null(x)) y else x
