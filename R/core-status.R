
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
