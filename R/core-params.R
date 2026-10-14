# Parameter and option specifications -----------------------------------------
#
# A method definition declares its inputs as parameters (numbers, integers,
# flags or labelled choices with a role, domain, unit and optional default)
# and its result-affecting choices as options (machine codes with optional
# default and applicability). Defaults are never hidden: analysis_spec()
# writes every applied default into the specification.

.ig_param_roles <- c("effect", "nuisance", "design", "error", "control")

.ig_name_pattern <- "^[a-z][a-z0-9_]*$"
.ig_code_pattern <- "^[a-z0-9][a-z0-9_]*$"

# Domains --------------------------------------------------------------------

#' Parameter domains
#'
#' Domain descriptors used by [ig_param()]. A domain checks the type, the
#' length and the admissible range of an input value.
#'
#' @param lower,upper Bounds of the admissible range. Use `-Inf`/`Inf` for an
#'   unbounded side.
#' @param lower_inclusive,upper_inclusive Whether the bound itself is
#'   admissible.
#' @param length Required vector length; `NA` accepts any length of at least
#'   one.
#' @param values Allowed character codes (for `ig_domain_choice()`).
#' @param scan Optional two positive finite distances, in increasing order,
#'   for the numerical effect scan of an unbounded real domain. Distances
#'   are measured from its finite bound, or from zero when both sides are
#'   unbounded. This controls the search, not admissible input values.
#' @return An object of class `ig_domain`.
#' @examples
#' ig_domain_real(0, 1, lower_inclusive = FALSE, upper_inclusive = FALSE)
#' ig_domain_integer(lower = 2)
#' @family method definition API
#' @name ig_domain
NULL

#' @rdname ig_domain
#' @export
ig_domain_real <- function(lower = -Inf, upper = Inf, lower_inclusive = TRUE,
                           upper_inclusive = TRUE, length = 1L, scan = NULL) {
  d <- .ig_new_domain("real", lower, upper, lower_inclusive, upper_inclusive, length)
  if (!is.null(scan)) {
    if (!is.numeric(scan) || length(scan) != 2L || any(!is.finite(scan)) ||
        any(scan <= 0) || scan[1] >= scan[2]) {
      stop("`scan` must contain two positive finite distances in increasing order", call. = FALSE)
    }
    if (is.finite(lower) && is.finite(upper)) {
      stop("`scan` is only supported for an unbounded real domain", call. = FALSE)
    }
    d$scan <- as.double(scan)
  }
  d
}

#' @rdname ig_domain
#' @export
ig_domain_integer <- function(lower = -Inf, upper = Inf, length = 1L) {
  .ig_new_domain("integer", lower, upper, TRUE, TRUE, length)
}

#' @rdname ig_domain
#' @export
ig_domain_logical <- function(length = 1L) {
  .ig_new_domain("logical", -Inf, Inf, TRUE, TRUE, length)
}

#' @rdname ig_domain
#' @export
ig_domain_choice <- function(values, length = 1L) {
  if (!is.character(values) || !length(values) || anyNA(values) ||
      anyDuplicated(values) || !all(grepl(.ig_code_pattern, values))) {
    stop("`values` must be unique lower-case codes (letters, digits, underscores)",
         call. = FALSE)
  }
  d <- .ig_new_domain("choice", -Inf, Inf, TRUE, TRUE, length)
  d$values <- values
  d
}

.ig_new_domain <- function(type, lower, upper, lower_inclusive, upper_inclusive, length) {
  if (!is.numeric(lower) || length(lower) != 1L || is.na(lower) ||
      !is.numeric(upper) || length(upper) != 1L || is.na(upper) || lower > upper) {
    stop("domain bounds must be single numbers with lower <= upper", call. = FALSE)
  }
  if (!.ig_is_flag(lower_inclusive) || !.ig_is_flag(upper_inclusive)) {
    stop("inclusivity flags must be TRUE or FALSE", call. = FALSE)
  }
  if (length(length) != 1L || (!is.na(length) && (length < 1 || length != round(length)))) {
    stop("`length` must be a positive whole number or NA", call. = FALSE)
  }
  structure(list(type = type, lower = as.double(lower), upper = as.double(upper),
                 lower_inclusive = lower_inclusive, upper_inclusive = upper_inclusive,
                 length = as.integer(length), values = NULL),
            class = "ig_domain")
}

# Check a value against a domain. Returns NULL when admissible, otherwise a
# short reason. `name` is used in the message.
.ig_domain_check <- function(domain, value, name) {
  if (is.null(value)) return(sprintf("'%s' is missing", name))
  n <- length(value)
  if (n == 0L) return(sprintf("'%s' is empty", name))
  if (!is.na(domain$length) && n != domain$length) {
    return(sprintf("'%s' must have length %d (got %d)", name, domain$length, n))
  }
  if (is.list(value) || !is.atomic(value)) {
    return(sprintf("'%s' must be an atomic value", name))
  }
  type <- domain$type
  if (type %in% c("real", "integer")) {
    if (!is.numeric(value)) return(sprintf("'%s' must be numeric", name))
    if (anyNA(value) || any(!is.finite(value))) {
      return(sprintf("'%s' must be finite (no NA, NaN or infinite values)", name))
    }
    if (type == "integer" && any(value != round(value))) {
      return(sprintf("'%s' must be a whole number", name))
    }
    if (type == "integer" && any(abs(value) > .Machine$integer.max)) {
      return(sprintf("'%s' exceeds the supported integer range", name))
    }
    lo_ok <- if (domain$lower_inclusive) value >= domain$lower else value > domain$lower
    hi_ok <- if (domain$upper_inclusive) value <= domain$upper else value < domain$upper
    if (!all(lo_ok & hi_ok)) {
      return(sprintf("'%s' must lie in %s (got %s)", name, .ig_format_domain(domain),
                     paste(.ig_format_double(as.double(value)), collapse = ", ")))
    }
    return(NULL)
  }
  if (type == "logical") {
    if (!is.logical(value) || anyNA(value)) return(sprintf("'%s' must be TRUE or FALSE", name))
    return(NULL)
  }
  if (type == "choice") {
    if (!is.character(value) || anyNA(value) || !all(value %in% domain$values)) {
      return(sprintf("'%s' must be one of: %s", name, paste(domain$values, collapse = ", ")))
    }
    return(NULL)
  }
  sprintf("'%s' has an unknown domain type '%s'", name, type)
}

# Coerce an admissible value to the domain's storage type.
.ig_domain_coerce <- function(domain, value) {
  switch(domain$type,
         real = as.double(value),
         integer = as.integer(round(value)),
         logical = as.logical(value),
         choice = as.character(value),
         value)
}

.ig_format_domain <- function(domain) {
  if (domain$type == "choice") return(paste0("{", paste(domain$values, collapse = ", "), "}"))
  lo <- if (is.infinite(domain$lower)) "-Inf" else .ig_format_double(domain$lower)
  hi <- if (is.infinite(domain$upper)) "Inf" else .ig_format_double(domain$upper)
  paste0(if (domain$lower_inclusive && is.finite(domain$lower)) "[" else "(",
         lo, ", ", hi,
         if (domain$upper_inclusive && is.finite(domain$upper)) "]" else ")")
}

#' @export
format.ig_domain <- function(x, ...) {
  len <- if (is.na(x$length)) "length >= 1" else paste("length", x$length)
  sprintf("%s %s, %s", x$type, .ig_format_domain(x), len)
}

#' @export
print.ig_domain <- function(x, ...) {
  cat("<ig_domain> ", format(x), "\n", sep = "")
  invisible(x)
}

# Parameters -----------------------------------------------------------------

#' Declare a method parameter
#'
#' Parameters are the numeric (or flag/choice) inputs of a method definition.
#' Each has a role, a domain, a unit and an optional default. A default is
#' written into every specification that does not supply the parameter, so
#' the serialised specification always shows the value that was used.
#'
#' @param name Parameter name (lower case, digits and underscores).
#' @param role One of `"effect"`, `"nuisance"`, `"design"`, `"error"` or
#'   `"control"` (computational settings such as an exact-computation
#'   threshold or solver bounds).
#' @param domain An [ig_domain] object.
#' @param unit Unit label (for example `"subjects"`); `NA` when unitless.
#' @param default Default value (validated against `domain`), or `NULL` when
#'   the parameter must be supplied.
#' @param label Display label, a list with elements `en` and optionally `de`.
#' @param description Optional longer description (English).
#' @param aliases Alternative input names accepted by [analysis_spec()].
#' @param modes Modes in which the parameter is an input; `NULL` means every
#'   mode in which it is not the solve target.
#' @param applies_when Optional applicability condition (see Details).
#' @details Applicability conditions use a small data grammar shared with the
#'   registry: `list(mode = c(...))`, `list(option = "<id>", "in" = c(...))`,
#'   `list(all = list(...))` and `list(any = list(...))`. A condition that
#'   refers to an option which is itself not applicable is false.
#' @return An object of class `ig_param`.
#' @examples
#' ig_param("d", role = "effect", domain = ig_domain_real(),
#'          label = list(en = "Effect size d"))
#' @family method definition API
#' @export
ig_param <- function(name, role, domain, unit = NA_character_, default = NULL,
                     label = NULL, description = NULL, aliases = character(),
                     modes = NULL, applies_when = NULL) {
  if (!.ig_is_string(name) || !grepl(.ig_name_pattern, name)) {
    stop("parameter `name` must match ", .ig_name_pattern, call. = FALSE)
  }
  if (!.ig_is_string(role) || !(role %in% .ig_param_roles)) {
    stop("parameter `role` must be one of: ", paste(.ig_param_roles, collapse = ", "),
         call. = FALSE)
  }
  if (!inherits(domain, "ig_domain")) stop("`domain` must be an ig_domain", call. = FALSE)
  if (length(unit) != 1L || !(is.na(unit) || is.character(unit))) {
    stop("`unit` must be a single string or NA", call. = FALSE)
  }
  if (!is.null(default)) {
    problem <- .ig_domain_check(domain, default, name)
    if (!is.null(problem)) stop("invalid default: ", problem, call. = FALSE)
    default <- .ig_domain_coerce(domain, default)
  }
  if (is.null(label)) label <- list(en = name)
  label <- .ig_check_label(label, "label")
  if (!is.null(description) && !.ig_is_string(description)) {
    stop("`description` must be a single string", call. = FALSE)
  }
  if (!is.character(aliases) || anyNA(aliases) ||
      !all(grepl(.ig_name_pattern, aliases)) || name %in% aliases) {
    stop("`aliases` must be valid names different from `name`", call. = FALSE)
  }
  if (!is.null(modes)) {
    if (!is.character(modes) || !length(modes) || !all(modes %in% .ig_modes)) {
      stop("`modes` must be a subset of: ", paste(.ig_modes, collapse = ", "), call. = FALSE)
    }
  }
  if (!is.null(applies_when)) .ig_check_condition(applies_when)
  structure(list(name = name, role = role, domain = domain,
                 unit = as.character(unit), default = default, label = label,
                 description = description, aliases = aliases, modes = modes,
                 applies_when = applies_when),
            class = "ig_param")
}

#' Collect method parameters
#'
#' @param ... [ig_param()] objects (or lists of them).
#' @return An object of class `ig_params`: a list named by parameter name.
#' @examples
#' ig_params(ig_param_alpha(), ig_param_power(), ig_param_n_total())
#' @family method definition API
#' @export
ig_params <- function(...) {
  items <- .ig_flatten_items(list(...), "ig_param")
  nms <- vapply(items, `[[`, "", "name")
  if (anyDuplicated(nms)) {
    stop("duplicated parameter names: ", paste(unique(nms[duplicated(nms)]), collapse = ", "),
         call. = FALSE)
  }
  all_names <- c(nms, unlist(lapply(items, `[[`, "aliases")))
  if (anyDuplicated(all_names)) {
    stop("parameter aliases collide with other names: ",
         paste(unique(all_names[duplicated(all_names)]), collapse = ", "), call. = FALSE)
  }
  names(items) <- nms
  structure(items, class = "ig_params")
}

#' Standard error-rate and sample-size parameters
#'
#' Convenience constructors for the inputs shared by all power analyses:
#' significance level `alpha`, power `power` (1 - beta), the compromise
#' ratio `error_ratio` (q = beta/alpha) and the total sample size `n_total`.
#' None has a default: they must be supplied whenever the mode needs them.
#'
#' @param lower Smallest admissible total sample size.
#' @return An [ig_param()] object.
#' @examples
#' ig_param_alpha()
#' ig_param_n_total(lower = 3)
#' @family method definition API
#' @name ig_standard_params
NULL

#' @rdname ig_standard_params
#' @export
ig_param_alpha <- function() {
  ig_param("alpha", role = "error",
           domain = ig_domain_real(0, 1, lower_inclusive = FALSE, upper_inclusive = FALSE),
           label = list(en = "\u03b1 err prob", de = "\u03b1-Fehlerwahrscheinlichkeit"))
}

#' @rdname ig_standard_params
#' @export
ig_param_power <- function() {
  ig_param("power", role = "error",
           domain = ig_domain_real(0, 1, lower_inclusive = FALSE, upper_inclusive = FALSE),
           label = list(en = "Power (1-\u03b2 err prob)", de = "Testst\u00e4rke (1-\u03b2)"))
}

#' @rdname ig_standard_params
#' @export
ig_param_error_ratio <- function() {
  ig_param("error_ratio", role = "error",
           domain = ig_domain_real(0, Inf, lower_inclusive = FALSE, upper_inclusive = FALSE),
           label = list(en = "\u03b2/\u03b1 ratio", de = "\u03b2/\u03b1-Verh\u00e4ltnis"),
           modes = "compromise")
}

#' @rdname ig_standard_params
#' @export
ig_param_n_total <- function(lower = 2) {
  ig_param("n_total", role = "design", domain = ig_domain_integer(lower = lower),
           unit = "subjects",
           label = list(en = "Total sample size", de = "Gesamtstichprobenumfang"))
}

# Options --------------------------------------------------------------------

#' Declare a result-affecting option
#'
#' Options are discrete, result-affecting choices identified by machine codes
#' (for example `tails = "two"`). Presentation settings are not options and
#' do not belong in a specification. An option without a default must be
#' supplied whenever it applies.
#'
#' @param name Option name.
#' @param values Allowed codes.
#' @param default Default code, or `NULL` when the value must be supplied.
#' @param labels Optional display labels, named by code.
#' @param applies_when Optional applicability condition (see [ig_param()]).
#' @param label Display label, a list with `en` and optionally `de`.
#' @return An object of class `ig_option`.
#' @examples
#' ig_option("tails", values = c("one", "two"))
#' @family method definition API
#' @export
ig_option <- function(name, values, default = NULL, labels = NULL,
                      applies_when = NULL, label = NULL) {
  if (!.ig_is_string(name) || !grepl(.ig_name_pattern, name)) {
    stop("option `name` must match ", .ig_name_pattern, call. = FALSE)
  }
  if (!is.character(values) || !length(values) || anyNA(values) ||
      anyDuplicated(values) || !all(grepl(.ig_code_pattern, values))) {
    stop("option `values` must be unique lower-case codes", call. = FALSE)
  }
  if (!is.null(default) && (!.ig_is_string(default) || !(default %in% values))) {
    stop("option `default` must be one of the values", call. = FALSE)
  }
  if (!is.null(labels)) {
    if (!is.character(labels) || is.null(names(labels)) || !all(names(labels) %in% values)) {
      stop("option `labels` must be a character vector named by value codes", call. = FALSE)
    }
  }
  if (!is.null(applies_when)) .ig_check_condition(applies_when)
  if (is.null(label)) label <- list(en = name)
  structure(list(name = name, values = values, default = default, labels = labels,
                 applies_when = applies_when, label = .ig_check_label(label, "label")),
            class = "ig_option")
}

#' Collect method options
#'
#' @param ... [ig_option()] objects (or lists of them). May be empty.
#' @return An object of class `ig_options`: a list named by option name.
#' @examples
#' ig_options(ig_option("tails", c("one", "two")))
#' @family method definition API
#' @export
ig_options <- function(...) {
  items <- .ig_flatten_items(list(...), "ig_option")
  nms <- vapply(items, `[[`, "", "name")
  if (anyDuplicated(nms)) {
    stop("duplicated option names: ", paste(unique(nms[duplicated(nms)]), collapse = ", "),
         call. = FALSE)
  }
  names(items) <- nms
  structure(items, class = "ig_options")
}

.ig_flatten_items <- function(items, cls) {
  out <- list()
  for (it in items) {
    if (inherits(it, cls)) {
      out[[length(out) + 1L]] <- it
    } else if (is.list(it) && all(vapply(it, inherits, logical(1), what = cls))) {
      out <- c(out, unclass(it))
    } else {
      stop("expected ", cls, " objects", call. = FALSE)
    }
  }
  unname(out)
}

.ig_check_label <- function(label, what) {
  if (.ig_is_string(label)) label <- list(en = label)
  if (!is.list(label) || !.ig_is_string(label$en) ||
      !all(names(label) %in% c("en", "de"))) {
    stop(sprintf("`%s` must be a list with a string `en` (and optional `de`)", what),
         call. = FALSE)
  }
  if (!is.null(label$de) && !(length(label$de) == 1L && (is.na(label$de) || is.character(label$de)))) {
    stop(sprintf("`%s$de` must be a single string", what), call. = FALSE)
  }
  label
}

# Applicability conditions -----------------------------------------------------

.ig_check_condition <- function(cond) {
  if (!is.list(cond) || length(cond) == 0L || is.null(names(cond))) {
    stop("an applicability condition must be a named list", call. = FALSE)
  }
  keys <- names(cond)
  if (identical(keys, "mode")) {
    modes <- .ig_flatten(cond$mode, "character")
    if (!length(modes) || !all(modes %in% .ig_modes)) stop("invalid modes in condition", call. = FALSE)
  } else if (setequal(keys, c("option", "in"))) {
    if (!.ig_is_string(cond$option)) stop("condition `option` must be a string", call. = FALSE)
    codes <- .ig_flatten(cond[["in"]], "character")
    if (!length(codes)) stop("condition `in` must list codes", call. = FALSE)
  } else if (identical(keys, "all") || identical(keys, "any")) {
    parts <- cond[[keys]]
    if (!is.list(parts) || !length(parts)) stop("`all`/`any` need a list of conditions", call. = FALSE)
    for (p in parts) .ig_check_condition(p)
  } else {
    stop("unknown condition form: ", paste(keys, collapse = ", "), call. = FALSE)
  }
  invisible(TRUE)
}

# Evaluate a condition. `options` holds the effective values of the options
# that are applicable (inapplicable options are absent).
.ig_eval_condition <- function(cond, mode, options) {
  if (is.null(cond)) return(TRUE)
  keys <- names(cond)
  if (identical(keys, "mode")) return(mode %in% .ig_flatten(cond$mode, "character"))
  if (setequal(keys, c("option", "in"))) {
    val <- options[[cond$option]]
    return(!is.null(val) && all(val %in% .ig_flatten(cond[["in"]], "character")))
  }
  if (identical(keys, "all")) {
    return(all(vapply(cond$all, .ig_eval_condition, logical(1), mode = mode, options = options)))
  }
  if (identical(keys, "any")) {
    return(any(vapply(cond$any, .ig_eval_condition, logical(1), mode = mode, options = options)))
  }
  FALSE
}

#' @export
format.ig_param <- function(x, ...) {
  def <- if (is.null(x$default)) "no default" else
    paste("default", paste(format(x$default), collapse = ", "))
  sprintf("%s [%s] %s; %s", x$name, x$role, format(x$domain), def)
}

#' @export
print.ig_param <- function(x, ...) {
  cat("<ig_param> ", format(x), "\n", sep = "")
  invisible(x)
}

#' @export
print.ig_params <- function(x, ...) {
  cat("<ig_params> ", length(x), " parameter(s)\n", sep = "")
  for (p in x) cat("  ", format(p), "\n", sep = "")
  invisible(x)
}

#' @export
print.ig_options <- function(x, ...) {
  cat("<ig_options> ", length(x), " option(s)\n", sep = "")
  for (o in x) {
    def <- if (is.null(o$default)) "no default" else paste("default", o$default)
    cat("  ", o$name, ": ", paste(o$values, collapse = " | "), " (", def, ")\n", sep = "")
  }
  invisible(x)
}
