# Configuration ------------------------------------------------------------
#
# inferguidR reads a small set of R options. Defaults are documented here and
# in ?inferguidR-package; .onLoad() never overwrites a value the user set.

# Known options and their defaults. A NULL default means "unset".
.ig_config_defaults <- function() {
  list(
    # Directory of the local G*Power reference archive (manual, binaries,
    # capture exports). Development-only; never required at run time.
    inferguidR.reference_archive_dir = NULL,
    # Validate imported JSON against the bundled JSON Schemas when the
    # 'jsonvalidate' package is available.
    inferguidR.validate_json_schema = TRUE
  )
}

#' Current inferguidR configuration
#'
#' Lists the R options that inferguidR reads, with their current values.
#'
#' * `inferguidR.reference_archive_dir`: path of a local reference archive
#'   (development use; default unset, i.e. `NULL`). It is never needed for
#'   calculations.
#' * `inferguidR.validate_json_schema`: if `TRUE` (default), [ig_from_json()]
#'   validates documents against the bundled JSON Schemas whenever the
#'   'jsonvalidate' package is installed.
#'
#' @return A named list of option values (`NULL` for unset options).
#' @examples
#' ig_config()
#' @export
ig_config <- function() {
  defaults <- .ig_config_defaults()
  out <- lapply(names(defaults), function(nm) getOption(nm, defaults[[nm]]))
  names(out) <- names(defaults)
  out
}

.ig_get_config <- function(name) {
  full <- paste0("inferguidR.", name)
  defaults <- .ig_config_defaults()
  if (!(full %in% names(defaults))) {
    stop("unknown inferguidR option: ", full, call. = FALSE)
  }
  getOption(full, defaults[[full]])
}

# Set defaults for options the user has not set. NULL defaults stay unset.
.ig_config_init <- function() {
  defaults <- .ig_config_defaults()
  current <- options()
  todo <- defaults[!(names(defaults) %in% names(current))]
  todo <- todo[!vapply(todo, is.null, logical(1))]
  if (length(todo)) options(todo)
  invisible(NULL)
}

# Path of the reference archive, or NULL if unset. A set but missing
# directory is an error (typed), because silently ignoring it would hide a
# misconfiguration.
.ig_reference_archive_dir <- function() {
  dir <- .ig_get_config("reference_archive_dir")
  if (is.null(dir)) return(NULL)
  if (!.ig_is_string(dir)) {
    .ig_abort("invalid_input",
              "option inferguidR.reference_archive_dir must be a single path string")
  }
  if (!dir.exists(dir)) {
    .ig_abort("invalid_input",
              sprintf("reference archive directory does not exist: %s", dir))
  }
  normalizePath(dir, winslash = "/", mustWork = TRUE)
}
