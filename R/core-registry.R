# Method registry ------------------------------------------------------------------
#
# Two layers:
#
# 1. Metadata: inst/registry/<family>.json, generated from the reference
#    inventory by project/tools/build_registry.R and loaded at .onLoad().
#    Every inventory procedure has an entry, implemented or not.
# 2. Definitions: ig_method() objects returned by provider functions named
#    `.ig_define_methods_<family>()` anywhere in the package namespace. They
#    are collected at load time and attached to their metadata entry. A
#    definition must use an inventory id, the same family, and only modes the
#    inventory lists.
#
# Mode solvers (phase 05) are collected from an optional hook
# `.ig_mode_solvers()` returning a list named by mode; each solver is a
# function(method, spec, entry) returning a power_result.

# Package state (created fresh at every load).
.ig_state <- new.env(parent = emptyenv())

.ig_registry_schema_version <- "inferguidR.registry/0.1.0"

.ig_registry_init <- function(ns = NULL) {
  .ig_state$metadata <- list()
  .ig_state$definitions <- list()
  .ig_state$solvers <- list()
  .ig_state$problems <- character()
  .ig_state$load_error <- NULL
  .ig_state$schema_validators <- new.env(parent = emptyenv())
  ok <- tryCatch({
    .ig_state$metadata <- .ig_load_registry_files(.ig_registry_dir())
    TRUE
  }, error = function(e) {
    .ig_state$load_error <- conditionMessage(e)
    FALSE
  })
  if (ok && !is.null(ns)) {
    .ig_collect_definitions(ns)
    .ig_collect_solvers(ns)
  }
  invisible(ok)
}

.ig_registry_dir <- function() {
  system.file("registry", package = "inferguidR")
}

# Read and normalise all registry files of a directory.
.ig_load_registry_files <- function(dir) {
  if (!nzchar(dir) || !dir.exists(dir)) {
    stop("the method registry directory (inst/registry) is missing from the installation",
         call. = FALSE)
  }
  files <- sort(list.files(dir, pattern = "\\.json$", full.names = TRUE))
  if (!length(files)) stop("no registry files found in ", dir, call. = FALSE)
  out <- list()
  for (f in files) {
    txt <- readChar(f, file.size(f), useBytes = TRUE)
    Encoding(txt) <- "UTF-8"
    doc <- tryCatch(.ig_json_parse(txt), error = function(e) {
      stop("cannot read registry file ", basename(f), ": ", conditionMessage(e), call. = FALSE)
    })
    if (!identical(doc$schema_version, .ig_registry_schema_version)) {
      stop("registry file ", basename(f), " has schema_version '",
           paste(format(doc$schema_version), collapse = ""), "', expected '",
           .ig_registry_schema_version, "'", call. = FALSE)
    }
    if (!.ig_is_string(doc$family) || !is.list(doc$methods)) {
      stop("registry file ", basename(f), " lacks `family` or `methods`", call. = FALSE)
    }
    for (m in doc$methods) {
      entry <- .ig_normalise_metadata(m, doc$family, basename(f))
      if (!is.null(out[[entry$method_id]])) {
        stop("method id '", entry$method_id, "' appears in more than one registry entry",
             call. = FALSE)
      }
      out[[entry$method_id]] <- entry
    }
  }
  idx <- vapply(out, function(e) if (is.null(e$inventory_index)) NA_real_ else
    as.double(e$inventory_index), double(1))
  out[order(idx, names(out), na.last = TRUE, method = "radix")]
}

# Turn one parsed registry entry into the R-side structure: arrays become
# atomic vectors, options become indexed lists.
.ig_normalise_metadata <- function(m, family, file) {
  id <- m$method_id
  if (!.ig_is_string(id)) stop("registry file ", file, " has an entry without method_id", call. = FALSE)
  if (!identical(m$family, family)) {
    stop("registry entry ", id, " has family '", format(m$family), "' in file for '",
         family, "'", call. = FALSE)
  }
  modes <- .ig_flatten(m$modes, "character")
  if (!length(modes) || !all(modes %in% .ig_modes)) {
    stop("registry entry ", id, " has invalid modes", call. = FALSE)
  }
  options <- lapply(m$options %||% list(), function(o) {
    values <- o$values %||% list()
    codes <- vapply(values, function(v) v$code, "")
    list(
      option_id = o$option_id,
      ui_label = o$ui_label,
      kind = o$kind,
      value_type = o$value_type,
      codes = codes,
      value_labels = vapply(values, function(v) v$label, ""),
      # Live control of each value (certified inventory), named by code.
      value_live = stats::setNames(lapply(values, function(v) v$live), codes),
      value_applies_when = stats::setNames(lapply(values, function(v) v$applies_when), codes),
      domain = o$domain,
      default = o$default,
      default_value = if (is.numeric(o$default_value)) as.double(o$default_value) else NULL,
      default_evidence = o$default_evidence,
      result_affecting = isTRUE(o$result_affecting),
      applies_when = o$applies_when,
      live = o$live
    )
  })
  names(options) <- vapply(options, `[[`, "", "option_id")
  m$modes <- modes
  m$modes_demonstrated <- .ig_flatten(m$modes_demonstrated, "character")
  m$aliases <- .ig_flatten(m$aliases, "character")
  m$documented_example_ids <- .ig_flatten(m$documented_example_ids, "character")
  m$gpower_inputs <- .ig_flatten(m$gpower_inputs, "character")
  m$outputs <- .ig_flatten(m$outputs, "character")
  m$notes <- .ig_flatten(m$notes, "character")
  m$options <- options
  # Live G*Power labels per analysis type (certified inventory only).
  if (is.list(m$outputs_by_mode_live)) {
    m$outputs_by_mode_live <- lapply(m$outputs_by_mode_live, .ig_flatten, "character")
  }
  if (is.list(m$inputs_by_mode_live)) {
    m$inputs_by_mode_live <- lapply(m$inputs_by_mode_live, function(rows) {
      lapply(rows, function(r) {
        r$choices <- .ig_flatten(r$choices, "character")
        r$buttons <- .ig_flatten(r$buttons, "character")
        r
      })
    })
  }
  m
}

# Live G*Power input and output labels of one method and analysis type, as
# certified against the running application (NULL elements when the
# inventory entry is not certified). Used by comparators and engines to map
# G*Power labels.
.ig_live_labels <- function(id, mode) {
  meta <- .ig_get_metadata(id)
  if (!.ig_is_string(mode) || !(mode %in% meta$modes)) {
    .ig_abort("invalid_input", sprintf("mode '%s' is not offered for %s", paste(format(mode), collapse = ""), id))
  }
  rows <- meta$inputs_by_mode_live[[mode]]
  list(
    inputs = if (is.null(rows)) NULL else vapply(rows, function(r) r$label, ""),
    input_defaults = if (is.null(rows)) NULL else
      stats::setNames(vapply(rows, function(r) if (is.null(r$default)) NA_character_ else r$default, ""),
                      vapply(rows, function(r) r$label, "")),
    outputs = meta$outputs_by_mode_live[[mode]]
  )
}

.ig_collect_definitions <- function(ns) {
  providers <- sort(ls(ns, all.names = TRUE, pattern = "^\\.ig_define_methods_[A-Za-z0-9_]+$"))
  for (p in providers) {
    fn <- get(p, envir = ns, inherits = FALSE)
    if (!is.function(fn)) next
    defs <- tryCatch(fn(), error = function(e) {
      .ig_state$problems <- c(.ig_state$problems,
                              sprintf("%s() failed: %s", p, conditionMessage(e)))
      list()
    })
    if (inherits(defs, "ig_method")) defs <- list(defs)
    for (d in defs) {
      tryCatch(.ig_register_method(d), error = function(e) {
        .ig_state$problems <- c(.ig_state$problems,
                                sprintf("%s(): %s", p, conditionMessage(e)))
      })
    }
  }
  invisible(NULL)
}

.ig_collect_solvers <- function(ns) {
  hook <- get0(".ig_mode_solvers", envir = ns, inherits = FALSE)
  if (!is.function(hook)) return(invisible(NULL))
  solvers <- tryCatch(hook(), error = function(e) {
    .ig_state$problems <- c(.ig_state$problems,
                            sprintf(".ig_mode_solvers() failed: %s", conditionMessage(e)))
    list()
  })
  ok <- is.list(solvers) && (length(solvers) == 0L || .ig_is_named(solvers)) &&
    all(names(solvers) %in% .ig_modes) && all(vapply(solvers, is.function, logical(1)))
  if (!ok) {
    .ig_state$problems <- c(.ig_state$problems,
                            ".ig_mode_solvers() must return a list of functions named by mode")
    return(invisible(NULL))
  }
  .ig_state$solvers <- solvers
  invisible(NULL)
}

# Registry access ------------------------------------------------------------------

.ig_registry_metadata <- function() {
  if (!is.null(.ig_state$load_error)) {
    # A broken installation is not a statistical outcome; plain error.
    stop("the inferguidR method registry could not be loaded: ", .ig_state$load_error,
         call. = FALSE)
  }
  if (is.null(.ig_state$metadata)) .ig_registry_init(NULL)
  .ig_state$metadata
}

.ig_get_metadata <- function(id) {
  meta <- .ig_registry_metadata()
  if (!.ig_is_string(id)) .ig_abort("invalid_input", "a method id must be a single string")
  entry <- meta[[id]]
  if (is.null(entry)) {
    .ig_abort("invalid_input", sprintf("unknown method '%s'; see list_methods()", id))
  }
  entry
}

.ig_get_method <- function(id) .ig_state$definitions[[id]]

.ig_list_method_ids <- function() names(.ig_registry_metadata())

.ig_registry_problems <- function() .ig_state$problems

.ig_register_method <- function(def, replace = FALSE) {
  if (!inherits(def, "ig_method")) stop("only ig_method objects can be registered", call. = FALSE)
  meta <- .ig_registry_metadata()[[def$id]]
  if (is.null(meta)) {
    stop("method '", def$id, "' has no registry metadata (not in the reference inventory)",
         call. = FALSE)
  }
  if (!identical(def$family, meta$family)) {
    stop("method '", def$id, "' declares family '", def$family, "' but the registry says '",
         meta$family, "'", call. = FALSE)
  }
  extra <- setdiff(def$modes, meta$modes)
  if (length(extra)) {
    stop("method '", def$id, "' offers modes absent from the reference inventory: ",
         paste(extra, collapse = ", "), call. = FALSE)
  }
  if (!replace && !is.null(.ig_state$definitions[[def$id]])) {
    stop("method '", def$id, "' is already registered", call. = FALSE)
  }
  .ig_state$definitions[[def$id]] <- def
  invisible(def)
}

.ig_unregister_method <- function(id) {
  .ig_state$definitions[[id]] <- NULL
  invisible(NULL)
}

.ig_solver_for <- function(mode) .ig_state$solvers[[mode]]

.ig_implementation_status <- function(id) {
  def <- .ig_get_method(id)
  if (is.null(def)) return("not_implemented")
  with_solver <- vapply(def$modes, function(m) is.function(.ig_solver_for(m)), logical(1))
  if (!any(with_solver)) return("definition_only")
  if (all(with_solver)) "implemented" else "partial"
}

.ig_method_validation_status <- function(id) {
  def <- .ig_get_method(id)
  if (is.null(def)) .ig_get_metadata(id)$reference_status else def$status
}

# State of the G*Power reference comparison: "validated" (authentic
# application captures), "not_applicable" (the registry entry is an
# inferguidR extension outside the G*Power inventory, so no G*Power
# reference can exist) or "reference_pending".
.ig_gpower_reference_status <- function(meta, status) {
  if (identical(status, "gpower_validated")) return("validated")
  if (identical(meta$reference_status, "not_applicable")) return("not_applicable")
  "reference_pending"
}

# Validation summary used inside results.
.ig_validation_record <- function(id) {
  meta <- .ig_get_metadata(id)
  def <- .ig_get_method(id)
  status <- .ig_method_validation_status(id)
  list(
    method_validation_status = status,
    gpower_reference_status = .ig_gpower_reference_status(meta, status),
    independent_checker = if (is.null(def)) NA_character_ else def$checker,
    documented_example_ids = I(meta$documented_example_ids),
    fixture_ids = I(if (is.null(def)) character() else as.character(def$reference$fixture_ids))
  )
}

# Public API ---------------------------------------------------------------------------

#' List the registered statistical methods
#'
#' Returns one row per registered method with its implementation and
#' validation state: the 49 procedures of the G*Power reference inventory,
#' followed by inferguidR extensions that are not G*Power procedures (family
#' `"boschloo"`: Boschloo's unconditional exact test). `implementation_status`
#' is `"not_implemented"` (no method definition), `"definition_only"`
#' (forward model registered, no mode solver yet), `"partial"` or
#' `"implemented"`. An implemented method has a registered engine and mode
#' solvers; individual parameter branches can still return a typed refusal.
#' `validation_status` is the method's validation status;
#' inventory methods without an engine report `"reference_pending"`, and
#' extensions without an engine report `"not_applicable"` (no G*Power
#' reference can exist for them). No method is G*Power-validated until
#' authentic application captures exist.
#'
#' @param family Optional family filter (`"exact"`, `"F"`, `"t"`, `"z"`,
#'   `"chi2"`, or `"boschloo"` for the extension).
#' @return A data frame with columns `id`, `family`, `label`, `modes`
#'   (comma-separated), `implementation_status` and `validation_status`.
#' @examples
#' head(list_methods())
#' list_methods(family = "chi2")
#' @seealso [validation_status()]
#' @export
list_methods <- function(family = NULL) {
  meta <- .ig_registry_metadata()
  if (!is.null(family)) {
    fams <- unique(vapply(meta, `[[`, "", "family"))
    if (!is.character(family) || !all(family %in% fams)) {
      .ig_abort("invalid_input", sprintf("unknown family; available: %s",
                                         paste(fams, collapse = ", ")))
    }
    meta <- meta[vapply(meta, function(e) e$family %in% family, logical(1))]
  }
  ids <- names(meta)
  data.frame(
    id = ids,
    family = vapply(meta, `[[`, "", "family", USE.NAMES = FALSE),
    label = vapply(meta, function(e) e$label$en, "", USE.NAMES = FALSE),
    modes = vapply(meta, function(e) paste(e$modes, collapse = ", "), "", USE.NAMES = FALSE),
    implementation_status = vapply(ids, .ig_implementation_status, "", USE.NAMES = FALSE),
    validation_status = vapply(ids, .ig_method_validation_status, "", USE.NAMES = FALSE),
    stringsAsFactors = FALSE
  )
}

#' Validation state of a method
#'
#' Reports what is known about the correctness evidence of one method:
#' whether an engine exists, its validation status, the state of the G*Power
#' reference comparison (`"validated"`, `"reference_pending"`, or
#' `"not_applicable"` for inferguidR extensions that are not G*Power
#' procedures), the independent checker and the documented manual examples
#' that can serve as (non-authoritative) evidence. Pending reference status
#' means that complete comparison is not certified, even when captured cases
#' and independent numerical checks already exist.
#'
#' @param method_id A method id from [list_methods()].
#' @return An object of class `ig_validation_status` (a list).
#' @examples
#' validation_status("t_means_two_independent")
#' @export
validation_status <- function(method_id) {
  meta <- .ig_get_metadata(method_id)
  def <- .ig_get_method(method_id)
  impl <- .ig_implementation_status(method_id)
  status <- .ig_method_validation_status(method_id)
  gp_validated <- identical(status, "gpower_validated")
  gp_status <- .ig_gpower_reference_status(meta, status)
  solver_modes <- if (is.null(def)) character() else
    def$modes[vapply(def$modes, function(m) is.function(.ig_solver_for(m)), logical(1))]
  explanation <- c(
    if (is.null(def)) "No engine is implemented for this method yet; calculations return status 'unsupported'.",
    if (!is.null(def) && !length(solver_modes)) "A method definition exists but no mode solver is implemented yet.",
    if (gp_status == "not_applicable") "This procedure is not part of the G*Power inventory; no G*Power reference exists.",
    if (gp_status == "reference_pending") "G*Power reference comparison is not fully certified for this method. Available captures and passing regression cases are limited evidence; complete agreement with G*Power is not claimed.",
    if (identical(status, "development")) "The engine has not been checked by an independent implementation yet.",
    if (length(meta$documented_example_ids)) "Documented manual examples are secondary evidence, not application captures."
  )
  structure(list(
    method_id = method_id,
    family = meta$family,
    label = meta$label$en,
    implementation_status = impl,
    validation_status = status,
    gpower_reference = list(
      status = gp_status,
      fixture_ids = if (is.null(def)) character() else as.character(def$reference$fixture_ids)
    ),
    independent_check = list(
      checker = if (is.null(def)) NA_character_ else def$checker,
      status = if (identical(status, "independently_checked") || gp_validated) "passed" else "none"
    ),
    modes = meta$modes,
    modes_with_solver = solver_modes,
    documented_example_ids = meta$documented_example_ids,
    inventory_status = meta$inventory_status,
    explanation = explanation
  ), class = "ig_validation_status")
}

#' @export
format.ig_validation_status <- function(x, ...) {
  c(sprintf("<ig_validation_status> %s", x$method_id),
    sprintf("  label:                 %s", x$label),
    sprintf("  implementation status: %s", x$implementation_status),
    sprintf("  validation status:     %s", x$validation_status),
    sprintf("  G*Power reference:     %s", x$gpower_reference$status),
    sprintf("  independent check:     %s", x$independent_check$status),
    sprintf("  documented examples:   %s",
            if (length(x$documented_example_ids)) paste(x$documented_example_ids, collapse = ", ") else "none"),
    paste0("  - ", x$explanation))
}

#' @export
print.ig_validation_status <- function(x, ...) {
  cat(format(x), sep = "\n")
  invisible(x)
}
