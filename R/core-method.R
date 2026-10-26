# Method definitions -----------------------------------------------------------
#
# One method definition per statistical procedure of the reference inventory.
# Definitions are created in R/methods-<family>.R by the owning phase and are
# collected at load time from provider functions named
# `.ig_define_methods_<family>()` (see R/core-registry.R). The definition
# holds the forward power model and design rule; the generic mode solvers
# (phase 05) derive every analysis type from them.

.ig_modes <- c("a_priori", "post_hoc", "sensitivity", "criterion", "compromise")

.ig_method_statuses <- c("development", "independently_checked",
                         "reference_pending", "gpower_validated")

.ig_effect_directions <- c("increasing", "decreasing", "unknown")

#' Define a statistical method
#'
#' Creates and validates a method definition. Each definition corresponds to
#' one procedure of the reference inventory (its `id` must be a registry
#' method id) and is registered when the package loads.
#'
#' @param id Method id, as in the registry (see [list_methods()]).
#' @param family Test family id (`"exact"`, `"F"`, `"t"`, `"z"`, `"chi2"`).
#' @param label Display label: a list with `en` and optionally `de`.
#' @param params An [ig_params()] object. It must declare `alpha` and
#'   `power`, every solve target, and `error_ratio` when the compromise mode
#'   is offered.
#' @param options An [ig_options()] object with the result-affecting options.
#' @param modes Offered analysis modes, a subset of `"a_priori"`,
#'   `"post_hoc"`, `"sensitivity"`, `"criterion"` and `"compromise"`.
#' @param power Forward model `function(p)`; `p` is the fully specified named
#'   list of inputs and options. It returns a named list with at least
#'   `power` and `alpha_attained`; undefined values are `NA` with a reason.
#'   Optional elements: `beta` (computed directly, not as 1 - power),
#'   `critical`, `df`, `ncp`, `n_total`, `group_sizes`,
#'   `statistic_distribution`, `outputs` (named list), `notes` (reasons for
#'   `NA` values) and `assumptions`; see the forward-model contract in
#'   `R/num-modes.R`.
#' @param design Allocation rule `function(p, n_total)` returning integer
#'   group sizes named by their parameters (for example
#'   `c(n1 = 10, n2 = 10)`), or `NULL` when `n_total` has no admissible
#'   design (for example an odd N under equal allocation). `NULL` as the
#'   argument value means a single-sample design (N is the design). The rule
#'   may accept an optional argument `continuous = TRUE` returning fractional
#'   group sizes, used only as the starting point of the a priori search, and
#'   may signal a typed condition (for example `reference_pending`) for an
#'   unsupported branch.
#' @param discrete `TRUE` for exact or otherwise discrete tests.
#' @param effect_direction Monotone direction of power in the effect:
#'   `"increasing"`, `"decreasing"` or `"unknown"`, or a function `function(p)`
#'   returning one of those strings from the validated inputs and options.
#' @param engine Backend provenance: a list with `package`, `fun` and
#'   `version` (use `NA` for `version` when it is taken at run time).
#' @param checker Identifier of the independent checker (a path under
#'   `validation/independent/`), or `NA`.
#' @param reference Reference metadata, for example
#'   `list(gpower_path = "...", manual_pages = "...")`. A status of
#'   `"gpower_validated"` requires `reference$fixture_ids` naming authentic
#'   application-capture fixture records.
#' @param status Validation status: `"development"`,
#'   `"independently_checked"` (requires a `checker`), `"reference_pending"`
#'   or `"gpower_validated"`.
#' @param targets Optional named list of solve targets per mode (character
#'   vectors; the first element is the default). Defaults: `n_total`
#'   (a priori), `power` (post hoc), the single effect-role parameter
#'   (sensitivity), `alpha` (criterion) and the joint pair `alpha`, `power`
#'   (compromise).
#' @param constraints Optional `function(p)` returning a character vector of
#'   cross-parameter inconsistencies (empty when consistent), for example
#'   group sizes that do not add up to `n_total`. A callback may additionally
#'   declare `mode` and/or `target` arguments to receive the validated analysis
#'   mode and selected solve target(s). Existing one-argument callbacks work unchanged.
#' @param a_priori_search Optional method-specific search of the a priori
#'   mode, `function(p, n_lo, n_max)`, used instead of the generic integer
#'   search (for example a discrete test that can exclude designs by a proven
#'   bound). `p` holds the specification's inputs and options, and
#'   `[n_lo, n_max]` is the admissible range of `n_total`. Return a list with
#'   `status = "ok"`, `n` (the smallest admissible total reaching the target),
#'   `power` (its power), and a `diagnostics` list. When there is no answer,
#'   return a list with `status` (the status code), `reason` (explanatory text)
#'   and `diagnostics`. The solver checks the returned design against
#'   `design` and recomputes its power with `power`.
#' @param sensitivity_search Optional `function(p, target_name)` replacing
#'   the generic sensitivity scan. Return a status/reason for failure, or a
#'   list with `status = "ok"`, `value` (the effect), a `diagnostics` list,
#'   and a character vector of `warnings`. The effect must be in its parameter domain
#'   and its recomputed power must reach the target within solver tolerance.
#'   A method-specific maximum below target is allowed only with
#'   `diagnostics$target_reached = FALSE`, `diagnostics$max_power` agreeing
#'   with recomputed power, and a non-empty explanatory warning.
#' @param sensitivity_domain Optional `function(p)` returning a scalar real
#'   [ig_domain_real()] domain for the generic sensitivity scan. It must be
#'   contained in the effect parameter domain. This restricts the search using
#'   validated inputs/options, for example to admissible correlations.
#' @return An object of class `ig_method`.
#' @family method definition API
#' @export
ig_method <- function(id, family, label, params, options = ig_options(), modes,
                      power, design = NULL, discrete = FALSE,
                      effect_direction = "unknown",
                      engine = list(package = NA_character_, fun = NA_character_,
                                    version = NA_character_),
                      checker = NA_character_, reference = list(),
                      status = "development", targets = NULL, constraints = NULL,
                      a_priori_search = NULL, sensitivity_search = NULL, sensitivity_domain = NULL) {
  problems <- character()
  add <- function(...) problems <<- c(problems, sprintf(...))

  if (!.ig_is_string(id) || !grepl("^[a-z0-9]+(_[a-z0-9]+)*$", id)) {
    add("`id` must be a lower-case identifier")
  }
  if (!.ig_is_string(family)) add("`family` must be a single string")
  label <- tryCatch(.ig_check_label(label, "label"), error = function(e) {
    add("%s", conditionMessage(e))
    NULL
  })
  if (!inherits(params, "ig_params")) add("`params` must be created with ig_params()")
  if (!inherits(options, "ig_options")) add("`options` must be created with ig_options()")
  if (!is.character(modes) || !length(modes) || anyNA(modes) ||
      anyDuplicated(modes) || !all(modes %in% .ig_modes)) {
    add("`modes` must be unique values from: %s", paste(.ig_modes, collapse = ", "))
  }
  if (!is.function(power) || length(formals(power)) < 1L) {
    add("`power` must be a function(p)")
  }
  if (!is.null(design) && (!is.function(design) || length(formals(design)) < 2L)) {
    add("`design` must be NULL or a function(p, n_total)")
  }
  if (!.ig_is_flag(discrete)) add("`discrete` must be TRUE or FALSE")
  if (is.function(effect_direction)) {
    if (length(formals(effect_direction)) < 1L) add("`effect_direction` must be a function(p)")
  } else if (!.ig_is_string(effect_direction) || !(effect_direction %in% .ig_effect_directions)) {
    add("`effect_direction` must be one of: %s, or a function(p)", paste(.ig_effect_directions, collapse = ", "))
  }
  if (!is.list(engine) || !all(c("package", "fun", "version") %in% names(engine))) {
    add("`engine` must be a list with package, fun and version")
  }
  if (length(checker) != 1L || !(is.na(checker) || .ig_is_string(checker))) {
    add("`checker` must be a single string or NA")
  }
  if (!is.list(reference) || (length(reference) && !.ig_is_named(reference))) {
    add("`reference` must be a named list")
  }
  if (!.ig_is_string(status) || !(status %in% .ig_method_statuses)) {
    add("`status` must be one of: %s", paste(.ig_method_statuses, collapse = ", "))
  } else {
    if (status == "independently_checked" && is.na(checker)) {
      add("status 'independently_checked' requires an independent `checker`")
    }
    if (status == "gpower_validated") {
      ids <- reference$fixture_ids
      if (!is.character(ids) || !length(ids) || anyNA(ids)) {
        add("status 'gpower_validated' requires reference$fixture_ids naming authentic application-capture fixtures")
      }
    }
  }
  if (!is.null(constraints) && !is.function(constraints)) {
    add("`constraints` must be NULL or a function(p)")
  }
  if (!is.null(a_priori_search) &&
      (!is.function(a_priori_search) ||
       !all(c("n_lo", "n_max") %in% names(formals(a_priori_search))))) {
    add("`a_priori_search` must be NULL or a function(p, n_lo, n_max)")
  }
  if (!is.null(a_priori_search) && is.character(modes) && !("a_priori" %in% modes)) {
    add("`a_priori_search` is given but the a priori mode is not offered")
  }
  if (!is.null(sensitivity_domain) &&
      (!is.function(sensitivity_domain) || length(formals(sensitivity_domain)) < 1L)) {
    add("`sensitivity_domain` must be NULL or a function(p)")
  }
  if (!is.null(sensitivity_domain) && is.character(modes) && !("sensitivity" %in% modes)) {
    add("`sensitivity_domain` is given but sensitivity mode is not offered")
  }
  if (!is.null(sensitivity_search) &&
      (!is.function(sensitivity_search) ||
       !all(c("p", "target_name") %in% names(formals(sensitivity_search))))) {
    add("`sensitivity_search` must be NULL or a function(p, target_name)")
  }
  if (!is.null(sensitivity_search) && is.character(modes) && !("sensitivity" %in% modes)) {
    add("`sensitivity_search` is given but sensitivity mode is not offered")
  }

  resolved_targets <- NULL
  if (!length(problems)) {
    pnames <- names(params)
    onames <- names(options)
    clash <- intersect(pnames, onames)
    if (length(clash)) add("names used for both a parameter and an option: %s", paste(clash, collapse = ", "))
    for (req in c("alpha", "power")) {
      if (!(req %in% pnames)) add("parameter '%s' must be declared", req)
    }
    if ("compromise" %in% modes && !("error_ratio" %in% pnames)) {
      add("the compromise mode needs an 'error_ratio' parameter")
    }
    resolved_targets <- tryCatch(.ig_resolve_targets(params, modes, targets),
                                 error = function(e) {
                                   add("%s", conditionMessage(e))
                                   NULL
                                 })
  }
  if (length(problems)) {
    stop(sprintf("invalid method definition '%s':\n%s",
                 if (.ig_is_string(id)) id else "?",
                 paste0("* ", problems, collapse = "\n")), call. = FALSE)
  }

  structure(list(id = id, family = family, label = label, params = params,
                 options = options, modes = .ig_modes[.ig_modes %in% modes],
                 targets = resolved_targets, power = power, design = design,
                 discrete = discrete, effect_direction = effect_direction,
                 engine = engine, checker = as.character(checker),
                 reference = reference, status = status, constraints = constraints,
                 a_priori_search = a_priori_search, sensitivity_search = sensitivity_search,
                 sensitivity_domain = sensitivity_domain),
            class = "ig_method")
}

# Resolve the allowed solve targets per mode.
.ig_resolve_targets <- function(params, modes, targets) {
  if (!is.null(targets)) {
    if (!is.list(targets) || !.ig_is_named(targets) || !all(names(targets) %in% .ig_modes)) {
      stop("`targets` must be a list named by mode", call. = FALSE)
    }
  }
  pnames <- names(params)
  roles <- vapply(params, `[[`, "", "role")
  out <- list()
  for (m in modes) {
    tg <- targets[[m]]
    if (is.null(tg)) {
      tg <- switch(m,
                   a_priori = "n_total",
                   post_hoc = "power",
                   sensitivity = {
                     eff <- pnames[roles == "effect"]
                     if (length(eff) != 1L) {
                       stop("the sensitivity target is ambiguous (", length(eff),
                            " effect parameters); declare `targets$sensitivity`", call. = FALSE)
                     }
                     eff
                   },
                   criterion = "alpha",
                   compromise = c("alpha", "power"))
    }
    if (!is.character(tg) || !length(tg) || anyNA(tg) || anyDuplicated(tg)) {
      stop("targets for mode '", m, "' must be unique parameter names", call. = FALSE)
    }
    missing <- setdiff(tg, pnames)
    if (length(missing)) {
      stop("solve target(s) for mode '", m, "' are not declared parameters: ",
           paste(missing, collapse = ", "), call. = FALSE)
    }
    if (m == "compromise" && !setequal(tg, c("alpha", "power"))) {
      stop("the compromise mode solves alpha and power jointly", call. = FALSE)
    }
    out[[m]] <- tg
  }
  out
}

#' @export
format.ig_method <- function(x, ...) {
  c(sprintf("<ig_method> %s (%s family)", x$id, x$family),
    sprintf("  label:   %s", x$label$en),
    sprintf("  modes:   %s", paste(x$modes, collapse = ", ")),
    sprintf("  params:  %s", paste(names(x$params), collapse = ", ")),
    sprintf("  options: %s", if (length(x$options)) paste(names(x$options), collapse = ", ") else "none"),
    sprintf("  status:  %s", x$status))
}

#' @export
print.ig_method <- function(x, ...) {
  cat(format(x), sep = "\n")
  invisible(x)
}
