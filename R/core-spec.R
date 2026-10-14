# Analysis specifications (power_spec) ----------------------------------------------
#
# A power_spec is the complete, immutable statement of one power-analysis
# question: method, analysis mode, solve target(s), every input and every
# result-affecting option, with all applied defaults written out. It is
# validated against the method definition when one is registered, otherwise
# against the registry metadata (and labelled as such). The spec hash is the
# SHA-256 of the canonical JSON of the analysis key (schema version, object
# type, method, mode, target, inputs, options).

# Generic solve targets used when only registry metadata is available. The
# sensitivity target is the method's effect-size parameter; its name is not
# known without a definition, except where the inventory enumerates
# alternative sensitivity targets.
.ig_registry_default_targets <- list(
  a_priori = "n_total", post_hoc = "power", sensitivity = "effect_size",
  criterion = "alpha", compromise = c("alpha", "power")
)

# Standard inputs validated for registry-only methods.
.ig_registry_standard_params <- function() {
  ig_params(ig_param_alpha(), ig_param_power(), ig_param_error_ratio(),
            ig_param_n_total(lower = 1))
}

# Standard inputs each mode needs (registry-only validation).
.ig_registry_required <- list(
  a_priori = c("alpha", "power"), post_hoc = "alpha",
  sensitivity = c("alpha", "power"), criterion = "power",
  compromise = "error_ratio"
)

#' Create an analysis specification
#'
#' Builds the immutable specification of one power analysis. The method
#' (see [list_methods()]), the analysis `mode` and the inputs are validated:
#' unknown or inapplicable inputs and options, domain violations, a supplied
#' solve target, missing required values and inconsistent duplicate inputs
#' are rejected with a typed `invalid_input` error. Every default that is
#' applied is written into the specification.
#'
#' Modes and their solve targets:
#'
#' * `"a_priori"`: total sample size (`n_total`) from effect, `alpha` and
#'   `power`;
#' * `"post_hoc"`: `power` from effect, `alpha` and the sample size;
#' * `"sensitivity"`: the effect size from `alpha`, `power` and the sample
#'   size;
#' * `"criterion"`: `alpha` (and critical values) from effect, `power` and
#'   the sample size;
#' * `"compromise"`: `alpha` and `power` jointly from effect, sample size and
#'   the ratio `error_ratio` = beta/alpha.
#'
#' For methods without a registered engine the specification is validated
#' against the registry metadata only (standard inputs, options and modes);
#' other inputs are recorded as `unvalidated_inputs`. Such specifications can
#' be stored and exported, and [calculate_power()] returns status
#' `"unsupported"` for them.
#'
#' Inputs may be given in `inputs`, in `...`, or both. The same name given
#' twice with identical values is accepted once; with different values it is
#' rejected.
#'
#' @param method Method id.
#' @param mode Analysis mode.
#' @param inputs Named list of inputs (numbers, integers, flags or codes).
#' @param options Named list (or named character vector) of option codes.
#' @param target Optional solve target. Defaults to the mode's target; needed
#'   only where a method offers alternative targets.
#' @param ... Further named inputs.
#' @return An immutable object of class `power_spec`.
#' @examples
#' spec <- analysis_spec("t_means_two_independent", "a_priori",
#'   inputs = list(alpha = 0.05, power = 0.8, d = 0.5),
#'   options = list(tails = "two"))
#' spec
#' # This method has an engine (phase 05 vertical slice), so a number is returned:
#' calculate_power(spec)
#' @seealso [update.power_spec()], [calculate_power()], [ig_to_json()]
#' @export
analysis_spec <- function(method, mode, inputs = list(), options = list(),
                          target = NULL, ...) {
  if (missing(method)) method <- NULL
  if (missing(mode)) mode <- NULL
  .ig_build_spec(method, mode, inputs, list(...), options, target)
}

.ig_build_spec <- function(method, mode, inputs, dots, options, target) {
  if (!.ig_is_string(method)) {
    .ig_abort("invalid_input", "`method` must be a single method id; see list_methods()")
  }
  meta <- .ig_get_metadata(method)
  def <- .ig_get_method(method)
  offered <- if (is.null(def)) meta$modes else def$modes
  if (is.null(mode)) {
    .ig_abort_invalid(sprintf("`mode` is missing; offered modes for %s: %s", method,
                              paste(offered, collapse = ", ")))
  }
  if (!.ig_is_string(mode) || !(mode %in% .ig_modes)) {
    .ig_abort_invalid(sprintf("unknown mode '%s'; modes are: %s",
                              paste(format(mode), collapse = ""), paste(.ig_modes, collapse = ", ")))
  }
  if (!(mode %in% offered)) {
    .ig_abort_invalid(sprintf("mode '%s' is not offered for %s; offered: %s", mode, method,
                              paste(offered, collapse = ", ")))
  }

  problems <- character()
  merged_inputs <- .ig_merge_named(list(inputs, dots), "input")
  problems <- c(problems, merged_inputs$problems)
  merged_options <- .ig_merge_named(list(options), "option")
  problems <- c(problems, merged_options$problems)

  # Validate options before choosing an option-dependent solve target.
  opt_specs <- if (is.null(def)) .ig_registry_option_specs(meta) else .ig_definition_option_specs(def)
  ov <- .ig_validate_options(merged_options$values, opt_specs, mode, method)
  problems <- c(problems, ov$problems)

  tg <- .ig_select_target(def, meta, mode, target, ov$values)
  problems <- c(problems, tg$problems)
  target <- tg$target

  # Inputs
  iv <- if (is.null(def)) {
    .ig_validate_inputs_registry(merged_inputs$values, meta, mode, target, ov)
  } else {
    .ig_validate_inputs_definition(merged_inputs$values, def, mode, target, ov$values)
  }
  problems <- c(problems, iv$problems)

  if (!length(problems) && !is.null(def) && is.function(def$constraints)) {
    context <- list(mode = mode, target = target)
    declared <- names(formals(def$constraints))
    args <- c(list(c(iv$values, ov$values)), context[intersect(names(context), declared)])
    extra <- tryCatch(do.call(def$constraints, args), error = function(e) {
      paste("constraint check failed:", conditionMessage(e))
    })
    if (length(extra)) problems <- c(problems, paste("inconsistent inputs:", as.character(extra)))
  }
  if (length(problems)) {
    .ig_abort_invalid(problems, sprintf("invalid specification for %s (%s)", method, mode))
  }

  .ig_new_spec(method_id = method, family = meta$family, mode = mode, target = target,
               inputs = iv$values, options = ov$values,
               defaults_applied = list(inputs = iv$defaults, options = ov$defaults),
               validation_level = if (is.null(def)) "registry_metadata" else "method_definition",
               unvalidated_inputs = iv$unvalidated)
}

.ig_new_spec <- function(method_id, family, mode, target, inputs, options,
                         defaults_applied, validation_level, unvalidated_inputs) {
  key <- .ig_spec_key(method_id, mode, target, inputs, options)
  structure(list(
    schema_version = .ig_schema_versions[["power_spec"]],
    object_type = "power_spec",
    method_id = method_id,
    family = family,
    mode = mode,
    target = target,
    joint_target = identical(mode, "compromise"),
    inputs = inputs,
    options = options,
    defaults_applied = list(inputs = as.character(defaults_applied$inputs),
                            options = as.character(defaults_applied$options)),
    validation_level = validation_level,
    unvalidated_inputs = as.character(unvalidated_inputs),
    spec_hash = .ig_hash_value(key)
  ), class = c("power_spec", "ig_frozen"))
}

# The hashed analysis key.
.ig_spec_key <- function(method_id, mode, target, inputs, options) {
  list(schema_version = .ig_schema_versions[["power_spec"]],
       object_type = "power_spec",
       method_id = method_id,
       mode = mode,
       target = I(as.character(target)),
       inputs = .ig_named_or_empty(inputs),
       options = .ig_named_or_empty(options))
}

# Merge several named lists. Repeated names with identical values collapse;
# repeated names with different values are reported.
.ig_merge_named <- function(sources, what) {
  problems <- character()
  pairs_n <- character()
  pairs_v <- list()
  for (src in sources) {
    if (is.null(src) || length(src) == 0L) next
    if (is.atomic(src) && !is.null(names(src))) src <- as.list(src)
    if (!is.list(src) || is.object(src)) {
      problems <- c(problems, sprintf("%ss must be given as a named list", what))
      next
    }
    nms <- names(src)
    if (is.null(nms) || anyNA(nms) || !all(nzchar(nms))) {
      problems <- c(problems, sprintf("every %s must be named", what))
      next
    }
    pairs_n <- c(pairs_n, nms)
    pairs_v <- c(pairs_v, unname(src))
  }
  values <- list()
  for (nm in unique(pairs_n)) {
    idx <- which(pairs_n == nm)
    vals <- pairs_v[idx]
    if (any(vapply(vals, is.null, logical(1)))) {
      problems <- c(problems, sprintf("%s '%s' is NULL", what, nm))
      next
    }
    first <- .ig_comparable(vals[[1]])
    same <- vapply(vals, function(v) identical(.ig_comparable(v), first), logical(1))
    if (!all(same)) {
      shown <- vapply(vals, .ig_show_value, character(1))
      problems <- c(problems, sprintf("inconsistent duplicate %s '%s': %s", what, nm,
                                      paste(unique(shown), collapse = " vs ")))
    }
    # The first value is kept so that later checks do not also report the
    # name as missing; the conflict above already invalidates the request.
    values[[nm]] <- vals[[1]]
  }
  list(values = values, problems = problems)
}

.ig_comparable <- function(v) {
  if (is.factor(v)) v <- as.character(v)
  if (is.integer(v)) v <- as.double(v)
  if (is.list(v) && !is.object(v)) v <- lapply(v, .ig_comparable)
  unname(v)
}

.ig_show_value <- function(v) {
  if (is.null(v)) return("NULL")
  if (is.list(v)) return("<list>")
  if (is.double(v)) v <- .ig_format_double(v)
  txt <- paste(as.character(v), collapse = ", ")
  if (length(v) > 1L) paste0("c(", txt, ")") else txt
}

.ig_select_target <- function(def, meta, mode, target, options = list()) {
  problems <- character()
  if (mode == "compromise") {
    if (!is.null(target) && !(is.character(target) && setequal(target, c("alpha", "power")))) {
      problems <- "the compromise mode solves 'alpha' and 'power' jointly; `target` must be c(\"alpha\", \"power\") or NULL"
    }
    return(list(target = c("alpha", "power"), problems = problems))
  }
  allowed <- if (!is.null(def)) def$targets[[mode]] else .ig_registry_targets(meta, mode)
  if (!is.null(def)) {
    active <- vapply(allowed, function(nm) {
      p <- def$params[[nm]]
      is.null(p) || .ig_eval_condition(p$applies_when, mode, options)
    }, logical(1))
    allowed <- allowed[active]
    if (!length(allowed)) return(list(target = character(),
      problems = "no solve target applies with these options"))
  }
  if (is.null(target)) {
    # A definition declares its default target first; the registry gives no
    # default where the inventory lists alternatives.
    if (length(allowed) == 1L || !is.null(def)) return(list(target = allowed[1], problems = problems))
    return(list(target = allowed[1], problems = sprintf(
      "missing solve target: mode '%s' of %s offers alternative targets (%s); supply `target`",
      mode, meta$method_id, paste(allowed, collapse = ", "))))
  }
  if (!.ig_is_string(target) || !(target %in% allowed)) {
    problems <- sprintf("target '%s' is not a solve target of mode '%s'; allowed: %s",
                        paste(format(target), collapse = ", "), mode, paste(allowed, collapse = ", "))
    return(list(target = allowed[1], problems = problems))
  }
  list(target = target, problems = problems)
}

# Registry-only targets: the inventory enumerates alternative sensitivity
# targets for a few methods (coverage dimension `sensitivity_target` or
# `sensitivity_target_index`); then no default is chosen.
.ig_registry_targets <- function(meta, mode) {
  if (mode == "sensitivity") {
    dims <- meta$coverage$branch_dimensions
    alt <- dims$sensitivity_target %||% dims$sensitivity_target_index
    if (!is.null(alt)) return(.ig_flatten(alt, "character"))
  }
  .ig_registry_default_targets[[mode]]
}

# Option specifications in a common shape --------------------------------------------

.ig_definition_option_specs <- function(def) {
  lapply(def$options, function(o) {
    list(name = o$name, codes = o$values, default = o$default,
         condition = o$applies_when, applicability_known = TRUE,
         result_affecting = TRUE, value_type = "choice", applies_text = NULL)
  })
}

.ig_registry_option_specs <- function(meta) {
  lapply(meta$options, function(o) {
    aw <- o$applies_when
    list(name = o$option_id, codes = o$codes, default = o$default,
         condition = if (is.null(aw)) NULL else aw$condition,
         applicability_known = is.null(aw) || !is.null(aw$condition),
         result_affecting = isTRUE(o$result_affecting),
         value_type = o$value_type,
         applies_text = if (is.null(aw)) NULL else aw$text,
         domain = o$domain)
  })
}

# Validate supplied options; returns effective values (applicable options
# only), the names whose value is a default, and problems.
.ig_validate_options <- function(supplied, specs, mode, method) {
  problems <- character()
  invalid <- character()
  choice_specs <- specs[vapply(specs, function(s) identical(s$value_type, "choice") &&
                                 isTRUE(s$result_affecting), logical(1))]
  for (nm in names(supplied)) {
    s <- specs[[nm]]
    if (is.null(s)) {
      problems <- c(problems, sprintf("unknown option '%s' for %s%s", nm, method,
                                      if (length(choice_specs)) paste0("; options: ",
                                        paste(names(choice_specs), collapse = ", ")) else
                                        " (the method has no options)"))
      next
    }
    if (identical(s$value_type, "continuous")) {
      problems <- c(problems, sprintf("'%s' is a numeric input; supply it in `inputs`", nm))
      next
    }
    if (!isTRUE(s$result_affecting)) {
      problems <- c(problems, sprintf(
        "'%s' is a presentation setting and not part of an analysis specification", nm))
      next
    }
    v <- supplied[[nm]]
    if (!.ig_is_string(v) || !(v %in% s$codes)) {
      problems <- c(problems, sprintf("option '%s' must be one of: %s (got %s)", nm,
                                      paste(s$codes, collapse = ", "), .ig_show_value(v)))
      invalid <- c(invalid, nm)
    }
  }
  valid_supplied <- supplied[names(supplied) %in% names(choice_specs)]
  valid_supplied <- valid_supplied[vapply(names(valid_supplied), function(nm) {
    .ig_is_string(valid_supplied[[nm]]) && valid_supplied[[nm]] %in% choice_specs[[nm]]$codes
  }, logical(1))]

  # Effective values: iterate to a fixed point because applicability of one
  # option can depend on the value of another.
  eff <- list()
  for (iter in seq_len(length(choice_specs) + 2L)) {
    new <- list()
    for (s in choice_specs) {
      applies <- if (s$applicability_known) .ig_eval_condition(s$condition, mode, eff) else
        !is.null(valid_supplied[[s$name]])
      val <- valid_supplied[[s$name]] %||% (if (s$applicability_known) s$default else NULL)
      if (applies && !is.null(val)) new[[s$name]] <- val
    }
    if (identical(new, eff)) break
    eff <- new
  }

  defaults <- character()
  for (s in choice_specs) {
    applies <- if (s$applicability_known) .ig_eval_condition(s$condition, mode, eff) else TRUE
    supplied_here <- !is.null(valid_supplied[[s$name]])
    if (supplied_here && !applies) {
      problems <- c(problems, sprintf("option '%s' does not apply here%s", s$name,
                                      if (!is.null(s$applies_text)) paste0(" (applies when: ", s$applies_text, ")") else
                                        .ig_condition_text(s$condition)))
    }
    if (applies && s$applicability_known && !supplied_here && !(s$name %in% invalid)) {
      if (is.null(s$default)) {
        problems <- c(problems, sprintf(
          "option '%s' must be supplied (no documented default); allowed: %s",
          s$name, paste(s$codes, collapse = ", ")))
      } else {
        defaults <- c(defaults, s$name)
      }
    }
  }
  ord <- names(choice_specs)[names(choice_specs) %in% names(eff)]
  list(values = eff[ord], defaults = intersect(ord, defaults), problems = problems)
}

.ig_condition_text <- function(cond) {
  if (is.null(cond)) return("")
  paste0(" (applies when: ", gsub("\\s+", " ", paste(deparse(cond), collapse = "")), ")")
}

# Inputs validated against a method definition -------------------------------------------

.ig_validate_inputs_definition <- function(supplied, def, mode, target, options) {
  problems <- character()
  params <- def$params
  alias_map <- character()
  for (p in params) for (a in p$aliases) alias_map[a] <- p$name

  # Resolve aliases, with the duplicate rule.
  canon <- list()
  for (nm in names(supplied)) {
    cn <- if (nm %in% names(params)) nm else if (nm %in% names(alias_map)) alias_map[[nm]] else NA
    if (is.na(cn)) {
      problems <- c(problems, sprintf("unknown input '%s' for %s; inputs: %s", nm, def$id,
                                      paste(names(params), collapse = ", ")))
      next
    }
    if (!is.null(canon[[cn]]) && !identical(.ig_comparable(canon[[cn]]), .ig_comparable(supplied[[nm]]))) {
      problems <- c(problems, sprintf("inconsistent duplicate input '%s' (given as '%s' with a different value)",
                                      cn, nm))
      next
    }
    canon[[cn]] <- supplied[[nm]]
  }

  values <- list()
  defaults <- character()
  for (p in params) {
    nm <- p$name
    is_target <- nm %in% target
    in_mode <- is.null(p$modes) || mode %in% p$modes
    applies <- in_mode && !is_target && .ig_eval_condition(p$applies_when, mode, options)
    given <- !is.null(canon[[nm]])
    if (given && is_target) {
      problems <- c(problems, sprintf(
        "input '%s' is the solve target of mode '%s' and must not be supplied", nm, mode))
      next
    }
    if (given && !applies) {
      problems <- c(problems, sprintf("input '%s' is not used in mode '%s' with these options", nm, mode))
      next
    }
    if (!applies) next
    if (!given) {
      if (is.null(p$default)) {
        problems <- c(problems, sprintf("missing input '%s' (%s)", nm, p$label$en))
        next
      }
      values[[nm]] <- p$default
      defaults <- c(defaults, nm)
      next
    }
    problem <- .ig_domain_check(p$domain, canon[[nm]], nm)
    if (!is.null(problem)) {
      problems <- c(problems, problem)
      next
    }
    values[[nm]] <- .ig_domain_coerce(p$domain, canon[[nm]])
  }
  list(values = values, defaults = defaults, unvalidated = character(), problems = problems)
}

# Inputs validated against registry metadata only ----------------------------------------

.ig_validate_inputs_registry <- function(supplied, meta, mode, target, ov) {
  problems <- character()
  std <- .ig_registry_standard_params()
  continuous <- Filter(function(o) identical(o$value_type, "continuous"), meta$options)
  values <- list()
  unvalidated <- character()

  for (nm in intersect(names(supplied), target)) {
    problems <- c(problems, sprintf(
      "input '%s' is the solve target of mode '%s' and must not be supplied", nm, mode))
  }
  for (p in std) {
    nm <- p$name
    if (nm %in% target || is.null(supplied[[nm]])) next
    if (!is.null(p$modes) && !(mode %in% p$modes)) {
      problems <- c(problems, sprintf("input '%s' is not used in mode '%s'", nm, mode))
      next
    }
    problem <- .ig_domain_check(p$domain, supplied[[nm]], nm)
    if (!is.null(problem)) problems <- c(problems, problem) else
      values[[nm]] <- .ig_domain_coerce(p$domain, supplied[[nm]])
  }
  for (nm in setdiff(.ig_registry_required[[mode]], c(target, names(supplied)))) {
    problems <- c(problems, sprintf("missing input '%s' (%s)", nm, std[[nm]]$label$en))
  }
  for (o in continuous) {
    nm <- o$option_id
    if (is.null(supplied[[nm]]) || nm %in% target) next
    aw <- o$applies_when
    if (!is.null(aw) && !is.null(aw$condition) &&
        !.ig_eval_condition(aw$condition, mode, ov$values)) {
      problems <- c(problems, sprintf("input '%s' does not apply here (applies when: %s)", nm, aw$text))
      next
    }
    problem <- .ig_registry_domain_check(o$domain, supplied[[nm]], nm)
    if (!is.null(problem)) problems <- c(problems, problem) else
      values[[nm]] <- as.double(supplied[[nm]])
  }
  known <- c(names(std), vapply(continuous, `[[`, "", "option_id"))
  choice_ids <- names(meta$options)[vapply(meta$options, function(o) {
    identical(o$value_type, "choice")
  }, logical(1))]
  others <- sort(setdiff(names(supplied), c(known, target)), method = "radix")
  for (nm in others) {
    v <- supplied[[nm]]
    if (nm %in% choice_ids) {
      problems <- c(problems, sprintf("'%s' is an option; supply it in `options`", nm))
      next
    }
    if (!grepl(.ig_name_pattern, nm)) {
      problems <- c(problems, sprintf("input name '%s' is not a valid identifier", nm))
      next
    }
    ok <- is.atomic(v) && length(v) >= 1L && !anyNA(v) &&
      ((is.numeric(v) && all(is.finite(v))) ||
         (length(v) == 1L && (is.character(v) || is.logical(v))))
    if (!ok) {
      problems <- c(problems, sprintf(
        "input '%s' must be finite numbers or a single string/flag", nm))
      next
    }
    values[[nm]] <- if (is.numeric(v)) as.double(v) else v
    unvalidated <- c(unvalidated, nm)
  }
  order_names <- c(names(std), vapply(continuous, `[[`, "", "option_id"), others)
  values <- values[order_names[order_names %in% names(values)]]
  list(values = values, defaults = character(), unvalidated = unvalidated, problems = problems)
}

# Domain check for a parsed registry domain (bounds may be NULL = unknown).
.ig_registry_domain_check <- function(domain, value, name) {
  if (!is.numeric(value) || length(value) != 1L || !is.finite(value)) {
    return(sprintf("'%s' must be a single finite number", name))
  }
  if (is.null(domain)) return(NULL)
  lo <- domain$lower
  hi <- domain$upper
  if (is.numeric(lo) && length(lo) == 1L) {
    ok <- if (isTRUE(domain$lower_inclusive)) value >= lo else value > lo
    if (!ok) return(sprintf("'%s' must lie in %s (got %s)", name, domain$text, .ig_format_double(value)))
  }
  if (is.numeric(hi) && length(hi) == 1L) {
    ok <- if (isTRUE(domain$upper_inclusive)) value <= hi else value < hi
    if (!ok) return(sprintf("'%s' must lie in %s (got %s)", name, domain$text, .ig_format_double(value)))
  }
  NULL
}

# Update ------------------------------------------------------------------------------------

#' Derive a new specification from an existing one
#'
#' Specifications are immutable; `update()` creates a new, fully validated
#' specification. Inputs and options given here replace the old ones; a value
#' of `NULL` removes an input or option. Applied defaults are re-derived.
#' When the mode changes, inputs that become solve targets of the new mode
#' are dropped.
#'
#' @param object A `power_spec`.
#' @param ... Further named inputs to set (or remove with `NULL`).
#' @param inputs Named list of inputs to set (or remove with `NULL`).
#' @param options Named list of options to set (or remove with `NULL`).
#' @param mode New analysis mode (optional).
#' @param target New solve target (optional).
#' @return A new `power_spec`.
#' @examples
#' spec <- analysis_spec("t_means_one_sample", "a_priori",
#'   inputs = list(alpha = 0.05, power = 0.8, d = 0.5), options = list(tails = "two"))
#' update(spec, inputs = list(alpha = 0.01))
#' update(spec, mode = "post_hoc", inputs = list(n_total = 40))
#' @importFrom stats update
#' @export
update.power_spec <- function(object, ..., inputs = list(), options = list(),
                              mode = NULL, target = NULL) {
  dots <- list(...)
  new_mode <- if (is.null(mode)) object$mode else mode
  base_inputs <- object$inputs[setdiff(names(object$inputs), object$defaults_applied$inputs)]
  base_options <- object$options[setdiff(names(object$options), object$defaults_applied$options)]
  apply_changes <- function(base, changes) {
    for (nm in names(changes)) {
      if (is.null(changes[[nm]])) base[[nm]] <- NULL else base[[nm]] <- changes[[nm]]
    }
    base
  }
  if (is.atomic(inputs) && !is.null(names(inputs))) inputs <- as.list(inputs)
  if (is.atomic(options) && !is.null(names(options))) options <- as.list(options)
  changes <- c(inputs, dots)
  nms <- names(changes)
  if (length(changes) && (is.null(nms) || anyNA(nms) || !all(nzchar(nms)))) {
    .ig_abort_invalid("every input must be named")
  }
  # The same name twice must carry the same value (or NULL twice).
  problems <- character()
  for (nm in unique(nms[duplicated(nms)])) {
    vals <- changes[nms == nm]
    same <- vapply(vals, function(v) identical(.ig_comparable(v), .ig_comparable(vals[[1]])),
                   logical(1))
    if (!all(same)) problems <- c(problems, sprintf("inconsistent duplicate input '%s'", nm))
  }
  if (length(problems)) .ig_abort_invalid(problems)
  changes <- changes[!duplicated(nms)]

  base_inputs <- apply_changes(base_inputs, changes)
  base_options <- apply_changes(base_options, options)
  if (!identical(new_mode, object$mode)) {
    # Inputs that become solve targets of the new mode are dropped.
    def <- .ig_get_method(object$method_id)
    drop <- if (identical(new_mode, "compromise")) {
      c("alpha", "power")
    } else if (!is.null(target)) {
      target
    } else if (!is.null(def)) {
      def$targets[[new_mode]][1]
    } else {
      .ig_registry_default_targets[[new_mode]]
    }
    base_inputs <- base_inputs[setdiff(names(base_inputs), drop)]
  } else if (is.null(target) && !identical(object$mode, "compromise")) {
    target <- object$target
  }
  .ig_build_spec(object$method_id, new_mode, base_inputs, list(), base_options, target)
}

# Serialisation ------------------------------------------------------------------------------

.ig_spec_to_list <- function(x) {
  list(
    schema_version = x$schema_version,
    object_type = x$object_type,
    method_id = x$method_id,
    family = x$family,
    mode = x$mode,
    target = I(x$target),
    joint_target = x$joint_target,
    inputs = .ig_named_or_empty(x$inputs),
    options = .ig_named_or_empty(x$options),
    defaults_applied = list(inputs = I(x$defaults_applied$inputs),
                            options = I(x$defaults_applied$options)),
    validation_level = x$validation_level,
    unvalidated_inputs = I(x$unvalidated_inputs),
    spec_hash = x$spec_hash
  )
}

# Rebuild a power_spec from parsed JSON by re-running the validators.
.ig_spec_from_list <- function(x) {
  need <- c("schema_version", "object_type", "method_id", "mode", "target", "inputs",
            "options", "defaults_applied", "spec_hash")
  miss <- setdiff(need, names(x))
  if (length(miss)) .ig_abort_invalid(sprintf("missing field '%s'", miss), "malformed power_spec")
  if (!identical(x$object_type, "power_spec") ||
      !identical(x$schema_version, .ig_schema_versions[["power_spec"]])) {
    .ig_abort("invalid_input", "not a supported power_spec document")
  }
  as_values <- function(lst, what) {
    if (!is.list(lst) || (length(lst) && is.null(names(lst)))) {
      .ig_abort("invalid_input", sprintf("`%s` must be an object", what))
    }
    out <- lapply(names(lst), function(nm) {
      v <- lst[[nm]]
      if (is.list(v)) {
        v <- tryCatch(.ig_flatten(v), error = function(e) {
          .ig_abort("invalid_input", sprintf("%s '%s' must be a scalar or an array of scalars", what, nm))
        })
      }
      if (is.null(v)) .ig_abort("invalid_input", sprintf("%s '%s' is null", what, nm))
      v
    })
    names(out) <- names(lst)
    out
  }
  inputs <- as_values(x$inputs, "input")
  options <- as_values(x$options, "option")
  d_in <- .ig_flatten(x$defaults_applied$inputs, "character") %||% character()
  d_opt <- .ig_flatten(x$defaults_applied$options, "character") %||% character()
  target <- .ig_flatten(x$target, "character")
  spec <- .ig_build_spec(x$method_id, x$mode,
                         inputs[setdiff(names(inputs), d_in)], list(),
                         options[setdiff(names(options), d_opt)],
                         if (identical(x$mode, "compromise")) NULL else target)
  if (!identical(spec$spec_hash, x$spec_hash)) {
    .ig_abort("invalid_input", paste(
      "spec_hash does not match the specification content (the document was modified,",
      "or the method's defaults differ from those used when it was written)"))
  }
  spec
}

# Display -------------------------------------------------------------------------------------

#' @export
format.power_spec <- function(x, ...) {
  meta <- tryCatch(.ig_get_metadata(x$method_id), error = function(e) NULL)
  show <- function(vals, defaults) {
    if (!length(vals)) return("none")
    paste(vapply(names(vals), function(nm) {
      paste0(nm, " = ", .ig_show_value(vals[[nm]]), if (nm %in% defaults) " (default)" else "")
    }, ""), collapse = ", ")
  }
  level <- if (identical(x$validation_level, "method_definition")) "method definition" else
    "registry metadata only (no engine yet)"
  c(sprintf("<power_spec> %s | %s", x$method_id, x$mode),
    if (!is.null(meta)) sprintf("  method:    %s", meta$label$en),
    sprintf("  solve for: %s%s", paste(x$target, collapse = " + "),
            if (isTRUE(x$joint_target)) " (joint)" else ""),
    sprintf("  inputs:    %s", show(x$inputs, x$defaults_applied$inputs)),
    sprintf("  options:   %s", show(x$options, x$defaults_applied$options)),
    if (length(x$unvalidated_inputs))
      sprintf("  not yet validated: %s", paste(x$unvalidated_inputs, collapse = ", ")),
    sprintf("  validated against: %s", level),
    sprintf("  spec hash: %s", x$spec_hash))
}

#' @export
print.power_spec <- function(x, ...) {
  cat(format(x), sep = "\n")
  invisible(x)
}

# Reproducible R call for a specification (an inert string).
.ig_spec_call <- function(spec) {
  lit_list <- function(vals) {
    if (!length(vals)) return("list()")
    paste0("list(", paste(vapply(names(vals), function(nm) {
      paste0(.ig_r_name(nm), " = ", .ig_r_literal(vals[[nm]]))
    }, ""), collapse = ", "), ")")
  }
  args <- c(sprintf("method = %s", .ig_r_literal(spec$method_id)),
            sprintf("mode = %s", .ig_r_literal(spec$mode)),
            sprintf("inputs = %s", lit_list(spec$inputs)),
            sprintf("options = %s", lit_list(spec$options)))
  if (!identical(spec$mode, "compromise")) {
    args <- c(args, sprintf("target = %s", .ig_r_literal(spec$target)))
  }
  paste0("inferguidR::analysis_spec(", paste(args, collapse = ", "), ")")
}

.ig_r_name <- function(nm) if (make.names(nm) == nm) nm else paste0("`", nm, "`")

.ig_r_literal <- function(v) {
  one <- if (is.character(v)) {
    vapply(v, function(s) encodeString(s, quote = "\""), "", USE.NAMES = FALSE)
  } else if (is.integer(v)) {
    paste0(v, "L")
  } else if (is.double(v)) {
    .ig_format_double(v)
  } else if (is.logical(v)) {
    ifelse(v, "TRUE", "FALSE")
  } else {
    stop("cannot write a literal for type ", typeof(v), call. = FALSE)
  }
  if (length(one) == 1L) one else paste0("c(", paste(one, collapse = ", "), ")")
}

# Immutability ---------------------------------------------------------------------------------

.ig_abort_frozen <- function(x) {
  .ig_abort("invalid_input", sprintf(
    "%s objects are immutable; create a new one (for specifications use update())",
    class(x)[1]))
}

#' Immutable inferguidR objects
#'
#' Specifications, results, envelopes, report snapshots and recommendations
#' cannot be modified in place. Assignment with `$<-`, `[[<-` or `[<-`
#' signals an `invalid_input` error; use [update.power_spec()] or the
#' constructors to create a new object.
#'
#' @param x An immutable inferguidR object.
#' @param name,i,... Ignored.
#' @param value Ignored.
#' @return Never returns; always signals an error.
#' @name ig_immutable
#' @keywords internal
NULL

#' @rdname ig_immutable
#' @export
`$<-.ig_frozen` <- function(x, name, value) .ig_abort_frozen(x)

#' @rdname ig_immutable
#' @export
`[[<-.ig_frozen` <- function(x, i, ..., value) .ig_abort_frozen(x)

#' @rdname ig_immutable
#' @export
`[<-.ig_frozen` <- function(x, i, ..., value) .ig_abort_frozen(x)

# Plain-list copy of a frozen object (for serialisation of auxiliary types).
.ig_unclass_frozen <- function(x) {
  out <- unclass(x)
  attr(out, "class") <- NULL
  out
}
