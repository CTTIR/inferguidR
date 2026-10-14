# Derive the editor from the same definitions used by analysis_spec().
# Context changes may retire old controls; every such reset is returned to UI.
.ig_workbench_form <- function(method, mode = NULL, options = list(), target = NULL) {
  meta <- .ig_get_metadata(method)
  def <- .ig_get_method(method)
  if (is.null(def)) .ig_abort("unsupported", "this method has no input definition")
  mode <- mode %||% if ("post_hoc" %in% def$modes) "post_hoc" else def$modes[1]
  if (!.ig_is_string(mode) || !mode %in% def$modes)
    .ig_abort("invalid_input", "the selected method does not offer this calculation mode")
  specs <- .ig_definition_option_specs(def)
  if (length(setdiff(names(options), names(specs))))
    .ig_abort("invalid_input", "unknown method option in the editor")
  for (nm in names(options)) {
    if (!.ig_is_string(options[[nm]]) || !options[[nm]] %in% specs[[nm]]$codes)
      .ig_abort("invalid_input", paste("invalid option value:", nm))
  }
  effective <- .ig_validate_options(options, specs, mode, method)$values
  dropped <- setdiff(names(options), names(effective))
  resets <- list()
  choices <- list()
  for (iteration in seq_len(length(specs) + 2L)) {
    old <- effective
    for (nm in names(effective)) {
      codes <- specs[[nm]]$codes
      conditions <- meta$options[[nm]]$value_applies_when
      keep <- vapply(codes, function(code) {
        cond <- conditions[[code]]
        is.null(cond) || .ig_eval_condition(cond$condition %||% cond, mode, effective)
      }, logical(1))
      codes <- codes[keep]
      if (!length(codes)) .ig_abort("unsupported", paste("no applicable option value:", nm))
      choices[[nm]] <- codes
      if (!effective[[nm]] %in% codes) {
        replacement <- if (isTRUE(specs[[nm]]$default %in% codes)) specs[[nm]]$default else codes[1]
        resets[[nm]] <- list(from = effective[[nm]], to = replacement)
        effective[[nm]] <- replacement
      }
    }
    effective <- .ig_validate_options(effective, specs, mode, method)$values
    if (identical(old, effective)) break
  }
  checked <- .ig_validate_options(effective, specs, mode, method)
  effective <- checked$values
  problems <- checked$problems
  choices <- list()
  for (nm in names(specs)) if (.ig_eval_condition(specs[[nm]]$condition, mode, effective)) {
    codes <- specs[[nm]]$codes
    conditions <- meta$options[[nm]]$value_applies_when
    choices[[nm]] <- codes[vapply(codes, function(code) {
      cond <- conditions[[code]]
      is.null(cond) || .ig_eval_condition(cond$condition %||% cond, mode, effective)
    }, logical(1))]
  }
  eligible <- if (mode == "compromise") c("alpha", "power") else def$targets[[mode]]
  if (mode != "compromise") eligible <- eligible[vapply(eligible, function(nm) {
    p <- def$params[[nm]]
    is.null(p) || .ig_eval_condition(p$applies_when, mode, effective)
  }, logical(1))]
  reset_target <- !is.null(target) && (if (mode == "compromise")
    !setequal(target, c("alpha", "power")) else !.ig_is_string(target) || !target %in% eligible)
  selected <- .ig_select_target(def, meta, mode, if (reset_target) NULL else target, effective)
  problems <- c(problems, selected$problems)
  params <- Filter(function(p) {
    !p$name %in% selected$target && (is.null(p$modes) || mode %in% p$modes) &&
      .ig_eval_condition(p$applies_when, mode, effective)
  }, def$params)
  dropped <- union(dropped, setdiff(names(options), names(effective)))
  list(method = method, definition = def, mode = mode, modes = def$modes,
       target = selected$target, targets = eligible, parameters = params,
       options = effective, option_choices = choices,
       dropped_options = dropped, reset_options = resets, reset_target = reset_target,
       ready = !length(problems), problems = problems)
}

.ig_workbench_field_id <- function(method, name, kind = "field") {
  paste(kind, method, name, sep = "_")
}

.ig_workbench_label <- function(labels, locale, fallback = "") {
  value <- labels[[locale]]
  if (!.ig_is_string(value)) value <- labels$en
  if (!.ig_is_string(value)) value <- fallback
  value
}

.ig_workbench_domain_hint <- function(p, locale) {
  de <- locale == "de"
  d <- p$domain
  parts <- if (d$type == "integer") {
    if (de) "Ganze Zahl" else "Whole number"
  } else if (de) "Endliche Zahl" else "Finite number"
  show <- function(x) {
    s <- .ig_format_double(x)
    if (de) sub(".", ",", s, fixed = TRUE) else s
  }
  if (is.finite(d$lower)) parts <- c(parts, paste(if (d$lower_inclusive) ">=" else ">", show(d$lower)))
  if (is.finite(d$upper)) parts <- c(parts, paste(if (d$upper_inclusive) "<=" else "<", show(d$upper)))
  if (!is.na(p$unit)) parts <- c(parts, p$unit)
  if (!is.null(p$default)) parts <- c(parts, paste(if (de) "Leer verwendet den Standardwert" else "Blank uses default", show(p$default)))
  paste(parts, collapse = "; ")
}

.ig_workbench_form_spec <- function(form, fields, locale = "en") {
  if (!isTRUE(form$ready)) .ig_abort_invalid(form$problems, "choose the required model settings")
  if (length(setdiff(names(fields), names(form$parameters))))
    .ig_abort("invalid_input", "the editor contains inactive inputs")
  if (any(vapply(fields, function(x) !is.null(x) &&
                !(is.character(x) && length(x) == 1L && !is.na(x)), logical(1))))
    .ig_abort("invalid_input", "editor entries must be scalar text")
  # Empty entries request a declared engine default; they never become zero.
  fields <- fields[vapply(fields, function(x) !is.null(x) && .ig_is_string(x) && nzchar(trimws(x)), logical(1))]
  values <- lapply(fields, .ig_workbench_parse_number, locale = locale)
  analysis_spec(form$method, form$mode, inputs = values, options = form$options, target = form$target)
}
