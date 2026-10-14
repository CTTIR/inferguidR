# X-Y plot data: G*Power's "X-Y plot for a range of values" ------------------------------
#
# Ownership: phase 10 utilities (group G12). Method notes:
# specs/method-notes/utilities.md, section "X-Y plot data".
#
# The 'G*Power' X-Y plot (specs/gpower-utilities.json util_xy_plot; manual
# 1.2.4) plots one of alpha, power, effect size and total sample size (y)
# against another (x), optionally as a family of up to ten curves over a
# third quantity, with the fourth held constant. Every number comes from
# calculate_power() on a complete specification:
#
#   y = power        -> post hoc        y = effect size -> sensitivity
#   y = sample size  -> a priori        y = alpha       -> criterion
#
# The plot inherits the method, the options and every other input from a
# successful base calculation, as G*Power does. Each point keeps its full
# specification and result, its status and reason; points that cannot be
# computed keep their status and have no number (no interpolation across
# gaps). Sample sizes are the smallest admissible integer designs of the
# a priori solver, so a sample-size curve is a staircase; G*Power's
# "interpolating points" values (the continuous relaxation of the a priori
# problem, manual Fig. 3) are reported separately as `n_total_continuous`
# and are never presented as achievable designs. On the x axis only
# admissible integer designs are evaluated.

.ig_xy_roles <- c("power", "sample_size", "effect_size", "alpha")
.ig_xy_mode_for_y <- c(power = "post_hoc", sample_size = "a_priori", effect_size = "sensitivity", alpha = "criterion")

# Effect-size parameter of a method: the default sensitivity target.
.ig_xy_effect_param <- function(def) {
  tg <- def$targets$sensitivity %||% character()
  if (length(tg)) tg[1] else NA_character_
}

.ig_xy_label <- function(def, role) {
  switch(role,
    power = "Power (1-\u03b2 err prob)",
    alpha = "\u03b1 err prob",
    sample_size = "Total sample size",
    effect_size = {
      ep <- .ig_xy_effect_param(def)
      p <- def$params[[ep]]
      if (!is.null(p) && .ig_is_string(p$label$en)) p$label$en else ep
    })
}

# The four quantities of a base result.
.ig_xy_base_values <- function(res, def) {
  spec <- res$spec
  ep <- .ig_xy_effect_param(def)
  pick <- function(nm, alt = NULL) {
    v <- spec$inputs[[nm]] %||% res$solved[[nm]] %||% alt
    if (is.null(v)) NA_real_ else as.double(v)
  }
  gs <- res$group_sizes
  n_total <- res$solved$n_total %||% res$quantities$n_total %||% (if (length(gs)) sum(gs) else NULL)
  ratio <- spec$inputs$allocation_ratio
  if (is.null(ratio) && length(gs) == 2L && all(gs > 0)) ratio <- as.double(gs[2]) / as.double(gs[1])
  list(
    power = pick("power", res$achieved$power),
    alpha = pick("alpha"),
    effect_size = if (is.na(ep)) NA_real_ else pick(ep),
    sample_size = if (is.null(n_total)) NA_real_ else as.double(n_total),
    group_sizes = gs,
    allocation_ratio = if (is.null(ratio)) NULL else as.double(ratio)
  )
}

# Parameters that are inputs of a mode.
.ig_xy_mode_inputs <- function(def, mode) {
  tg <- if (identical(mode, "compromise")) c("alpha", "power") else def$targets[[mode]][1]
  nms <- character()
  for (p in def$params) {
    if (p$name %in% tg) next
    if (is.null(p$modes) || mode %in% p$modes) nms <- c(nms, p$name)
  }
  nms
}

# Inputs of one point: the role values, the sample size (total or group
# sizes through the method's allocation rule) and every other input of the
# base specification. Returns list(inputs) or list(status, reason).
.ig_xy_point_inputs <- function(def, mode, vals, n_from_axis, base_vals, base_spec) {
  want <- .ig_xy_mode_inputs(def, mode)
  ep <- .ig_xy_effect_param(def)
  inputs <- list()
  for (nm in want) {
    if (nm == "alpha") inputs$alpha <- vals$alpha
    else if (nm == "power") inputs$power <- vals$power
    else if (identical(nm, ep)) inputs[[nm]] <- vals$effect_size
    else if (nm == "allocation_ratio" && !is.null(base_vals$allocation_ratio)) inputs$allocation_ratio <- base_vals$allocation_ratio
    else if (!is.null(base_spec$inputs[[nm]])) inputs[[nm]] <- base_spec$inputs[[nm]]
  }
  design_pars <- vapply(def$params, function(p) identical(p$role, "design"), logical(1))
  group_pars <- intersect(names(def$params)[design_pars], setdiff(want, c("allocation_ratio", "n_total")))
  if (mode != "a_priori") {
    if ("n_total" %in% want) {
      inputs$n_total <- as.integer(vals$sample_size)
    } else if (length(group_pars)) {
      if (!n_from_axis && length(base_vals$group_sizes) && all(group_pars %in% names(base_vals$group_sizes))) {
        for (g in group_pars) inputs[[g]] <- as.integer(base_vals$group_sizes[[g]])
      } else {
        if (!is.function(def$design)) {
          return(list(status = "unsupported", reason = "the method has no allocation rule to split a total sample size"))
        }
        p <- c(base_spec$inputs, list(allocation_ratio = base_vals$allocation_ratio %||% 1))
        split <- tryCatch(def$design(p, vals$sample_size), inferguidR_error = function(e) e)
        if (inherits(split, "inferguidR_error")) return(list(status = split$status$code, reason = conditionMessage(split)))
        if (is.null(split)) {
          return(list(status = "invalid_input", reason = sprintf(
            "total sample size %s is not an admissible design under allocation ratio N2/N1 = %s (the method's allocation rule has no integer split for it)",
            .ig_format_double(vals$sample_size), .ig_format_double(base_vals$allocation_ratio %||% 1))))
        }
        for (g in names(split)) inputs[[g]] <- as.integer(split[[g]])
      }
    }
  }
  missing <- setdiff(want, names(inputs))
  missing <- missing[vapply(missing, function(m) is.null(def$params[[m]]$default), logical(1))]
  if (length(missing)) {
    return(list(status = "unsupported", reason = sprintf("the X-Y plot cannot set input(s) %s of mode %s from the base analysis",
                                                         paste(missing, collapse = ", "), mode)))
  }
  list(inputs = inputs)
}

.ig_xy_seq <- function(from, to, by) {
  k <- floor((to - from) / by * (1 + 1e-12) + 1e-9)
  from + by * seq.int(0, k)
}

#' Data of an X-Y plot for a range of values
#'
#' Computes the data of the 'G*Power' "X-Y plot for a range of values": one
#' of power, total sample size, effect size and alpha (`y`) as a function of
#' another (`x`), optionally as a family of curves over a third quantity,
#' with the fourth held constant at its value in the base analysis (or at
#' `fixed`). Every point is a complete specification solved with
#' [calculate_power()]: `y = "power"` uses post hoc mode, `"sample_size"`
#' a priori, `"effect_size"` sensitivity and `"alpha"` criterion mode. The
#' options and all other inputs come from the base analysis.
#'
#' Sample sizes are integer designs. A sample-size curve is therefore a
#' staircase (`render = "step"`); the fractional totals of the continuous
#' relaxation (what 'G*Power' tabulates with "interpolating points") are in
#' column `n_total_continuous` and are not achievable designs. With the total
#' sample size on the x axis only whole numbers are allowed, and totals
#' without an admissible split into groups (for example an odd total under
#' equal allocation) are kept as points with status `invalid_input`.
#'
#' @param base A `power_result` with status `"ok"` or a `power_spec` (it is
#'   calculated first).
#' @param y,x Quantities on the axes: `"power"`, `"sample_size"`,
#'   `"effect_size"` or `"alpha"`; they must differ.
#' @param from,to,by Range of `x` ("from", "through to", "in steps of").
#' @param curves Optional family of curves: a list with `parameter` (one of
#'   the two remaining quantities) and either `values` (1 to 10 values) or
#'   `from`, `by` and `n` (1 to 10 curves).
#' @param fixed Optional named list with a value for the quantity that is
#'   held constant (default: its value in the base analysis).
#' @param max_points Largest number of points per curve.
#' @return An object of class `ig_xy_plot_data`: a list with `points` (a
#'   data frame with one row per point: `curve`, `curve_value`, `x`, `y`,
#'   `status`, `reason`, `n_total`, `n_total_continuous`, `achieved_power`,
#'   `spec_hash`, `result_hash` and the list columns `group_sizes` and
#'   `result`), `axes`, `curves`, `fixed` (every input and option held
#'   constant, with labels), `method_id`, `mode` (mode used for y),
#'   `title`, `subtitle`, `render` (`"step"` for sample-size curves,
#'   `"line"` otherwise), `validation` and `notes`. Invalid ranges signal an
#'   `inferguidR_invalid_input` condition.
#' @examples
#' base <- calculate_power(analysis_spec("t_means_two_independent", "post_hoc",
#'   inputs = list(d = 0.5, alpha = 0.05, n1 = 50, n2 = 50), options = list(tails = "two")))
#' xy <- xy_plot_data(base, y = "power", x = "sample_size", from = 20, to = 200, by = 20)
#' xy$points[, c("x", "y", "status")]
#' xy2 <- xy_plot_data(base, y = "sample_size", x = "power", from = 0.6, to = 0.95, by = 0.05,
#'   curves = list(parameter = "effect_size", values = c(0.4, 0.5, 0.6)))
#' head(xy2$points[, c("curve_value", "x", "y", "n_total_continuous")])
#' @seealso [sweep_design()], [compare_designs()]
#' @export
xy_plot_data <- function(base, y, x, from, to, by, curves = NULL, fixed = list(), max_points = 2000L) {
  if (inherits(base, "power_spec")) base <- calculate_power(base)
  if (!inherits(base, "power_result")) .ig_abort_invalid("`base` must be a power_result or a power_spec")
  if (!identical(base$status$code, "ok")) {
    .ig_abort_invalid(sprintf("the base analysis has no result (%s); the X-Y plot needs a successful calculation",
                              format(base$status)))
  }
  if (!.ig_is_string(y) || !(y %in% .ig_xy_roles)) .ig_abort_invalid(sprintf("`y` must be one of %s", paste(.ig_xy_roles, collapse = ", ")))
  if (!.ig_is_string(x) || !(x %in% .ig_xy_roles)) .ig_abort_invalid(sprintf("`x` must be one of %s", paste(.ig_xy_roles, collapse = ", ")))
  if (identical(x, y)) .ig_abort_invalid("`x` and `y` must be different quantities")
  def <- .ig_get_method(base$method_id)
  if (is.null(def)) .ig_abort("unsupported", sprintf("no engine for %s", base$method_id))
  meta <- .ig_get_metadata(base$method_id)
  mode <- unname(.ig_xy_mode_for_y[y])
  if (!(mode %in% def$modes)) {
    .ig_abort("unsupported", sprintf("%s offers no %s mode, which a plot of %s needs", base$method_id, mode, y))
  }
  if (is.na(.ig_xy_effect_param(def))) .ig_abort("unsupported", "the method declares no effect-size parameter")
  num1 <- function(v, nm) {
    if (!is.numeric(v) || length(v) != 1L || !is.finite(v)) .ig_abort_invalid(sprintf("`%s` must be one finite number", nm))
    as.double(v)
  }
  from <- num1(from, "from")
  to <- num1(to, "to")
  by <- num1(by, "by")
  problems <- character()
  # Range rules of the G*Power X-Y plot (messages captured live:
  # GPW-UT-XYPLOT-UTL-U026-I01..I04, SUP-G12-UT-XYPLOT-X03, -X05, -X08).
  if (!(to > from)) problems <- c(problems, "the end value of the x range must be greater than the start value")
  if (!(by > 0)) problems <- c(problems, "the step size of the x variable must be a positive value")
  if (x == "sample_size" && (from != round(from) || by != round(by) || by < 1)) {
    problems <- c(problems, "the x variable total sample size needs a whole-number start value and a whole-number step of at least 1")
  }
  if (length(problems)) .ig_abort_invalid(problems, "invalid X-Y plot range")
  xs <- .ig_xy_seq(from, to, by)
  if (length(xs) > max_points) .ig_abort_invalid(sprintf("the x range has %d points (more than max_points = %d)", length(xs), max_points))

  base_vals <- .ig_xy_base_values(base, def)
  remaining <- setdiff(.ig_xy_roles, c(x, y))
  curve_par <- NULL
  curve_vals <- NA_real_
  if (!is.null(curves)) {
    if (!is.list(curves) || !.ig_is_string(curves$parameter) || !(curves$parameter %in% remaining)) {
      .ig_abort_invalid(sprintf("`curves$parameter` must be one of %s", paste(remaining, collapse = ", ")))
    }
    curve_par <- curves$parameter
    curve_vals <- if (!is.null(curves$values)) curves$values else {
      n <- curves$n
      if (!is.numeric(n) || length(n) != 1L || !is.finite(n) || n != round(n) || n < 1 || n > 10) {
        .ig_abort_invalid("`curves$n` must be a whole number from 1 to 10")
      }
      num1(curves$from, "curves$from") + num1(curves$by, "curves$by") * (seq_len(n) - 1)
    }
    if (!is.numeric(curve_vals) || !length(curve_vals) || length(curve_vals) > 10L || any(!is.finite(curve_vals))) {
      .ig_abort_invalid("a family of curves has 1 to 10 finite values")
    }
    if (curve_par == "sample_size" && any(curve_vals != round(curve_vals))) {
      .ig_abort_invalid("curve values of the total sample size must be whole numbers")
    }
  }
  const_par <- setdiff(remaining, curve_par)
  if (is.null(curve_par)) const_par <- remaining
  vals0 <- base_vals[.ig_xy_roles]
  if (length(fixed)) {
    if (!.ig_is_named(fixed) || !all(names(fixed) %in% const_par)) {
      .ig_abort_invalid(sprintf("`fixed` may only set the constant quantit%s %s", if (length(const_par) > 1L) "ies" else "y",
                                paste(const_par, collapse = ", ")))
    }
    for (nm in names(fixed)) vals0[[nm]] <- num1(fixed[[nm]], paste0("fixed$", nm))
  }
  # Probability ranges and the power >= alpha rule of G*Power.
  check_prob <- function(v, nm) {
    if (any(v <= 0 | v >= 1)) problems <<- c(problems, sprintf("%s values must lie strictly between 0 and 1", nm))
  }
  if (x %in% c("alpha", "power")) check_prob(xs, x)
  if (!is.null(curve_par) && curve_par %in% c("alpha", "power")) check_prob(curve_vals, curve_par)
  for (nm in intersect(const_par, c("alpha", "power"))) check_prob(vals0[[nm]], nm)
  if (x == "power" && "alpha" %in% const_par && from < vals0$alpha) {
    problems <- c(problems, sprintf("the start value of the x variable power (%s) must be at least the constant alpha (%s)",
                                    .ig_format_double(from), .ig_format_double(vals0$alpha)))
  }
  if (x == "sample_size" && from < 1) problems <- c(problems, "total sample sizes must be at least 1")
  if (length(problems)) .ig_abort_invalid(unique(problems), "invalid X-Y plot range")

  rows <- list()
  for (ci in seq_along(curve_vals)) {
    for (xv in xs) {
      vals <- vals0
      vals[[x]] <- xv
      if (!is.null(curve_par)) vals[[curve_par]] <- curve_vals[ci]
      n_axis <- identical(x, "sample_size") || identical(curve_par, "sample_size") || ("sample_size" %in% names(fixed))
      pin <- .ig_xy_point_inputs(def, mode, vals, n_axis, base_vals, base$spec)
      row <- list(curve = ci, curve_value = if (is.null(curve_par)) NA_real_ else curve_vals[ci], x = xv,
                  y = NA_real_, status = NA_character_, reason = "", n_total = NA_real_, n_total_continuous = NA_real_,
                  achieved_power = NA_real_, spec_hash = NA_character_, result_hash = NA_character_,
                  group_sizes = NULL, result = NULL)
      if (is.null(pin$inputs)) {
        row$status <- pin$status
        row$reason <- pin$reason
        rows[[length(rows) + 1L]] <- row
        next
      }
      spec <- tryCatch(analysis_spec(base$method_id, mode, inputs = pin$inputs, options = base$spec$options),
                       inferguidR_error = function(e) e)
      if (inherits(spec, "inferguidR_error")) {
        row$status <- spec$status$code
        row$reason <- conditionMessage(spec)
        rows[[length(rows) + 1L]] <- row
        next
      }
      res <- calculate_power(spec)
      row$status <- res$status$code
      row$reason <- res$status$reason
      row$spec_hash <- res$spec_hash
      row$result_hash <- res$result_hash
      row$result <- res
      if (identical(res$status$code, "ok")) {
        ep <- .ig_xy_effect_param(def)
        row$y <- switch(y,
          power = res$achieved$power,
          sample_size = res$solved$n_total,
          effect_size = res$solved[[ep]],
          alpha = res$solved$alpha)
        row$n_total <- as.double(res$solved$n_total %||% res$quantities$n_total %||% NA_real_)
        cr <- res$diagnostics$continuous_relaxation$n_total
        if (identical(y, "sample_size") && is.numeric(cr) && length(cr) == 1L) row$n_total_continuous <- cr
        row$achieved_power <- res$achieved$power %||% NA_real_
        row$group_sizes <- res$group_sizes
      }
      rows[[length(rows) + 1L]] <- row
    }
  }
  points <- data.frame(
    curve = vapply(rows, function(r) as.integer(r$curve), integer(1)),
    curve_value = vapply(rows, `[[`, numeric(1), "curve_value"),
    x = vapply(rows, `[[`, numeric(1), "x"),
    y = vapply(rows, function(r) as.double(r$y), numeric(1)),
    status = vapply(rows, `[[`, "", "status"),
    reason = vapply(rows, `[[`, "", "reason"),
    n_total = vapply(rows, `[[`, numeric(1), "n_total"),
    n_total_continuous = vapply(rows, `[[`, numeric(1), "n_total_continuous"),
    achieved_power = vapply(rows, `[[`, numeric(1), "achieved_power"),
    spec_hash = vapply(rows, `[[`, "", "spec_hash"),
    result_hash = vapply(rows, `[[`, "", "result_hash"),
    stringsAsFactors = FALSE
  )
  points$group_sizes <- I(lapply(rows, `[[`, "group_sizes"))
  points$result <- I(lapply(rows, `[[`, "result"))

  # Constant inputs and options, with labels (the subtitle of G*Power's plot).
  fixed_out <- list()
  for (nm in const_par) fixed_out[[.ig_xy_label(def, nm)]] <- if (nm == "sample_size" && !("sample_size" %in% names(fixed)) &&
                                                                     length(base_vals$group_sizes)) {
    paste(sprintf("%s = %d", names(base_vals$group_sizes), base_vals$group_sizes), collapse = ", ")
  } else vals0[[nm]]
  for (o in names(base$spec$options)) {
    lab <- meta$options[[o]]$ui_label %||% o
    code <- base$spec$options[[o]]
    fixed_out[[lab]] <- meta$options[[o]]$value_live[[as.character(code)]]$text %||% code
  }
  if (!is.null(base_vals$allocation_ratio) && "allocation_ratio" %in% .ig_xy_mode_inputs(def, mode)) {
    fixed_out[["Allocation ratio N2/N1"]] <- base_vals$allocation_ratio
  }
  show <- function(v) if (is.numeric(v)) .ig_format_double(v) else as.character(v)
  subtitle <- paste(sprintf("%s = %s", names(fixed_out), vapply(fixed_out, show, "")), collapse = ", ")
  family_text <- unname(.ig_protocol_family_text[meta$family])
  title <- paste0(if (is.na(family_text)) meta$family else family_text, " - ", meta$label_raw$windows_3197_live %||% meta$label$en)
  notes <- c(
    if (y == "sample_size") paste("y is the smallest admissible integer total sample size of each a priori problem (a staircase);",
                                   "n_total_continuous is the continuous relaxation (G*Power 'interpolating points') and is not an achievable design"),
    if (x == "sample_size") "only whole-number totals with an admissible split into groups are evaluated on the x axis",
    if (any(points$status != "ok")) sprintf("%d point(s) without a result keep their status and reason (no interpolation)", sum(points$status != "ok"))
  )
  structure(list(
    schema_version = "1.0.0",
    object_type = "xy_plot_data",
    method_id = base$method_id,
    base_mode = base$mode,
    mode = mode,
    points = points,
    axes = list(x = list(quantity = x, label = .ig_xy_label(def, x), from = from, to = to, by = by),
                y = list(quantity = y, label = .ig_xy_label(def, y))),
    curves = if (is.null(curve_par)) NULL else list(quantity = curve_par, label = .ig_xy_label(def, curve_par), values = curve_vals),
    fixed = fixed_out,
    title = title,
    subtitle = subtitle,
    render = if (y == "sample_size") "step" else "line",
    validation = base$validation,
    base_spec_hash = base$spec_hash,
    base_result_hash = base$result_hash,
    notes = notes
  ), class = "ig_xy_plot_data")
}

#' @export
format.ig_xy_plot_data <- function(x, ...) {
  ok <- x$points$status == "ok"
  c(sprintf("<ig_xy_plot_data> %s", x$title),
    sprintf("  %s as a function of %s (%s mode), %d point(s), %d without a result", x$axes$y$label, x$axes$x$label,
            x$mode, nrow(x$points), sum(!ok)),
    if (!is.null(x$curves)) sprintf("  curves: %s = %s", x$curves$label, paste(.ig_format_double(x$curves$values), collapse = ", ")),
    sprintf("  constant: %s", x$subtitle),
    if (length(x$notes)) paste0("  note: ", x$notes))
}

#' @export
print.ig_xy_plot_data <- function(x, ...) {
  cat(format(x), sep = "\n")
  invisible(x)
}
