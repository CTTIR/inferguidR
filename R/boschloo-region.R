
#' Power of Boschloo's test for a fixed design
#'
#' Computes the exact power at the alternative `(p1, p2)` as the double
#' binomial sum over the rejection region of [boschloo_region()]: the
#' probability that a future study with group sizes `n1` and `n2` produces a
#' table that Boschloo's test rejects. No asymptotic approximation is used.
#'
#' @inheritParams boschloo_region
#' @param p1,p2 Event probabilities of group 1 and group 2 under the
#'   alternative (in `[0, 1]`). With `p1 = p2` the result is a null
#'   rejection probability, at most the attained size.
#' @return An immutable object of class `boschloo_power` with the power, the
#'   attained size, the number of rejected tables, the region hash and the
#'   inputs.
#' @examples
#' boschloo_power(20, 20, p1 = 0.2, p2 = 0.6, alpha = 0.05, variant = "central")
#' @family Boschloo
#' @export
boschloo_power <- function(n1, n2, p1, p2, alpha = 0.05, variant, control = boschloo_control()) {
  info <- .bl_variant(if (missing(variant)) NULL else variant)
  .bl_check_design(n1, n2)
  alpha <- .bl_check_alpha(alpha)
  p1 <- .bl_check_prob(p1, "p1")
  p2 <- .bl_check_prob(p2, "p2")
  control <- .bl_check_control(control)
  reg <- .bl_region(n1, n2, alpha, info, control)
  .bl_power_object(reg, p1, p2, control)
}

.bl_power_object <- function(reg, p1, p2, control) {
  region <- .bl_region_object(reg, control)
  ok <- identical(region$status$code, "ok")
  pw <- if (ok) .bl_power_idx(reg$idx, reg$n1, reg$n2, p1, p2) else NA_real_
  body <- list(
    object_type = "boschloo_power",
    status = region$status,
    variant = region$variant,
    design = region$design,
    alpha = reg$alpha,
    alternative_probabilities = list(p1 = p1, p2 = p2, difference = p2 - p1),
    power = pw,
    beta = if (ok) 1 - pw else NA_real_,
    size = region$size,
    n_rejected = region$n_rejected,
    region_hash = region$region_hash,
    warnings = region$warnings,
    undefined = if (ok) .ig_empty_object() else list(power = region$status$reason),
    settings = region$settings,
    scope = "hypothetical future study with fixed group sizes n1 and n2",
    provenance = region$provenance,
    r_call = sprintf("inferguidR::boschloo_power(n1 = %dL, n2 = %dL, p1 = %s, p2 = %s, alpha = %s, variant = %s%s)",
                     reg$n1, reg$n2, .ig_r_literal(p1), .ig_r_literal(p2), .ig_r_literal(reg$alpha),
                     .ig_r_literal(reg$info$variant),
                     if (is.null(.bl_control_call(control))) "" else paste0(", ", .bl_control_call(control)))
  )
  structure(body, class = c("boschloo_power", "ig_frozen"))
}

#' @export
format.boschloo_power <- function(x, ...) {
  c(sprintf("<boschloo_power> %s", x$variant$label),
    sprintf("  design:      n1 = %d, n2 = %d; alpha = %s; p1 = %s, p2 = %s",
            as.integer(x$design$n1), as.integer(x$design$n2), format(x$alpha),
            format(x$alternative_probabilities$p1), format(x$alternative_probabilities$p2)),
    sprintf("  status:      %s", format(.ig_status_from_list(x$status))),
    if (identical(x$status$code, "ok")) c(
      sprintf("  power:       %s", .bl_fmt(x$power, 7)),
      sprintf("  attained size: %s (upper bound %s); %d tables rejected", .bl_fmt(x$size$value, 7),
              .bl_fmt(x$size$upper, 7), as.integer(x$n_rejected))),
    if (length(x$warnings)) paste0("  warning: ", x$warnings))
}

#' @export
print.boschloo_power <- function(x, ...) {
  cat(format(x), sep = "\n")
  invisible(x)
}
