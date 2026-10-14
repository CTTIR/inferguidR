# Observed Boschloo test ---------------------------------------------------------------
#
# boschloo_test() analyses one observed 2 x 2 table of two independent
# binomial samples. The result is an `observed_test_result` that pairs with
# its `observed_test_spec` in an analysis_envelope(kind = "observed_test").
# It refers to observed data only; planning quantities (power, sample size)
# come from boschloo_power() and boschloo_sample_size(), which use the same
# ordering, tie rule and nuisance maximisation.

.bl_result_schema <- "inferguidR.observed_test_result/0.1.0-draft"

#' Boschloo's unconditional exact test for two independent binomial samples
#'
#' Computes the Boschloo p-value of an observed table. Group 1 is `(x1, n1)`,
#' group 2 is `(x2, n2)`; `x1` and `x2` count the events (the coding is
#' recorded through `groups` and `event`). The null hypothesis is
#' `p1 = p2`. Tables are ordered by a Fisher conditional p-value, and the
#' p-value is the largest null probability, over the common event
#' probability, of all tables that are at least as extreme as the observed
#' one (ties included).
#'
#' Variants (always stated explicitly):
#'
#' * `"minlike"`: two-sided, ordered by the two-sided Fisher p-value of
#'   [stats::fisher.test()] (sum of conditional probabilities not larger than
#'   the observed one). Same definition as `Exact::exact.test(...,
#'   method = "boschloo", tsmethod = "square")`.
#' * `"central"`: two-sided, `min(1, 2 * min(p_less, p_greater))`, each
#'   one-sided p-value maximised separately.
#' * `"less"`: one-sided, alternative `p1 < p2`.
#' * `"greater"`: one-sided, alternative `p1 > p2`.
#'
#' Mid-p values, randomisation, Berger-Boos confidence-set restrictions and
#' non-zero null margins are not part of this test.
#'
#' The p-value is reported with a rigorous upper bound (`p_upper`, see
#' [boschloo_control()]). With `alpha`, the decision is `"reject"` when
#' `p_upper <= alpha`, `"do_not_reject"` when `p_value > alpha` and
#' `"unresolved"` otherwise (which requires the supremum to lie within the
#' floating-point allowance of `alpha`). The same rule defines the rejection
#' regions of [boschloo_region()], so an observed table is rejected exactly
#' when it lies in the region of its design.
#'
#' @param x1,x2 Numbers of events in group 1 and group 2.
#' @param n1,n2 Group sizes (fixed by design).
#' @param variant `"minlike"`, `"central"`, `"less"` or `"greater"`.
#' @param alpha Optional significance level for a decision.
#' @param groups Labels of group 1 and group 2 (recorded, not used in the
#'   computation).
#' @param event Label of the counted outcome level.
#' @param control Numerical settings from [boschloo_control()].
#' @return An immutable object of classes `boschloo_test` and
#'   `observed_test_result` with the p-value, its upper bound, the
#'   maximising nuisance probability, the ordering statistic, the Fisher
#'   p-value of the same sidedness convention, the observed specification
#'   (`$spec`), diagnostics and a reproducible call. Wrap it with
#'   `analysis_envelope(res$spec, res)`.
#' @examples
#' res <- boschloo_test(3, 15, 9, 15, variant = "minlike", alpha = 0.05)
#' res
#' res$p_value
#' boschloo_test(1, 20, 8, 20, variant = "less")$p_value
#' @family Boschloo
#' @export
boschloo_test <- function(x1, n1, x2, n2, variant, alpha = NULL,
                          groups = c("group 1", "group 2"), event = "event",
                          control = boschloo_control()) {
  info <- .bl_variant(if (missing(variant)) NULL else variant)
  .bl_check_design(n1, n2)
  p <- character()
  if (!.bl_is_count(x1) || x1 > n1) p <- c(p, "`x1` must be a whole number between 0 and n1")
  if (!.bl_is_count(x2) || x2 > n2) p <- c(p, "`x2` must be a whole number between 0 and n2")
  if (!is.character(groups) || length(groups) != 2L || anyNA(groups) || !all(nzchar(groups)) ||
      groups[1] == groups[2]) {
    p <- c(p, "`groups` must be two different non-empty labels")
  }
  if (!.ig_is_string(event)) p <- c(p, "`event` must be a non-empty label")
  if (length(p)) .ig_abort_invalid(p, "invalid input for boschloo_test()")
  if (!is.null(alpha)) alpha <- .bl_check_alpha(alpha)
  control <- .bl_check_control(control)
  x1 <- as.integer(x1); x2 <- as.integer(x2); n1 <- as.integer(n1); n2 <- as.integer(n2)

  inputs <- list(x1 = x1, n1 = n1, x2 = x2, n2 = n2, group_labels = I(groups),
                 event_label = event)
  if (!is.null(alpha)) inputs$alpha <- alpha
  spec <- .ig_new_observed_test_spec(info$variant_id, inputs = inputs,
                                     options = c(list(variant = info$variant),
                                                 .bl_control_record(control)))
  i <- x1 + 1L + x2 * (n1 + 1L)
  alpha_side <- if (is.null(alpha)) NULL else if (info$variant == "central") alpha / 2 else alpha
  sides <- lapply(info$sides, function(side) .bl_side_pvalue(n1, n2, side, i, control, alpha_side))
  names(sides) <- info$sides

  # Combine the sides.
  if (info$variant == "central") {
    p_value <- min(1, 2 * min(sides$lower$sup$value, sides$upper$sup$value))
    p_upper <- min(1, 2 * min(sides$lower$sup$upper, sides$upper$sup$upper))
    fisher_p <- min(1, 2 * min(sides$lower$ordering$value, sides$upper$ordering$value))
    smaller <- if (sides$lower$sup$value <= sides$upper$sup$value) "less" else "greater"
  } else {
    p_value <- sides[[1]]$sup$value
    p_upper <- sides[[1]]$sup$upper
    fisher_p <- sides[[1]]$ordering$value
    smaller <- NULL
  }

  problems <- character()
  for (s in sides) {
    if (isTRUE(s$ambiguity_touches)) {
      problems <- c(problems, sprintf(
        "the %s ordering has an ambiguous tie at the observed table (relative gap in the ambiguity band) and n1 + n2 = %d exceeds %d, the limit of exact integer arithmetic",
        s$side, n1 + n2, .bl_exact_max_n))
    }
    if (!isTRUE(s$sup$converged)) {
      problems <- c(problems, sprintf(
        "the nuisance supremum of the %s tail did not converge within %d bisection levels",
        s$side, control$max_levels))
    }
  }
  status <- if (length(problems)) {
    ig_status("numerical_failure", paste(problems, collapse = "; "))
  } else {
    ig_status("ok")
  }
  ok <- status$code == "ok"
  decision <- NULL
  if (!is.null(alpha) && ok) {
    decision <- if (p_upper <= alpha) "reject" else if (p_value > alpha) "do_not_reject" else "unresolved"
  }
  undefined <- list()
  if (!ok) {
    # No finite number is reported without an ok status.
    undefined <- list(p_value = status$reason, p_upper = status$reason,
                      fisher_p_value = status$reason)
    p_value <- NA_real_
    p_upper <- NA_real_
  }

  res <- list(
    schema_version = .bl_result_schema,
    object_type = "observed_test_result",
    method = "boschloo",
    status = .ig_status_to_list(status),
    spec = spec,
    spec_hash = spec$spec_hash,
    variant = .bl_variant_record(info),
    data = list(x1 = x1, n1 = n1, x2 = x2, n2 = n2,
                proportion_1 = x1 / n1, proportion_2 = x2 / n2,
                group_labels = I(groups), event_label = event,
                orientation = sprintf("group 1 = %s (x1 = %d of n1 = %d), group 2 = %s (x2 = %d of n2 = %d); events = %s",
                                      groups[1], x1, n1, groups[2], x2, n2, event)),
    p_value = p_value,
    p_upper = p_upper,
    alpha = alpha,
    decision = decision,
    smaller_one_sided = smaller,
    fisher_p_value = if (ok) fisher_p else NA_real_,
    fisher_note = paste("Fisher conditional p-value with the same sidedness convention",
                        "(the Boschloo p-value never exceeds it); not the tested statistic."),
    sides = lapply(sides, .bl_side_record, ok = ok),
    undefined = .ig_named_or_empty(undefined),
    settings = .bl_control_record(control),
    scope = "observed data: p-value of the observed table; no planning quantity",
    provenance = .bl_provenance(),
    r_call = .bl_test_call(x1, n1, x2, n2, info$variant, alpha, groups, event, control)
  )
  structure(res, class = c("boschloo_test", "observed_test_result", "ig_frozen"))
}

# p-value of one side for table index i.
.bl_side_pvalue <- function(n1, n2, side, i, control, alpha = NULL) {
  ord <- .bl_ordering(n1, n2, side, control)
  j <- ord$level[i]
  idx <- .bl_level_tables(ord, j)
  s <- .bl_sup(.bl_tail_c(ord, idx), control, alpha)
  list(side = side, sup = s, level = j, tail_tables = length(idx),
       ordering = .bl_ordering_value(ord, i), arithmetic = ord$arithmetic,
       fallback_reason = ord$fallback_reason, ambiguous_gaps = ord$ambiguous_gaps,
       ambiguity_touches = .bl_ambiguity_touches(ord, j, "tail"))
}

.bl_side_record <- function(s, ok = TRUE) {
  sup <- .bl_sup_record(s$sup)
  num <- function(v) if (ok) v else NA_real_
  if (!ok) {
    sup$value <- NA_real_
    sup$upper <- NA_real_
    sup$gap <- NA_real_
    sup$argmax <- NA_real_
  }
  list(side = s$side,
       direction = if (s$side == "lower") "p1 < p2" else if (s$side == "upper") "p1 > p2" else "p1 != p2",
       p_value = sup$value, p_upper = sup$upper, nuisance_argmax = sup$argmax,
       ordering_statistic = num(s$ordering$value), ordering_complement = num(s$ordering$complement),
       level = s$level, n_levels = s$ordering$n_levels, tail_tables = s$tail_tables,
       ordering_arithmetic = s$arithmetic, ambiguous_gaps = s$ambiguous_gaps,
       ambiguity_touches_result = isTRUE(s$ambiguity_touches),
       fallback_reason = s$fallback_reason, supremum = sup)
}

.bl_variant_record <- function(info) {
  list(variant = info$variant, variant_id = info$variant_id,
       advisor_variant_id = info$advisor_variant_id, label = info$label,
       alternative = info$alternative, null = "p1 = p2",
       design = paste("two independent binomial samples with fixed group sizes;",
                      "not for paired or matched binary data (McNemar-type tests)"),
       ordering_statistic = info$ordering_statistic, construction = info$construction,
       ties = "tables whose ordering statistic ties the observed one (relative tolerance) are in the tail",
       not_included = "no mid-p, no randomisation, no Berger-Boos restriction, no null margin",
       validation_status = .bl_evidence$validation_status,
       planning_compatibility = .bl_evidence$planning_compatibility,
       gpower_relation = sprintf("%s (%s)", .bl_evidence$gpower_relation,
                                 .bl_evidence$gpower_classification_status))
}

# Non-default control settings as an R expression (inert string).
.bl_control_call <- function(control) {
  def <- boschloo_control()
  parts <- character()
  for (nm in c("tie_tolerance", "grid_points", "nuisance_tol", "polish_tol", "max_levels")) {
    if (!identical(control[[nm]], def[[nm]])) {
      parts <- c(parts, paste0(nm, " = ", .ig_r_literal(control[[nm]])))
    }
  }
  if (!identical(control$ambiguity_band, def$ambiguity_band)) {
    parts <- c(parts, paste0("ambiguity_band = ", .ig_r_literal(control$ambiguity_band)))
  }
  if (!length(parts)) return(NULL)
  paste0("control = inferguidR::boschloo_control(", paste(parts, collapse = ", "), ")")
}

.bl_test_call <- function(x1, n1, x2, n2, variant, alpha, groups, event, control) {
  args <- c(sprintf("x1 = %dL, n1 = %dL, x2 = %dL, n2 = %dL", x1, n1, x2, n2),
            sprintf("variant = %s", .ig_r_literal(variant)),
            if (!is.null(alpha)) sprintf("alpha = %s", .ig_r_literal(alpha)),
            if (!identical(groups, c("group 1", "group 2"))) sprintf("groups = %s", .ig_r_literal(groups)),
            if (!identical(event, "event")) sprintf("event = %s", .ig_r_literal(event)),
            .bl_control_call(control))
  paste0("inferguidR::boschloo_test(", paste(args, collapse = ", "), ")")
}

.bl_fmt <- function(x, digits = 6) {
  if (is.null(x) || length(x) == 0L) return("NA")
  ifelse(is.na(x), "NA", formatC(x, digits = digits, format = "g"))
}

#' @export
format.boschloo_test <- function(x, ...) {
  d <- x$data
  lines <- c(sprintf("<boschloo_test> %s", x$variant$label),
             sprintf("  data:        %s", d$orientation),
             sprintf("  H0: p1 = p2;  H1: %s", x$variant$alternative),
             sprintf("  status:      %s", format(.ig_status_from_list(x$status))))
  if (identical(x$status$code, "ok")) {
    lines <- c(lines,
               sprintf("  p-value:     %s (rigorous upper bound %s)", .bl_fmt(x$p_value, 7),
                       .bl_fmt(x$p_upper, 7)),
               if (!is.null(x$decision)) sprintf("  decision at alpha = %s: %s", format(x$alpha), x$decision),
               sprintf("  Fisher p (same convention, comparison only): %s", .bl_fmt(x$fisher_p_value, 7)),
               vapply(x$sides, function(s) sprintf(
                 "  %s tail: sup at pi = %s; ordering statistic %s; tail of %d tables; %s ordering",
                 s$side, .bl_fmt(s$nuisance_argmax, 5), .bl_fmt(s$ordering_statistic, 7),
                 as.integer(s$tail_tables), s$ordering_arithmetic), character(1)))
  }
  c(lines,
    sprintf("  settings:    tie tolerance %s; %d grid points + bisection; bound tolerance %s",
            format(x$settings$tie_tolerance), as.integer(x$settings$grid_points),
            format(x$settings$nuisance_tol)),
    sprintf("  scope:       %s", x$scope))
}

#' @export
print.boschloo_test <- function(x, ...) {
  cat(format(x), sep = "\n")
  invisible(x)
}
