# Boschloo rejection regions, attained size and power ------------------------------
#
# For a design (n1, n2), level alpha and variant, the rejection region is
# R = {x : p_upper(x) <= alpha}. Tail sets are nested in the ordering, so R
# is the union of the levels 1..j* where j* is the last level whose
# supremum bound is at most alpha; j* is found by binary search over the
# levels (each step one rigorous supremum). The central variant rejects when
# 2 * min(p_less, p_greater) <= alpha, i.e. R = R_less(alpha / 2) union
# R_greater(alpha / 2).
#
# The attained size is the supremum over the nuisance probability of the
# null probability of R (with its own rigorous bound), and the power at
# (p1, p2) is the double binomial sum over R. The observed test, the region,
# the plots and the planning numbers all use this one region object.

# Region of one side at level alpha (cached).
.bl_region_side <- function(n1, n2, side, alpha, control) {
  key <- paste("region", n1, n2, side, .ig_format_double(alpha), .bl_control_key(control), sep = "|")
  cached <- .bl_cache_get(key, control)
  if (!is.null(cached)) return(cached)
  ord <- .bl_ordering(n1, n2, side, control)
  L <- ord$n_levels
  sups <- list()
  evaluate <- function(j) {
    s <- .bl_sup(.bl_tail_c(ord, .bl_level_tables(ord, j)), control, alpha)
    sups[[as.character(j)]] <<- s
    s
  }
  # Invariant: level lo is admissible (lo = 0 is the empty region), level hi
  # is not (the full sample space has null probability 1 > alpha).
  lo <- 0L
  hi <- L
  while (hi - lo > 1L) {
    mid <- (lo + hi) %/% 2L
    s <- evaluate(mid)
    if (s$upper <= alpha) lo <- mid else hi <- mid
  }
  cut_sup <- if (lo == 0L) NULL else sups[[as.character(lo)]]
  next_sup <- if (hi >= L) NULL else sups[[as.character(hi)]] %||% evaluate(hi)
  idx <- .bl_level_tables(ord, lo)
  trace <- data.frame(level = as.integer(names(sups)),
                      value = vapply(sups, `[[`, 0, "value"),
                      upper = vapply(sups, `[[`, 0, "upper"),
                      row.names = NULL)
  out <- list(side = side, alpha = alpha, level = lo, n_levels = L, idx = idx,
              cut_sup = cut_sup, next_sup = next_sup,
              next_resolved = is.null(next_sup) || isTRUE(next_sup$resolved),
              converged = all(vapply(sups, function(s) isTRUE(s$converged), logical(1))),
              ambiguity_touches = .bl_ambiguity_touches(ord, lo, "cut"),
              arithmetic = ord$arithmetic, ambiguous_gaps = ord$ambiguous_gaps,
              threshold = if (lo > 0L) exp(max(ord$lt[idx])) else NA_real_,
              next_threshold = if (hi <= L && hi > lo) {
                exp(min(ord$lt[ord$perm[ord$level_start[hi]:ord$level_end[hi]]]))
              } else NA_real_,
              trace = trace[order(trace$level), , drop = FALSE])
  .bl_cache_set(key, out, control)
}

# Full region for a variant: membership, size and checks.
.bl_region <- function(n1, n2, alpha, info, control) {
  n1 <- as.integer(n1)
  n2 <- as.integer(n2)
  key <- paste("fullregion", n1, n2, info$variant, .ig_format_double(alpha), .bl_control_key(control),
               sep = "|")
  cached <- .bl_cache_get(key, control)
  if (!is.null(cached)) return(cached)
  alpha_side <- if (info$variant == "central") alpha / 2 else alpha
  sides <- lapply(info$sides, function(s) .bl_region_side(n1, n2, s, alpha_side, control))
  names(sides) <- info$sides
  n_tab <- (n1 + 1L) * (n2 + 1L)
  member <- logical(n_tab)
  for (s in sides) member[s$idx] <- TRUE
  idx <- which(member)
  ord1 <- .bl_ordering(n1, n2, info$sides[1], control)
  size <- if (length(sides) == 1L) {
    sides[[1]]$cut_sup %||% .bl_sup(numeric(n1 + n2 + 1L), control)
  } else {
    .bl_sup(.bl_tail_c(ord1, idx), control)
  }
  problems <- character()
  for (s in sides) {
    if (isTRUE(s$ambiguity_touches)) {
      problems <- c(problems, sprintf(
        "the %s ordering has an ambiguous tie at the region boundary and n1 + n2 = %d exceeds %d, the limit of exact integer arithmetic",
        s$side, n1 + n2, .bl_exact_max_n))
    }
    if (!isTRUE(s$converged)) {
      problems <- c(problems, sprintf("a nuisance supremum of the %s ordering did not converge", s$side))
    }
  }
  if (!isTRUE(size$converged)) problems <- c(problems, "the attained-size supremum did not converge")
  warnings <- character()
  for (s in sides) {
    if (!isTRUE(s$next_resolved)) {
      warnings <- c(warnings, sprintf(
        "the first excluded %s level has a supremum within the floating-point allowance of alpha; it is excluded (conservative)",
        s$side))
    }
  }
  out <- list(n1 = n1, n2 = n2, alpha = alpha, info = info, sides = sides, member = member,
              idx = idx, size = size, problems = problems, warnings = warnings,
              x1 = ord1$x1, x2 = ord1$x2, t = ord1$t)
  .bl_cache_set(key, out, control)
}

# Power: sum over the region of dbinom(x1; n1, p1) * dbinom(x2; n2, p2).
.bl_power_idx <- function(idx, n1, n2, p1, p2) {
  if (!length(idx)) return(0)
  P <- outer(stats::dbinom(0:n1, n1, p1), stats::dbinom(0:n2, n2, p2))
  min(1, sum(P[idx]))
}

.bl_region_status <- function(reg) {
  if (length(reg$problems)) ig_status("numerical_failure", paste(reg$problems, collapse = "; ")) else ig_status("ok")
}

.bl_side_region_record <- function(s, ok) {
  num <- function(v) if (ok) v else NA_real_
  list(side = s$side,
       direction = if (s$side == "lower") "p1 < p2" else if (s$side == "upper") "p1 > p2" else "p1 != p2",
       alpha_side = s$alpha, last_level_in_region = s$level, n_levels = s$n_levels,
       tables_in_region = length(s$idx),
       ordering_threshold = num(s$threshold),
       next_ordering_value = num(s$next_threshold),
       p_at_boundary = if (is.null(s$cut_sup)) NULL else
         list(value = num(s$cut_sup$value), upper = num(s$cut_sup$upper)),
       p_first_excluded = if (is.null(s$next_sup)) NULL else
         list(value = num(s$next_sup$value), upper = num(s$next_sup$upper),
              resolved = isTRUE(s$next_sup$resolved)),
       ordering_arithmetic = s$arithmetic, ambiguous_gaps = s$ambiguous_gaps,
       binary_search = list(levels = s$trace$level, value = if (ok) s$trace$value else NA_real_,
                            upper = if (ok) s$trace$upper else NA_real_))
}

#' Rejection region of Boschloo's test for a fixed design
#'
#' Enumerates all `(n1 + 1) * (n2 + 1)` possible tables of a planned study
#' and returns the tables that Boschloo's test rejects at level `alpha`,
#' together with the attained size (the largest null rejection probability
#' over the common event probability). The region is exactly the set of
#' tables for which [boschloo_test()] decides `"reject"`.
#'
#' @param n1,n2 Planned group sizes.
#' @param alpha Nominal significance level.
#' @param variant `"minlike"`, `"central"`, `"less"` or `"greater"` (see
#'   [boschloo_test()]).
#' @param control Numerical settings from [boschloo_control()].
#' @return An immutable object of class `boschloo_region` with the rejected
#'   tables (`rejected`, columns `x1`, `x2`), the attained size (`size`,
#'   with value, rigorous upper bound and maximising nuisance probability),
#'   the boundary of each ordering side, the status, the settings and a
#'   region hash.
#' @examples
#' reg <- boschloo_region(10, 12, alpha = 0.05, variant = "central")
#' reg
#' reg$size$value
#' @family Boschloo
#' @export
boschloo_region <- function(n1, n2, alpha = 0.05, variant, control = boschloo_control()) {
  info <- .bl_variant(if (missing(variant)) NULL else variant)
  .bl_check_design(n1, n2)
  alpha <- .bl_check_alpha(alpha)
  control <- .bl_check_control(control)
  reg <- .bl_region(n1, n2, alpha, info, control)
  .bl_region_object(reg, control)
}

.bl_region_object <- function(reg, control) {
  status <- .bl_region_status(reg)
  ok <- status$code == "ok"
  num <- function(v) if (ok) v else NA_real_
  rejected <- data.frame(x1 = reg$x1[reg$idx], x2 = reg$x2[reg$idx])
  if (!ok) rejected <- rejected[0, , drop = FALSE]
  body <- list(
    object_type = "boschloo_region",
    status = .ig_status_to_list(status),
    variant = .bl_variant_record(reg$info),
    design = list(n1 = reg$n1, n2 = reg$n2, n_total = reg$n1 + reg$n2,
                  n_tables = length(reg$member)),
    alpha = reg$alpha,
    rejected = rejected,
    n_rejected = if (ok) length(reg$idx) else NA_integer_,
    size = list(value = num(reg$size$value), upper = num(reg$size$upper),
                argmax = num(reg$size$argmax), converged = isTRUE(reg$size$converged),
                note = "attained size = sup over the common event probability of the null rejection probability"),
    sides = lapply(reg$sides, .bl_side_region_record, ok = ok),
    warnings = reg$warnings,
    undefined = if (ok) .ig_empty_object() else list(region = status$reason),
    settings = .bl_control_record(control),
    scope = "hypothetical future study with fixed group sizes n1 and n2",
    provenance = .bl_provenance()
  )
  body$region_hash <- .ig_hash_value(list(design = body$design, alpha = body$alpha,
                                          variant_id = reg$info$variant_id,
                                          settings = body$settings,
                                          rejected = list(x1 = I(rejected$x1), x2 = I(rejected$x2))))
  structure(body, class = c("boschloo_region", "ig_frozen"))
}

#' @export
format.boschloo_region <- function(x, ...) {
  lines <- c(sprintf("<boschloo_region> %s", x$variant$label),
             sprintf("  design:      n1 = %d, n2 = %d (%d tables); alpha = %s",
                     as.integer(x$design$n1), as.integer(x$design$n2),
                     as.integer(x$design$n_tables), format(x$alpha)),
             sprintf("  status:      %s", format(.ig_status_from_list(x$status))))
  if (identical(x$status$code, "ok")) {
    lines <- c(lines,
               sprintf("  rejected:    %d tables", as.integer(x$n_rejected)),
               sprintf("  attained size: %s (upper bound %s) at pi = %s", .bl_fmt(x$size$value, 7),
                       .bl_fmt(x$size$upper, 7), .bl_fmt(x$size$argmax, 5)))
  }
  c(lines, if (length(x$warnings)) paste0("  warning: ", x$warnings),
    sprintf("  region hash: %s", x$region_hash))
}

#' @export
print.boschloo_region <- function(x, ...) {
  cat(format(x), sep = "\n")
  invisible(x)
}

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

# Fisher comparison ----------------------------------------------------------------------

# Fisher conditional p-values of every table for a variant.
.bl_fisher_p <- function(n1, n2, info, control, convention = "same") {
  if (convention == "minlike" || info$variant == "minlike") {
    return(exp(.bl_ordering(n1, n2, "minlike", control)$lt))
  }
  if (info$variant == "central") {
    tl <- exp(.bl_ordering(n1, n2, "lower", control)$lt)
    tu <- exp(.bl_ordering(n1, n2, "upper", control)$lt)
    return(pmin(1, 2 * pmin(tl, tu)))
  }
  exp(.bl_ordering(n1, n2, info$sides, control)$lt)
}

#' Compare Boschloo's and Fisher's rejection regions
#'
#' Classifies every possible table of a design as rejected by both tests,
#' by Boschloo's test only, by Fisher's conditional test only, or by
#' neither, and reports the attained size of both regions and (optionally)
#' their power at `(p1, p2)`. Both tests are evaluated on the same design,
#' level and sidedness convention unless `fisher_convention = "minlike"`
#' requests the two-sided Fisher p-value of [stats::fisher.test()] (which
#' is a different two-sided convention from the central variant and is
#' labelled as such).
#'
#' With the same convention, Boschloo's p-value never exceeds Fisher's. The
#' result checks this on the p-values of every table (`dominance_holds`,
#' computed for designs with at most 5000 tables). Fisher p-values within the
#' relative tie tolerance of `alpha` count as equal to `alpha`.
#'
#' @inheritParams boschloo_power
#' @param p1,p2 Optional alternative event probabilities for a power
#'   comparison.
#' @param fisher_convention `"same"` (default) or `"minlike"`.
#' @return An object of class `boschloo_fisher_comparison` with a table of
#'   all `(x1, x2)` with both p-values (the Boschloo p-value is `NA` above
#'   5000 tables), both decisions and the category, counts per category, the
#'   dominance check, attained sizes and power of both tests.
#' @examples
#' cmp <- boschloo_compare_fisher(12, 12, alpha = 0.05, variant = "central",
#'                                p1 = 0.2, p2 = 0.7)
#' cmp$counts
#' @family Boschloo
#' @export
boschloo_compare_fisher <- function(n1, n2, alpha = 0.05, variant, p1 = NULL, p2 = NULL,
                                    fisher_convention = c("same", "minlike"),
                                    control = boschloo_control()) {
  info <- .bl_variant(if (missing(variant)) NULL else variant)
  .bl_check_design(n1, n2)
  alpha <- .bl_check_alpha(alpha)
  fisher_convention <- match.arg(fisher_convention)
  if (is.null(p1) != is.null(p2)) .ig_abort_invalid("give both `p1` and `p2`, or neither")
  if (!is.null(p1)) {
    p1 <- .bl_check_prob(p1, "p1")
    p2 <- .bl_check_prob(p2, "p2")
  }
  control <- .bl_check_control(control)
  n1 <- as.integer(n1)
  n2 <- as.integer(n2)
  reg <- .bl_region(n1, n2, alpha, info, control)
  status <- .bl_region_status(reg)
  fp <- .bl_fisher_p(n1, n2, info, control, fisher_convention)
  fisher_member <- fp <= alpha * (1 + control$tie_tolerance)
  ord1 <- .bl_ordering(n1, n2, info$sides[1], control)
  fisher_size <- .bl_sup(.bl_tail_c(ord1, which(fisher_member)), control)
  cat_of <- ifelse(reg$member & fisher_member, "both",
                   ifelse(reg$member, "boschloo_only", ifelse(fisher_member, "fisher_only", "neither")))
  same <- fisher_convention == "same" || info$variant == "minlike"
  # Boschloo p-values of every table (one supremum per ordering level). The
  # dominance property p_Boschloo <= p_Fisher is checked on the p-values:
  # region decisions alone could show a "Fisher only" table when a p-value
  # lies within the numerical tolerance of alpha.
  n_tab <- length(reg$member)
  pv <- if (n_tab <= 5000L && status$code == "ok") .bl_pvalues_all(n1, n2, info, control) else NULL
  dominance <- NA
  if (same && !is.null(pv)) {
    dominance <- all(pv$value <= fp * (1 + control$tie_tolerance) + .bl_fp_allowance(n1 + n2))
  }
  ok <- status$code == "ok"
  if (!ok) cat_of <- rep(NA_character_, n_tab)    # no Boschloo decision without a certified region
  tables <- data.frame(x1 = reg$x1, x2 = reg$x2, fisher_p = fp,
                       boschloo_p = if (is.null(pv)) NA_real_ else pv$value,
                       boschloo_reject = if (ok) reg$member else NA,
                       fisher_reject = fisher_member, category = cat_of, stringsAsFactors = FALSE)
  counts <- vapply(c("both", "boschloo_only", "fisher_only", "neither"),
                   function(k) if (ok) sum(cat_of == k) else NA_integer_, integer(1))
  num <- function(v) if (ok) v else NA_real_
  out <- list(
    object_type = "boschloo_fisher_comparison",
    status = .ig_status_to_list(status),
    variant = .bl_variant_record(info),
    fisher_test = list(
      convention = if (same) "same sidedness convention as the Boschloo variant" else
        "two-sided Fisher p-value of stats::fisher.test (minlike); a different convention from the central Boschloo variant",
      conditional_on = "both margins (total number of events fixed)"),
    design = list(n1 = n1, n2 = n2), alpha = alpha,
    tables = tables, counts = counts,
    dominance_holds = dominance,
    dominance_note = if (!same) "not applicable across different two-sided conventions" else
      if (is.null(pv)) "not checked (more than 5000 tables, or the region has no ok status)" else
        paste("checked on the p-values of every table: Boschloo p <= Fisher p for the same ordering",
              "statistic (a 'fisher_only' table can only arise within numerical tolerance of alpha)"),
    size = list(boschloo = list(value = num(reg$size$value), upper = num(reg$size$upper)),
                fisher = list(value = fisher_size$value, upper = fisher_size$upper)),
    power = if (is.null(p1)) NULL else list(
      p1 = p1, p2 = p2,
      boschloo = num(.bl_power_idx(reg$idx, n1, n2, p1, p2)),
      fisher = .bl_power_idx(which(fisher_member), n1, n2, p1, p2)),
    undefined = if (ok) .ig_empty_object() else list(boschloo = status$reason),
    settings = .bl_control_record(control),
    scope = "hypothetical future study; comparison of two different tests on the same design"
  )
  structure(out, class = "boschloo_fisher_comparison")
}

#' @export
print.boschloo_fisher_comparison <- function(x, ...) {
  cat(sprintf("<boschloo_fisher_comparison> %s vs Fisher (%s)\n", x$variant$label, x$fisher_test$convention),
      sprintf("  design n1 = %d, n2 = %d, alpha = %s\n", x$design$n1, x$design$n2, format(x$alpha)),
      sprintf("  rejected by both: %d; Boschloo only: %d; Fisher only: %d; neither: %d\n",
              x$counts[["both"]], x$counts[["boschloo_only"]], x$counts[["fisher_only"]],
              x$counts[["neither"]]),
      sprintf("  attained size: Boschloo %s, Fisher %s\n", .bl_fmt(x$size$boschloo$value),
              .bl_fmt(x$size$fisher$value)), sep = "")
  if (!is.null(x$power)) {
    cat(sprintf("  power at p1 = %s, p2 = %s: Boschloo %s, Fisher %s\n", format(x$power$p1),
                format(x$power$p2), .bl_fmt(x$power$boschloo), .bl_fmt(x$power$fisher)))
  }
  invisible(x)
}
