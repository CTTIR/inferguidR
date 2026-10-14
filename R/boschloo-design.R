# Boschloo planning: allocation, effect conversion and sample-size search ---------------
#
# Power of an exact unconditional test is a saw-toothed function of the
# sample size: adding a participant can lower the power. The search therefore
# never brackets or interpolates. It walks the declared allocation sequence
# of designs upward and evaluates the exact power of every design it reports
# on, and it reports two answers:
#
# (a) n_first: the first design whose power reaches the target. Designs
#     before it are excluded either by direct evaluation or by the
#     Neyman-Pearson bound below, so (a) is the minimum over the declared
#     admissible designs up to `n_max`.
# (b) n_stable: the first design from which the next `stability_window`
#     designs also reach the target. Nothing is claimed beyond that window.
#
# Neyman-Pearson bound. A Boschloo region has null rejection probability at
# most alpha at every common probability pi0. Its power at (p1, p2) can
# therefore not exceed the power of the most powerful (randomised) level-alpha
# test of the simple null (pi0, pi0) against the simple alternative (p1, p2),
# for any pi0 in (0, 1). That bound does not decrease along an allocation
# sequence in which n1 and n2 never decrease (a larger study can ignore its
# extra observations), so a binary search finds the first design whose bound
# reaches the target; all earlier designs fail for every level-alpha test.

.bl_allocations <- c("n2_ceiling", "total_nearest")

# Group sizes of the m-th design of an allocation sequence.
#
# * n2_ceiling:    n1 = m, n2 = max(1, ceiling(r * n1)) (r = n2 / n1); the
#                  ceiling is guarded against the floating representation of r.
# * total_nearest: N = m + 1, n1 = round(N / (1 + r)) (halves rounded up),
#                  kept in [1, N - 1], n2 = N - n1.
.bl_design_at <- function(m, r, allocation) {
  m <- as.integer(m)
  if (allocation == "n2_ceiling") {
    raw <- r * m
    n2 <- pmax(1L, as.integer(ceiling(raw - 1e-9 * pmax(1, raw))))
    return(list(n1 = m, n2 = n2))
  }
  N <- m + 1L
  n1 <- as.integer(floor(N / (1 + r) + 0.5 + 1e-9))
  n1 <- pmin(pmax(n1, 1L), N - 1L)
  list(n1 = n1, n2 = N - n1)
}

# Index of the largest design with n1 + n2 <= n_max (0 when even the first
# design is larger).
.bl_max_index <- function(n_max, r, allocation) {
  if (allocation == "total_nearest") return(as.integer(n_max) - 1L)
  total <- function(m) sum(unlist(.bl_design_at(m, r, allocation)))
  m <- max(1L, floor(n_max / (1 + r)))
  while (m > 1L && total(m) > n_max) m <- m - 1L
  if (total(m) > n_max) return(0L)
  while (total(m + 1L) <= n_max) m <- m + 1L
  as.integer(m)
}

# Power margin of the Neyman-Pearson exclusion: a design is excluded only if
# its bound is below the target by more than this (floating-point safety).
.bl_np_margin <- 1e-6

# Power of the most powerful randomised level-alpha test of (pi0, pi0)
# against (p1, p2) for the design (n1, n2).
.bl_np_power <- function(n1, n2, p1, p2, alpha, pi0) {
  tab <- .bl_tables(n1, n2)
  l0 <- stats::dbinom(tab$x1, n1, pi0, log = TRUE) + stats::dbinom(tab$x2, n2, pi0, log = TRUE)
  l1 <- stats::dbinom(tab$x1, n1, p1, log = TRUE) + stats::dbinom(tab$x2, n2, p2, log = TRUE)
  llr <- l1 - l0                                # l0 is finite for pi0 in (0, 1)
  o <- order(llr, decreasing = TRUE)
  llr <- llr[o]
  p0 <- exp(l0[o])
  q1 <- exp(l1[o])
  # Blocks of (numerically) tied likelihood ratios are randomised together.
  brk <- c(TRUE, abs(diff(llr)) > 1e-9 * pmax(1, abs(llr[-1L])) | !is.finite(diff(llr)))
  block <- cumsum(brk)
  b0 <- as.vector(tapply(p0, block, sum))
  b1 <- as.vector(tapply(q1, block, sum))
  c0 <- cumsum(b0)
  full <- c0 <= alpha
  pw <- sum(b1[full])
  k <- which(!full)[1L]
  if (!is.na(k)) {
    room <- alpha - (if (k > 1L) c0[k - 1L] else 0)
    pw <- pw + b1[k] * room / b0[k]
  }
  min(1, pw)
}

.bl_np_bound <- function(n1, n2, p1, p2, alpha) {
  cand <- unique(c(p1, p2, (p1 + p2) / 2))
  cand <- cand[cand > 0 & cand < 1]
  if (!length(cand)) cand <- 0.5
  min(vapply(cand, function(pi0) .bl_np_power(n1, n2, p1, p2, alpha, pi0), double(1)))
}

#' Event probabilities from a baseline and an effect
#'
#' Converts a baseline event probability of group 1 (the reference group)
#' and one effect measure into the probability of group 2, checks that it is
#' a valid probability and returns all effect scales. Effects are oriented
#' as group 2 relative to group 1: difference `p2 - p1`, ratio `p2 / p1`,
#' odds ratio `odds(p2) / odds(p1)`.
#'
#' @param p1 Event probability of group 1, in `[0, 1]`.
#' @param p2,difference,ratio,odds_ratio Exactly one of: the probability of
#'   group 2, the risk difference, the risk ratio or the odds ratio.
#' @return A list with `p1`, `p2`, `difference`, `ratio` and `odds_ratio`
#'   (`NA` with a reason in `undefined` when a scale is not defined, for
#'   example a ratio with `p1 = 0`).
#' @examples
#' boschloo_probabilities(0.2, difference = 0.15)
#' boschloo_probabilities(0.2, odds_ratio = 3)
#' @family Boschloo
#' @export
boschloo_probabilities <- function(p1, p2 = NULL, difference = NULL, ratio = NULL,
                                   odds_ratio = NULL) {
  p1 <- .bl_check_prob(p1, "p1")
  given <- c(p2 = !is.null(p2), difference = !is.null(difference), ratio = !is.null(ratio),
             odds_ratio = !is.null(odds_ratio))
  if (sum(given) != 1L) {
    .ig_abort_invalid("give exactly one of `p2`, `difference`, `ratio` or `odds_ratio`")
  }
  num1 <- function(x) is.numeric(x) && length(x) == 1L && is.finite(x)
  if (given[["p2"]]) {
    p2 <- .bl_check_prob(p2, "p2")
  } else if (given[["difference"]]) {
    if (!num1(difference)) .ig_abort_invalid("`difference` must be a finite number")
    p2 <- p1 + difference
  } else if (given[["ratio"]]) {
    if (!num1(ratio) || ratio < 0) .ig_abort_invalid("`ratio` must be a finite number >= 0")
    if (p1 == 0) .ig_abort_invalid("a risk ratio cannot define p2 when p1 = 0")
    p2 <- p1 * ratio
  } else {
    if (!num1(odds_ratio) || odds_ratio <= 0) .ig_abort_invalid("`odds_ratio` must be a finite number > 0")
    if (p1 %in% c(0, 1)) .ig_abort_invalid("an odds ratio cannot define p2 when p1 is 0 or 1")
    o2 <- odds_ratio * p1 / (1 - p1)
    p2 <- o2 / (1 + o2)
  }
  if (!is.finite(p2) || p2 < -1e-12 || p2 > 1 + 1e-12) {
    .ig_abort_invalid(sprintf("the effect implies p2 = %s, which is not a probability in [0, 1]",
                              format(p2)))
  }
  p2 <- min(1, max(0, p2))
  undefined <- list()
  ratio_v <- if (p1 > 0) p2 / p1 else {
    undefined$ratio <- "risk ratio undefined for p1 = 0"
    NA_real_
  }
  or_v <- if (p1 > 0 && p1 < 1 && p2 > 0 && p2 < 1) (p2 / (1 - p2)) / (p1 / (1 - p1)) else {
    undefined$odds_ratio <- "odds ratio undefined when a probability is 0 or 1"
    NA_real_
  }
  list(p1 = p1, p2 = p2, difference = p2 - p1, ratio = ratio_v, odds_ratio = or_v,
       orientation = "group 2 relative to group 1 (reference)",
       undefined = .ig_named_or_empty(undefined))
}

#' Sample size for Boschloo's test
#'
#' Finds the smallest design of a declared allocation sequence whose exact
#' Boschloo power at `(p1, p2)` reaches `power`, recomputing the achieved
#' power and attained size of every reported design from its rejection
#' region (the same region as [boschloo_region()] and [boschloo_test()]).
#'
#' Because power is not monotone in the sample size, two designs are
#' reported:
#'
#' * `first`: the first design in the sequence whose power reaches the
#'   target. Every earlier design up to it was either evaluated and failed,
#'   or is excluded by a Neyman-Pearson bound (the power of the most
#'   powerful level-`alpha` test of a simple null `(pi0, pi0)` against
#'   `(p1, p2)`, which no level-`alpha` test can exceed). It is therefore the
#'   minimum over the admissible designs with `n1 + n2 <= n_max`.
#' * `stable`: the first design from which the next `stability_window`
#'   designs of the sequence also reach the target (the window `[N, N + K]`
#'   in design steps). Nothing is claimed for designs beyond the window.
#'
#' Allocation policies (allocation ratio r = `allocation_ratio` = n2 / n1):
#'
#' * `"n2_ceiling"`: `n1 = 1, 2, ...` and `n2 = max(1, ceiling(r * n1))`;
#' * `"total_nearest"`: `N = 2, 3, ...`, `n1 = round(N / (1 + r))` (halves
#'   rounded up, kept between 1 and `N - 1`) and `n2 = N - n1`.
#'
#' @param p1,p2 Event probabilities of group 1 (reference) and group 2 under
#'   the alternative; see [boschloo_probabilities()] for baseline-plus-effect
#'   input.
#' @param power Target power.
#' @inheritParams boschloo_region
#' @param allocation_ratio The allocation ratio r = n2 / n1 (default 1).
#' @param allocation `"n2_ceiling"` (default) or `"total_nearest"`.
#' @param stability_window `K`, the number of following designs that must
#'   also reach the target for the `stable` design (default 10).
#' @param n_max Largest total sample size searched (default 500).
#' @return An immutable object of class `boschloo_sample_size` with status
#'   (`"ok"`; `"infeasible"` when no design up to `n_max` reaches the target
#'   or when `p1 = p2`; `"invalid_input"` for a one-sided variant whose
#'   direction contradicts the alternative or a target power not above
#'   `alpha`; `"numerical_failure"` when a region cannot be certified), the
#'   `first` and `stable` designs (group sizes, total, achieved power,
#'   attained size), the evaluated staircase, the Neyman-Pearson exclusion,
#'   the minimality evidence and the search policy. Malformed arguments (for
#'   example a probability outside `[0, 1]`) signal an
#'   `inferguidR_invalid_input` condition instead.
#' @examples
#' ss <- boschloo_sample_size(p1 = 0.2, p2 = 0.7, power = 0.8, alpha = 0.05,
#'                            variant = "central")
#' ss
#' ss$first
#' @family Boschloo
#' @export
boschloo_sample_size <- function(p1, p2, power = 0.8, alpha = 0.05, variant,
                                 allocation_ratio = 1, allocation = c("n2_ceiling", "total_nearest"),
                                 stability_window = 10L, n_max = 500L,
                                 control = boschloo_control()) {
  info <- .bl_variant(if (missing(variant)) NULL else variant)
  p1 <- .bl_check_prob(p1, "p1")
  p2 <- .bl_check_prob(p2, "p2")
  target <- .bl_check_alpha(power, "power")
  alpha <- .bl_check_alpha(alpha)
  allocation <- match.arg(allocation)
  pr <- character()
  if (!is.numeric(allocation_ratio) || length(allocation_ratio) != 1L ||
      !is.finite(allocation_ratio) || allocation_ratio <= 0 || allocation_ratio > 100) {
    pr <- c(pr, "`allocation_ratio` must be a number in (0, 100]")
  }
  if (!.bl_is_count(stability_window, 0) || stability_window > 1000) {
    pr <- c(pr, "`stability_window` must be a whole number between 0 and 1000")
  }
  if (!.bl_is_count(n_max, 2) || n_max > 3000) pr <- c(pr, "`n_max` must be a whole number between 2 and 3000")
  if (length(pr)) .ig_abort_invalid(pr, "invalid input for boschloo_sample_size()")
  control <- .bl_check_control(control)
  r <- as.double(allocation_ratio)
  K <- as.integer(stability_window)
  n_max <- as.integer(n_max)

  policy <- list(allocation = allocation, allocation_ratio = r,
                 rule = if (allocation == "n2_ceiling") "n1 = 1, 2, ...; n2 = max(1, ceiling(r * n1))" else
                   "N = 2, 3, ...; n1 = round(N / (1 + r)) in [1, N - 1]; n2 = N - n1",
                 stability_window = K, n_max = n_max,
                 minimality = "first = minimum over the allocation sequence with n1 + n2 <= n_max; stable = first design followed by K designs that all reach the target (no claim beyond)")
  inputs <- list(p1 = p1, p2 = p2, power = target, alpha = alpha, variant = info$variant)
  call <- sprintf("inferguidR::boschloo_sample_size(p1 = %s, p2 = %s, power = %s, alpha = %s, variant = %s, allocation_ratio = %s, allocation = %s, stability_window = %dL, n_max = %dL%s)",
                  .ig_r_literal(p1), .ig_r_literal(p2), .ig_r_literal(target), .ig_r_literal(alpha),
                  .ig_r_literal(info$variant), .ig_r_literal(r), .ig_r_literal(allocation), K, n_max,
                  if (is.null(.bl_control_call(control))) "" else paste0(", ", .bl_control_call(control)))
  finish <- function(status, first = NULL, stable = NULL, staircase = NULL, np = NULL, notes = character(),
                     minimality = NULL) {
    body <- list(object_type = "boschloo_sample_size", status = .ig_status_to_list(status),
                 variant = .bl_variant_record(info), inputs = inputs, policy = policy,
                 first = first, stable = stable, staircase = staircase,
                 neyman_pearson = np, minimality = minimality, notes = notes,
                 settings = .bl_control_record(control),
                 scope = "planning of a future study (hypothetical data)",
                 provenance = .bl_provenance(), r_call = call)
    structure(body, class = c("boschloo_sample_size", "ig_frozen"))
  }

  # Requests that are well formed but cannot be answered return a result with
  # a typed status (no numbers) instead of stopping.
  if (target <= alpha) {
    return(finish(ig_status("invalid_input", sprintf(
      "the target power (%s) must exceed alpha (%s)", format(target), format(alpha)))))
  }
  if ((info$variant == "less" && p1 >= p2) || (info$variant == "greater" && p1 <= p2)) {
    return(finish(ig_status("invalid_input", sprintf(
      "the one-sided variant '%s' (H1: %s) does not match the alternative p1 = %s, p2 = %s",
      info$variant, info$alternative, format(p1), format(p2)))))
  }
  if (p1 == p2) {
    return(finish(ig_status("infeasible", sprintf(
      "with p1 = p2 the rejection probability is a null probability, at most alpha = %s < power = %s",
      format(alpha), format(target)))))
  }

  m_max <- .bl_max_index(n_max, r, allocation)
  if (m_max < 1L) {
    return(finish(ig_status("infeasible", "no design of the allocation sequence fits n_max")))
  }
  # Neyman-Pearson exclusion (binary search on a non-decreasing bound).
  margin <- .bl_np_margin
  np_at <- function(m) {
    d <- .bl_design_at(m, r, allocation)
    .bl_np_bound(d$n1, d$n2, p1, p2, alpha)
  }
  if (np_at(m_max) < target - margin) {
    d <- .bl_design_at(m_max, r, allocation)
    return(finish(ig_status("infeasible", sprintf(
      "no level-%s test can reach power %s with n1 + n2 <= %d under this allocation (Neyman-Pearson bound %s at n1 = %d, n2 = %d)",
      format(alpha), format(target), n_max, format(np_at(m_max), digits = 6), d$n1, d$n2)),
      np = list(first_index = NA_integer_, bound_at_n_max = np_at(m_max))))
  }
  lo <- 0L
  hi <- m_max
  while (hi - lo > 1L) {
    mid <- (lo + hi) %/% 2L
    if (np_at(mid) >= target - margin) hi <- mid else lo <- mid
  }
  m_lb <- hi
  d_lb <- .bl_design_at(m_lb, r, allocation)
  np <- list(first_index = m_lb, n1 = d_lb$n1, n2 = d_lb$n2, n_total = d_lb$n1 + d_lb$n2,
             bound = np_at(m_lb), margin = margin,
             statement = "every design before this one has a Neyman-Pearson power bound below the target and is excluded for any level-alpha test")

  rows <- list()
  first <- NULL
  stable <- NULL
  run <- 0L
  failure <- NULL
  for (m in m_lb:m_max) {
    d <- .bl_design_at(m, r, allocation)
    reg <- .bl_region(d$n1, d$n2, alpha, info, control)
    if (length(reg$problems)) {
      failure <- sprintf("design n1 = %d, n2 = %d: %s", d$n1, d$n2, paste(reg$problems, collapse = "; "))
      break
    }
    pw <- .bl_power_idx(reg$idx, d$n1, d$n2, p1, p2)
    meets <- pw >= target
    row <- list(index = m, n1 = d$n1, n2 = d$n2, n_total = d$n1 + d$n2, power = pw,
                size = reg$size$value, size_upper = reg$size$upper, n_rejected = length(reg$idx),
                meets_target = meets)
    rows[[length(rows) + 1L]] <- row
    if (meets && is.null(first)) first <- row
    run <- if (meets) run + 1L else 0L
    if (run >= K + 1L) {
      stable <- rows[[length(rows) - K]]
      break
    }
  }
  staircase <- do.call(rbind, lapply(rows, as.data.frame))
  if (!is.null(failure)) {
    return(finish(ig_status("numerical_failure", failure), staircase = NULL, np = np))
  }
  notes <- character()
  if (is.null(first)) {
    return(finish(ig_status("infeasible", sprintf(
      "no design with n1 + n2 <= %d reaches power %s (all designs from the Neyman-Pearson bound on were evaluated)",
      n_max, format(target))), staircase = staircase, np = np))
  }
  if (is.null(stable)) {
    notes <- c(notes, sprintf(
      "no stable design: a run of %d consecutive designs reaching the target was not found up to n1 + n2 <= %d",
      K + 1L, n_max))
  }
  if (!is.null(stable) && stable$index > first$index) {
    notes <- c(notes, "power is not monotone here: some designs after the first one fall below the target")
  }
  below_after_first <- staircase$index > first$index & !staircase$meets_target
  minimality <- list(
    first_is_minimum = TRUE,
    scope = sprintf("all designs of the allocation sequence with n1 + n2 <= %d", n_max),
    designs_excluded_by_bound = m_lb - 1L,
    designs_evaluated = nrow(staircase),
    designs_evaluated_below_first = sum(staircase$index < first$index),
    designs_below_target_after_first = sum(below_after_first),
    stability_window = K,
    stable_found = !is.null(stable),
    rule = paste("first: every earlier design was evaluated and failed or is excluded by the",
                 "Neyman-Pearson bound; stable: the first design whose next K designs all reach",
                 "the target (window [first stable design, K designs later]); no claim beyond the window"))
  finish(ig_status("ok"), first = first, stable = stable, staircase = staircase, np = np, notes = notes,
         minimality = minimality)
}

#' @export
format.boschloo_sample_size <- function(x, ...) {
  des <- function(d) {
    if (is.null(d)) return("not established")
    sprintf("n1 = %d, n2 = %d (N = %d); power %s; attained size %s", as.integer(d$n1),
            as.integer(d$n2), as.integer(d$n_total), .bl_fmt(d$power, 6), .bl_fmt(d$size, 6))
  }
  c(sprintf("<boschloo_sample_size> %s", x$variant$label),
    sprintf("  inputs:      p1 = %s, p2 = %s, alpha = %s, target power = %s", format(x$inputs$p1),
            format(x$inputs$p2), format(x$inputs$alpha), format(x$inputs$power)),
    sprintf("  allocation:  %s (r = %s)", x$policy$rule, format(x$policy$allocation_ratio)),
    sprintf("  status:      %s", format(.ig_status_from_list(x$status))),
    if (identical(x$status$code, "ok")) c(
      sprintf("  first design reaching the target: %s", des(x$first)),
      sprintf("  stable design (next %d designs also reach it): %s", as.integer(x$policy$stability_window),
              des(x$stable))),
    if (length(x$notes)) paste0("  note: ", x$notes))
}

#' @export
print.boschloo_sample_size <- function(x, ...) {
  cat(format(x), sep = "\n")
  invisible(x)
}
