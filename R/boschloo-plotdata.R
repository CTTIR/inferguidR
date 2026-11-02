
#' Plot data: Boschloo power versus sample size
#'
#' Exact power and attained size of every design of an allocation sequence
#' between two total sample sizes, with the target line. The staircase is not
#' smoothed; designs are the integer designs of the declared allocation.
#'
#' @inheritParams boschloo_sample_size
#' @param n_range Range of total sample sizes `c(from, to)`.
#' @param fisher Also compute the power of Fisher's conditional test with the
#'   same convention (comparison only).
#' @return An object of class `boschloo_plotdata` (plot id
#'   `"boschloo_power_n"`).
#' @examples
#' pd <- boschloo_plotdata_power_n(0.2, 0.7, alpha = 0.05, variant = "central",
#'                                 n_range = c(10, 40))
#' head(pd$data)
#' @family Boschloo
#' @export
boschloo_plotdata_power_n <- function(p1, p2, alpha = 0.05, variant, n_range, power = NULL,
                                      allocation_ratio = 1,
                                      allocation = c("n2_ceiling", "total_nearest"),
                                      fisher = FALSE, control = boschloo_control()) {
  info <- .bl_variant(if (missing(variant)) NULL else variant)
  p1 <- .bl_check_prob(p1, "p1")
  p2 <- .bl_check_prob(p2, "p2")
  alpha <- .bl_check_alpha(alpha)
  if (!is.null(power)) power <- .bl_check_alpha(power, "power")
  allocation <- match.arg(allocation)
  if (!is.numeric(n_range) || length(n_range) != 2L || !all(vapply(n_range, .bl_is_count, logical(1), lower = 2)) ||
      n_range[1] > n_range[2] || n_range[2] > 3000) {
    .ig_abort_invalid("`n_range` must be two whole numbers 2 <= from <= to <= 3000")
  }
  if (!is.numeric(allocation_ratio) || length(allocation_ratio) != 1L || !is.finite(allocation_ratio) ||
      allocation_ratio <= 0) {
    .ig_abort_invalid("`allocation_ratio` must be a positive number")
  }
  control <- .bl_check_control(control)
  rows <- list()
  m <- 1L
  repeat {
    d <- .bl_design_at(m, allocation_ratio, allocation)
    N <- d$n1 + d$n2
    if (N > n_range[2]) break
    if (N >= n_range[1]) {
      reg <- .bl_region(d$n1, d$n2, alpha, info, control)
      ok <- !length(reg$problems)
      pw <- if (ok) .bl_power_idx(reg$idx, d$n1, d$n2, p1, p2) else NA_real_
      fpw <- NA_real_
      if (fisher && ok) {
        fm <- which(.bl_fisher_p(d$n1, d$n2, info, control) <= alpha * (1 + control$tie_tolerance))
        fpw <- .bl_power_idx(fm, d$n1, d$n2, p1, p2)
      }
      rows[[length(rows) + 1L]] <- data.frame(
        n1 = d$n1, n2 = d$n2, n_total = N, power = pw,
        size = if (ok) reg$size$value else NA_real_,
        size_upper = if (ok) reg$size$upper else NA_real_,
        n_rejected = if (ok) length(reg$idx) else NA_integer_,
        meets_target = if (is.null(power) || !ok) NA else pw >= power,
        fisher_power = fpw, status = if (ok) "ok" else "numerical_failure")
    }
    m <- m + 1L
  }
  data <- if (length(rows)) do.call(rbind, rows) else data.frame()
  .bl_plotdata(
    "boschloo_power_n", sprintf("Power versus total sample size, %s", info$label), data,
    fixed = list(p1 = p1, p2 = p2, alpha = alpha, variant = info$variant,
                 variant_id = info$variant_id, allocation = allocation,
                 allocation_ratio = allocation_ratio),
    scope = "hypothetical future study",
    reference = list(target_power = power, alpha = alpha),
    notes = c("integer designs of the declared allocation; power is exact and not smoothed",
              if (fisher) "fisher_power: Fisher's conditional test, same convention, comparison only"))
}

#' Plot data: Boschloo power versus the alternative
#'
#' Exact power of one fixed design (one rejection region) over a set of
#' alternatives `(p1, p2)`. `p1` and `p2` are recycled to a common length,
#' so a scalar `p1` with a vector `p2` gives power versus the group-2
#' probability, and a vector `p1` with `p2 = p1 + d` gives power versus the
#' baseline at a fixed difference `d`.
#'
#' @inheritParams boschloo_power
#' @param p1,p2 Event probabilities of group 1 and group 2 (vectors,
#'   recycled to a common length).
#' @return An object of class `boschloo_plotdata` (plot id
#'   `"boschloo_power_effect"`).
#' @examples
#' pd <- boschloo_plotdata_power_effect(20, 20, p1 = 0.2, p2 = seq(0, 1, 0.1),
#'                                      variant = "central")
#' pd$data
#' base <- seq(0.05, 0.6, by = 0.05)
#' boschloo_plotdata_power_effect(20, 20, p1 = base, p2 = base + 0.3,
#'                                variant = "central")$data
#' @family Boschloo
#' @export
boschloo_plotdata_power_effect <- function(n1, n2, p1, p2, alpha = 0.05, variant,
                                           control = boschloo_control()) {
  info <- .bl_variant(if (missing(variant)) NULL else variant)
  .bl_check_design(n1, n2)
  probs <- function(x) is.numeric(x) && length(x) && !anyNA(x) && all(x >= 0 & x <= 1)
  if (!probs(p1)) .ig_abort_invalid("`p1` must be a vector of probabilities in [0, 1]")
  if (!probs(p2)) .ig_abort_invalid("`p2` must be a vector of probabilities in [0, 1]")
  len <- max(length(p1), length(p2))
  if (len %% length(p1) != 0L || len %% length(p2) != 0L) {
    .ig_abort_invalid("the lengths of `p1` and `p2` must be multiples of each other")
  }
  p1 <- rep_len(as.double(p1), len)
  p2 <- rep_len(as.double(p2), len)
  alpha <- .bl_check_alpha(alpha)
  control <- .bl_check_control(control)
  n1 <- as.integer(n1)
  n2 <- as.integer(n2)
  reg <- .bl_region(n1, n2, alpha, info, control)
  status <- .bl_region_status(reg)
  ok <- status$code == "ok"
  pw <- if (ok) mapply(function(a, b) .bl_power_idx(reg$idx, n1, n2, a, b), p1, p2) else NA_real_
  title <- if (length(unique(p1)) == 1L) {
    sprintf("Power versus p2 (p1 = %s), %s", format(p1[1]), info$label)
  } else {
    sprintf("Power versus the alternative (p1, p2), %s", info$label)
  }
  .bl_plotdata(
    "boschloo_power_effect", title,
    data.frame(p1 = p1, p2 = p2, difference = p2 - p1, power = pw),
    fixed = list(n1 = n1, n2 = n2, alpha = alpha, variant = info$variant,
                 variant_id = info$variant_id, p1 = if (length(unique(p1)) == 1L) p1[1]),
    scope = "hypothetical future study with fixed group sizes",
    reference = list(attained_size = if (ok) reg$size$value else NA_real_, alpha = alpha),
    notes = "where p2 = p1 the value is a null rejection probability (at most the attained size)",
    status = status)
}

#' Plot data: Boschloo power versus the allocation at a fixed total
#'
#' Exact power and attained size of every split `n1 + n2 = n_total` (from
#' `n1 = 1` to `n1 = n_total - 1`, or the splits in `n1`), each with its own
#' rejection region. Shows how unequal allocation changes power and how
#' discreteness makes neighbouring splits differ.
#'
#' @inheritParams boschloo_power
#' @param n_total Fixed total sample size (at least 2).
#' @param n1 Optional subset of group-1 sizes to evaluate.
#' @return An object of class `boschloo_plotdata` (plot id
#'   `"boschloo_power_allocation"`).
#' @examples
#' pd <- boschloo_plotdata_allocation(0.2, 0.6, n_total = 30, variant = "central",
#'                                    n1 = 8:22)
#' pd$data[which.max(pd$data$power), ]
#' @family Boschloo
#' @export
boschloo_plotdata_allocation <- function(p1, p2, n_total, alpha = 0.05, variant, n1 = NULL,
                                         control = boschloo_control()) {
  info <- .bl_variant(if (missing(variant)) NULL else variant)
  p1 <- .bl_check_prob(p1, "p1")
  p2 <- .bl_check_prob(p2, "p2")
  if (!.bl_is_count(n_total, 2) || n_total > 3000) {
    .ig_abort_invalid("`n_total` must be a whole number between 2 and 3000")
  }
  n_total <- as.integer(n_total)
  if (is.null(n1)) n1 <- seq_len(n_total - 1L)
  if (!is.numeric(n1) || !length(n1) || anyNA(n1) || any(n1 != round(n1)) || any(n1 < 1 | n1 > n_total - 1)) {
    .ig_abort_invalid("`n1` must be whole numbers between 1 and n_total - 1")
  }
  alpha <- .bl_check_alpha(alpha)
  control <- .bl_check_control(control)
  rows <- lapply(sort(unique(as.integer(n1))), function(a) {
    b <- n_total - a
    .bl_check_design(a, b)
    reg <- .bl_region(a, b, alpha, info, control)
    ok <- !length(reg$problems)
    data.frame(n1 = a, n2 = b, allocation_ratio = b / a,
               power = if (ok) .bl_power_idx(reg$idx, a, b, p1, p2) else NA_real_,
               size = if (ok) reg$size$value else NA_real_,
               n_rejected = if (ok) length(reg$idx) else NA_integer_,
               status = if (ok) "ok" else "numerical_failure", stringsAsFactors = FALSE)
  })
  data <- do.call(rbind, rows)
  .bl_plotdata(
    "boschloo_power_allocation", sprintf("Power versus allocation at N = %d, %s", n_total, info$label), data,
    fixed = list(p1 = p1, p2 = p2, n_total = n_total, alpha = alpha, variant = info$variant,
                 variant_id = info$variant_id),
    scope = "hypothetical future study with a fixed total sample size",
    reference = list(alpha = alpha, balanced_n1 = n_total %/% 2L),
    notes = c("each split has its own rejection region; power is exact and not smoothed",
              if (any(data$status != "ok")) "rows with status numerical_failure have no numbers"))
}
