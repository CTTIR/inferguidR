# Boschloo plot data ---------------------------------------------------------------------
#
# Plot-data helpers return the numbers behind the Boschloo panels (table
# lattice with rejection membership, nuisance envelope, power versus sample
# size and versus the alternative). Every number comes from the same
# ordering, supremum and region code as the observed test and the planning
# functions; nothing is smoothed. Each object states whether it concerns the
# observed data or a hypothetical future study, and lists its fixed inputs.

.bl_plotdata <- function(plot_id, title, data, fixed, scope, reference = list(), notes = character(),
                         status = ig_status("ok")) {
  structure(list(object_type = "boschloo_plotdata", plot_id = plot_id, title = title,
                 status = .ig_status_to_list(status), scope = scope, fixed_inputs = fixed,
                 data = data, reference = reference, notes = notes),
            class = "boschloo_plotdata")
}

#' @export
print.boschloo_plotdata <- function(x, ...) {
  cat(sprintf("<boschloo_plotdata> %s: %s\n", x$plot_id, x$title),
      sprintf("  scope:  %s\n", x$scope),
      sprintf("  status: %s\n", format(.ig_status_from_list(x$status))),
      sprintf("  rows:   %d\n", if (is.data.frame(x$data)) nrow(x$data) else 0L), sep = "")
  invisible(x)
}

# Boschloo p-value (value and upper bound) of every table: one supremum per
# level of each ordering side.
.bl_pvalues_all <- function(n1, n2, info, control) {
  per_side <- lapply(info$sides, function(side) {
    ord <- .bl_ordering(n1, n2, side, control)
    val <- numeric(ord$n_levels)
    upp <- numeric(ord$n_levels)
    conv <- logical(ord$n_levels)
    for (j in seq_len(ord$n_levels)) {
      s <- .bl_sup(.bl_tail_c(ord, .bl_level_tables(ord, j)), control)
      val[j] <- s$value
      upp[j] <- s$upper
      conv[j] <- isTRUE(s$converged)
    }
    list(value = val[ord$level], upper = upp[ord$level], converged = all(conv),
         ordering = exp(ord$lt), level = ord$level, ambiguity = !is.null(ord$ambiguity))
  })
  if (length(per_side) == 1L) {
    s <- per_side[[1]]
    return(list(value = s$value, upper = s$upper, ordering = s$ordering, level = s$level,
                converged = s$converged, ambiguity = s$ambiguity))
  }
  lo <- per_side[[1]]
  up <- per_side[[2]]
  list(value = pmin(1, 2 * pmin(lo$value, up$value)), upper = pmin(1, 2 * pmin(lo$upper, up$upper)),
       ordering = pmin(1, 2 * pmin(lo$ordering, up$ordering)), level = NULL,
       level_less = lo$level, level_greater = up$level,
       converged = lo$converged && up$converged, ambiguity = lo$ambiguity || up$ambiguity)
}

#' Plot data: lattice of possible tables with the Boschloo rejection region
#'
#' One row per possible table `(x1, x2)` of a planned design: the sample
#' proportions, the Fisher ordering statistic, the Boschloo p-value with its
#' upper bound, rejection membership at `alpha`, and the probability of the
#' table under the null at the size-maximising common probability and under
#' the alternative `(p1, p2)` when given. An `observed` table can be marked;
#' marking it does not change the planning inputs.
#'
#' @inheritParams boschloo_region
#' @param p1,p2 Optional alternative event probabilities.
#' @param observed Optional `c(x1, x2)` to mark.
#' @param pvalues Compute the Boschloo p-value of every table (one supremum
#'   per ordering level; skipped above 5000 tables).
#' @return An object of class `boschloo_plotdata` (plot id
#'   `"boschloo_lattice"`).
#' @examples
#' pd <- boschloo_plotdata_lattice(8, 10, alpha = 0.05, variant = "minlike",
#'                                 p1 = 0.2, p2 = 0.7)
#' head(pd$data)
#' @family Boschloo
#' @export
boschloo_plotdata_lattice <- function(n1, n2, alpha = 0.05, variant, p1 = NULL, p2 = NULL,
                                      observed = NULL, pvalues = TRUE,
                                      control = boschloo_control()) {
  info <- .bl_variant(if (missing(variant)) NULL else variant)
  .bl_check_design(n1, n2)
  alpha <- .bl_check_alpha(alpha)
  if (is.null(p1) != is.null(p2)) .ig_abort_invalid("give both `p1` and `p2`, or neither")
  if (!is.null(p1)) {
    p1 <- .bl_check_prob(p1, "p1")
    p2 <- .bl_check_prob(p2, "p2")
  }
  if (!is.null(observed) && (!is.numeric(observed) || length(observed) != 2L ||
                             !.bl_is_count(observed[1]) || !.bl_is_count(observed[2]) ||
                             observed[1] > n1 || observed[2] > n2)) {
    .ig_abort_invalid("`observed` must be c(x1, x2) with 0 <= x1 <= n1 and 0 <= x2 <= n2")
  }
  control <- .bl_check_control(control)
  n1 <- as.integer(n1)
  n2 <- as.integer(n2)
  reg <- .bl_region(n1, n2, alpha, info, control)
  status <- .bl_region_status(reg)
  ok <- status$code == "ok"
  n_tab <- length(reg$member)
  notes <- character()
  pv <- NULL
  if (pvalues && n_tab <= 5000L && ok) {
    pv <- .bl_pvalues_all(n1, n2, info, control)
  } else if (pvalues) {
    notes <- c(notes, if (!ok) "p-values not computed: the region has no ok status" else
      sprintf("p-values not computed for %d tables (limit 5000); membership is still exact", n_tab))
  }
  pi_star <- if (ok && is.finite(reg$size$argmax)) reg$size$argmax else 0.5
  data <- data.frame(
    x1 = reg$x1, x2 = reg$x2, total = reg$t,
    proportion_1 = reg$x1 / n1, proportion_2 = reg$x2 / n2,
    fisher_ordering = .bl_fisher_p(n1, n2, info, control),
    boschloo_p = if (is.null(pv)) NA_real_ else pv$value,
    boschloo_p_upper = if (is.null(pv)) NA_real_ else pv$upper,
    in_region = if (ok) reg$member else NA,
    prob_null = stats::dbinom(reg$x1, n1, pi_star) * stats::dbinom(reg$x2, n2, pi_star),
    prob_alternative = if (is.null(p1)) NA_real_ else
      stats::dbinom(reg$x1, n1, p1) * stats::dbinom(reg$x2, n2, p2),
    observed = if (is.null(observed)) FALSE else reg$x1 == observed[1] & reg$x2 == observed[2]
  )
  if (!is.null(pv) && ok) {
    consistent <- identical(data$in_region, data$boschloo_p_upper <= alpha)
    if (!consistent) {
      notes <- c(notes, "per-table bounds and region membership differ at the boundary (tables within the floating-point allowance of alpha)")
    }
  }
  .bl_plotdata(
    "boschloo_lattice", sprintf("Possible tables and %s rejection region", info$label), data,
    fixed = list(n1 = n1, n2 = n2, alpha = alpha, variant = info$variant,
                 variant_id = info$variant_id, p1 = p1, p2 = p2,
                 null_probability_at = pi_star),
    scope = if (is.null(observed)) "hypothetical future study (all possible tables)" else
      "hypothetical future study; the marked table is the observed one (marking does not change planning inputs)",
    reference = list(attained_size = if (ok) reg$size$value else NA_real_,
                     attained_size_upper = if (ok) reg$size$upper else NA_real_,
                     power = if (!is.null(p1) && ok) .bl_power_idx(reg$idx, n1, n2, p1, p2) else NA_real_,
                     n_rejected = if (ok) length(reg$idx) else NA_integer_),
    notes = notes, status = status)
}

#' Plot data: nuisance envelope
#'
#' The null probability of a tail set as a function of the common event
#' probability `pi`. For an observed test (a [boschloo_test()] result) the
#' curve is `P_pi(tables at least as extreme as the observed one)`, whose
#' maximum is the p-value (two curves for the central variant). For a
#' rejection region (a [boschloo_region()] result) it is the null rejection
#' probability, whose maximum is the attained size.
#'
#' @param x A `boschloo_test` or `boschloo_region` object.
#' @param grid_points Number of equally spaced values of `pi` (the
#'   maximiser is added).
#' @param control Numerical settings. By default they are taken from `x`; a
#'   supplied value must match the settings recorded in `x`.
#' @return An object of class `boschloo_plotdata` (plot id
#'   `"boschloo_nuisance"`).
#' @examples
#' res <- boschloo_test(3, 15, 9, 15, variant = "minlike")
#' env <- boschloo_plotdata_nuisance(res)
#' max(env$data$probability)
#' @family Boschloo
#' @export
boschloo_plotdata_nuisance <- function(x, grid_points = 401L, control = NULL) {
  if (!.bl_is_count(grid_points, 3) || grid_points > 1e5) {
    .ig_abort_invalid("`grid_points` must be a whole number between 3 and 100000")
  }
  if (!inherits(x, "boschloo_test") && !inherits(x, "boschloo_region")) {
    .ig_abort("invalid_input", "`x` must be a boschloo_test or boschloo_region object")
  }
  recorded <- .bl_control_from_record(x$settings)
  if (is.null(control)) {
    control <- recorded
  } else {
    control <- .bl_check_control(control)
    if (!identical(.bl_control_key(control), .bl_control_key(recorded))) {
      .ig_abort("invalid_input", "`control` differs from the settings recorded in `x`; omit it to use those")
    }
  }
  if (inherits(x, "boschloo_test")) {
    d <- x$data
    info <- .bl_variant_info[[x$variant$variant]]
    n1 <- as.integer(d$n1); n2 <- as.integer(d$n2)
    fixed <- list(x1 = d$x1, n1 = n1, x2 = d$x2, n2 = n2, variant = info$variant,
                  variant_id = info$variant_id)
    title <- "Null probability of the tail set versus the common event probability"
    scope <- "observed data: tail set of the observed table"
    status <- .ig_status_from_list(unclass(x$status))
    if (status$code != "ok") {
      # No curve for a p-value that could not be certified.
      return(.bl_plotdata("boschloo_nuisance", title,
                          data.frame(side = character(), pi = numeric(), probability = numeric()),
                          fixed = fixed, scope = scope, notes = status$reason, status = status))
    }
    i <- as.integer(d$x1) + 1L + as.integer(d$x2) * (n1 + 1L)
    rows <- list()
    marks <- list()
    for (side in info$sides) {
      ord <- .bl_ordering(n1, n2, side, control)
      cvec <- .bl_tail_c(ord, .bl_level_tables(ord, ord$level[i]))
      s <- .bl_sup(cvec, control)
      pis <- sort(unique(c(seq(0, 1, length.out = grid_points), if (is.finite(s$argmax)) s$argmax)))
      rows[[side]] <- data.frame(side = side, pi = pis, probability = .bl_null_curve(cvec, pis))
      marks[[side]] <- list(argmax = s$argmax, supremum = s$value, upper = s$upper)
    }
    data <- do.call(rbind, rows)
    rownames(data) <- NULL
    return(.bl_plotdata("boschloo_nuisance", title, data, fixed = fixed, scope = scope,
                        reference = list(p_value = x$p_value, p_upper = x$p_upper, sides = marks,
                                         construction = info$construction),
                        status = status))
  }
  info <- .bl_variant_info[[x$variant$variant]]
  n1 <- as.integer(x$design$n1); n2 <- as.integer(x$design$n2)
  reg <- .bl_region(n1, n2, x$alpha, info, control)
  status <- .bl_region_status(reg)
  fixed <- list(n1 = n1, n2 = n2, alpha = x$alpha, variant = info$variant,
                variant_id = info$variant_id)
  title <- "Null rejection probability versus the common event probability"
  scope <- "hypothetical future study: rejection region of the design"
  if (status$code != "ok") {
    # No curve for a region that could not be certified.
    return(.bl_plotdata("boschloo_nuisance", title,
                        data.frame(side = character(), pi = numeric(), probability = numeric()),
                        fixed = fixed, scope = scope, reference = list(alpha = x$alpha),
                        notes = status$reason, status = status))
  }
  ord1 <- .bl_ordering(n1, n2, info$sides[1], control)
  cvec <- .bl_tail_c(ord1, reg$idx)
  pis <- sort(unique(c(seq(0, 1, length.out = grid_points),
                       if (is.finite(reg$size$argmax)) reg$size$argmax)))
  data <- data.frame(side = "region", pi = pis, probability = .bl_null_curve(cvec, pis))
  .bl_plotdata("boschloo_nuisance", title, data, fixed = fixed, scope = scope,
               reference = list(alpha = x$alpha, attained_size = reg$size$value,
                                attained_size_upper = reg$size$upper, argmax = reg$size$argmax),
               status = status)
}

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
