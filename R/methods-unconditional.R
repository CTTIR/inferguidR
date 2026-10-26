# Method definitions: unconditional exact tests of two proportions (group G9) --------------
#
# Ownership: R/methods-unconditional.R belongs to task group G9 (phases 06-10 for
# exact_proportions_unconditional and exact_proportions_unconditional_offset).
# Method notes: specs/method-notes/exact_proportions_unconditional.md and
# specs/method-notes/exact_proportions_unconditional_offset.md.
#
# G*Power 3.1.9.7: Exact > "Proportions: Inequality, two independent groups
# (unconditional)" and "Proportions: Inequality (offset), two independent groups
# (unconditional)". The manual has no chapter on either test; the semantics below
# were determined from the certified live inventory and the emulator captures
# (validation/fixtures/emulator_capture, ..._supplement/g9) and are documented in
# the method notes. In short:
#
# * The test is the large-sample test of the chosen statistic (z pooled/unpooled
#   with or without continuity correction, Mantel-Haenszel, likelihood ratio, t,
#   and for the offset test the Gart-Nam, Farrington-Manning and
#   Miettinen-Nurminen score statistics) with the large-sample critical value
#   (normal quantile; t quantile with N - 2 df for the t test). "Exact" means
#   that its power and size are computed exactly by enumerating all
#   (n1 + 1)(n2 + 1) tables of two independent binomials; nothing is adjusted to
#   keep the size at or below alpha (G*Power's "Actual alpha" exceeds alpha in
#   many captures).
# * "Actual alpha" is the rejection probability at the null point with the
#   given p2 (p1 = p2 for the inequality test, p1 = p1|H0 for the offset test);
#   no supremum over the nuisance parameter is taken.
# * One-sided tests reject in the direction of the effect (sign of p1 - p2,
#   resp. of p1|H1 - p1|H0); in sensitivity mode the direction is the option
#   effect_direction. Two-sided tests reject |T| >= the alpha/2 quantile.
# * Large-sample computation: normal approximations of the same statistics
#   (noncentral t for the t test), with actual alpha = alpha.
# * The exact computation is used when N = n1 + n2 does not exceed the
#   threshold of the Options dialog ("Use exact distribution if N <", default
#   10000; G*Power computes N equal to the threshold exactly).

.ig_define_methods_exact_unconditional <- function() {
  list(.ig_method_exact_proportions_unconditional(),
       .ig_method_exact_proportions_unconditional_offset())
}

# Option codes and labels (registry codes; labels as in the certified inventory) ------------

.ig_unc_stat_labels <- c(
  z_pooled = "z-test pooled var",
  z_pooled_cc = "z-test pooled var with cont. corr.",
  z_unpooled = "z-test unpooled var",
  z_unpooled_cc = "z-test unpooled var with cont. corr.",
  mantel_haenszel = "Mantel-Haenszel test",
  likelihood_ratio = "Likelihood ratio",
  t_test = "t-test",
  score_gart_nam = "Likelihood score: (Gart & Nam)",
  score_farrington_manning = "Likelihood score: (Farrington & Manning)",
  score_miettinen_nurminen = "Likelihood score: (Miettinen & Nurminen)"
)

.ig_unc_index_labels <- c(proportion_p1 = "Proportion p1", difference = "Difference (p1-p2)",
                          risk_ratio = "Risk ratio (p1/p2)", odds_ratio = "Odds ratio")

# Statistics offered per effect-size index of the offset test (live list probe,
# validation/raw/emulator/cases/_layout_survey/statistic-lists.json): the list
# depends on the index; G*Power preselects the first item.
.ig_unc_offset_statistics <- list(
  proportion_p1 = c("z_pooled", "z_pooled_cc", "z_unpooled", "z_unpooled_cc", "t_test",
                    "score_gart_nam", "score_farrington_manning", "score_miettinen_nurminen"),
  difference = c("z_pooled", "z_pooled_cc", "z_unpooled", "z_unpooled_cc", "t_test",
                 "score_gart_nam", "score_farrington_manning", "score_miettinen_nurminen"),
  risk_ratio = c("score_gart_nam", "score_farrington_manning", "score_miettinen_nurminen"),
  odds_ratio = c("score_farrington_manning", "score_miettinen_nurminen")
)

# Statistics whose value is non-decreasing in x1 for every fixed x2 (proofs in
# the method notes): their rejection regions are found by a threshold search
# per x2 instead of evaluating every table. Only the inequality test (null
# difference 0) qualifies; the offset statistics are not monotone near the
# table boundaries and are always enumerated in full.
.ig_unc_monotone <- c("z_pooled", "z_pooled_cc", "z_unpooled", "z_unpooled_cc", "mantel_haenszel",
                      "likelihood_ratio", "t_test")

# Effect-size indices ---------------------------------------------------------------------

# p1 from an effect-size index value and p2 (NA when the index value has no
# proportion in (0, 1)).
.ig_unc_p1_from_index <- function(index, value, p2) {
  p1 <- switch(index,
               proportion_p1 = value,
               difference = p2 + value,
               risk_ratio = value * p2,
               odds_ratio = value * p2 / (1 - p2 + value * p2),
               stop("unknown effect-size index ", index, call. = FALSE))
  if (!is.finite(p1)) NA_real_ else p1
}

# All four indices of the pair (p1, p2).
.ig_unc_indices <- function(p1, p2) {
  list(proportion_p1 = p1, difference = p1 - p2, risk_ratio = p1 / p2,
       odds_ratio = (p1 / (1 - p1)) / (p2 / (1 - p2)))
}

# Restricted maximum likelihood estimates under a null on the difference,
# ratio or odds-ratio scale (vectorised over tables) -------------------------------------

# Difference p1 - p2 = d0: Farrington & Manning (1990), closed-form root of the
# likelihood cubic (Miettinen & Nurminen 1985, Appendix). Returns a two-column
# matrix (p1, p2).
.ig_unc_rmle_diff <- function(x1, x2, n1, n2, d0) {
  if (d0 == 0) {
    pp <- (x1 + x2) / (n1 + n2)
    return(unname(cbind(pp, pp)))
  }
  p1h <- x1 / n1
  p2h <- x2 / n2
  th <- n2 / n1
  a <- 1 + th
  b <- -(1 + th + p1h + th * p2h + d0 * (th + 2))
  c <- d0^2 + d0 * (2 * p1h + th + 1) + p1h + th * p2h
  d <- -p1h * d0 * (1 + d0)
  v <- b^3 / (3 * a)^3 - b * c / (6 * a^2) + d / (2 * a)
  u <- sign(v) * sqrt(pmax(b^2 / (3 * a)^2 - c / (3 * a), 0))
  arg <- ifelse(u == 0, 0, v / u^3)
  w <- (pi + acos(pmin(pmax(arg, -1), 1))) / 3
  r1 <- 2 * u * cos(w) - b / (3 * a)
  r1 <- pmin(pmax(r1, max(0, d0)), min(1, 1 + d0))
  unname(cbind(r1, r1 - d0))
}

# Ratio p1 / p2 = phi: smaller root of N phi p^2 - (phi n1 + x1 + n2 + phi x2) p
# + (x1 + x2) = 0 for p = p2 (Farrington & Manning 1990).
.ig_unc_rmle_ratio <- function(x1, x2, n1, n2, phi) {
  A <- (n1 + n2) * phi
  B <- phi * n1 + x1 + n2 + phi * x2
  C <- x1 + x2
  disc <- sqrt(pmax(B^2 - 4 * A * C, 0))
  # Rationalised smaller root (no cancellation): 2C / (B + sqrt(B^2 - 4AC)).
  p2 <- ifelse(C == 0, 0, 2 * C / (B + disc))
  unname(cbind(phi * p2, p2))
}

# Odds ratio psi: the score equation n1 p1 + n2 p2 = x1 + x2 with
# p1 = psi p2 / (1 - p2 + psi p2) (Miettinen & Nurminen 1985), i.e.
# n2 (psi - 1) p^2 + (n1 psi + n2 - m (psi - 1)) p - m = 0 for p = p2.
.ig_unc_rmle_or <- function(x1, x2, n1, n2, psi) {
  m <- x1 + x2
  if (psi == 1) {
    pp <- m / (n1 + n2)
    return(unname(cbind(pp, pp)))
  }
  A <- n2 * (psi - 1)
  B <- n1 * psi + n2 - m * (psi - 1)
  disc <- sqrt(pmax(B^2 + 4 * A * m, 0))
  p2 <- ifelse(m == 0, 0, ifelse(B > 0, 2 * m / (B + disc), (disc - B) / (2 * A)))
  p2 <- pmin(pmax(p2, 0), 1)
  unname(cbind(psi * p2 / (1 - p2 + psi * p2), p2))
}

# Test statistics ---------------------------------------------------------------------------
#
# .ig_unc_stat(code, x1, x2, n1, n2, h0) returns the statistic for the tables
# (x1, x2) (vectors of equal length). `h0` is list(scale, value) with scale
# "difference" (value d0 = p1 - p2 under H0; 0 for the inequality test),
# "ratio" (phi0) or "odds_ratio" (psi0). Conventions (all confirmed by the
# captures, see the method notes): a statistic whose numerator is 0 is 0; a
# non-zero numerator over a zero standard error is +-Inf (the table is in the
# rejection region of its side).

.ig_unc_finish <- function(num, se) {
  z <- num / se
  z[num == 0] <- 0
  inf <- num != 0 & se == 0
  if (any(inf)) z[inf] <- sign(num[inf]) * Inf
  z
}

.ig_unc_count_square <- function(numerator, denominator, multiplier) {
  # Reducing the fraction before multiplying by N or N-1 gives identical
  # floating-point operands for equal rational statistics. This also avoids
  # rounding a large product before a cancellation (notably for Mantel).
  a <- numerator^2
  if (any(a > 2^53 | denominator > 2^53)) {
    .ig_abort("numerical_failure", paste(
      "exact pooled count products exceed the exactly representable integer range;",
      "reduce the design or explicitly select the large_sample computation"))
  }
  # Corrected numerators may be half integers; multiplying numerator and
  # denominator by four would exceed the exact range earlier. Instead scale
  # the numerator itself by two for the gcd, retaining its fractional part.
  half <- numerator != floor(numerator)
  a[half] <- (2 * numerator[half])^2
  denominator[half] <- 4 * denominator[half]
  if (any(a > 2^53 | denominator > 2^53)) {
    .ig_abort("numerical_failure", "exact corrected pooled count products exceed the exactly representable integer range")
  }
  gcd <- function(a, b) {
    repeat {
      open <- b != 0
      if (!any(open)) break
      remainder <- a[open] %% b[open]
      a[open] <- b[open]
      b[open] <- remainder
    }
    a
  }
  zero <- a == 0
  denominator[zero] <- 1
  divisor <- gcd(a, denominator)
  a <- a / divisor; b <- denominator / divisor
  factor <- gcd(rep(multiplier, length(b)), b)
  out <- (a * (multiplier / factor)) / (b / factor)
  out[zero] <- 0
  out
}

.ig_unc_stat <- function(code, x1, x2, n1, n2, h0) {
  n1 <- as.double(n1)
  n2 <- as.double(n2)
  x1 <- as.double(x1)
  x2 <- as.double(x2)
  N <- n1 + n2
  k <- 1 / n1 + 1 / n2
  p1 <- x1 / n1
  p2 <- x2 / n2
  scale <- h0$scale
  h <- h0$value
  if (scale == "difference") {
    if (h == 0 && code %in% c("z_pooled", "z_pooled_cc", "mantel_haenszel")) {
      # Keep algebraically tied tables tied at a discrete critical boundary.
      # Forming proportions and their standard errors separately can split
      # one equality class by several ulps, including only some tied tables.
      cross <- x1 * n2 - x2 * n1
      total <- x1 + x2
      numerator <- if (code == "z_pooled_cc") pmax(abs(cross) - N / 2, 0) else abs(cross)
      multiplier <- if (code == "mantel_haenszel") N - 1 else N
      squared <- .ig_unc_count_square(numerator, n1 * n2 * total * (N - total), multiplier)
      return(sign(cross) * sqrt(squared))
    }
    # The difference is formed in integer arithmetic so that equal proportions
    # give an exact 0 (x1 n2 - x2 n1 is an exact integer in double precision).
    num <- (x1 * n2 - x2 * n1) / (n1 * n2)
    if (h != 0) num <- num - h
    cc_num <- function() {
      if (h == 0) {
        a <- abs(x1 * n2 - x2 * n1) - (n1 + n2) / 2
        sign(num) * pmax(a, 0) / (n1 * n2)
      } else {
        sign(num) * pmax(abs(num) - k / 2, 0)
      }
    }
    switch(code,
           z_pooled = {
             pb <- (x1 + x2) / N
             .ig_unc_finish(num, sqrt(pb * (1 - pb) * k))
           },
           z_pooled_cc = {
             pb <- (x1 + x2) / N
             .ig_unc_finish(cc_num(), sqrt(pb * (1 - pb) * k))
           },
           z_unpooled = .ig_unc_finish(num, sqrt(p1 * (1 - p1) / n1 + p2 * (1 - p2) / n2)),
           z_unpooled_cc = .ig_unc_finish(cc_num(), sqrt(p1 * (1 - p1) / n1 + p2 * (1 - p2) / n2)),
           mantel_haenszel = {
             pb <- (x1 + x2) / N
             .ig_unc_finish(num, sqrt(pb * (1 - pb) * k)) * sqrt((N - 1) / N)
           },
           likelihood_ratio = {
             pb <- (x1 + x2) / N
             g <- .ig_unc_xlogx(x1, n1 * pb) + .ig_unc_xlogx(n1 - x1, n1 * (1 - pb)) +
               .ig_unc_xlogx(x2, n2 * pb) + .ig_unc_xlogx(n2 - x2, n2 * (1 - pb))
             sign(num) * sqrt(pmax(2 * g, 0))
           },
           t_test = {
             s2 <- (x1 * (n1 - x1) / n1 + x2 * (n2 - x2) / n2) / (N - 2)
             .ig_unc_finish(num, sqrt(s2 * k))
           },
           score_farrington_manning = ,
           score_miettinen_nurminen = {
             r <- .ig_unc_rmle_diff(x1, x2, n1, n2, h)
             v <- r[, 1] * (1 - r[, 1]) / n1 + r[, 2] * (1 - r[, 2]) / n2
             if (code == "score_miettinen_nurminen") v <- v * N / (N - 1)
             .ig_unc_finish(num, sqrt(v))
           },
           score_gart_nam = .ig_unc_stat_gart_nam(num, x1, x2, n1, n2, h0),
           stop("statistic ", code, " is not defined on the difference scale", call. = FALSE))
  } else if (scale == "ratio") {
    num <- p1 - h * p2
    switch(code,
           score_farrington_manning = ,
           score_miettinen_nurminen = {
             r <- .ig_unc_rmle_ratio(x1, x2, n1, n2, h)
             v <- r[, 1] * (1 - r[, 1]) / n1 + h^2 * r[, 2] * (1 - r[, 2]) / n2
             if (code == "score_miettinen_nurminen") v <- v * N / (N - 1)
             .ig_unc_finish(num, sqrt(v))
           },
           score_gart_nam = .ig_unc_stat_gart_nam(num, x1, x2, n1, n2, h0),
           stop("statistic ", code, " is not defined on the ratio scale", call. = FALSE))
  } else if (scale == "odds_ratio") {
    r <- .ig_unc_rmle_or(x1, x2, n1, n2, h)
    num <- x1 - n1 * r[, 1]
    # The score is formed from rounded MLEs; values within rounding of zero
    # (m = 0 or m = N, or psi = 1 with x1 / n1 = x2 / n2) are exact zeros.
    num[abs(num) <= 64 * .Machine$double.eps * pmax(1, abs(x1))] <- 0
    switch(code,
           score_farrington_manning = ,
           score_miettinen_nurminen = {
             a1 <- n1 * r[, 1] * (1 - r[, 1])
             a2 <- n2 * r[, 2] * (1 - r[, 2])
             v <- ifelse(a1 > 0 & a2 > 0, a1 * a2 / (a1 + a2), 0)
             if (code == "score_miettinen_nurminen") v <- v * N / (N - 1)
             .ig_unc_finish(num, sqrt(v))
           },
           stop("statistic ", code, " is not defined on the odds-ratio scale", call. = FALSE))
  } else {
    stop("unknown null scale ", scale, call. = FALSE)
  }
}

.ig_unc_xlogx <- function(o, e) ifelse(o > 0, o * log(o / e), 0)

# Gart & Nam (1988, ratio; 1990, difference): the skewness-corrected score
# statistic. With the Farrington-Manning score z (restricted MLE variance v)
# and the skewness gamma = mu3 / v^(3/2) of the numerator at the restricted
# MLE, mu3 = p1 q1 (q1 - p1) / n1^2 - phi^3 p2 q2 (q2 - p2) / n2^2 (phi = 1 for
# the difference, phi0 for the ratio), the corrected statistic S solves
# z = S + gamma (S^2 - 1) / 6, i.e. S = 2 (z + k) / (1 + sqrt(1 + 4 k (z + k)))
# with k = gamma / 6 (identified from the boundaries of the supplement
# criterion captures SUP-G9-OFF-CRI-B02..B13 to about 1e-6; method notes).
# A negative discriminant leaves S undefined (NaN; the table is not
# rejected); a zero variance follows the conventions of .ig_unc_finish().
.ig_unc_stat_gart_nam <- function(num, x1, x2, n1, n2, h0) {
  if (identical(h0$scale, "difference")) {
    r <- .ig_unc_rmle_diff(x1, x2, n1, n2, h0$value)
    phi <- 1
  } else {
    r <- .ig_unc_rmle_ratio(x1, x2, n1, n2, h0$value)
    phi <- h0$value
  }
  a <- r[, 1]
  b <- r[, 2]
  v <- a * (1 - a) / n1 + phi^2 * b * (1 - b) / n2
  z <- .ig_unc_finish(num, sqrt(v))
  mu3 <- a * (1 - a) * (1 - 2 * a) / n1^2 - phi^3 * b * (1 - b) * (1 - 2 * b) / n2^2
  k <- ifelse(v > 0, mu3 / v^1.5 / 6, 0)
  disc <- 1 + 4 * k * (z + k)
  S <- ifelse(disc >= 0, 2 * (z + k) / (1 + sqrt(pmax(disc, 0))), NaN)
  fin <- !is.finite(z)
  S[fin] <- z[fin]
  S[v == 0 & num == 0] <- 0
  S
}

# Critical value of the chosen statistic.
.ig_unc_critical <- function(code, alpha, tails, N) {
  a <- if (identical(tails, "two")) alpha / 2 else alpha
  if (identical(code, "t_test")) stats::qt(a, N - 2, lower.tail = FALSE) else stats::qnorm(a, lower.tail = FALSE)
}

# Exact rejection and acceptance probabilities --------------------------------------------
#
# The rejection region of one design is a set of tables (x1, x2): {T >= c}
# (one-sided, upper), {T <= -c} (one-sided, lower) or {T >= c} u {T <= -c}
# (two-sided). It is stored per x2 as runs of x1 values (from a threshold
# search for statistics that are monotone in x1, otherwise from the full table
# of statistic values), so that the probability of the region and of its
# complement under any (p1, p2) costs one pbinom() evaluation per run. For
# every probability pair (p1, p2) in `pairs` the rejection probability and the
# acceptance probability are each summed over their own tables (a small beta
# keeps its relative accuracy).

# Fixed numerical settings of the exact engine (reported in diagnostics).
.ig_unc_settings <- list(
  # Statistic values of designs with at most this many tables are kept (one
  # design at a time) and reused when only the critical value changes
  # (criterion and compromise searches).
  cache_tables = 2^22,
  # Budget of the exhaustive a priori search: table evaluations of the full
  # enumeration plus statistic evaluations of the threshold search, summed
  # over all designs evaluated.
  a_priori_budget = 1e8,
  # Designs evaluated one by one in the large-sample part of the a priori
  # search before the search jumps to a window below the continuous solution.
  large_sample_scan = 2e5,
  large_sample_window = 200L
)

.ig_unc_cache <- new.env(parent = emptyenv())

# Budget of the exhaustive a priori search (a function, so that tests can
# lower it).
.ig_unc_budget <- function() .ig_unc_settings$a_priori_budget

.ig_unc_exact_probs <- function(code, n1, n2, h0, crit, side, pairs) {
  reg <- .ig_unc_region(code, as.integer(n1), as.integer(n2), h0, crit, side)
  pr <- .ig_unc_region_probs(reg, as.integer(n1), as.integer(n2), pairs)
  c(pr, list(method = reg$method, tables_rejected = reg$tables_rejected))
}

# Whether the region of a statistic is found by the threshold search.
.ig_unc_is_monotone <- function(code, h0) {
  code %in% .ig_unc_monotone && identical(h0$scale, "difference") && h0$value == 0
}

# Cost of evaluating one design exactly (table or statistic evaluations).
.ig_unc_exact_cost <- function(code, n1, n2, h0) {
  n1 <- as.double(n1)
  n2 <- as.double(n2)
  if (.ig_unc_is_monotone(code, h0)) (n2 + 1) * (ceiling(log2(n1 + 2)) + 1) else (n1 + 1) * (n2 + 1)
}

.ig_unc_region <- function(code, n1, n2, h0, crit, side) {
  key <- paste(code, n1, n2, h0$scale, sprintf("%.17g", h0$value), sprintf("%.17g", crit), side)
  hit <- .ig_unc_cache$region
  if (!is.null(hit) && identical(hit$key, key)) return(hit)
  reg <- if (.ig_unc_is_monotone(code, h0)) {
    .ig_unc_region_threshold(code, n1, n2, h0, crit, side)
  } else {
    .ig_unc_region_full(code, n1, n2, h0, crit, side)
  }
  reg$key <- key
  .ig_unc_cache$region <- reg
  reg
}

# Statistic values of all tables, as an (n1 + 1) x (n2 + 1) matrix; kept for
# the last design when it has at most cache_tables tables.
.ig_unc_stat_matrix <- function(code, n1, n2, h0, idx = 0:n2) {
  m1 <- n1 + 1L
  matrix(.ig_unc_stat(code, rep(0:n1, length(idx)), rep(idx, each = m1), n1, n2, h0), nrow = m1)
}

.ig_unc_region_full <- function(code, n1, n2, h0, crit, side) {
  m1 <- n1 + 1L
  cells <- as.double(m1) * (n2 + 1)
  key <- paste(code, n1, n2, h0$scale, sprintf("%.17g", h0$value))
  test <- function(z) switch(side, upper = z >= crit, lower = z <= -crit, both = z >= crit | z <= -crit)
  runs <- list()
  nrej <- 0
  if (cells <= .ig_unc_settings$cache_tables) {
    hit <- .ig_unc_cache$stat
    Tm <- if (!is.null(hit) && identical(hit$key, key)) hit$T else {
      v <- .ig_unc_stat_matrix(code, n1, n2, h0)
      .ig_unc_cache$stat <- list(key = key, T = v)
      v
    }
    rm <- test(Tm)
    runs <- list(.ig_unc_runs(rm, 0L))
    nrej <- sum(rm)
  } else {
    block <- max(1L, as.integer(floor(2^20 / m1)))
    start <- 0L
    while (start <= n2) {
      idx <- start:min(n2, start + block - 1L)
      rm <- test(.ig_unc_stat_matrix(code, n1, n2, h0, idx))
      runs[[length(runs) + 1L]] <- .ig_unc_runs(rm, start)
      nrej <- nrej + sum(rm)
      start <- start + block
    }
  }
  runs <- do.call(rbind, runs)
  list(type = "runs", runs = runs, tables_rejected = nrej, method = "full enumeration of all tables")
}

# Runs of TRUE (rejected) and FALSE (accepted) values in the columns of a
# logical matrix whose first column is x2 = offset. Returns a data frame with
# columns x2, from, to (x1 values) and rejected.
.ig_unc_runs <- function(rm, offset) {
  m1 <- nrow(rm)
  m <- ncol(rm)
  v <- as.vector(rm)
  x1 <- rep.int(seq_len(m1) - 1L, m)
  x2 <- rep(seq_len(m) - 1L + offset, each = m1)
  first <- c(TRUE, v[-1L] != v[-length(v)] | x1[-1L] == 0L)
  last <- c(first[-1L], TRUE)
  data.frame(x2 = x2[first], from = x1[first], to = x1[last], rejected = v[first])
}

# Threshold search (statistics non-decreasing in x1 for fixed x2). For each
# x2, k_up = min{x1: T >= c} (n1 + 1 when empty) and k_lo = max{x1: T <= -c}
# (-1 when empty) are found by a vectorised bisection over x1.
.ig_unc_region_threshold <- function(code, n1, n2, h0, crit, side) {
  x2 <- 0:n2
  first_ge <- function(cond) {
    # smallest x1 in 0..n1 with cond(T(x1, x2)) TRUE, where cond is monotone
    # FALSE...TRUE in x1; n1 + 1 when never TRUE.
    lo <- rep(-1L, n2 + 1L)
    hi <- rep(n1 + 1L, n2 + 1L)
    repeat {
      open <- hi - lo > 1L
      if (!any(open)) break
      mid <- (lo[open] + hi[open]) %/% 2L
      ok <- cond(.ig_unc_stat(code, mid, x2[open], n1, n2, h0))
      hi[open][ok] <- mid[ok]
      lo[open][!ok] <- mid[!ok]
    }
    hi
  }
  k_up <- if (side %in% c("upper", "both")) first_ge(function(z) z >= crit) else rep(n1 + 1L, n2 + 1L)
  # max{x1: T <= -c} = (min{x1: T > -c}) - 1
  k_lo <- if (side %in% c("lower", "both")) first_ge(function(z) z > -crit) - 1L else rep(-1L, n2 + 1L)
  # A two-sided region with c <= 0 cannot occur (alpha / 2 < 0.5); a
  # one-sided region uses one of the two thresholds only.
  rows <- list(
    data.frame(x2 = x2, from = 0L, to = k_lo, rejected = TRUE),
    data.frame(x2 = x2, from = k_lo + 1L, to = k_up - 1L, rejected = FALSE),
    data.frame(x2 = x2, from = k_up, to = n1, rejected = TRUE))
  runs <- do.call(rbind, rows)
  runs <- runs[runs$from <= runs$to, , drop = FALSE]
  list(type = "runs", runs = runs,
       tables_rejected = sum(pmax(k_lo + 1L, 0L)) + sum(pmax(n1 - k_up + 1L, 0L)),
       method = "threshold search per x2 (statistic non-decreasing in x1; see the method notes)")
}

# P(a <= X <= b) for X ~ Bin(n, p), from whichever tail representation avoids
# cancellation.
.ig_unc_binom_interval <- function(a, b, n, p) {
  Fa <- stats::pbinom(a - 1L, n, p)
  out <- ifelse(Fa < 0.5,
                stats::pbinom(b, n, p) - Fa,
                stats::pbinom(a - 1L, n, p, lower.tail = FALSE) - stats::pbinom(b, n, p, lower.tail = FALSE))
  out[b < a] <- 0
  pmax(out, 0)
}

.ig_unc_region_probs <- function(reg, n1, n2, pairs) {
  r <- reg$runs
  rej <- numeric(length(pairs))
  acc <- numeric(length(pairs))
  for (j in seq_along(pairs)) {
    p1 <- pairs[[j]][1]
    p2 <- pairs[[j]][2]
    pr <- .ig_unc_binom_interval(r$from, r$to, n1, p1) * stats::dbinom(r$x2, n2, p2)
    rej[j] <- sum(pr[r$rejected])
    acc[j] <- sum(pr[!r$rejected])
  }
  list(reject = pmin(rej, 1), accept = pmin(acc, 1))
}

# Large-sample power ------------------------------------------------------------------------
#
# Normal approximations of the chosen statistic under H1 (noncentral t for the
# t test), as G*Power computes them (forms identified from the captures; see
# the method notes). Two forms occur:
#
#   form A  power = Phi((d - cc - c sd0) / sd1)      z tests, LR, FM/MN (difference, ratio)
#   form B  power = Phi(d / sd0 - c)                 Mantel-Haenszel, FM/MN (odds ratio)
#
# The Gart-Nam statistic uses form A with sd0 = delta / S_E, S_E the statistic
# at the expected H1 table (see gart_nam() below).
#
# with d = dir * (expected numerator under H1), sd0 the null standard error of
# the numerator (from the H1 proportions, or the restricted MLE evaluated at
# the expected H1 table), sd1 its standard error under H1 and cc the
# continuity correction (1/n1 + 1/n2) / 2. Two-sided tests add the opposite
# rejection side (-d instead of d). Returns power and beta (beta computed as
# the non-rejection probability of the main side minus the opposite side).

.ig_unc_large_sample <- function(code, n1, n2, p1, p2, h0, alpha, tails, dir) {
  n1 <- as.double(n1)
  n2 <- as.double(n2)
  N <- n1 + n2
  k <- 1 / n1 + 1 / n2
  q1 <- 1 - p1
  q2 <- 1 - p2
  two <- identical(tails, "two")
  c <- .ig_unc_critical(code, alpha, tails, N)
  scale <- h0$scale
  h <- h0$value
  sides <- function(z_main, z_other) {
    z_main <- unname(z_main)
    z_other <- unname(z_other)
    pw_main <- stats::pnorm(z_main)
    beta_main <- stats::pnorm(z_main, lower.tail = FALSE)
    if (!two) return(c(power = pw_main, beta = beta_main))
    pw_other <- stats::pnorm(z_other)
    c(power = pw_main + pw_other, beta = max(0, beta_main - pw_other))
  }
  form_a <- function(delta, sd0, sd1, cc = 0) {
    d <- dir * delta
    sides((d - cc - c * sd0) / sd1, (-d - cc - c * sd0) / sd1)
  }
  form_b <- function(delta, sd0) {
    d <- dir * delta
    sides(d / sd0 - c, -d / sd0 - c)
  }
  mn <- if (identical(code, "score_miettinen_nurminen")) sqrt(N / (N - 1)) else 1
  # Gart-Nam: G*Power evaluates the statistic at the expected H1 table
  # (x1 = n1 p1, x2 = n2 p2), S_E, and uses form A with the null standard
  # error replaced by delta / S_E, the value for which the statistic of form A
  # at the expected table is S_E: power = Phi((d - c delta / S_E) / sd1).
  # Identified from the captures GPW-EX-UNCOFF-POS-B157/B160/B163/B349/B352/B355
  # (O01, E01) and SUP-G9-OFF-POS-B18..B23, which it reproduces to every
  # displayed digit. The form has no value when S_E is undefined or zero, when
  # S_E and delta have opposite signs, or when delta = 0 (the ratio delta / S_E
  # is then not a standard error); it is unreliable when the effect is within
  # about |skewness| / 6 null standard errors of the null (method notes).
  gart_nam <- function(delta, sd1) {
    S <- .ig_unc_stat("score_gart_nam", n1 * p1, n2 * p2, n1, n2, h0)
    if (delta == 0 || !is.finite(S) || S == 0 || sign(S) != sign(delta)) {
      return(structure(c(power = NA_real_, beta = NA_real_), reason = sprintf(paste(
        "G*Power's large-sample form of the Gart-Nam statistic has no value here: it divides the effect",
        "(%s) by the statistic at the expected H1 table (%s), which must be finite, non-zero and on the",
        "side of the effect"), format(delta), format(S))))
    }
    form_a(delta, delta / S, sd1)
  }
  if (scale == "difference") {
    delta <- p1 - p2 - h
    sd1 <- sqrt(p1 * q1 / n1 + p2 * q2 / n2)
    pb <- (n1 * p1 + n2 * p2) / N
    sd_pooled <- sqrt(pb * (1 - pb) * k)
    switch(code,
           z_pooled = ,
           likelihood_ratio = form_a(delta, sd_pooled, sd1),
           z_pooled_cc = form_a(delta, sd_pooled, sd1, cc = k / 2),
           z_unpooled = form_a(delta, sd1, sd1),
           z_unpooled_cc = form_a(delta, sd1, sd1, cc = k / 2),
           mantel_haenszel = form_b(delta, sd_pooled * sqrt(N / (N - 1))),
           t_test = {
             st <- sqrt((n1 * p1 * q1 + n2 * p2 * q2) / (N - 2) * k)
             ncp <- dir * delta / st
             df <- N - 2
             pw <- stats::pt(c, df, ncp = ncp, lower.tail = FALSE)
             beta <- stats::pt(c, df, ncp = ncp)
             if (two) {
               other <- stats::pt(-c, df, ncp = ncp)
               pw <- pw + other
               beta <- max(0, beta - other)
             }
             c(power = pw, beta = beta)
           },
           score_farrington_manning = ,
           score_miettinen_nurminen = {
             r <- .ig_unc_rmle_diff(n1 * p1, n2 * p2, n1, n2, h)
             sd0 <- sqrt(r[1, 1] * (1 - r[1, 1]) / n1 + r[1, 2] * (1 - r[1, 2]) / n2) * mn
             form_a(delta, sd0, sd1)
           },
           score_gart_nam = gart_nam(delta, sd1),
           stop("statistic ", code, " has no large-sample form on the difference scale", call. = FALSE))
  } else if (scale == "ratio") {
    delta <- p1 - h * p2
    sd1 <- sqrt(p1 * q1 / n1 + h^2 * p2 * q2 / n2)
    switch(code,
           score_farrington_manning = ,
           score_miettinen_nurminen = {
             r <- .ig_unc_rmle_ratio(n1 * p1, n2 * p2, n1, n2, h)
             sd0 <- sqrt(r[1, 1] * (1 - r[1, 1]) / n1 + h^2 * r[1, 2] * (1 - r[1, 2]) / n2) * mn
             form_a(delta, sd0, sd1)
           },
           score_gart_nam = gart_nam(delta, sd1),
           stop("statistic ", code, " has no large-sample form on the risk-ratio scale", call. = FALSE))
  } else {
    switch(code,
           score_farrington_manning = ,
           score_miettinen_nurminen = {
             r <- .ig_unc_rmle_or(n1 * p1, n2 * p2, n1, n2, h)
             a1 <- n1 * r[1, 1] * (1 - r[1, 1])
             a2 <- n2 * r[1, 2] * (1 - r[1, 2])
             form_b(n1 * (p1 - r[1, 1]), sqrt(a1 * a2 / (a1 + a2)) * mn)
           },
           stop("statistic ", code, " has no large-sample form on the odds-ratio scale", call. = FALSE))
  }
}

# Method definitions ------------------------------------------------------------------------

.ig_unc_common_params <- function() {
  list(
    ig_param_alpha(),
    ig_param_power(),
    ig_param_error_ratio(),
    ig_param("n_total", role = "design", domain = ig_domain_integer(lower = 2), unit = "subjects",
             label = list(en = "Total sample size", de = "Gesamtstichprobenumfang"),
             description = "Solved in a priori mode (n1 + n2).",
             modes = "a_priori"),
    ig_param("n1", role = "design", domain = ig_domain_integer(lower = 1), unit = "subjects",
             label = list(en = "Sample size group 1", de = "Stichprobenumfang Gruppe 1"),
             aliases = "sample_size_group_1",
             modes = c("post_hoc", "sensitivity", "criterion", "compromise")),
    ig_param("n2", role = "design", domain = ig_domain_integer(lower = 1), unit = "subjects",
             label = list(en = "Sample size group 2", de = "Stichprobenumfang Gruppe 2"),
             aliases = "sample_size_group_2",
             modes = c("post_hoc", "sensitivity", "criterion", "compromise")),
    ig_param("allocation_ratio", role = "design",
             domain = ig_domain_real(0, Inf, lower_inclusive = FALSE, upper_inclusive = FALSE),
             default = 1,
             label = list(en = "Allocation ratio N2/N1", de = "Allokationsverh\u00e4ltnis N2/N1"),
             description = paste("Group-size ratio n2/n1 of the a priori design: n1 = 1, 2, ... and",
                                 "n2 = floor(r n1 + 0.5) (G*Power's rule, determined from the captures);",
                                 "default 1 as in G*Power's main window."),
             aliases = "allocation_ratio_n2_n1",
             modes = "a_priori"),
    ig_param("exact_n_limit", role = "control", domain = ig_domain_integer(lower = 1),
             default = 10000L, unit = "subjects",
             label = list(en = "Use exact distribution if N <", de = "Exakte Verteilung verwenden, wenn N <"),
             description = paste("Options dialog threshold 'Use exact distribution if N <': the exact enumeration",
                                 "is used for designs with N = n1 + n2 up to and including this value (G*Power's",
                                 "behaviour, captures SUP-G9-UNC-POS-A05/A06/A07), the large-sample approximation",
                                 "otherwise (G*Power default 10000)."),
             applies_when = list(option = "computation", "in" = "exact"))
  )
}

.ig_unc_common_options <- function(statistics, index_default, direction_labels) {
  ig_options(
    ig_option("tails", values = c("one", "two"), default = "one", labels = c(one = "One", two = "Two"),
              label = list(en = "Tail(s)", de = "Seiten")),
    ig_option("computation", values = c("exact", "large_sample"), default = "exact",
              labels = c(exact = "Use exact distribution if N < x", large_sample = "Use large sample approximation"),
              label = list(en = "Computation", de = "Berechnung")),
    ig_option("statistic", values = statistics, default = "z_pooled",
              labels = .ig_unc_stat_labels[statistics],
              label = list(en = "Statistic", de = "Teststatistik")),
    ig_option("effect_size_index", values = names(.ig_unc_index_labels), default = index_default,
              labels = .ig_unc_index_labels,
              label = list(en = "Effect size index", de = "Effektgr\u00f6\u00dfenindex")),
    ig_option("effect_direction", values = c("p_h1_le_p_h0", "p_h1_ge_p_h0"), default = "p_h1_ge_p_h0",
              labels = direction_labels, applies_when = list(mode = "sensitivity"),
              label = list(en = "Effect direction", de = "Effektrichtung"))
  )
}

.ig_method_exact_proportions_unconditional <- function() {
  # In sensitivity mode the solve target is p1 for every index (the index
  # only selects which transformation of the solved p1 is displayed), so the
  # other index parameters are inputs of the remaining modes only.
  index_param <- function(name, index, domain, label, description, modes = NULL) {
    ig_param(name, role = "effect", domain = domain, label = label, description = description,
             applies_when = if (name == "p1") list(any = list(list(mode = "sensitivity"),
               list(option = "effect_size_index", "in" = index))) else
               list(option = "effect_size_index", "in" = index), modes = modes)
  }
  not_sens <- c("a_priori", "post_hoc", "criterion", "compromise")
  ig_method(
    id = "exact_proportions_unconditional",
    family = "exact",
    label = list(en = "Proportions: Inequality, two independent groups (unconditional)",
                 de = "Anteile: Ungleichheit, zwei unabh\u00e4ngige Gruppen (unbedingt)"),
    params = ig_params(
      index_param("p1", "proportion_p1", ig_domain_real(0, 1, FALSE, FALSE),
                  list(en = "Proportion p1", de = "Anteil p1"),
                  "Event probability of group 1 under H1 (effect size index 'Proportion p1')."),
      index_param("difference", "difference", ig_domain_real(-1, 1, FALSE, FALSE),
                  list(en = "Difference (p1-p2)", de = "Differenz (p1-p2)"),
                  "p1 - p2 under H1 (effect size index 'Difference (p1-p2)').", not_sens),
      index_param("risk_ratio", "risk_ratio", ig_domain_real(0, Inf, FALSE, FALSE),
                  list(en = "Risk ratio (p1/p2)", de = "Risikoverh\u00e4ltnis (p1/p2)"),
                  "p1 / p2 under H1 (effect size index 'Risk ratio (p1/p2)').", not_sens),
      index_param("odds_ratio", "odds_ratio", ig_domain_real(0, Inf, FALSE, FALSE),
                  list(en = "Odds ratio", de = "Odds Ratio"),
                  "(p1 / (1 - p1)) / (p2 / (1 - p2)) under H1 (effect size index 'Odds ratio', G*Power's default).",
                  not_sens),
      ig_param("p2", role = "nuisance", domain = ig_domain_real(0, 1, FALSE, FALSE),
               label = list(en = "Proportion p2", de = "Anteil p2"),
               description = "Event probability of group 2 (the reference group); H0: p1 = p2.",
               aliases = "proportion_p2"),
      .ig_unc_common_params()
    ),
    options = .ig_unc_common_options(
      statistics = c("z_pooled", "z_pooled_cc", "z_unpooled", "z_unpooled_cc", "mantel_haenszel",
                     "likelihood_ratio", "t_test"),
      index_default = "odds_ratio",
      direction_labels = c(p_h1_le_p_h0 = "p1\u2264p2", p_h1_ge_p_h0 = "p1\u2265p2")),
    modes = c("a_priori", "post_hoc", "sensitivity", "criterion", "compromise"),
    targets = list(sensitivity = "p1"),
    power = .ig_power_exact_proportions_unconditional,
    design = .ig_design_unconditional,
    discrete = TRUE,
    effect_direction = "unknown",
    engine = list(package = "inferguidR",
                  fun = paste("exact enumeration of the two-binomial sample space (stats::dbinom/pbinom;",
                              "threshold search per x2 for statistics monotone in x1); large-sample",
                              "normal approximations (stats::pnorm, noncentral stats::pt)"),
                  version = NA_character_),
    checker = "validation/independent/g9/check_unconditional.R",
    reference = list(
      gpower_path = "Tests > Proportions > Two independent groups: Inequality, unconditional exact",
      gpower_label = "Exact - Proportions: Inequality, two independent groups (unconditional)",
      manual_pages = "none (the manual has no chapter on this test)",
      allocation_rule = paste("a priori: n1 = 1, 2, ..., n2 = floor(r n1 + 0.5) (r = allocation ratio N2/N1);",
                              "the smallest design whose power reaches the target (every earlier design evaluated)"),
      method_notes = "specs/method-notes/exact_proportions_unconditional.md"
    ),
    status = "independently_checked",
    constraints = .ig_unc_constraints,
    a_priori_search = .ig_unc_a_priori_search_inequality
  )
}

.ig_method_exact_proportions_unconditional_offset <- function() {
  pair <- function(stem, index, domain, en, de, what) {
    list(
      ig_param(paste0(stem, "_h1"), role = "effect", domain = domain,
               label = list(en = paste0(en, "|H1"), de = paste0(de, "|H1")),
               description = sprintf("%s under H1 (effect size index '%s').", what, .ig_unc_index_labels[[index]]),
               applies_when = if (stem == "p1") list(any = list(list(mode = "sensitivity"),
                 list(option = "effect_size_index", "in" = index))) else
                 list(option = "effect_size_index", "in" = index),
               # sensitivity mode solves p1|H1 for every index
               modes = if (stem == "p1") NULL else c("a_priori", "post_hoc", "criterion", "compromise")),
      ig_param(paste0(stem, "_h0"), role = "effect", domain = domain,
               label = list(en = paste0(en, "|H0"), de = paste0(de, "|H0")),
               description = sprintf("%s under H0, the null offset (effect size index '%s').", what,
                                     .ig_unc_index_labels[[index]]),
               applies_when = list(option = "effect_size_index", "in" = index),
               modes = c("a_priori", "post_hoc", "sensitivity", "criterion", "compromise"))
    )
  }
  ig_method(
    id = "exact_proportions_unconditional_offset",
    family = "exact",
    label = list(en = "Proportions: Inequality (offset), two independent groups (unconditional)",
                 de = "Anteile: Ungleichheit (mit Offset), zwei unabh\u00e4ngige Gruppen (unbedingt)"),
    params = ig_params(
      pair("p1", "proportion_p1", ig_domain_real(0, 1, FALSE, FALSE), "Prop. p1", "Anteil p1",
           "Event probability p1 of group 1"),
      pair("difference", "difference", ig_domain_real(-1, 1, FALSE, FALSE), "(p1-p2)", "(p1-p2)",
           "Difference p1 - p2"),
      pair("risk_ratio", "risk_ratio", ig_domain_real(0, Inf, FALSE, FALSE), "(p1/p2)", "(p1/p2)",
           "Risk ratio p1 / p2"),
      pair("odds_ratio", "odds_ratio", ig_domain_real(0, Inf, FALSE, FALSE), "Odds ratio", "Odds Ratio",
           "Odds ratio (p1 / (1 - p1)) / (p2 / (1 - p2))"),
      ig_param("p2", role = "nuisance", domain = ig_domain_real(0, 1, FALSE, FALSE),
               label = list(en = "Proportion p2", de = "Anteil p2"),
               description = "Event probability of group 2 (the reference group), the same under H0 and H1.",
               aliases = "proportion_p2"),
      .ig_unc_common_params()
    ),
    options = .ig_unc_common_options(
      statistics = c("z_pooled", "z_pooled_cc", "z_unpooled", "z_unpooled_cc", "t_test", "score_gart_nam",
                     "score_farrington_manning", "score_miettinen_nurminen"),
      index_default = "proportion_p1",
      direction_labels = c(p_h1_le_p_h0 = "p|H1\u2264p|H0", p_h1_ge_p_h0 = "p|H1\u2265p|H0")),
    modes = c("a_priori", "post_hoc", "sensitivity", "criterion", "compromise"),
    targets = list(sensitivity = "p1_h1"),
    power = .ig_power_exact_proportions_unconditional_offset,
    design = .ig_design_unconditional,
    discrete = TRUE,
    effect_direction = "unknown",
    engine = list(package = "inferguidR",
                  fun = paste("exact enumeration of all (n1 + 1)(n2 + 1) tables (stats::dbinom); restricted",
                              "MLEs in closed form (Farrington-Manning cubic, ratio and odds-ratio quadratics);",
                              "large-sample normal approximations (stats::pnorm, noncentral stats::pt)"),
                  version = NA_character_),
    checker = "validation/independent/g9/check_unconditional.R",
    reference = list(
      gpower_path = "Tests > Proportions > Two independent groups: Inequality with offset, unconditional exact",
      gpower_label = "Exact - Proportions: Inequality (offset), two independent groups (unconditional)",
      manual_pages = "none (the manual has no chapter on this test)",
      allocation_rule = paste("a priori: n1 = 1, 2, ..., n2 = floor(r n1 + 0.5) (r = allocation ratio N2/N1);",
                              "the smallest design whose power reaches the target (every earlier design evaluated)"),
      method_notes = "specs/method-notes/exact_proportions_unconditional_offset.md"
    ),
    status = "independently_checked",
    constraints = .ig_unc_constraints_offset,
    a_priori_search = .ig_unc_a_priori_search_offset
  )
}

# Effect resolution -------------------------------------------------------------------------

.ig_unc_index_params <- c(proportion_p1 = "p1", difference = "difference", risk_ratio = "risk_ratio",
                          odds_ratio = "odds_ratio")

# p1 under H1 from a fully specified p: in sensitivity mode (the only mode with
# the option effect_direction) the solver sets p1 (p1_h1); otherwise exactly
# the effect parameter of the chosen index is given.
.ig_unc_effect_p1 <- function(p, suffix = "") {
  p2 <- as.double(p$p2)
  if (!is.null(p$effect_direction)) {
    nm <- paste0("p1", suffix)
    return(list(p1 = as.double(p[[nm]]), name = nm, value = as.double(p[[nm]])))
  }
  idx <- p$effect_size_index
  want <- paste0(.ig_unc_index_params[[idx]], suffix)
  given <- intersect(paste0(.ig_unc_index_params, suffix), names(p))
  if (!identical(given, want)) {
    .ig_abort("invalid_input", sprintf("the effect parameter must be '%s' for effect_size_index = '%s' (got: %s)",
                                       want, idx, if (length(given)) paste(given, collapse = ", ") else "none"))
  }
  v <- as.double(p[[want]])
  list(p1 = .ig_unc_p1_from_index(idx, v, p2), name = want, value = v)
}

# Cross-parameter checks at specification time.
.ig_unc_constraints <- function(p) {
  out <- character()
  idx <- p$effect_size_index
  par <- .ig_unc_index_params[[idx]]
  if (!is.null(p[[par]]) && !is.null(p$p2)) {
    p1 <- .ig_unc_p1_from_index(idx, as.double(p[[par]]), as.double(p$p2))
    if (is.na(p1) || p1 <= 0 || p1 >= 1) {
      out <- c(out, sprintf("the effect size %s = %s with p2 = %s gives p1 = %s outside (0, 1)", par,
                            format(p[[par]]), format(p$p2), format(p1)))
    }
  }
  c(out, .ig_unc_constraints_design(p))
}

.ig_unc_constraints_offset <- function(p) {
  out <- character()
  idx <- p$effect_size_index
  if (!(p$statistic %in% .ig_unc_offset_statistics[[idx]])) {
    out <- c(out, sprintf(paste(
      "statistic '%s' is not offered with effect_size_index '%s'; G*Power offers: %s (it preselects '%s')"),
      p$statistic, idx, paste(.ig_unc_offset_statistics[[idx]], collapse = ", "),
      .ig_unc_offset_statistics[[idx]][1]))
  }
  par <- .ig_unc_index_params[[idx]]
  for (suffix in c("_h1", "_h0")) {
    nm <- paste0(par, suffix)
    if (!is.null(p[[nm]]) && !is.null(p$p2)) {
      p1 <- .ig_unc_p1_from_index(idx, as.double(p[[nm]]), as.double(p$p2))
      if (is.na(p1) || p1 <= 0 || p1 >= 1) {
        out <- c(out, sprintf("the effect size %s = %s with p2 = %s gives p1 = %s outside (0, 1)", nm,
                              format(p[[nm]]), format(p$p2), format(p1)))
      }
    }
  }
  if (!is.null(p$n1) && !is.null(p$n2) && as.double(p$n1) + as.double(p$n2) < 3) {
    out <- c(out, "the offset test needs n1 + n2 >= 3 (G*Power's admissible range of the total sample size)")
  }
  c(out, .ig_unc_constraints_design(p))
}

.ig_unc_constraints_design <- function(p) {
  out <- character()
  # a priori mode is the only mode with an allocation ratio
  if (!is.null(p$allocation_ratio) && !is.null(p$power) && !is.null(p$alpha) && p$power <= p$alpha) {
    out <- c(out, sprintf("the target power (%s) must exceed alpha (%s) in a priori mode", format(p$power),
                          format(p$alpha)))
  }
  if (identical(p$statistic, "t_test") && !is.null(p$n1) && !is.null(p$n2) &&
      as.double(p$n1) + as.double(p$n2) < 3) {
    out <- c(out, "the t test needs n1 + n2 >= 3 (df = n1 + n2 - 2 >= 1)")
  }
  out
}

# Allocation rule (a priori). G*Power searches n1 = 1, 2, ... with
# n2 = floor(r n1 + 0.5); every N = n1 + n2 is admitted at most once because
# N increases strictly with n1. Designs with n2 < 1, and designs with
# n1 + n2 < 3 for the t test, are not admissible.
.ig_design_unconditional <- function(p, n_total, continuous = FALSE) {
  r <- as.double(p$allocation_ratio)
  N <- as.double(n_total)
  if (continuous) return(c(n1 = N / (1 + r), n2 = r * N / (1 + r)))
  n2_of <- function(n1) floor(r * n1 + 0.5)
  n1 <- max(1, floor(N / (1 + r)))
  while (n1 > 1 && n1 + n2_of(n1) > N) n1 <- n1 - 1
  while (n1 + n2_of(n1) < N) n1 <- n1 + 1
  n2 <- n2_of(n1)
  if (n1 + n2 != N || n2 < 1) return(NULL)
  if (identical(p$statistic, "t_test") && N < 3) return(NULL)
  c(n1 = n1, n2 = n2)
}

# Forward models ----------------------------------------------------------------------------

.ig_power_exact_proportions_unconditional <- function(p) {
  p2 <- as.double(p$p2)
  e <- .ig_unc_effect_p1(p)
  .ig_unc_forward(p, p1 = e$p1, p10 = p2, p2 = p2, h0 = list(scale = "difference", value = 0),
                  offset = FALSE, effect = e)
}

.ig_power_exact_proportions_unconditional_offset <- function(p) {
  p2 <- as.double(p$p2)
  idx <- p$effect_size_index
  if (!(p$statistic %in% .ig_unc_offset_statistics[[idx]])) {
    .ig_abort("invalid_input", sprintf("statistic '%s' is not offered with effect_size_index '%s'",
                                       p$statistic, idx))
  }
  e <- .ig_unc_effect_p1(p, "_h1")
  v0 <- as.double(p[[paste0(.ig_unc_index_params[[idx]], "_h0")]])
  p10 <- .ig_unc_p1_from_index(idx, v0, p2)
  if (is.na(p10) || p10 <= 0 || p10 >= 1) {
    .ig_abort("invalid_input", sprintf("the null value %s gives p1|H0 = %s outside (0, 1)", format(v0),
                                       format(p10)))
  }
  h0 <- switch(idx,
               proportion_p1 = ,
               difference = list(scale = "difference", value = p10 - p2),
               risk_ratio = list(scale = "ratio", value = v0),
               odds_ratio = list(scale = "odds_ratio", value = v0))
  .ig_unc_forward(p, p1 = e$p1, p10 = p10, p2 = p2, h0 = h0, offset = TRUE, effect = e, value0 = v0)
}

# Shared forward model. p1, p2: H1 proportions; (p10, p2): the null point at
# which the actual alpha is evaluated; h0: null hypothesis of the statistic.
.ig_unc_forward <- function(p, p1, p10, p2, h0, offset, effect, value0 = NULL) {
  code <- p$statistic
  tails <- p$tails
  n1 <- as.double(p$n1)
  n2 <- as.double(p$n2)
  N <- n1 + n2
  alpha <- as.double(p$alpha)
  effect_name <- effect$name
  value <- effect$value
  undefined <- function(reason) {
    list(power = NA_real_, beta = NA_real_, alpha_attained = NA_real_,
         notes = list(power = reason, beta = reason, alpha_attained = reason))
  }
  if (is.na(p1) || p1 <= 0 || p1 >= 1) {
    return(undefined(sprintf("%s = %s gives p1 = %s outside (0, 1)", effect_name, format(value), format(p1))))
  }
  # Rejection side of a one-sided test: the effect direction (sensitivity) or
  # the sign of the effect (p1 = p1|H0: upper tail).
  if (!is.null(p$effect_direction)) {
    dir <- if (identical(p$effect_direction, "p_h1_ge_p_h0")) 1 else -1
    if (dir * (p1 - p10) < 0) {
      return(undefined(sprintf(paste(
        "%s = %s lies on the other side of the null than the effect direction '%s' of the sensitivity",
        "search; the search is restricted to that side"), effect_name, format(value), p$effect_direction)))
    }
  } else {
    dir <- if (p1 >= p10) 1 else -1
  }
  side <- if (identical(tails, "two")) "both" else if (dir > 0) "upper" else "lower"
  # G*Power uses the exact distribution for N <= the threshold (capture
  # SUP-G9-UNC-POS-A05: N = 80 with threshold 80 is computed exactly).
  exact <- identical(p$computation, "exact") && N <= as.double(p$exact_n_limit)
  if (exact && !is.null(p$error_ratio)) {
    # Compromise mode (the only mode with an error ratio).
    .ig_abort("reference_pending", paste(
      "compromise mode with the exact computation is reference_pending: the target is beta = q * actual alpha",
      "(the attained size), while the generic compromise solver uses nominal alpha. Deterministic rejection",
      "regions give discrete attained sizes, so the requested equality may have no solution. A qualified",
      "attained-size solver and an explicit rule for unattainable ratios are still required."))
  }
  crit <- .ig_unc_critical(code, alpha, tails, N)
  if (exact) {
    ex <- .ig_unc_exact_probs(code, n1, n2, h0, crit, side, list(c(p1, p2), c(p10, p2)))
    power <- ex$reject[1]
    beta <- ex$accept[1]
    alpha_attained <- ex$reject[2]
    region_method <- ex$method
  } else {
    ls <- .ig_unc_large_sample(code, n1, n2, p1, p2, h0, alpha, tails, dir)
    if (is.na(ls[["power"]])) {
      # The approximation has no value at this point. Inside a search (the
      # sensitivity scan over p1, the a priori scan over designs) the point is
      # reported as undefined and skipped; a calculation at a fixed effect and
      # design is declined with a typed status.
      if (!is.null(p$effect_direction) || isTRUE(p$.ig_unc_search)) return(undefined(attr(ls, "reason")))
      .ig_abort("unsupported", attr(ls, "reason"))
    }
    power <- ls[["power"]]
    beta <- ls[["beta"]]
    alpha_attained <- alpha
    region_method <- "large-sample approximation"
  }
  idx1 <- .ig_unc_indices(p1, p2)
  outputs <- list(
    computation_used = if (exact) "exact" else "large_sample",
    statistic = .ig_unc_stat_labels[[code]],
    rejection_side = side,
    p1 = p1, p2 = p2,
    difference = idx1$difference, risk_ratio = idx1$risk_ratio, odds_ratio = idx1$odds_ratio,
    p1_h0 = p10,
    region_method = region_method
  )
  if (offset) {
    idx0 <- .ig_unc_indices(p10, p2)
    outputs$difference_h0 <- idx0$difference
    outputs$risk_ratio_h0 <- idx0$risk_ratio
    outputs$odds_ratio_h0 <- idx0$odds_ratio
  }
  stat_scale <- if (identical(code, "t_test")) "t" else "z"
  dist <- if (exact) {
    list(h0 = list(family = "unconditional_two_binomial", n1 = n1, n2 = n2, p1 = p10, p2 = p2),
         h1 = list(family = "unconditional_two_binomial", n1 = n1, n2 = n2, p1 = p1, p2 = p2),
         statistic = code, statistic_scale = stat_scale, rejection_side = side)
  } else {
    list(h0 = list(family = if (stat_scale == "t") "t" else "normal"),
         h1 = list(family = if (stat_scale == "t") "noncentral_t" else "normal_approximation"),
         statistic = code, statistic_scale = stat_scale, rejection_side = side)
  }
  null_text <- if (offset) {
    sprintf("H0: %s = %s (p1|H0 = %s, p2 = %s)", .ig_unc_index_labels[[p$effect_size_index]],
            format(value0), format(p10), format(p2))
  } else {
    "H0: p1 = p2"
  }
  list(
    power = power,
    beta = beta,
    alpha_attained = alpha_attained,
    critical = if (identical(tails, "two")) c(-crit, crit) else dir * crit,
    n_total = N,
    group_sizes = c(n1 = n1, n2 = n2),
    statistic_distribution = dist,
    outputs = outputs,
    assumptions = c(
      "two independent binomial samples X1 ~ Bin(n1, p1), X2 ~ Bin(n2, p2) with fixed group sizes",
      null_text,
      sprintf(paste("test: %s with the large-sample critical value (%s); %s"),
              .ig_unc_stat_labels[[code]], if (stat_scale == "t") "t quantile, df = N - 2" else "normal quantile",
              if (identical(tails, "two")) "two-sided: |T| >= the alpha/2 quantile" else
                "one-sided in the direction of the effect"),
      if (exact) paste("power and actual alpha computed exactly over all tables; actual alpha is the",
                       "rejection probability at the null point with the given p2 (no supremum over the",
                       "nuisance parameter; it can exceed alpha)") else
        "large-sample approximation (actual alpha = alpha)",
      if (!exact && identical(code, "score_gart_nam")) paste(
        "Gart-Nam large-sample form as G*Power computes it: the null standard error is the effect divided by",
        "the statistic at the expected H1 table; not reliable for effects within about |skewness| / 6 null",
        "standard errors of the null")
    )
  )
}

# A priori search -----------------------------------------------------------------------------
#
# Hook ig_method(a_priori_search =). The designs of the allocation sequence
# (n1 = 1, 2, ..., n2 = floor(r n1 + 0.5)) are evaluated in order and the first
# one whose power reaches the target is returned, so minimality over the
# sequence is proven by enumeration (exact power is not monotone in N). Exact
# evaluations are charged to a fixed budget (.ig_unc_settings$a_priori_budget);
# a search that would exceed it stops with status "unsupported" instead of
# returning an unverified design. In the large-sample part of the sequence
# (N > exact_n_limit, or computation = "large_sample"), where every
# evaluation is cheap, at most large_sample_scan designs are evaluated one by
# one; beyond that the search moves to a window of large_sample_window designs
# below the root of the power with fractional group sizes and says so in its
# diagnostics.

.ig_unc_a_priori_search_inequality <- function(p, n_lo, n_max) {
  .ig_unc_a_priori_search(p, n_lo, n_max, .ig_power_exact_proportions_unconditional, offset = FALSE)
}

.ig_unc_a_priori_search_offset <- function(p, n_lo, n_max) {
  .ig_unc_a_priori_search(p, n_lo, n_max, .ig_power_exact_proportions_unconditional_offset, offset = TRUE)
}

.ig_unc_a_priori_search <- function(p, n_lo, n_max, power_fun, offset) {
  r <- as.double(p$allocation_ratio)
  target <- as.double(p$power)
  n_min <- if (offset || identical(p$statistic, "t_test")) 3 else 2
  n2_of <- function(n1) floor(r * n1 + 0.5)
  at <- function(n1) {
    q <- p
    q$n1 <- n1
    q$n2 <- n2_of(n1)
    # a design at which the approximation has no value is skipped, not an error
    q$.ig_unc_search <- TRUE
    q
  }
  # Effect, null point and null hypothesis of the statistic.
  exact_mode <- identical(p$computation, "exact")
  limit <- if (exact_mode) as.double(p$exact_n_limit) else 0
  idx <- p$effect_size_index
  p2 <- as.double(p$p2)
  e <- .ig_unc_effect_p1(p, if (offset) "_h1" else "")
  if (!offset) {
    p10 <- p2
    h0 <- list(scale = "difference", value = 0)
  } else {
    v0 <- as.double(p[[paste0(.ig_unc_index_params[[idx]], "_h0")]])
    p10 <- .ig_unc_p1_from_index(idx, v0, p2)
    h0 <- switch(idx, proportion_p1 = , difference = list(scale = "difference", value = p10 - p2),
                 risk_ratio = list(scale = "ratio", value = v0), odds_ratio = list(scale = "odds_ratio", value = v0))
  }
  if (is.na(e$p1) || e$p1 <= 0 || e$p1 >= 1 || is.na(p10) || p10 <= 0 || p10 >= 1) {
    return(list(status = "invalid_input", reason = "the effect or null value gives a proportion outside (0, 1)"))
  }
  if (e$p1 == p10) {
    return(list(status = "infeasible", reason = paste(
      "with no effect (p1 equal to its null value) the rejection probability is the actual alpha of each",
      "design, not a power that grows with N; no sample size reaches the target power")))
  }
  dir <- if (e$p1 >= p10) 1 else -1
  budget <- .ig_unc_budget()
  cost <- 0
  n_exact <- 0
  n_large <- 0
  seq_n <- numeric()
  seq_p <- numeric()
  minimality <- "proven: every earlier design of the allocation sequence was evaluated and failed"
  jumped <- NULL
  n1 <- 1
  large_run <- 0
  repeat {
    n2 <- n2_of(n1)
    N <- n1 + n2
    if (N > n_max) {
      return(list(status = "infeasible", reason = sprintf(
        "no design of the allocation sequence with N <= %s reaches power %s", format(n_max, scientific = FALSE),
        format(target)), diagnostics = list(designs_evaluated = n_exact + n_large, cost = cost)))
    }
    if (n2 < 1 || N < n_min || N < n_lo) {
      n1 <- n1 + 1
      next
    }
    exact <- exact_mode && N <= limit
    if (exact) {
      cost <- cost + .ig_unc_exact_cost(p$statistic, n1, n2, h0)
      if (cost > budget) {
        return(list(status = "unsupported", reason = sprintf(paste(
          "the exhaustive exact search would exceed inferguidR's computation budget (%s table evaluations) at",
          "n1 = %s, n2 = %s before reaching power %s; use computation = 'large_sample' or a lower",
          "exact_n_limit, or evaluate chosen designs in post hoc mode"), format(budget), format(n1), format(n2),
          format(target)),
          diagnostics = list(designs_evaluated = n_exact + n_large, last_n_total = N, cost = cost, budget = budget)))
      }
      n_exact <- n_exact + 1
      large_run <- 0
    } else {
      n_large <- n_large + 1
      large_run <- large_run + 1
      if (large_run > .ig_unc_settings$large_sample_scan && is.null(jumped)) {
        # Jump to a window below the root of the power with fractional group
        # sizes (smooth and increasing in N for these approximations).
        frac <- function(Nf) {
          q <- p
          q$n1 <- Nf / (1 + r)
          q$n2 <- r * Nf / (1 + r)
          v <- .ig_unc_large_sample(p$statistic, q$n1, q$n2, e$p1, p2, h0, as.double(p$alpha), p$tails,
                                    dir)[["power"]]
          # an undefined power (Gart-Nam form without a value) reaches no target
          if (is.na(v)) -target else v - target
        }
        hi <- N
        while (frac(hi) < 0 && hi < n_max) hi <- min(n_max, hi * 2)
        if (frac(hi) >= 0) {
          root <- stats::uniroot(frac, c(N, hi), tol = 1e-6)$root
          n1_new <- max(n1, floor(root / (1 + r)) - .ig_unc_settings$large_sample_window)
          jumped <- c(from_n1 = n1, to_n1 = n1_new, continuous_root_n_total = root)
          minimality <- sprintf(paste(
            "proven up to n1 = %s by enumeration; designs from n1 = %s to %s (large-sample approximation) were",
            "skipped and minimality above them is checked on a window of %d designs below the root of the",
            "power with fractional group sizes"), format(n1 - 1), format(n1), format(n1_new - 1),
            .ig_unc_settings$large_sample_window)
          if (n1_new > n1) {
            n1 <- n1_new
            next
          }
        }
      }
    }
    pw <- power_fun(at(n1))$power
    seq_n <- utils::tail(c(seq_n, N), 10)
    seq_p <- utils::tail(c(seq_p, pw), 10)
    if (!is.na(pw) && pw >= target) {
      return(list(status = "ok", n = N, power = pw, diagnostics = c(list(
        search = "allocation sequence n1 = 1, 2, ..., n2 = floor(r n1 + 0.5), evaluated in order",
        minimality = minimality,
        monotonicity = "not assumed (exact power is not monotone in N)",
        designs_evaluated = n_exact + n_large, exact_designs = n_exact, large_sample_designs = n_large,
        exact_cost = cost, exact_budget = budget,
        power_sequence = list(n_total = seq_n, power = seq_p)),
        if (!is.null(jumped)) list(large_sample_jump = as.list(jumped)))))
    }
    n1 <- n1 + 1
  }
}
