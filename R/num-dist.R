
# Rejection regions ---------------------------------------------------------------------

# Critical values of a test whose statistic is symmetric about zero under H0
# (t, z). `qupper(p)` is the upper-tail quantile function of the H0
# distribution. One-sided tests reject on `side` ("upper" or "lower"); a
# two-sided test puts alpha/2 into each tail. The two-sided values are
# exactly symmetric (-c, c).
.ig_critical_symmetric <- function(alpha, tails, side = "upper", qupper) {
  if (identical(tails, "two")) {
    c <- qupper(alpha / 2)
    return(c(-c, c))
  }
  c <- qupper(alpha)
  if (identical(side, "lower")) -c else c
}

# Critical values for an H0 distribution that need not be symmetric (for
# example chi-square or F). `qtail(p, lower_tail)` is its quantile function
# for either tail, so both critical values are taken from the tail they
# bound (an upper critical value from qtail(alpha / 2, FALSE), never from
# qtail(1 - alpha / 2, TRUE)). Two-sided: c(lower, upper) with alpha/2 in
# each tail; one-sided: the single value on `side`.
.ig_critical_values <- function(alpha, tails, side = "upper", qtail) {
  if (identical(tails, "two")) {
    return(c(qtail(alpha / 2, TRUE), qtail(alpha / 2, FALSE)))
  }
  if (identical(side, "lower")) qtail(alpha, TRUE) else qtail(alpha, FALSE)
}

# Probabilities of a two-sided region. The test rejects when X <= lo or
# X > hi (lo <= hi); `ptail(q, lower_tail)` is the tail function of the
# distribution of X under the hypothesis of interest (for example
# function(q, lt) .ig_p_t(q, df, ncp, lt)). Returns
#   reject  P(X <= lo) + P(X > hi)
#   accept  P(lo < X <= hi), a difference of two tails taken on the side
#           where both are the smaller ones, so that a small acceptance
#           probability (a small beta) keeps its relative accuracy; it is
#           never formed as 1 - reject;
#   lower, upper  the two rejection tails separately.
# For a continuous distribution the boundary convention (<= vs <) does not
# matter; discrete tests use their own region functions.
.ig_two_sided_prob <- function(ptail, lo, hi) {
  if (is.na(lo) || is.na(hi) || lo > hi) stop("a two-sided region needs lo <= hi", call. = FALSE)
  lower <- ptail(lo, TRUE)
  upper <- ptail(hi, FALSE)
  below_hi <- ptail(hi, TRUE)
  accept <- if (below_hi <= 0.5) {
    below_hi - lower                  # P(X <= hi) - P(X <= lo)
  } else {
    ptail(lo, FALSE) - upper          # P(X > lo) - P(X > hi)
  }
  list(reject = lower + upper, accept = max(accept, 0), lower = lower, upper = upper)
}

# Power of a t test with H0 distribution t(df) and H1 distribution
# t(df, ncp), following the convention of G*Power's t tests (manual 25.5):
#
# * two-sided: critical values -t_c and t_c with alpha/2 in each tail; the
#   power counts both rejection tails (including the opposite one);
# * one-sided: the critical value has the sign of the noncentrality
#   parameter (rejection on the side of the effect); ncp = 0 uses the upper
#   tail, where the power equals alpha on either side.
#
# beta is computed directly as the non-rejection probability (for two tails
# a difference of two tails of the same side, see .ig_two_sided_prob()),
# not as 1 - power.
.ig_power_t_test <- function(df, ncp, alpha, tails) {
  qupper <- function(p) .ig_q_t(p, df, lower_tail = FALSE)
  if (identical(tails, "two")) {
    crit <- .ig_critical_symmetric(alpha, "two", qupper = qupper)
    region <- .ig_two_sided_prob(function(q, lt) .ig_p_t(q, df, ncp, lower_tail = lt),
                                 crit[1], crit[2])
    return(list(power = region$reject, beta = region$accept, critical = crit, side = "both"))
  }
  side <- if (ncp >= 0) "upper" else "lower"
  crit <- .ig_critical_symmetric(alpha, "one", side, qupper = qupper)
  if (side == "upper") {
    power <- .ig_p_t(crit, df, ncp, lower_tail = FALSE)
    beta <- .ig_p_t(crit, df, ncp, lower_tail = TRUE)
  } else {
    power <- .ig_p_t(crit, df, ncp, lower_tail = TRUE)
    beta <- .ig_p_t(crit, df, ncp, lower_tail = FALSE)
  }
  list(power = power, beta = beta, critical = crit, side = side)
}
