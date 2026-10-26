# Copyright (c) 2026 Raban Heller. MIT licence; see ../../../LICENSE.md.
# Selected independent reference primitives from tests/testthat/helper-reference.R.
# Noncentral t uses a normal-CDF integral over the chi-square density;
# noncentral chi-square uses convolution. These are different representations
# from the production algorithms, but share central base-R primitives and
# QUADPACK integration. They do not call inferguidR or noncentral R routines.

# Integrate exp(logf(u)) over the real line, split around the mode of logf.
# Returns the log of the integral.
ref_log_integrate <- function(logf, lower = -Inf, upper = Inf, grid) {
  v <- logf(grid)
  if (!any(is.finite(v))) return(-Inf)
  u0 <- grid[which.max(v)]
  step <- if (length(grid) > 1L) diff(range(grid)) / (length(grid) - 1L) else 1
  opt <- stats::optimize(function(u) {
    w <- logf(u)
    if (is.finite(w)) -w else 1e300
  }, c(max(lower, u0 - 2 * step), min(upper, u0 + 2 * step)), tol = 1e-12)
  um <- if (logf(opt$minimum) > max(v)) opt$minimum else u0
  lm <- logf(um)
  # Width of the peak from the curvature of logf at the mode; the pieces are
  # cut at multiples of it so that a narrow peak is never missed.
  h <- 1e-4 * max(1, abs(um))
  curv <- (logf(um + h) - 2 * lm + logf(um - h)) / h^2
  width <- if (is.finite(curv) && curv < 0) 1 / sqrt(-curv) else 1
  g <- function(u) exp(logf(u) - lm)
  offsets <- width * c(1, 4, 16, 64, 256)
  cuts <- sort(unique(c(lower, um - offsets, um, um + offsets, upper)))
  cuts <- cuts[cuts >= lower & cuts <= upper]
  # The scaled integrand inherits the rounding error of logf, about |lm| eps,
  # so no quadrature can deliver a smaller relative error than that; asking
  # for 1e-12 where |lm| is ~1e6 only produces roundoff errors.
  rel_tol <- max(1e-12, 64 * .Machine$double.eps * abs(lm))
  total <- 0
  for (i in seq_len(length(cuts) - 1L)) {
    if (cuts[i] < cuts[i + 1L]) {
      total <- total + stats::integrate(g, cuts[i], cuts[i + 1L], rel.tol = rel_tol, abs.tol = 0,
                                        subdivisions = 5000L)$value
    }
  }
  log(total) + lm
}

# log(a log a - a - lgamma(a)) without cancellation for large a: by
# Stirling's formula it equals 0.5 log(a) - 0.5 log(2 pi) - e(a), where the
# Stirling error e(a) = lgamma(a) - ((a - 0.5) log a - a + 0.5 log(2 pi)) is
# taken from its asymptotic series for a >= 15.
ref_chisq_const <- function(a) {
  if (a < 15) return(a * log(a) - a - lgamma(a))
  stirling_error <- 1 / (12 * a) - 1 / (360 * a^3) + 1 / (1260 * a^5) - 1 / (1680 * a^7)
  0.5 * log(a) - 0.5 * log(2 * pi) - stirling_error
}

# Noncentral t tail probability, P(T <= t) (lower) or P(T > t) (upper):
#   P(T <= t) = int_0^Inf Phi(t sqrt(v / df) - ncp) f_chi2(v; df) dv.
# The integral is taken over x = log(v / df), where
#   f_chi2(v; df) dv = exp(a (x - expm1(x)) + C(a)) dx,  a = df / 2,
# which stays accurate for very large df and where exp() underflows.
ref_pt <- function(t, df, ncp, lower_tail = TRUE, log_p = FALSE) {
  if (is.infinite(t)) {
    one <- (t > 0) == lower_tail
    return(if (log_p) (if (one) 0 else -Inf) else as.double(one))
  }
  a <- df / 2
  const <- ref_chisq_const(a)
  logf <- function(x) {
    stats::pnorm(t * exp(x / 2) - ncp, lower.tail = lower_tail, log.p = TRUE) +
      a * (x - expm1(x)) + const
  }
  grid <- seq(log(stats::qchisq(1e-12, df) / df) - 30,
              log(stats::qchisq(1e-12, df, lower.tail = FALSE) / df) + 5, length.out = 400L)
  lp <- ref_log_integrate(logf, grid = grid)
  if (log_p) lp else exp(lp)
}

# Upper-tail quantile of the central t distribution from the integral CDF:
# the x with P(T > x) = p. f below is decreasing with f(0) = log(0.5 / p), so
# the root is positive for p < 1/2 and negative for p > 1/2 (a one-sided
# criterion alpha above one half).
ref_qt_upper <- function(p, df) {
  if (p == 0.5) return(0)
  f <- function(x) log(ref_pt(x, df, 0, lower_tail = FALSE)) - log(p)
  if (p < 0.5) {
    lo <- 0
    hi <- 1
    while (f(hi) > 0) hi <- hi * 2
  } else {
    hi <- 0
    lo <- -1
    while (f(lo) < 0) lo <- lo * 2
  }
  stats::uniroot(f, c(lo, hi), tol = 1e-14, maxiter = 1000L)$root
}

# Power of the two-sample t test (pooled variance) from its definition:
# df = n1 + n2 - 2, ncp = d sqrt(n1 n2 / (n1 + n2)); two-sided: alpha/2 per
# tail, both tails counted; one-sided: rejection in the direction of d.
ref_power_t2 <- function(d, n1, n2, alpha, tails) {
  df <- n1 + n2 - 2
  ncp <- d * sqrt(n1 * n2 / (n1 + n2))
  if (tails == "two") {
    tc <- ref_qt_upper(alpha / 2, df)
    power <- ref_pt(tc, df, ncp, lower_tail = FALSE) + ref_pt(-tc, df, ncp, lower_tail = TRUE)
    critical <- c(-tc, tc)
  } else {
    tc <- ref_qt_upper(alpha, df)
    if (d >= 0) {
      power <- ref_pt(tc, df, ncp, lower_tail = FALSE)
      critical <- tc
    } else {
      power <- ref_pt(-tc, df, ncp, lower_tail = TRUE)
      critical <- -tc
    }
  }
  list(power = power, critical = critical, df = df, ncp = ncp)
}

# Noncentral chi-square tail by convolution (k > 1 uses chi2(k - 1); k = 1
# is closed form). Only central pnorm/pchisq/dnorm are used.
ref_pchisq <- function(y, k, lambda, lower_tail = TRUE) {
  if (y <= 0) return(if (lower_tail) 0 else 1)
  s <- sqrt(y)
  mu <- sqrt(lambda)
  if (k == 1) {
    inside <- stats::pnorm(s - mu) - stats::pnorm(-s - mu)
    outside <- stats::pnorm(s - mu, lower.tail = FALSE) + stats::pnorm(-s - mu)
    return(if (lower_tail) inside else outside)
  }
  inner <- function(u) {
    stats::dnorm(u - mu) * stats::pchisq(pmax(y - u^2, 0), k - 1, lower.tail = lower_tail)
  }
  val <- stats::integrate(inner, -s, s, rel.tol = 1e-13, abs.tol = 0, subdivisions = 2000L)$value
  if (lower_tail) {
    val
  } else {
    val + stats::pnorm(s - mu, lower.tail = FALSE) + stats::pnorm(-s - mu)
  }
}
