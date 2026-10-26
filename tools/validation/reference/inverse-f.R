# Copyright (c) 2026 Raban Heller. MIT licence; see ../../../LICENSE.md.
# Stable coordinates for the independent chi-square/F convolution. Raw-scale
# adaptive integration can fail at otherwise ordinary nodes in the integer
# design scans under different R versions; neither tolerance nor design grid
# is relaxed here.

# (Z + sqrt(lambda))^2 + chi-square(k-1). Set u = sqrt(y)*sin(theta)
# in the original normal/chi-square convolution. The sine coordinate removes
# endpoint behavior, and explicit normal-center panels resolve narrow peaks.
ref_inverse_pchisq <- function(y, k, lambda, lower_tail = TRUE) {
  if (y <= 0) return(if (lower_tail) 0 else 1)
  s <- sqrt(y); mu <- sqrt(lambda)
  if (k == 1) {
    if (lower_tail) return(stats::pnorm(s - mu) - stats::pnorm(-s - mu))
    return(stats::pnorm(s - mu, lower.tail = FALSE) + stats::pnorm(-s - mu))
  }
  if (k < 1) stop("inverse F reference requires numerator degrees of freedom >= 1")
  inner <- function(theta) {
    u <- s * sin(theta)
    stats::dnorm(u - mu) * stats::pchisq((s * cos(theta))^2, k - 1, lower.tail = lower_tail) * s * cos(theta)
  }
  cuts <- sort(unique(asin(pmax(-1, pmin(1, c(-s, mu - 8, mu, mu + 8, s) / s)))))
  value <- sum(vapply(seq_len(length(cuts) - 1L), function(i) {
    stats::integrate(inner, cuts[i], cuts[i + 1L], rel.tol = 1e-13, abs.tol = 0, subdivisions = 2000L)$value
  }, 0.0))
  if (lower_tail) value else value + stats::pnorm(s - mu, lower.tail = FALSE) + stats::pnorm(-s - mu)
}

# F = (X1/df1)/(X2/df2), integrated over s=log(X2). The exp(s)
# Jacobian is retained in log form. Only the denominator's two 1e-14 tails
# are omitted, bounding the resulting absolute error by 2e-14.
ref_inverse_pf <- function(x, df1, df2, lambda, lower_tail = TRUE) {
  outer <- function(s) {
    v <- exp(s)
    vapply(v, function(vi) ref_inverse_pchisq(x * df1 * vi / df2, df1, lambda, lower_tail), numeric(1)) *
      exp(stats::dchisq(v, df2, log = TRUE) + s)
  }
  lo <- log(stats::qchisq(1e-14, df2))
  hi <- log(stats::qchisq(1e-14, df2, lower.tail = FALSE))
  md <- log(stats::qchisq(0.5, df2))
  stats::integrate(outer, lo, md, rel.tol = 1e-11, abs.tol = 0)$value +
    stats::integrate(outer, md, hi, rel.tol = 1e-11, abs.tol = 0)$value
}
