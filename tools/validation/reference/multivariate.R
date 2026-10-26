# Copyright (c) 2026 Raban Heller. MIT licence; see ../../../LICENSE.md.
# F distribution from independent noncentral chi-square convolution.
ref_pf <- function(x, df1, df2, lambda, lower_tail = TRUE) {
  outer <- function(v) {
    vapply(v, function(vi) ref_pchisq(x * df1 * vi / df2, df1, lambda, lower_tail), numeric(1)) *
      stats::dchisq(v, df2)
  }
  lo <- stats::qchisq(1e-14, df2)
  hi <- stats::qchisq(1e-14, df2, lower.tail = FALSE)
  md <- stats::qchisq(0.5, df2)
  stats::integrate(outer, lo, md, rel.tol = 1e-11, abs.tol = 0)$value +
    stats::integrate(outer, md, hi, rel.tol = 1e-11, abs.tol = 0)$value
}
