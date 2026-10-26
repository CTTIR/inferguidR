# Copyright (c) 2026 Raban Heller. MIT licence; see ../../../LICENSE.md.
# Independent inversion of the Fleiss-Tytun-Ury relation by bisection.
ref_ftu_effective <- function(n1, n2, D) {
  r <- n2 / n1
  a <- (r + 1) / (2 * r * D)
  if (n1 <= a) return(c(1, max(1, ceiling(r - 1e-12))))
  f <- function(n) (n / 4) * (1 + sqrt(1 + 2 * (r + 1) / (n * r * D)))^2 - n1
  lo <- 1e-12
  hi <- n1
  for (i in 1:200) {
    mid <- (lo + hi) / 2
    if (f(mid) < 0) lo <- mid else hi <- mid
  }
  e1 <- (lo + hi) / 2
  c(max(1, ceiling(e1 - 1e-9 * e1)), max(1, ceiling(r * e1 - 1e-9 * r * e1)))
}

ref_z2p_power <- function(p1, p2, n1, n2, alpha, tails, cc = FALSE, arcsine = FALSE) {
  if (cc) {
    e <- ref_ftu_effective(n1, n2, abs(p2 - p1))
    n1 <- e[1]
    n2 <- e[2]
  }
  z <- ck_zq(if (tails == "two") alpha / 2 else alpha)
  if (arcsine) {
    eff <- abs(2 * asin(sqrt(p2)) - 2 * asin(sqrt(p1)))
    s0 <- sqrt(1 / n1 + 1 / n2)
    s1 <- s0
  } else {
    eff <- abs(p2 - p1)
    pbar <- (n1 * p1 + n2 * p2) / (n1 + n2)
    s0 <- sqrt(pbar * (1 - pbar) * (1 / n1 + 1 / n2))
    s1 <- sqrt(p1 * (1 - p1) / n1 + p2 * (1 - p2) / n2)
  }
  pw <- ck_Phi((eff - z * s0) / s1)
  if (tails == "two") pw <- pw + ck_Phi((-eff - z * s0) / s1)
  list(power = pw, z = z, s0 = s0)
}
