# Copyright (c) 2026 Raban Heller. MIT licence; see ../../../LICENSE.md.
# Alternative SAS/Muller-Peterson parameterisation.
ref_terms <- function(approx, ne, p, q, f2 = NA_real_) {
  s <- min(p, q)
  m <- (abs(p - q) - 1) / 2
  n <- (ne - p - 1) / 2
  df1 <- s * (2 * m + s + 1)                       # = p q
  if (approx == "pillai_v") {
    df2 <- s * (2 * n + s + 1)
    g <- s
    crit <- if (is.finite(f2)) f2 * s / (1 + f2) else NA_real_
  } else if (approx == "hotelling_lawley_t1") {      # Pillai-Samson form (SAS: n <= 0 branch)
    df2 <- 2 * (s * n + 1)
    g <- s
    crit <- f2 * s
  } else if (approx == "wilks_u") {                  # Rao (1951): F = (1 - U^(1/t)) / U^(1/t) * df2 / df1
    t <- if (p^2 + q^2 - 5 > 0) sqrt((p^2 * q^2 - 4) / (p^2 + q^2 - 5)) else 1
    r <- ne - (p - q + 1) / 2
    u <- (p * q - 2) / 4
    df2 <- r * t - 2 * u
    g <- t
    crit <- (1 + f2)^(-t)
  } else if (approx == "hotelling_lawley_t2") {      # McKeon (1974), SAS n > 0 form
    b <- (p + 2 * n) * (q + 2 * n) / (2 * (2 * n + 1) * (n - 1))
    df2 <- 4 + (p * q + 2) / (b - 1)
    g <- (2 + (p * q + 2) / (b - 1)) / (2 * n)       # McKeon's c
    crit <- g * f2
  } else stop("unknown approximation")
  list(df1 = df1, df2 = df2, g = g, s = s, crit = crit)
}

ref_lambda <- function(algorithm, N, terms, f2) {
  if (algorithm == "obrien_shieh_1999") N * terms$g * f2 else terms$df2 * f2
}
