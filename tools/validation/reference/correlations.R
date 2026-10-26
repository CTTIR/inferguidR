# Copyright (c) 2026 Raban Heller. MIT licence; see ../../../LICENSE.md.
# Orthant quadrature, Pearson variance and Isserlis covariance references.
g7_gl <- function(n = 60) {
  k <- seq_len(n - 1)
  b <- k / sqrt(4 * k^2 - 1)
  J <- matrix(0, n, n)
  J[cbind(k, k + 1)] <- b
  J[cbind(k + 1, k)] <- b
  e <- eigen(J, symmetric = TRUE)
  list(x = e$values, w = 2 * e$vectors[1, ]^2)
}
g7_nodes <- g7_gl(60)

# P(X > a, Y > b) by composite Gauss-Legendre quadrature over x in [a, a + 30].
g7_orthant <- function(a, b, rho) {
  s <- sqrt(1 - rho^2)
  f <- function(x) dnorm(x) * pnorm((b - rho * x) / s, lower.tail = FALSE)
  edges <- seq(max(a, -30), max(a, -30) + 30, length.out = 31)
  sum(vapply(seq_len(30), function(i) {
    lo <- edges[i]; hi <- edges[i + 1]
    (hi - lo) / 2 * sum(g7_nodes$w * f((hi - lo) / 2 * g7_nodes$x + (hi + lo) / 2))
  }, numeric(1)))
}

g7_cells <- function(rho, px, py) {
  h <- qnorm(1 - px); k <- qnorm(1 - py)
  c(a = g7_orthant(-h, -k, rho), b = g7_orthant(h, -k, -rho), c = g7_orthant(-h, k, -rho),
    d = g7_orthant(h, k, rho), h = h, k = k)
}

g7_se_pearson <- function(cl, rho, n) {
  a <- cl[["a"]]; b <- cl[["b"]]; c <- cl[["c"]]; d <- cl[["d"]]; h <- cl[["h"]]; k <- cl[["k"]]
  s <- sqrt(1 - rho^2)
  F1 <- pnorm((h - rho * k) / s) - 0.5
  F2 <- pnorm((k - rho * h) / s) - 0.5
  v <- (a + d) * (b + c) / 4 + (a + c) * (b + d) * F2^2 + (a + b) * (c + d) * F1^2 +
    2 * (a * d - b * c) * F1 * F2 - (a * b - c * d) * F2 - (a * c - b * d) * F1
  sqrt(v / n) / (exp(-(h^2 - 2 * rho * h * k + k^2) / (2 * s^2)) / (2 * pi * s))
}

# Steiger covariance (N - 3) Cov(z_ij, z_kl) from Isserlis' theorem and the delta method.
g7_cov_z <- function(R, ij, kl) {
  i <- ij[1]; j <- ij[2]; k <- kl[1]; l <- kl[2]
  cs <- function(a, b, c, d) R[a, c] * R[b, d] + R[a, d] * R[b, c]   # N Cov(s_ab, s_cd)
  r1 <- R[i, j]; r2 <- R[k, l]
  # r_ij = s_ij - r_ij (s_ii + s_jj)/2 to first order
  terms <- list(list(1, i, j), list(-r1 / 2, i, i), list(-r1 / 2, j, j))
  terms2 <- list(list(1, k, l), list(-r2 / 2, k, k), list(-r2 / 2, l, l))
  v <- 0
  for (t1 in terms) for (t2 in terms2) v <- v + t1[[1]] * t2[[1]] * cs(t1[[2]], t1[[3]], t2[[2]], t2[[3]])
  v / ((1 - r1^2) * (1 - r2^2))
}

g7_zpower <- function(mu, sd1, alpha, tails, side) {
  if (tails == "two") {
    zc <- qnorm(1 - alpha / 2)
    return(pnorm((zc - mu) / sd1, lower.tail = FALSE) + pnorm((-zc - mu) / sd1))
  }
  zc <- qnorm(1 - alpha)
  if (side > 0) pnorm((zc - mu) / sd1, lower.tail = FALSE) else pnorm((-zc - mu) / sd1)
}


# Bonett-Price delta-method SE via Richardson numerical differentiation.
ref_bp <- function(cells, px, py, n) {
  a <- cells[["a"]]; b <- cells[["b"]]; c <- cells[["c"]]; d <- cells[["d"]]
  pm <- min(px, py, 1 - px, 1 - py)
  cc <- (1 - abs((a + b) - (a + c)) / 5 - (0.5 - pm)^2) / 2
  lo <- log(a) + log(d) - log(b) - log(c)
  theta <- function(l) pi / (1 + exp(cc * l))
  rstar <- cos(theta(lo))
  g <- if (rstar >= 0) function(l) -2 * sin(theta(l) / 2)^2 else function(l) 2 * cos(theta(l) / 2)^2
  D <- function(h) (g(lo + h) - g(lo - h)) / (2 * h)
  h <- 1e-3
  deriv <- (4 * D(h / 2) - D(h)) / 3
  sl <- sqrt((1 / a + 1 / b + 1 / c + 1 / d) / n)
  list(r = rstar, se = abs(deriv) * sl, sl = sl, lo = lo, c = cc)
}
