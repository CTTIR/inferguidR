# Copyright (c) 2026 Raban Heller. MIT licence; see ../../../LICENSE.md.
# Independent parent moments and Lehmann power; no engine adapter.
gl_nodes <- function(n) {
  k <- seq_len(n - 1L)
  b <- k / sqrt(4 * k^2 - 1)
  J <- matrix(0, n, n)
  J[cbind(k, k + 1L)] <- b
  J[cbind(k + 1L, k)] <- b
  e <- eigen(J, symmetric = TRUE)
  o <- order(e$values)
  list(x = e$values[o], w = 2 * e$vectors[1L, o]^2)
}
GL <- gl_nodes(30L)

# Composite Gauss-Legendre over consecutive panel edges.
gl_panels <- function(f, edges) {
  a <- edges[-length(edges)]
  h <- diff(edges) / 2
  X <- outer(h, GL$x) + (a + h)
  W <- outer(h, GL$w)
  sum(W * f(X))
}

# Panel edges over [-L, L]: fine panels near the centre and near the points
# `brk`, geometric growth outwards; every point of `brk` is an edge.
ref_edges <- function(brk, L) {
  core <- seq(-14, 14, by = 0.1)
  tail <- 14 * 1.15^(1:60)
  tail <- tail[tail < L]
  e <- sort(unique(c(-L, -rev(tail), core, tail, L, brk)))
  e[e >= -L & e <= L]
}

ref_parent <- function(parent) {
  switch(parent,
    normal = list(f = function(x) exp(-x^2 / 2) / sqrt(2 * pi), F = function(x) stats::pnorm(x), L = 40),
    laplace = {
      b <- 1 / sqrt(2)
      list(f = function(x) exp(-abs(x) / b) / (2 * b),
           F = function(x) ifelse(x < 0, 0.5 * exp(x / b), 1 - 0.5 * exp(-x / b)), L = 520, b = b)
    },
    logistic = {
      s <- sqrt(3) / pi
      list(f = function(x) { z <- exp(-abs(x) / s); z / (s * (1 + z)^2) },
           F = function(x) ifelse(x >= 0, 1 / (1 + exp(-x / s)), exp(x / s) / (1 + exp(x / s))), L = 420, s = s)
    })
}

# Owen's T(h, a) by Gauss-Legendre on [0, a] (smooth integrand).
owen_t <- function(h, a) {
  x <- a / 2 * (GL$x + 1)
  sum(a / 2 * GL$w * exp(-h^2 * (1 + x^2) / 2) / (1 + x^2)) / (2 * pi)
}

# G1(t) = P(X1 - X2 <= t) for two independent parent variables (closed forms).
ref_g1 <- function(t, parent) {
  if (parent == "normal") return(stats::pnorm(t / sqrt(2)))
  if (parent == "laplace") {
    b <- 1 / sqrt(2)
    u <- abs(t)
    up <- 1 - 0.5 * exp(-u / b) * (1 + u / (2 * b))
    return(if (t >= 0) up else 1 - up)
  }
  s <- sqrt(3) / pi
  x <- abs(t) / s
  # CDF of the difference of two standard logistic variables:
  # C(x) = e^x (e^x - 1 - x) / (e^x - 1)^2 = 1 + h(x) + x h'(x) with
  # h(x) = 1 / (e^x - 1); near 0 its series 1/2 + x/6 - x^3/180 + x^5/5040
  # - x^7/151200 (from the Bernoulli series of h; truncation below 3e-16 for
  # x < 0.1) avoids the cancellation of the closed form; otherwise the closed
  # form with q = e^-x.
  cx <- if (x < 0.1) 0.5 + x / 6 - x^3 / 180 + x^5 / 5040 - x^7 / 151200 else {
    q <- exp(-x)
    (-expm1(-x) - x * q) / expm1(-x)^2
  }
  if (t >= 0) cx else 1 - cx
}

# G2(t) = int f(u) F(u + t)^2 du: normal by Owen's T, otherwise by the own
# Gauss-Legendre rule with the kinks of the Laplace parent as panel ends.
ref_g2 <- function(t, parent) {
  if (parent == "normal") {
    h <- t / sqrt(2)
    return(stats::pnorm(h) - 2 * owen_t(h, 1 / sqrt(3)))
  }
  P <- ref_parent(parent)
  gl_panels(function(u) P$f(u) * P$F(u + t)^2, ref_edges(c(0, -t), P$L))
}

# Same quadrature applied to G1 (to show its accuracy against the closed forms).
ref_g1_quad <- function(t, parent) {
  P <- ref_parent(parent)
  gl_panels(function(u) P$f(u) * P$F(u + t), ref_edges(c(0, -t), P$L))
}

ref_moments <- function(kind, parent, delta) {
  P <- ref_parent(parent)
  if (kind == "signed_rank") c(p1 = P$F(delta), p2 = ref_g1(2 * delta, parent), p3 = ref_g2(2 * delta, parent)) else {
    g2 <- ref_g2(delta, parent)
    c(p1 = ref_g1(delta, parent), p2 = g2, p3 = g2)
  }
}

# Shift with a given p1 (bisection on the closed forms).
ref_shift <- function(kind, parent, p1) {
  g <- if (kind == "signed_rank") ref_parent(parent)$F else function(x) ref_g1(x, parent)
  lo <- -60
  hi <- 60
  for (i in 1:200) {
    m <- (lo + hi) / 2
    if (m == lo || m == hi) break
    if (g(m) < p1) lo <- m else hi <- m
  }
  (lo + hi) / 2
}

ref_are_k <- function(parent) {
  if (parent == "min_are") return(108 / 125)
  P <- ref_parent(parent)
  # variance 1 by construction: k = 12 (int f^2)^2
  i2 <- gl_panels(function(u) P$f(u)^2, ref_edges(0, P$L))
  12 * i2^2
}

# Lehmann normal approximation (written from Lehmann 1975 and Shieh, Jan &
# Randles 2007 with the three corrections of specs/engine-decisions.md G3).
ref_lehmann <- function(kind, n, m, alpha, tails, cc) {
  p1 <- m[[1]]; p2 <- m[[2]]; p3 <- m[[3]]
  if (kind == "signed_rank") {
    N <- n
    mu0 <- N * (N + 1) / 4
    var0 <- N * (N + 1) * (2 * N + 1) / 24
    mu1 <- N * p1 + choose(N, 2) * p2
    var1 <- N * p1 * (1 - p1) + choose(N, 2) * (2 * (p1 - p2)^2 + 3 * p2 * (1 - p2)) +
      6 * choose(N, 3) * (p3 - p2^2)
  } else {
    a <- n[1]; b <- n[2]
    mu0 <- a * b / 2
    var0 <- a * b * (a + b + 1) / 12
    mu1 <- a * b * p1
    var1 <- a * b * p1 * (1 - p1) + a * b * (b - 1) * (p2 - p1^2) + a * b * (a - 1) * (p3 - p1^2)
  }
  z <- stats::qnorm(if (tails == "two") alpha / 2 else alpha, lower.tail = FALSE)
  corr <- if (cc) 0.5 else 0
  hi <- mu0 + z * sqrt(var0) + corr
  lo <- mu0 - z * sqrt(var0) - corr
  s1 <- sqrt(var1)
  up <- stats::pnorm((hi - mu1) / s1, lower.tail = FALSE)
  if (tails == "two") {
    low <- stats::pnorm((lo - mu1) / s1)
    power <- up + low
    # beta as P(lo < X <= hi) from the smaller tails
    beta <- if (hi <= mu1) stats::pnorm((hi - mu1) / s1) - stats::pnorm((lo - mu1) / s1) else
      stats::pnorm((lo - mu1) / s1, lower.tail = FALSE) - up
    list(power = power, beta = beta, z = z, T_lower = lo, T_upper = hi)
  } else {
    list(power = up, beta = stats::pnorm((hi - mu1) / s1), z = z, T_lower = hi, T_upper = hi)
  }
}
