# Copyright (c) 2026 Raban Heller. MIT licence; see ../../../LICENSE.md.
# Independent density integrals for correlation and random-model R-squared.
g5r_gl <- local({
  n <- 32L
  k <- seq_len(n - 1L)
  b <- k / sqrt(4 * k^2 - 1)
  J <- matrix(0, n, n)
  J[cbind(k, k + 1L)] <- b
  J[cbind(k + 1L, k)] <- b
  e <- eigen(J, symmetric = TRUE)
  o <- order(e$values)
  list(x = e$values[o], w = 2 * e$vectors[1, o]^2)
})

# log of sum(exp(v)) (finite entries).
g5r_lse <- function(v) {
  v <- v[is.finite(v)]
  if (!length(v)) return(-Inf)
  m <- max(v)
  m + log(sum(exp(v - m)))
}

# log integral of exp(logf(u)) over [lo, hi] with `panels` composite GL panels.
g5r_log_int <- function(logf, lo, hi, panels = 200L) {
  if (!(hi > lo)) return(-Inf)
  edges <- seq(lo, hi, length.out = panels + 1L)
  h <- diff(edges) / 2
  mid <- (edges[-1L] + edges[-length(edges)]) / 2
  u <- as.vector(outer(g5r_gl$x, h) + rep(mid, each = length(g5r_gl$x)))
  w <- as.vector(outer(g5r_gl$w, h))
  g5r_lse(logf(u) + log(w))
}

# Standard normal: own erfc, vectorised. erfc(y), y >= 0: Taylor series of
# erf for y < 2 (81 terms), Lentz continued fraction (400 steps) above.
g5r_log_erfc <- function(y) {
  out <- numeric(length(y))
  s <- y < 2
  if (any(s)) {
    k <- 0:80
    coef <- (-1)^k / (factorial(k) * (2 * k + 1))
    M <- outer(y[s], 2 * k + 1, "^")
    erf <- 2 / sqrt(pi) * as.vector(M %*% coef)
    out[s] <- log1p(-erf)
  }
  if (any(!s)) {
    yl <- y[!s]
    tiny <- 1e-300
    f <- yl
    C <- yl
    D <- rep(0, length(yl))
    for (j in 1:400) {
      a <- j / 2
      D <- yl + a * D
      D[abs(D) < tiny] <- tiny
      C <- yl + a / C
      C[abs(C) < tiny] <- tiny
      D <- 1 / D
      f <- f * C * D
    }
    out[!s] <- -yl^2 - 0.5 * log(pi) - log(f)
  }
  out
}

# log Phi(x), vectorised.
g5r_log_pnorm <- function(x) {
  out <- numeric(length(x))
  out[x == Inf] <- 0
  out[x == -Inf] <- -Inf
  neg <- is.finite(x) & x <= 0
  pos <- is.finite(x) & x > 0
  if (any(neg)) out[neg] <- log(0.5) + g5r_log_erfc(-x[neg] / sqrt(2))
  if (any(pos)) out[pos] <- log1p(-exp(log(0.5) + g5r_log_erfc(x[pos] / sqrt(2))))
  out
}
g5r_pnorm <- function(x, lower = TRUE) exp(g5r_log_pnorm(if (lower) x else -x))

# Root of an increasing function: bracket expansion from x0 (step doubling,
# within [lo, hi]) and the Illinois variant of regula falsi; x tolerance
# 1e-14 relative. A start value only seeds the bracket; the root is the one
# of f itself.
g5r_root <- function(f, x0, step, lo = -Inf, hi = Inf, tol = 1e-14, maxit = 200L) {
  a <- x0
  fa <- f(a)
  if (fa == 0) return(a)
  dir <- if (fa < 0) 1 else -1
  b <- a
  fb <- fa
  for (i in 1:200) {
    b <- min(max(a + dir * step, lo), hi)
    fb <- f(b)
    if (sign(fb) != sign(fa)) break
    if (b == lo || b == hi) stop("g5r_root: no sign change in the range", call. = FALSE)
    a <- b
    fa <- fb
    step <- step * 2
  }
  if (fb == 0) return(b)
  side <- 0
  for (i in seq_len(maxit)) {
    c0 <- (a * fb - b * fa) / (fb - fa)
    if (!is.finite(c0) || c0 <= min(a, b) || c0 >= max(a, b)) c0 <- (a + b) / 2
    fc <- f(c0)
    if (fc == 0 || abs(b - a) <= tol * max(1, abs(c0))) return(c0)
    if (sign(fc) == sign(fb)) {
      b <- c0; fb <- fc
      if (side == -1) fa <- fa / 2
      side <- -1
    } else {
      a <- c0; fa <- fc
      if (side == 1) fb <- fb / 2
      side <- 1
    }
  }
  (a + b) / 2
}

# Bisection (kept for robust brackets with known end points).
g5r_bisect <- function(f, lo, hi, tol = 1e-14, maxit = 400L) {
  flo <- f(lo)
  fhi <- f(hi)
  if (flo > 0 || fhi < 0) stop("g5r_bisect: no sign change", call. = FALSE)
  for (i in seq_len(maxit)) {
    mid <- (lo + hi) / 2
    if (mid == lo || mid == hi || hi - lo <= tol * max(1, abs(mid))) break
    if (f(mid) >= 0) hi <- mid else lo <- mid
  }
  (lo + hi) / 2
}

g5r_qnorm_upper <- function(p) {
  g5r_root(function(x) log(p) - g5r_log_pnorm(-x), 0, 1, -40, 40)
}

# Chi-square log density.
g5r_ldchisq <- function(v, k) (k / 2 - 1) * log(v) - v / 2 - (k / 2) * log(2) - lgamma(k / 2)

# Range of u = log(v) carrying the chi-square(k) mass (log density within 70 of its maximum).
g5r_chisq_logrange <- function(k) {
  g <- function(u) g5r_ldchisq(exp(u), k) + u
  mode <- log(k)
  gm <- g(mode)
  sd <- sqrt(2 / k)
  lo <- mode - sd
  while (g(lo) > gm - 70) lo <- lo - sd
  hi <- mode + sd
  while (g(hi) > gm - 70) hi <- hi + sd
  c(lo, hi)
}

# Noncentral t: P(T <= t) or P(T > t), T = (Z + delta) / sqrt(V / nu). The
# integral over u = log(v) is split where the normal argument t sqrt(v/nu) - delta
# changes sign (the integrand's steep part for large |delta|), with dense
# panels on both sides.
g5r_pt <- function(t, nu, delta, lower = TRUE, log_p = FALSE) {
  rg <- g5r_chisq_logrange(nu)
  logf <- function(u) {
    v <- exp(u)
    arg <- t * sqrt(v / nu) - delta
    g5r_log_pnorm(if (lower) arg else -arg) + g5r_ldchisq(v, nu) + u
  }
  cuts <- rg
  if (t != 0 && delta / t > 0) {
    us <- log(nu * (delta / t)^2)
    if (us > rg[1] && us < rg[2]) cuts <- c(rg[1], us, rg[2])
  }
  parts <- vapply(seq_len(length(cuts) - 1L), function(i) g5r_log_int(logf, cuts[i], cuts[i + 1L], panels = 80L), numeric(1))
  lp <- g5r_lse(parts)
  if (log_p) lp else exp(lp)
}

# Central t upper quantile by bisection on the reference CDF.
g5r_qt_upper <- function(p, nu, start = 0) {
  g5r_root(function(x) log(p) - g5r_pt(x, nu, 0, lower = FALSE, log_p = TRUE), start, 0.01 * max(1, abs(start)), -1e4, 1e4)
}

# t-test power in G*Power's convention (two tails: alpha/2 each tail, both
# tails counted; one tail: the side of delta); beta computed directly.
g5r_t_power <- function(nu, delta, alpha, tails) {
  if (tails == "two") {
    tc <- g5r_qt_upper(alpha / 2, nu)
    up <- g5r_pt(tc, nu, delta, lower = FALSE)
    lo <- g5r_pt(-tc, nu, delta, lower = TRUE)
    below_hi <- g5r_pt(tc, nu, delta, lower = TRUE)
    beta <- if (below_hi <= 0.5) below_hi - lo else g5r_pt(-tc, nu, delta, lower = FALSE) - up
    return(list(power = up + lo, beta = beta, critical = c(-tc, tc)))
  }
  tc <- g5r_qt_upper(alpha, nu)
  if (delta >= 0) list(power = g5r_pt(tc, nu, delta, lower = FALSE), beta = g5r_pt(tc, nu, delta, lower = TRUE), critical = tc) else
    list(power = g5r_pt(-tc, nu, delta, lower = TRUE), beta = g5r_pt(-tc, nu, delta, lower = FALSE), critical = -tc)
}

# Noncentral chi-square by convolution (Appendix A): (Z + sqrt(lambda))^2 +
# chi^2(k - 1); vectorised over y (lambda > 0, k > 1).
g5r_pchisq_nc <- function(y, k, lambda, lower = TRUE) {
  m <- sqrt(lambda)
  s <- sqrt(pmax(y, 0))
  if (abs(k - 1) < 1e-12) {
    lo <- g5r_pnorm(s - m) - g5r_pnorm(-s - m)
    return(if (lower) lo else g5r_pnorm(m - s) + g5r_pnorm(-m - s))
  }
  out <- numeric(length(y))
  for (i in seq_along(y)) {
    a <- max(-s[i], m - 40)
    b <- min(s[i], m + 40)
    inner <- 0
    if (b > a) {
      # u = s sin(theta) removes the square-root behaviour of P(chi^2 <= y - u^2)
      # at u = +-s (y - u^2 = s^2 cos^2(theta), du = s cos(theta) dtheta).
      ta <- asin(max(-1, min(1, a / s[i])))
      tb <- asin(max(-1, min(1, b / s[i])))
      edges <- seq(ta, tb, length.out = 25L)
      h <- diff(edges) / 2
      mid <- (edges[-1L] + edges[-25L]) / 2
      th <- as.vector(outer(g5r_gl$x, h) + rep(mid, each = 32L))
      w <- as.vector(outer(g5r_gl$w, h))
      u <- s[i] * sin(th)
      inner <- sum(w * s[i] * cos(th) * exp(-(u - m)^2 / 2) / sqrt(2 * pi) *
                     stats::pchisq((s[i] * cos(th))^2, k - 1, lower.tail = lower))
    }
    out[i] <- if (lower) inner else inner + g5r_pnorm(m - s[i]) + g5r_pnorm(-m - s[i])
  }
  out
}

# Noncentral F: P(F <= x) or P(F > x), F = (X1/d1)/(X2/d2), X1 ~ chi^2(d1, lambda), X2 ~ chi^2(d2).
g5r_pf <- function(x, d1, d2, lambda, lower = TRUE) {
  rg <- g5r_chisq_logrange(d2)
  edges <- seq(rg[1], rg[2], length.out = 41L)
  h <- diff(edges) / 2
  mid <- (edges[-1L] + edges[-41L]) / 2
  u <- as.vector(outer(g5r_gl$x, h) + rep(mid, each = 32L))
  w <- as.vector(outer(g5r_gl$w, h))
  v <- exp(u)
  wt <- w * exp(g5r_ldchisq(v, d2) + u)
  y <- x * d1 * v / d2
  inner <- if (lambda == 0) stats::pchisq(y, d1, lower.tail = lower) else g5r_pchisq_nc(y, d1, lambda, lower)
  sum(wt * inner)
}

g5r_qf_upper <- function(p, d1, d2, start = 1) {
  lx <- g5r_root(function(z) log(p) - log(g5r_pf(exp(z), d1, d2, 0, lower = FALSE)), log(start), 0.01, -60, 60, tol = 1e-14)
  exp(lx)
}

# ---- exact distribution of r (Hotelling 1953 density) -------------------------------

# log 2F1(1/2, 1/2; c; w) by its own series, vectorised over w in [0, 1).
g5r_log_2f1_half <- function(cc, w) {
  wmax <- max(w)
  K <- if (wmax <= 0) 1L else as.integer(min(3e5, ceiling(log(1e-18) / log(wmax)) + 200))
  if (cc > 200) K <- min(K, 400L)
  k <- 0:(K - 1)
  base <- 2 * (lgamma(k + 0.5) - lgamma(0.5)) - (lgamma(cc + k) - lgamma(cc)) - lgamma(k + 1)
  vapply(w, function(ww) if (ww <= 0) 0 else g5r_lse(base + k * log(ww)), numeric(1))
}

# log density of r at x (sample size n, population rho), vectorised.
g5r_ldr <- function(x, rho, n) {
  ok <- abs(x) < 1
  out <- rep(-Inf, length(x))
  r <- x[ok]
  out[ok] <- log(n - 2) + lgamma(n - 1) - 0.5 * log(2 * pi) - lgamma(n - 0.5) +
    ((n - 1) / 2) * log1p(-rho^2) + ((n - 4) / 2) * log1p(-r^2) - (n - 1.5) * log1p(-rho * r) +
    g5r_log_2f1_half(n - 0.5, (1 + rho * r) / 2)
  out
}

# n = 3: the density of r behaves like (1 - r^2)^(-1/2) at r = +-1, and in
# z = atanh(r) its tails decay only like exp(-|z|), so the z range [-19, 19]
# (beyond which tanh(z) is 1 in double precision) would cut off up to 1e-7 of
# the mass. With r = sin(theta) the factor (1 - r^2)^((n - 4)/2) dr becomes
# cos(theta)^(n - 3) dtheta, a bounded integrand on [-pi/2, pi/2]; the panels
# are graded geometrically towards both ends (the factor
# (1 - rho r)^(-(n - 3/2)) peaks there for |rho| near 1).
g5r_pr_theta <- function(x, rho, n, lower = TRUE, log_p = FALSE) {
  tx <- asin(x)
  a <- if (lower) -pi / 2 else tx
  b <- if (lower) tx else pi / 2
  cuts <- sort(unique(c(-pi / 2 + 2^-(0:45), 0, pi / 2 - 2^-(0:45))))
  edges <- c(a, cuts[cuts > a & cuts < b], b)
  const <- log(n - 2) + lgamma(n - 1) - 0.5 * log(2 * pi) - lgamma(n - 0.5) + ((n - 1) / 2) * log1p(-rho^2)
  logf <- function(th) {
    s <- sin(th)
    const + (if (n == 3) 0 else (n - 3) * log(cos(th))) - (n - 1.5) * log1p(-rho * s) +
      g5r_log_2f1_half(n - 0.5, (1 + rho * s) / 2)
  }
  parts <- vapply(seq_len(length(edges) - 1L), function(i) g5r_log_int(logf, edges[i], edges[i + 1L], panels = 1L), numeric(1))
  lp <- g5r_lse(parts)
  if (log_p) lp else exp(lp)
}

# P(r <= x) or P(r > x): density integrated in z = atanh(r) (n >= 4) or in
# theta = asin(r) (n = 3, see g5r_pr_theta).
g5r_pr <- function(x, rho, n, lower = TRUE, log_p = FALSE) {
  if (x <= -1) return(if (lower) (if (log_p) -Inf else 0) else (if (log_p) 0 else 1))
  if (x >= 1) return(if (lower) (if (log_p) 0 else 1) else (if (log_p) -Inf else 0))
  if (n < 4) return(g5r_pr_theta(x, rho, n, lower, log_p))
  z0 <- atanh(rho)
  sd <- 1 / sqrt(max(n - 3, 0.25))
  zx <- atanh(x)
  lo <- min(z0 - 45 * sd, zx - 1) ; hi <- max(z0 + 45 * sd, zx + 1)
  lo <- max(lo, -19); hi <- min(hi, 19)
  logf <- function(z) g5r_ldr(tanh(z), rho, n) + log1p(-tanh(z)^2)
  lp <- if (lower) g5r_log_int(logf, lo, zx, panels = 40L) else g5r_log_int(logf, zx, hi, panels = 40L)
  if (log_p) lp else exp(lp)
}

g5r_qr <- function(p, rho, n, lower = TRUE, start = NULL) {
  f <- if (lower) function(z) g5r_pr(tanh(z), rho, n, TRUE, TRUE) - log(p) else
    function(z) log(p) - g5r_pr(tanh(z), rho, n, FALSE, TRUE)
  z0 <- if (is.null(start)) atanh(rho) else atanh(start)
  tanh(g5r_root(f, z0, 0.02 / sqrt(max(n - 3, 1)), -18, 18, tol = 1e-14))
}

# ---- exact distribution of R^2 (Fisher 1928 density) ---------------------------------

# log 2F1(a, a; c; z), a = (N-1)/2, by its own log-scale series (0 <= z < 1).
g5r_log_2f1_aa <- function(a, cc, z) {
  if (z == 0) return(0)
  lz <- log(z)
  # terms peak near k* with (a + k)^2 z ~ (c + k)(k + 1)
  kpk <- max(0, (2 * a * z - cc + sqrt((2 * a * z - cc)^2 + 4 * (1 - z) * a^2 * z)) / (2 * (1 - z)))
  half <- ceiling(25 * sqrt(kpk + 1) / (1 - sqrt(z)) + 100)
  k <- seq(max(0, floor(kpk) - half), floor(kpk) + half)
  lt <- 2 * (lgamma(a + k) - lgamma(a)) - (lgamma(cc + k) - lgamma(cc)) - lgamma(k + 1) + k * lz
  g5r_lse(lt)
}

g5r_ldr2 <- function(x, rho2, m, N) {
  a <- (N - 1) / 2
  b <- (N - m - 1) / 2
  vapply(x, function(v) {
    if (v <= 0 || v >= 1) return(-Inf)
    a * log1p(-rho2) + (m / 2 - 1) * log(v) + (b - 1) * log1p(-v) -
      (lgamma(m / 2) + lgamma(b) - lgamma(m / 2 + b)) + g5r_log_2f1_aa(a, m / 2, rho2 * v)
  }, numeric(1))
}

# P(R^2 <= x) or P(R^2 > x): density integrated in logit(x).
g5r_pr2 <- function(x, rho2, m, N, lower = TRUE, log_p = FALSE) {
  if (x <= 0) return(if (lower) (if (log_p) -Inf else 0) else (if (log_p) 0 else 1))
  if (x >= 1) return(if (lower) (if (log_p) 0 else 1) else (if (log_p) -Inf else 0))
  logit <- function(p) log(p) - log1p(-p)
  inv <- function(u) 1 / (1 + exp(-u))
  mu <- rho2 + (1 - rho2) * m / (N - 1)
  sdx <- sqrt(max((4 * rho2 * (1 - rho2)^2 + 2 * m * (1 - rho2)^2 / (N - 1)) / (N - 1), 1e-16))
  mode_u <- logit(min(max(mu, 1e-12), 1 - 1e-12))
  span <- 60 * sdx / (min(max(mu, 1e-12), 1 - 1e-12) * (1 - min(max(mu, 1e-12), 1 - 1e-12)))
  lo <- min(mode_u - span, logit(x) - 2, -60)
  hi <- max(mode_u + span, logit(x) + 2)
  hi <- min(hi, 40)
  logf <- function(u) { v <- inv(u); g5r_ldr2(v, rho2, m, N) + log(v) + log1p(-v) }
  ux <- logit(x)
  lp <- if (lower) g5r_log_int(logf, lo, ux, panels = 60L) else g5r_log_int(logf, ux, hi, panels = 60L)
  if (log_p) lp else exp(lp)
}

g5r_qr2 <- function(p, rho2, m, N, lower = TRUE, start = NULL, cdf = g5r_pr2) {
  f <- if (lower) function(u) log(cdf(1 / (1 + exp(-u)), rho2, m, N, TRUE)) - log(p) else
    function(u) log(p) - log(cdf(1 / (1 + exp(-u)), rho2, m, N, FALSE))
  x0 <- if (is.null(start)) min(max(rho2 + (1 - rho2) * m / (N - 1), 1e-6), 1 - 1e-6) else start
  1 / (1 + exp(-g5r_root(f, log(x0) - log1p(-x0), 0.01, -60, 36, tol = 1e-14)))
}

# Lee (1971) three-moment approximation: parameters derived from the cumulants
# of chi^2_m(phi W), W ~ chi^2(N - 1): kappa_r = 2^(r-1)(r-1)! (m + (N-1)((1+phi)^r - 1)),
# matched by g chi^2_nu(lambda) (independent derivation and code).
g5r_lee <- function(rho2, m, N) {
  phi <- rho2 / (1 - rho2)
  n <- N - 1
  k1 <- m + n * phi
  k2 <- 2 * (m + n * ((1 + phi)^2 - 1))
  k3 <- 8 * (m + n * ((1 + phi)^3 - 1))
  # scaled noncentral chi-square: k1 = g(nu + l), k2 = 2 g^2 (nu + 2l), k3 = 8 g^3 (nu + 3l)
  A1 <- k1; A2 <- k2 / 2; A3 <- k3 / 8
  g <- (A2 - sqrt(max(A2^2 - A1 * A3, 0))) / A1
  l <- max((A2 - g * A1) / g^2, 0)
  list(g = g, lambda = l, nu = A1 / g - l)
}

g5r_pr2_lee <- function(x, rho2, m, N, lower = TRUE) {
  if (rho2 == 0) return(g5r_pr2(x, 0, m, N, lower))
  lp <- g5r_lee(rho2, m, N)
  d2 <- N - m - 1
  g5r_pf(x / (1 - x) * d2 / (lp$g * lp$nu), lp$nu, d2, lp$lambda, lower)
}
