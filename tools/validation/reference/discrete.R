# Copyright (c) 2026 Raban Heller. MIT licence; see ../../../LICENSE.md.
# Enumerated binomial, McNemar and Fisher reference regions.
ref_lse <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) return(-Inf)
  m <- max(x)
  m + log(sum(exp(x - m)))
}

ref_ldbinom <- function(k, n, p) lchoose(n, k) + k * log(p) + (n - k) * log1p(-p)

# P(X <= k) and P(X >= k) for X ~ Bin(n, p) (exact 0/1 outside the support).
ref_bcdf <- function(k, n, p) {
  if (k < 0) return(0)
  if (k >= n) return(1)
  min(1, exp(ref_lse(ref_ldbinom(0:k, n, p))))
}
ref_bsf <- function(k, n, p) {
  if (k <= 0) return(1)
  if (k > n) return(0)
  min(1, exp(ref_lse(ref_ldbinom(k:n, n, p))))
}

# Rejection regions of the binomial test (method notes of
# exact_proportion_one_sample_binomial): reject X <= lo or X >= up.
ref_lower_crit <- function(n, p, a) {
  k <- -1
  while (k + 1 <= n && ref_bcdf(k + 1, n, p) <= a) k <- k + 1
  k
}
ref_upper_crit <- function(n, p, a) {
  k <- n + 1
  while (k - 1 >= 0 && ref_bsf(k - 1, n, p) <= a) k <- k - 1
  k
}
ref_region_prob <- function(r, n, p) ref_bcdf(r[1], n, p) + ref_bsf(r[2], n, p)
# P(lo < X < up), summed directly
ref_region_accept <- function(r, n, p) {
  a <- max(r[1] + 1, 0)
  b <- min(r[2] - 1, n)
  if (b < a) return(0)
  min(1, exp(ref_lse(ref_ldbinom(a:b, n, p))))
}

ref_binom_region <- function(n, p0, alpha, tails, balancing, p1) {
  up <- p1 >= p0
  if (tails == "one") return(if (up) c(-1, ref_upper_crit(n, p0, alpha)) else c(ref_lower_crit(n, p0, alpha), n + 1))
  kl <- ref_lower_crit(n, p0, alpha / 2)
  ku <- ref_upper_crit(n, p0, alpha / 2)
  if (balancing == "half_each_side") return(c(kl, ku))
  if (balancing == "minor_tail_half_rest_major") {
    A <- c(kl, ref_upper_crit(n, p0, alpha - ref_bcdf(kl, n, p0)))
    B <- c(ref_lower_crit(n, p0, alpha - ref_bsf(ku, n, p0)), ku)
    return(if (ref_region_accept(A, n, p1) <= ref_region_accept(B, n, p1)) A else B)
  }
  r <- c(kl, ku)
  repeat {
    tot <- ref_region_prob(r, n, p0)
    il <- if (r[1] + 1 < r[2]) exp(ref_ldbinom(r[1] + 1, n, p0)) else Inf
    iu <- if (r[2] - 1 > r[1]) exp(ref_ldbinom(r[2] - 1, n, p0)) else Inf
    st <- min(il, iu)
    if (!is.finite(st) || tot + st > alpha) break
    tie <- is.finite(il) && is.finite(iu) && abs(il - iu) <= 1e-10 * max(il, iu)
    if ((!tie && il < iu) || (tie && !up)) r[1] <- r[1] + 1 else r[2] <- r[2] - 1
  }
  r
}

ref_binom_test <- function(n, p0, p1, alpha, tails, balancing) {
  r <- ref_binom_region(n, p0, alpha, tails, balancing, p1)
  list(region = r, power = min(1, ref_region_prob(r, n, p1)), alpha = min(1, ref_region_prob(r, n, p0)))
}

# Normal distribution without pnorm/qnorm.
ref_pgamma_lower <- function(a, x) {
  if (x <= 0) return(0)
  if (x < a + 1) {
    term <- 1 / a
    s <- term
    k <- 1
    repeat {
      term <- term * x / (a + k)
      s <- s + term
      if (abs(term) < abs(s) * 1e-17 || k > 100000) break
      k <- k + 1
    }
    return(min(1, exp(-x + a * log(x) - lgamma(a)) * s))
  }
  # continued fraction for Q (modified Lentz)
  tiny <- 1e-300
  b <- x + 1 - a
  cc <- 1 / tiny
  d <- 1 / b
  h <- d
  for (i in 1:100000) {
    an <- -i * (i - a)
    b <- b + 2
    d <- an * d + b
    if (abs(d) < tiny) d <- tiny
    cc <- b + an / cc
    if (abs(cc) < tiny) cc <- tiny
    d <- 1 / d
    del <- d * cc
    h <- h * del
    if (abs(del - 1) < 1e-16) break
  }
  1 - exp(-x + a * log(x) - lgamma(a)) * h
}
ref_upper_gamma <- function(a, x) {
  # Q(a, x) directly for large x (continued fraction), avoiding 1 - P cancellation
  if (x < a + 1) return(1 - ref_pgamma_lower(a, x))
  tiny <- 1e-300
  b <- x + 1 - a
  cc <- 1 / tiny
  d <- 1 / b
  h <- d
  for (i in 1:100000) {
    an <- -i * (i - a)
    b <- b + 2
    d <- an * d + b
    if (abs(d) < tiny) d <- tiny
    cc <- b + an / cc
    if (abs(cc) < tiny) cc <- tiny
    d <- 1 / d
    del <- d * cc
    h <- h * del
    if (abs(del - 1) < 1e-16) break
  }
  exp(-x + a * log(x) - lgamma(a)) * h
}
ref_Phi <- function(x) {
  if (!is.finite(x)) return(if (x > 0) 1 else 0)
  if (x >= 0) 1 - 0.5 * ref_upper_gamma(0.5, x * x / 2) else 0.5 * ref_upper_gamma(0.5, x * x / 2)
}
ref_qnorm_upper <- function(p) {
  # z with 1 - Phi(z) = p, by bisection
  lo <- -40
  hi <- 40
  for (i in 1:200) {
    mid <- (lo + hi) / 2
    if (1 - ref_Phi(mid) > p) lo <- mid else hi <- mid
    if (hi - lo < 1e-15 * max(1, abs(mid))) break
  }
  (lo + hi) / 2
}

ref_bisect <- function(f, lo, hi, tol = 1e-13, maxit = 300) {
  flo <- f(lo)
  for (i in seq_len(maxit)) {
    mid <- (lo + hi) / 2
    fm <- f(mid)
    if (sign(fm) == sign(flo)) {
      lo <- mid
      flo <- fm
    } else hi <- mid
    if (abs(hi - lo) <= tol * max(1, abs(mid))) break
  }
  (lo + hi) / 2
}

# ---------------------------------------------------------------------------
# McNemar's test (method notes of exact_proportions_mcnemar)
# ---------------------------------------------------------------------------

# Exact unconditional power and size: sum over the number of discordant pairs
# i ~ Bin(N, pd) of the conditional binomial test (theta0 = 1/2, theta1 =
# OR / (1 + OR)); i = 0 never rejects.
ref_mcnemar_exact <- function(N, or, pd, alpha, tails, balancing) {
  theta <- or / (1 + or)
  pw <- 0
  al <- 0
  for (i in seq_len(N)) {
    w <- exp(ref_ldbinom(i, N, pd))
    if (w < 1e-22) next
    t <- ref_binom_test(i, 0.5, theta, alpha, tails, balancing)
    pw <- pw + w * t$power
    al <- al + w * t$alpha
  }
  list(power = pw, alpha = al)
}

# Fast approximation: one binomial test with N_D = ceiling(N pd) discordant
# pairs (an integer product is that integer).
ref_mcnemar_fast <- function(N, or, pd, alpha, tails, balancing) {
  x <- N * pd
  nd <- max(1, if (abs(x - round(x)) < 1e-9) round(x) else ceiling(x))
  t <- ref_binom_test(nd, 0.5, or / (1 + or), alpha, tails, balancing)
  c(t, list(nd = nd))
}

# ---------------------------------------------------------------------------
# Fisher's exact test: unconditional power of the conditional test (method
# notes of exact_proportions_fisher). Hypergeometric and binomial
# probabilities from lchoose(); all tables enumerated; p-values by direct sums.
# ---------------------------------------------------------------------------

ref_lr_stat <- function(x, m, n1, n2) {
  N <- n1 + n2
  o <- cbind(x, n1 - x, m - x, n2 - m + x)
  e <- c(n1 * m / N, n1 * (N - m) / N, n2 * m / N, n2 * (N - m) / N)
  t <- numeric(length(x))
  for (k in 1:4) t <- t + ifelse(o[, k] > 0, o[, k] * log(o[, k] / e[k]), 0)
  2 * t
}

# Rejection indicator within one margin. `tie` resolves equal increments of
# the "maximize alpha" rule: "far" (the tail away from the effect first, the
# engine's rule) or "effect" (the tail on the side of the effect first); the
# two bracket G*Power's floating-point tie breaking.
ref_fisher_region <- function(x, h, m, n1, n2, alpha, tails, rule, up, tie = "far") {
  K <- length(x)
  lim <- alpha * (1 + 1e-12)
  if (tails == "one") {
    pv <- vapply(seq_len(K), function(j) if (up) sum(h[j:K]) else sum(h[1:j]), 0)
    return(pv <= lim)
  }
  if (rule == "maximize_alpha") {
    # two tails of the Fisher-statistic region (scanned inwards from both ends)
    pf <- vapply(seq_len(K), function(j) sum(h[h <= h[j] * (1 + 1e-7)]), 0) <= lim
    if (all(pf)) return(pf)
    a <- 0
    while (a + 1 <= K && pf[a + 1]) a <- a + 1
    b <- K + 1
    while (b - 1 >= 1 && pf[b - 1]) b <- b - 1
    repeat {
      tot <- (if (a > 0) sum(h[1:a]) else 0) + (if (b <= K) sum(h[b:K]) else 0)
      il <- if (a + 1 < b) h[a + 1] else Inf
      iu <- if (b - 1 > a) h[b - 1] else Inf
      st <- min(il, iu)
      if (!is.finite(st) || tot + st > lim) break
      eq <- is.finite(il) && is.finite(iu) && abs(il - iu) <= 1e-10 * max(il, iu)
      lower <- if (!eq) il < iu else if (tie == "effect") !up else up
      if (lower) a <- a + 1 else b <- b - 1
    }
    return(seq_len(K) <= a | seq_len(K) >= b)
  }
  pv <- switch(rule,
    symmetric_fisher_statistic = vapply(seq_len(K), function(j) sum(h[h <= h[j] * (1 + 1e-7)]), 0),
    symmetric_pearson_chi2_statistic = {
      s <- abs(x - n1 * m / (n1 + n2))
      vapply(seq_len(K), function(j) sum(h[s >= s[j] - 1e-7 * abs(s[j])]), 0)
    },
    symmetric_lr_statistic = {
      s <- ref_lr_stat(x, m, n1, n2)
      vapply(seq_len(K), function(j) sum(h[s >= s[j] - 1e-7 * abs(s[j])]), 0)
    })
  pv <= lim
}

ref_fisher_exact <- function(n1, n2, p1, p2, alpha, tails, rule = "symmetric_fisher_statistic", tie = "far") {
  N <- n1 + n2
  up <- p1 > p2
  l1 <- lchoose(n1, 0:n1) + (0:n1) * log(p1) + (n1 - 0:n1) * log1p(-p1)
  l2 <- lchoose(n2, 0:n2) + (0:n2) * log(p2) + (n2 - 0:n2) * log1p(-p2)
  pw <- 0
  al <- 0
  for (m in 0:N) {
    x <- max(0, m - n2):min(n1, m)
    h <- exp(lchoose(n1, x) + lchoose(n2, m - x) - lchoose(N, m))
    rej <- ref_fisher_region(x, h, m, n1, n2, alpha, tails, rule, up, tie)
    pw <- pw + sum(exp(l1[x + 1] + l2[m - x + 1])[rej])
    al <- al + exp(lchoose(N, m) + m * log(p2) + (N - m) * log1p(-p2)) * sum(h[rej])
  }
  list(power = pw, alpha = al)
}

# Large-sample approximation: continuity-corrected z test, pooled H0 SD,
# unpooled H1 SD; both rejection boundaries move outwards by c.
ref_fisher_large <- function(n1, n2, p1, p2, alpha, tails) {
  D <- abs(p2 - p1)
  pbar <- (n1 * p1 + n2 * p2) / (n1 + n2)
  s0 <- sqrt(pbar * (1 - pbar) * (1 / n1 + 1 / n2))
  s1 <- sqrt(p1 * (1 - p1) / n1 + p2 * (1 - p2) / n2)
  cc <- (1 / n1 + 1 / n2) / 2
  z <- ref_qnorm_upper(if (tails == "two") alpha / 2 else alpha)
  pw <- ref_Phi((D - cc - z * s0) / s1)
  if (tails == "two") pw <- pw + ref_Phi((-D - cc - z * s0) / s1)
  pw
}

# ---------------------------------------------------------------------------
# z test of two proportions (method notes of z_proportions_two_independent)
# ---------------------------------------------------------------------------

# Effective group sizes under the continuity correction: the uncorrected
# sizes n with n' = (n / 4) (1 + sqrt(1 + 2 (r + 1) / (n r D)))^2 (Fleiss,
# Tytun & Ury 1980), found by bisection on that relation (the engine uses the
# closed-form inverse), each rounded up; group 1 at or below the minimum
# (r + 1) / (2 r D) of the relation gives the smallest sizes.
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
  z <- ref_qnorm_upper(if (tails == "two") alpha / 2 else alpha)
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
  pw <- ref_Phi((eff - z * s0) / s1)
  if (tails == "two") pw <- pw + ref_Phi((-eff - z * s0) / s1)
  list(power = pw, z = z, s0 = s0)
}
