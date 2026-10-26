# Copyright (c) 2026 Raban Heller. MIT licence; see ../../../LICENSE.md.
# Independent checker for the inferguidR Boschloo engine.
#
# Shares no numerical routine with the engine. The engine uses
# stats::dhyper/phyper (ordering), stats::dbinom (Bernstein basis, power),
# a rigorous branch-and-bound bound and stats::optimize. This checker uses
#
# * exact tie classes of the ordering statistics from exact_classes.py
#   (Python integers and fractions, no floating point);
# * binomial coefficients from Pascal's rule (exact for n1 + n2 <= 56,
#   ordinary floating point above);
# * null and alternative probabilities as polynomials evaluated with `^`;
# * a dense uniform grid of the nuisance probability (default 2001 points,
#   both endpoints included) followed by a golden-section search around every
#   grid local maximum within 1e-3 of the grid maximum.
#
# Only base R is used. The public driver sources these functions.

ck_choose <- function(n) {
  row <- 1
  for (i in seq_len(n)) row <- c(row, 0) + c(0, row)
  row
}

ck_design <- function(n1, n2) {
  N <- n1 + n2
  x1 <- rep(0:n1, n2 + 1)                       # x1 varies fastest
  x2 <- rep(0:n2, each = n1 + 1)
  c1 <- ck_choose(n1)
  c2 <- ck_choose(n2)
  list(n1 = n1, n2 = n2, N = N, x1 = x1, x2 = x2, t = x1 + x2,
       w = c1[x1 + 1] * c2[x2 + 1], c1 = c1, c2 = c2)
}

# Exact ranks written by exact_classes.py (x1 varies fastest, as here).
ck_read_ranks <- function(dir, n1, n2) {
  f <- file.path(dir, sprintf("ranks_%d_%d.csv", n1, n2))
  r <- utils::read.csv(f)
  stopifnot(nrow(r) == (n1 + 1) * (n2 + 1),
            identical(r$x1, rep(0:n1, n2 + 1)), identical(r$x2, rep(0:n2, each = n1 + 1)))
  r
}

ck_golden_max <- function(f, a, b, tol = 1e-13) {
  g <- (sqrt(5) - 1) / 2
  c1 <- b - g * (b - a)
  c2 <- a + g * (b - a)
  f1 <- f(c1)
  f2 <- f(c2)
  while (b - a > tol) {
    if (f1 >= f2) {
      b <- c2; c2 <- c1; f2 <- f1; c1 <- b - g * (b - a); f1 <- f(c1)
    } else {
      a <- c1; c1 <- c2; f1 <- f2; c2 <- a + g * (b - a); f2 <- f(c2)
    }
  }
  max(f1, f2, f(a), f(b))
}

ck_basis <- function(N, q) {
  tt <- 0:N
  outer(tt, q, function(t, p) p^t * (1 - p)^(N - t))
}

# Supremum over pi in [0, 1] of sum_t W[t] pi^t (1 - pi)^(N - t); one row of W
# per set of tables. Returns value and maximiser.
ck_sup <- function(W, N, G = 2001L) {
  W <- matrix(W, ncol = N + 1L)
  grid <- seq(0, 1, length.out = G)
  V <- W %*% ck_basis(N, grid)
  tt <- 0:N
  out <- matrix(NA_real_, nrow(W), 2L, dimnames = list(NULL, c("value", "argmax")))
  for (r in seq_len(nrow(W))) {
    v <- V[r, ]
    best <- max(v)
    arg <- grid[which.max(v)]
    f <- function(q) sum(W[r, ] * q^tt * (1 - q)^(N - tt))
    loc <- which(v >= best - 1e-3 & v >= c(-Inf, v[-G]) & v >= c(v[-1L], -Inf))
    for (j in loc) {
      a <- grid[max(1L, j - 1L)]
      b <- grid[min(G, j + 1L)]
      m <- ck_golden_max(f, a, b)
      if (m > best) {
        best <- m
        arg <- (a + b) / 2
      }
    }
    out[r, ] <- c(min(1, best), arg)
  }
  out
}

ck_W <- function(d, member) {
  W <- numeric(d$N + 1L)
  if (any(member)) {
    agg <- tapply(d$w[member], d$t[member], sum)
    W[as.integer(names(agg)) + 1L] <- agg
  }
  W
}

# One-sided or minlike p-values of every table from exact ranks.
ck_side_p <- function(d, rank) {
  R <- max(rank)
  Wc <- matrix(0, R, d$N + 1L)
  for (i in seq_along(d$w)) Wc[rank[i], d$t[i] + 1L] <- Wc[rank[i], d$t[i] + 1L] + d$w[i]
  Wc <- matrix(apply(Wc, 2L, cumsum), nrow = R)
  ck_sup(Wc, d$N)[, "value"][rank]
}

ck_pvalues <- function(d, ranks, variant) {
  switch(variant,
         less = ck_side_p(d, ranks$rank_lower),
         greater = ck_side_p(d, ranks$rank_upper),
         minlike = ck_side_p(d, ranks$rank_minlike),
         central = pmin(1, 2 * pmin(ck_side_p(d, ranks$rank_lower), ck_side_p(d, ranks$rank_upper))))
}

ck_power <- function(d, member, p1, p2) {
  pr <- d$c1[d$x1 + 1] * p1^d$x1 * (1 - p1)^(d$n1 - d$x1) *
    d$c2[d$x2 + 1] * p2^d$x2 * (1 - p2)^(d$n2 - d$x2)
  sum(pr[member])
}

ck_size <- function(d, member) ck_sup(ck_W(d, member), d$N)[1, "value"]
