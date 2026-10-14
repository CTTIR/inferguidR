# Boschloo engine versus an independent brute-force enumeration.
#
# The reference below shares no numerical routine with the engine (which uses
# stats::dhyper, stats::phyper, stats::dbinom and stats::optimize):
#
# * binomial coefficients from Pascal's rule (exact integers in doubles for
#   n1 + n2 <= 56);
# * ordering statistics as exact fractions a / C(N, t), ranked with an exact
#   comparison of cross products held in base-2^18 limbs, so ties are exact
#   ties (no tolerance);
# * null and alternative probabilities as polynomials evaluated with `^`;
# * the nuisance supremum from a dense grid (2001 points) followed by a
#   golden-section search around every grid local maximum near the top.

rf_choose <- function(n) {
  row <- 1
  for (i in seq_len(n)) row <- c(row, 0) + c(0, row)
  row
}

# Exact product of two non-negative integers below 2^53 as five base-2^18
# limbs, most significant first.
rf_mul <- function(a, b) {
  B <- 2^18
  la <- c(a %/% B^2, (a %/% B) %% B, a %% B)
  lb <- c(b %/% B^2, (b %/% B) %% B, b %% B)
  r <- numeric(5)
  for (i in 1:3) for (j in 1:3) r[i + j - 1] <- r[i + j - 1] + la[i] * lb[j]
  for (k in 5:2) {
    carry <- r[k] %/% B
    r[k] <- r[k] - carry * B
    r[k - 1] <- r[k - 1] + carry
  }
  r
}

# Exact sign of a1 / d1 - a2 / d2.
rf_cmp <- function(a1, d1, a2, d2) {
  u <- rf_mul(a1, d2)
  v <- rf_mul(a2, d1)
  k <- which(u != v)
  if (length(k)) sign(u[k[1]] - v[k[1]]) else 0
}

rf_design <- function(n1, n2) {
  N <- n1 + n2
  x1 <- rep(0:n1, n2 + 1)                      # x1 varies fastest, as in the engine
  x2 <- rep(0:n2, each = n1 + 1)
  c1 <- rf_choose(n1)
  c2 <- rf_choose(n2)
  cN <- rf_choose(N)
  stopifnot(max(cN) < 2^53)
  d <- list(n1 = n1, n2 = n2, N = N, x1 = x1, x2 = x2, t = x1 + x2,
            w = c1[x1 + 1] * c2[x2 + 1], D = cN[x1 + x2 + 1], c1 = c1, c2 = c2)
  # Vandermonde: the weights of a stratum add up to C(N, t).
  stopifnot(all(tapply(d$w, d$t, sum) == cN))
  d
}

rf_numerators <- function(d, side) {
  vapply(seq_along(d$w), function(i) {
    s <- d$t == d$t[i]
    switch(side,
           lower = sum(d$w[s & d$x1 <= d$x1[i]]),
           upper = sum(d$w[s & d$x1 >= d$x1[i]]),
           minlike = sum(d$w[s & d$w <= d$w[i]]))
  }, numeric(1))
}

# Dense ranks of the exact ordering statistic (1 = most extreme, ties share a rank).
rf_ranks <- function(d, side) {
  a <- rf_numerators(d, side)
  o <- order(a / d$D)
  M <- length(o)
  # s[k] = sign(value of k-th - value of (k+1)-th): -1 for a strict increase,
  # 0 for an exact tie; +1 would mean the floating order is wrong.
  s <- vapply(seq_len(M - 1L), function(k) rf_cmp(a[o[k]], d$D[o[k]], a[o[k + 1L]], d$D[o[k + 1L]]),
              numeric(1))
  if (any(s > 0)) stop("floating order disagrees with the exact order")
  r <- integer(M)
  r[o] <- cumsum(c(1L, as.integer(s < 0)))
  list(rank = r, value = a / d$D)
}

rf_golden_max <- function(f, a, b, tol = 1e-13) {
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

# Supremum over pi in [0, 1] of sum_t W[t] pi^t (1 - pi)^(N - t), one row of
# W per set of tables.
rf_sup <- function(W, N, G = 2001L) {
  W <- matrix(W, ncol = N + 1L)
  grid <- seq(0, 1, length.out = G)
  tt <- 0:N
  B <- outer(tt, grid, function(t, q) q^t * (1 - q)^(N - t))
  V <- W %*% B
  vapply(seq_len(nrow(W)), function(r) {
    # Empty and complete sample spaces have exact probabilities 0 and 1.
    # Pascal coefficients give the latter by the binomial theorem.
    if (all(W[r, ] == 0)) return(0)
    if (identical(as.numeric(W[r, ]), rf_choose(N))) return(1)
    v <- V[r, ]
    best <- max(v)
    f <- function(q) sum(W[r, ] * q^tt * (1 - q)^(N - tt))
    loc <- which(v >= best - 1e-3 & v >= c(-Inf, v[-G]) & v >= c(v[-1L], -Inf))
    for (j in loc) best <- max(best, rf_golden_max(f, grid[max(1L, j - 1L)], grid[min(G, j + 1L)]))
    min(1, best)
  }, numeric(1))
}

# Coefficients W_t of a set of tables.
rf_W <- function(d, member) {
  W <- numeric(d$N + 1L)
  if (any(member)) {
    agg <- tapply(d$w[member], d$t[member], sum)
    W[as.integer(names(agg)) + 1L] <- agg
  }
  W
}

# Boschloo p-values of all tables for one side: p = sup of the null
# probability of all tables ranked at most as extreme (ties included).
rf_side_p <- function(d, side) {
  rk <- rf_ranks(d, side)
  R <- max(rk$rank)
  Wc <- matrix(0, R, d$N + 1L)
  for (i in seq_along(d$w)) Wc[rk$rank[i], d$t[i] + 1L] <- Wc[rk$rank[i], d$t[i] + 1L] + d$w[i]
  Wc <- matrix(apply(Wc, 2L, cumsum), nrow = R)
  list(p = rf_sup(Wc, d$N)[rk$rank], rank = rk$rank, value = rk$value)
}

rf_pvalues <- function(d, variant) {
  if (variant == "central") {
    return(pmin(1, 2 * pmin(rf_side_p(d, "lower")$p, rf_side_p(d, "upper")$p)))
  }
  rf_side_p(d, c(less = "lower", greater = "upper", minlike = "minlike")[[variant]])$p
}

rf_power <- function(d, member, p1, p2) {
  pr <- d$c1[d$x1 + 1] * p1^d$x1 * (1 - p1)^(d$n1 - d$x1) *
    d$c2[d$x2 + 1] * p2^d$x2 * (1 - p2)^(d$n2 - d$x2)
  sum(pr[member])
}

# Engine p-values (value and rigorous upper bound) of every table.
bl_all <- function(n1, n2, variant, control = boschloo_control()) {
  .bl_pvalues_all(as.integer(n1), as.integer(n2), .bl_variant(variant), control)
}

variants <- c("less", "greater", "minlike", "central")
designs <- list(c(1, 1), c(1, 2), c(2, 1), c(2, 2), c(1, 6), c(3, 4), c(4, 4), c(5, 7), c(7, 5),
                c(6, 6), c(3, 10), c(10, 3), c(8, 9), c(15, 15))

test_that("the reference reproduces the phase 02 brute-force values", {
  # Values published with 8 decimals in specs/engine-decisions-aux.md (section
  # 2.3); compared to half a unit of the last decimal.
  d <- rf_design(5, 7)
  i <- which(d$x1 == 1 & d$x2 == 6)
  expect_lt(abs(rf_pvalues(d, "less")[i] - 0.02317922), 5e-9)
  expect_lt(abs(rf_pvalues(d, "minlike")[i] - 0.03857422), 5e-9)
  expect_lt(abs(rf_pvalues(d, "central")[i] - 0.04635844), 5e-9)
  d <- rf_design(15, 15)
  i <- which(d$x1 == 3 & d$x2 == 9)
  expect_lt(abs(rf_pvalues(d, "minlike")[i] - 0.03006801), 5e-9)
})

test_that("p-values of every table agree with the enumeration reference", {
  for (nn in designs) {
    d <- rf_design(nn[1], nn[2])
    for (v in variants) {
      eng <- bl_all(nn[1], nn[2], v)
      ref <- rf_pvalues(d, v)
      info <- sprintf("design (%d, %d), variant %s", nn[1], nn[2], v)
      expect_true(isTRUE(eng$converged), info = info)
      expect_lt(max(abs(eng$value - ref)), 1e-9)
      # The reference value (a lower bound of the supremum) never exceeds the
      # engine's rigorous upper bound.
      expect_true(all(ref <= eng$upper + 1e-12), info = info)
      expect_true(all(eng$value <= eng$upper), info = info)
    }
  }
})

test_that("the engine's tie levels equal the exact tie classes", {
  for (nn in designs) {
    d <- rf_design(nn[1], nn[2])
    for (side in c("lower", "upper", "minlike")) {
      ord <- .bl_ordering(as.integer(nn[1]), as.integer(nn[2]), side, boschloo_control())
      rk <- rf_ranks(d, side)$rank
      expect_identical(as.integer(ord$level), rk,
                       info = sprintf("design (%d, %d), side %s", nn[1], nn[2], side))
    }
  }
})

test_that("rejection regions, attained size and power agree with the enumeration", {
  alphas <- list(less = c(0.025, 0.05), greater = c(0.025, 0.05), minlike = c(0.01, 0.05),
                 central = c(0.01, 0.05))
  probs <- rbind(c(0.1, 0.6), c(0.5, 0.5), c(0.7, 0.2), c(0, 1), c(1, 0), c(0, 0), c(0.3, 0.9))
  boundary <- 0L
  for (nn in designs[-c(1, 2, 3)]) {
    d <- rf_design(nn[1], nn[2])
    for (v in variants) {
      ref_p <- rf_pvalues(d, v)
      for (a in alphas[[v]]) {
        info <- sprintf("design (%d, %d), %s, alpha %s", nn[1], nn[2], v, a)
        reg <- boschloo_region(nn[1], nn[2], alpha = a, variant = v)
        expect_identical(reg$status$code, "ok")
        member <- logical(length(d$w))
        member[reg$rejected$x1 + 1 + reg$rejected$x2 * (nn[1] + 1)] <- TRUE
        near <- abs(ref_p - a) < 1e-8
        boundary <- boundary + sum(near)
        expect_identical(member[!near], (ref_p <= a)[!near], info = info)
        size_ref <- rf_sup(rf_W(d, member), d$N)
        expect_lt(abs(reg$size$value - size_ref), 1e-9)
        expect_lte(reg$size$upper, a)
        expect_lte(size_ref, a + 1e-12)
        for (k in seq_len(nrow(probs))) {
          pw <- boschloo_power(nn[1], nn[2], probs[k, 1], probs[k, 2], alpha = a, variant = v)
          expect_lt(abs(pw$power - rf_power(d, member, probs[k, 1], probs[k, 2])), 1e-12)
        }
      }
    }
  }
  expect_identical(boundary, 0L)
})

test_that("the public observed test equals the per-table enumeration", {
  set.seed(20261003)
  for (nn in list(c(5, 7), c(10, 3), c(8, 9))) {
    d <- rf_design(nn[1], nn[2])
    pick <- sample(seq_along(d$w), 8)
    for (v in variants) {
      ref <- rf_pvalues(d, v)
      for (i in pick) {
        res <- boschloo_test(d$x1[i], nn[1], d$x2[i], nn[2], variant = v)
        expect_lt(abs(res$p_value - ref[i]), 1e-9)
        expect_true(ref[i] <= res$p_upper + 1e-12)
      }
    }
  }
})

test_that("the exact integer ordering equals the exact tie classes and replaces ambiguous floating ties", {
  # Direct check of the exact-arithmetic ordering on every design above.
  for (nn in c(designs, list(c(15, 20), c(20, 30)))) {
    d <- rf_design(nn[1], nn[2])
    for (side in c("lower", "upper", "minlike")) {
      ex <- .bl_exact_order(as.integer(nn[1]), as.integer(nn[2]), side,
                            .bl_tables(as.integer(nn[1]), as.integer(nn[2])))
      rank <- integer(length(ex$perm))
      rank[ex$perm] <- cumsum(c(1L, as.integer(!ex$tie)))
      expect_identical(rank, rf_ranks(d, side)$rank,
                       info = sprintf("design (%d, %d), side %s", nn[1], nn[2], side))
    }
  }
  # With a wide ambiguity band the fallback is triggered (n1 + n2 <= 56) and
  # gives the same levels and p-values as the default floating ordering.
  wide <- boschloo_control(ambiguity_band = c(1e-12, 1e-3))
  used <- 0L
  for (nn in list(c(15, 20), c(20, 30))) {
    d <- rf_design(nn[1], nn[2])
    for (side in c("lower", "upper", "minlike")) {
      ord <- .bl_ordering(as.integer(nn[1]), as.integer(nn[2]), side, wide)
      if (ord$arithmetic == "exact_integer") used <- used + 1L
      expect_identical(as.integer(ord$level), rf_ranks(d, side)$rank)
      expect_null(ord$ambiguity)
    }
  }
  for (v in variants) {
    a <- bl_all(15, 20, v, wide)$value
    b <- bl_all(15, 20, v)$value
    expect_lt(max(abs(a - b)), 1e-12)
  }
  expect_gte(used, 4L)
})


test_that("independent polynomial supremum preserves sample-space identities", {
  for (N in c(2L, 10L, 30L)) {
    expect_equal(rf_sup(rep(0, N + 1L), N), 0)
    expect_equal(rf_sup(rf_choose(N), N), 1)
  }
  # pi(1-pi) has its unique maximum 1/4 at pi=1/2.
  expect_equal(rf_sup(c(0, 1, 0), 2L), .25, tolerance = 1e-14)
})
