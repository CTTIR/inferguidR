# Copyright (c) 2026 Raban Heller. MIT licence; see ../../../LICENSE.md.
# Count-based independent table enumeration, including boundary probabilities.
g9_pmf <- function(n, p) {
  x <- 0:n
  choose(n, x) * p^x * (1 - p)^(n - x)
}

g9_ref_stat <- function(stat, x1, x2, n1, n2) {
  N <- n1 + n2
  m <- x1 + x2
  cr <- x1 * n2 - x2 * n1
  fin <- function(num, den) {
    out <- num / den
    out[num == 0] <- 0
    i <- num != 0 & den == 0
    out[i] <- sign(num[i]) * Inf
    out
  }
  pooled <- sqrt(n1 * n2 * m * (N - m) / N)
  unp <- sqrt(n2^2 * x1 * (n1 - x1) / n1 + n1^2 * x2 * (n2 - x2) / n2)
  cc <- sign(cr) * pmax(abs(cr) - N / 2, 0)
  switch(stat,
         z_pooled = sign(cr) * sqrt(ifelse(m == 0 | m == N, 0, N * cr^2 / (n1 * n2 * m * (N - m)))),
         z_pooled_cc = sign(cc) * sqrt(ifelse(m == 0 | m == N, 0, N * cc^2 / (n1 * n2 * m * (N - m)))),
         z_unpooled = fin(cr, unp),
         z_unpooled_cc = fin(cc, unp),
         mantel_haenszel = sign(cr) * sqrt(ifelse(m == 0 | m == N, 0, (N - 1) * cr^2 / (n1 * n2 * m * (N - m)))),
         likelihood_ratio = {
           e <- function(o, ex) ifelse(o == 0, 0, o * log(o / ex))
           sign(cr) * sqrt(pmax(2 * (e(x1, n1 * m / N) + e(n1 - x1, n1 * (N - m) / N) + e(x2, n2 * m / N) +
                                       e(n2 - x2, n2 * (N - m) / N)), 0))
         },
         t_test = fin(cr / (n1 * n2), sqrt((x1 * (n1 - x1) / n1 + x2 * (n2 - x2) / n2) / (N - 2) * (1 / n1 + 1 / n2))))
}

g9_ref_power <- function(stat, n1, n2, p1, p2, alpha, tails) {
  g <- expand.grid(x1 = 0:n1, x2 = 0:n2)
  z <- g9_ref_stat(stat, g$x1, g$x2, n1, n2)
  a <- if (tails == "two") alpha / 2 else alpha
  c <- if (stat == "t_test") stats::qt(a, n1 + n2 - 2, lower.tail = FALSE) else stats::qnorm(a, lower.tail = FALSE)
  dir <- if (p1 >= p2) 1 else -1
  rej <- if (tails == "two") z >= c | z <= -c else if (dir > 0) z >= c else z <= -c
  P1 <- as.vector(outer(g9_pmf(n1, p1), g9_pmf(n2, p2)))
  P0 <- as.vector(outer(g9_pmf(n1, p2), g9_pmf(n2, p2)))
  c(power = sum(P1[rej]), beta = sum(P1[!rej]), alpha = sum(P0[rej]))
}

g9_offset_power <- function(stat, n1, n2, p1_h0, p1_h1, p2, alpha, tails) {
  tab <- expand.grid(x1 = 0:n1, x2 = 0:n2)
  delta <- p1_h0 - p2
  difference <- tab$x1 / n1 - tab$x2 / n2 - delta
  corrected <- sign(difference) * pmax(abs(difference) - (1 / n1 + 1 / n2) / 2, 0)
  m <- tab$x1 + tab$x2; N <- n1 + n2
  pooled_variance <- (m / N) * (1 - m / N) * (1 / n1 + 1 / n2)
  unpooled_variance <- tab$x1 * (n1 - tab$x1) / n1^3 + tab$x2 * (n2 - tab$x2) / n2^3
  t_variance <- (tab$x1 * (n1 - tab$x1) / n1 + tab$x2 * (n2 - tab$x2) / n2) / (N - 2) * (1 / n1 + 1 / n2)
  numerator <- if (grepl("_cc$", stat)) corrected else difference
  variance <- if (stat == "t_test") t_variance else if (grepl("unpooled", stat)) unpooled_variance else pooled_variance
  z <- numerator / sqrt(variance)
  z[numerator == 0] <- 0
  critical <- if (stat == "t_test") qt(if (tails == "two") alpha / 2 else alpha, N - 2, lower.tail = FALSE) else
    qnorm(if (tails == "two") alpha / 2 else alpha, lower.tail = FALSE)
  member <- if (tails == "two") abs(z) >= critical else if (p1_h1 >= p1_h0) z >= critical else z <= -critical
  null <- as.vector(outer(g9_pmf(n1, p1_h0), g9_pmf(n2, p2)))
  alt <- as.vector(outer(g9_pmf(n1, p1_h1), g9_pmf(n2, p2)))
  c(power = sum(alt[member]), beta = sum(alt[!member]), alpha = sum(null[member]))
}
