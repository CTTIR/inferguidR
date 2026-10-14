# Boschloo ordering statistics and tie classes ------------------------------------
#
# Every table (x1, x2) of a design (n1, n2) gets the Fisher conditional
# p-value T(x) of one side given its total t = x1 + x2:
#
# * lower:   T = P(X1 <= x1 | t)  (evidence for p1 < p2 when small)
# * upper:   T = P(X1 >= x1 | t)  (evidence for p1 > p2 when small)
# * minlike: T = sum of P(X1 = y | t) over y with P(X1 = y | t) <= P(X1 = x1 | t),
#            compared with the relative tolerance (as stats::fisher.test).
#
# Conditional probabilities are hypergeometric and computed in log space. The
# complementary value Q = 1 - T is computed directly as well, so that the
# ordering of values near 1 is decided on Q, where they are well separated.
#
# Tables are sorted from most to least extreme and grouped into tie classes
# ("levels"): consecutive sorted values whose relative gap is at most the tie
# tolerance belong to the same level. The tail set of a table is the union of
# all levels up to its own, so tail sets are nested; this is what makes the
# binary search for a rejection region valid.
#
# Gaps between distinct values that fall inside the ambiguity band are
# recorded. For n1 + n2 <= 56 such an ordering is recomputed in exact integer
# arithmetic: binomial coefficients from Pascal's recurrence (exact in
# doubles up to N = 56) and comparisons of fractions by error-free products
# (Dekker's two-product). Larger designs keep the floating ordering and the
# callers check whether the ambiguity touches their result.

# Table index of (x1, x2): x1 varies fastest, i.e. index = x1 + 1 + x2 * (n1 + 1).
.bl_tables <- function(n1, n2) {
  x1 <- rep.int(0:n1, n2 + 1L)
  x2 <- rep(0:n2, each = n1 + 1L)
  list(x1 = x1, x2 = x2, t = x1 + x2)
}

# log(cumsum(exp(l))) for a sequence sorted in ascending order, without
# underflow: blocks whose values span less than 600 nats are summed relative
# to their largest value and combined in log space.
.bl_log_cumsum_asc <- function(l) {
  n <- length(l)
  out <- numeric(n)
  acc <- -Inf
  start <- 1L
  while (start <= n) {
    if (l[start] == -Inf) {
      out[start] <- acc
      start <- start + 1L
      next
    }
    end <- start - 1L + max(which(l[start:n] - l[start] <= 600))
    ref <- l[end]
    loc <- log(cumsum(exp(l[start:end] - ref))) + ref
    if (is.finite(acc)) {
      hi <- pmax(acc, loc)
      loc <- hi + log1p(exp(pmin(acc, loc) - hi))
    }
    out[start:end] <- loc
    acc <- out[end]
    start <- end + 1L
  }
  out
}

# Ordering statistic on the log scale: list(lt = log T, lq = log Q).
.bl_stat_log <- function(n1, n2, side, tab, ld, tol) {
  x1 <- tab$x1
  t <- tab$t
  if (side == "lower") {
    return(list(lt = stats::phyper(x1, n1, n2, t, log.p = TRUE),
                lq = stats::phyper(x1, n1, n2, t, lower.tail = FALSE, log.p = TRUE)))
  }
  if (side == "upper") {
    return(list(lt = stats::phyper(x1 - 1L, n1, n2, t, lower.tail = FALSE, log.p = TRUE),
                lq = stats::phyper(x1 - 1L, n1, n2, t, log.p = TRUE)))
  }
  # minlike: within each stratum, sum the probabilities not larger than the
  # table's own (relative tolerance on the probability, in log space).
  lt <- numeric(length(t))
  lq <- numeric(length(t))
  ltol <- log1p(tol)
  for (ix in split(seq_along(t), t)) {
    if (length(ix) == 1L) {
      lt[ix] <- 0
      lq[ix] <- -Inf
      next
    }
    o <- order(ld[ix])
    ls <- ld[ix][o]
    k <- findInterval(ls + ltol, ls)            # number of values <= own + tolerance
    lcs <- .bl_log_cumsum_asc(ls)
    m <- ls[length(ls)]
    rc <- rev(cumsum(rev(exp(ls - m))))          # rc[j] = sum_{i >= j} exp(ls[i] - m)
    lqs <- rep(-Inf, length(ls))
    rest <- k < length(ls)
    lqs[rest] <- log(rc[k[rest] + 1L]) + m
    lt[ix[o]] <- lcs[k]
    lq[ix[o]] <- lqs
  }
  list(lt = lt, lq = lq)
}

# Build the ordering of all tables of a design for one side.
.bl_ordering <- function(n1, n2, side, control) {
  key <- paste("ordering", n1, n2, side, .bl_control_key(control), sep = "|")
  cached <- .bl_cache_get(key, control)
  if (!is.null(cached)) return(cached)

  n1 <- as.integer(n1)
  n2 <- as.integer(n2)
  N <- n1 + n2
  tab <- .bl_tables(n1, n2)
  ld <- stats::dhyper(tab$x1, n1, n2, tab$t, log = TRUE)
  st <- .bl_stat_log(n1, n2, side, tab, ld, control$tie_tolerance)
  lt <- pmin(st$lt, 0)
  lq <- pmin(st$lq, 0)
  low <- lt <= log(0.5)

  # Sort: low values by T ascending, then high values by Q descending.
  idx_low <- which(low)
  idx_high <- which(!low)
  perm <- c(idx_low[order(lt[idx_low], method = "radix")],
            idx_high[order(-lq[idx_high], method = "radix")])
  gaps <- .bl_sorted_gaps(lt[perm], lq[perm], low[perm])
  band <- control$ambiguity_band
  ambiguous <- which(gaps > band[1] & gaps <= band[2])
  arithmetic <- "floating_tolerance"
  fallback_reason <- NULL
  if (length(ambiguous) && N <= .bl_exact_max_n) {
    ex <- .bl_exact_order(n1, n2, side, tab)
    perm <- ex$perm
    tie <- ex$tie
    arithmetic <- "exact_integer"
    fallback_reason <- sprintf(
      "%d gap(s) between ordering values fell in the ambiguity band; ordering recomputed in exact integer arithmetic",
      length(ambiguous))
  } else {
    tie <- gaps <= control$tie_tolerance
  }
  # Levels: a new level starts after every gap that is not a tie.
  M <- length(perm)
  new_level <- c(TRUE, !tie)
  level_sorted <- cumsum(new_level)
  level <- integer(M)
  level[perm] <- level_sorted
  n_levels <- level_sorted[M]
  level_end <- c(which(!tie), M)
  level_start <- c(1L, level_end[-n_levels] + 1L)

  amb <- NULL
  if (length(ambiguous) && arithmetic == "floating_tolerance") {
    amb <- list(position = ambiguous, gap = gaps[ambiguous], merged = tie[ambiguous],
                level_low = level_sorted[ambiguous], level_high = level_sorted[ambiguous + 1L],
                value_low = exp(lt[perm][ambiguous]))
  }
  stratum_size <- tabulate(tab$t + 1L, N + 1L)
  ord <- list(n1 = n1, n2 = n2, N = N, side = side,
              x1 = tab$x1, x2 = tab$x2, t = tab$t, d = exp(ld),
              lt = lt, lq = lq, perm = perm, level = level, n_levels = n_levels,
              level_start = level_start, level_end = level_end,
              stratum_size = stratum_size, arithmetic = arithmetic,
              fallback_reason = fallback_reason,
              ambiguous_gaps = length(ambiguous), ambiguity = amb,
              underflow = sum(!is.finite(lt)))
  class(ord) <- "bl_ordering"
  .bl_cache_set(key, ord, control)
}

# Relative gaps between consecutive sorted ordering values (on their accurate
# scale). Gap i lies between sorted positions i and i + 1.
.bl_sorted_gaps <- function(lt, lq, low) {
  M <- length(lt)
  if (M < 2L) return(numeric())
  a <- seq_len(M - 1L)
  b <- a + 1L
  g <- numeric(M - 1L)
  both_low <- low[a] & low[b]
  both_high <- !low[a] & !low[b]
  cross <- !both_low & !both_high
  # Low scale: (T_b - T_a) / T_a.
  g[both_low] <- expm1(lt[b][both_low] - lt[a][both_low])
  # High scale: (Q_a - Q_b) / Q_a.
  g[both_high] <- -expm1(lq[b][both_high] - lq[a][both_high])
  g[cross] <- expm1(lt[b][cross] - lt[a][cross])
  g[is.nan(g)] <- 0                             # both values 0 on their scale: a tie
  pmax(g, 0)
}

# Exact integer ordering (n1 + n2 <= 56) ----------------------------------------------

# Row n of Pascal's triangle by repeated addition (exact in doubles for n <= 56;
# base::choose() is not exact there).
.bl_pascal <- function(n) {
  row <- 1
  for (i in seq_len(n)) row <- c(row, 0) + c(0, row)
  row
}

# Error-free product a * b = hi + lo for doubles (Dekker), valid for integers
# below 2^53 under round-to-nearest arithmetic.
.bl_two_prod <- function(a, b) {
  split <- function(x) {
    cc <- 134217729 * x                         # 2^27 + 1
    hi <- cc - (cc - x)
    list(hi = hi, lo = x - hi)
  }
  p <- a * b
  sa <- split(a)
  sb <- split(b)
  err <- ((sa$hi * sb$hi - p) + sa$hi * sb$lo + sa$lo * sb$hi) + sa$lo * sb$lo
  list(hi = p, lo = err)
}

# Exact sign of a1 / d1 - a2 / d2 for non-negative integers below 2^53.
.bl_cmp_frac <- function(a1, d1, a2, d2) {
  u <- .bl_two_prod(a1, d2)
  v <- .bl_two_prod(a2, d1)
  ifelse(u$hi != v$hi, sign(u$hi - v$hi), sign(u$lo - v$lo))
}

# Exact numerators a(x) and denominators D(t) of the ordering statistic.
.bl_exact_values <- function(n1, n2, side, tab) {
  c1 <- .bl_pascal(n1)
  c2 <- .bl_pascal(n2)
  w <- c1[tab$x1 + 1L] * c2[tab$x2 + 1L]       # each product <= C(N, t) < 2^53
  a <- numeric(length(w))
  D <- numeric(length(w))
  for (ix in split(seq_along(w), tab$t)) {
    ww <- w[ix]
    o <- order(tab$x1[ix])
    den <- sum(ww)
    D[ix] <- den
    if (side == "lower") {
      a[ix[o]] <- cumsum(ww[o])
    } else if (side == "upper") {
      a[ix[o]] <- rev(cumsum(rev(ww[o])))
    } else {
      sw <- sort(ww)
      a[ix] <- cumsum(sw)[findInterval(ww, sw)]  # sum of weights <= own weight (exact)
    }
  }
  list(a = a, D = D)
}

.bl_exact_order <- function(n1, n2, side, tab) {
  ev <- .bl_exact_values(n1, n2, side, tab)
  a <- ev$a
  D <- ev$D
  approx <- a / D
  perm <- order(approx, method = "radix")
  M <- length(perm)
  # Values whose floating ratios differ by more than 1e-9 (relative) are
  # ordered correctly; runs of closer values are sorted with the exact
  # comparison and their exact ties identified.
  rel <- if (M > 1L) diff(approx[perm]) / pmax(approx[perm][-1L], .Machine$double.xmin) else numeric()
  close <- rel < 1e-9
  i <- 1L
  while (i < M) {
    if (!close[i]) {
      i <- i + 1L
      next
    }
    j <- i
    while (j < M && close[j]) j <- j + 1L
    run <- perm[i:j]
    for (u in seq_along(run)[-1L]) {            # insertion sort with the exact comparator
      v <- u
      while (v > 1L && .bl_cmp_frac(a[run[v - 1L]], D[run[v - 1L]], a[run[v]], D[run[v]]) > 0) {
        run[c(v - 1L, v)] <- run[c(v, v - 1L)]
        v <- v - 1L
      }
    }
    perm[i:j] <- run
    i <- j + 1L
  }
  tie <- if (M > 1L) {
    .bl_cmp_frac(a[perm[-M]], D[perm[-M]], a[perm[-1L]], D[perm[-1L]]) == 0
  } else logical()
  list(perm = perm, tie = tie, a = a, D = D)
}

# Tail sets ----------------------------------------------------------------------------

# Conditional probabilities c_k = P_H0(tail | X1 + X2 = k), k = 0..N, of the
# tail consisting of a set of table indices. Strata that are completely in
# the tail get exactly 1.
.bl_tail_c <- function(ord, idx) {
  N <- ord$N
  cvec <- numeric(N + 1L)
  if (!length(idx)) return(cvec)
  tt <- ord$t[idx]
  s <- rowsum(ord$d[idx], tt, reorder = TRUE)
  cvec[as.integer(rownames(s)) + 1L] <- s[, 1L]
  full <- tabulate(tt + 1L, N + 1L) == ord$stratum_size
  cvec[full] <- 1
  pmin(cvec, 1)
}

# Table indices of the tail made of levels 1..j.
.bl_level_tables <- function(ord, j) {
  if (j <= 0L) return(integer())
  ord$perm[seq_len(ord$level_end[j])]
}

# Does an ambiguous gap (floating ordering, N > 56) change a result?
#
# * kind "tail": the tail of a table at level j (levels 1..j). A merged gap
#   inside level j could split it; a separating gap right after level j
#   could merge level j + 1 into it.
# * kind "cut": a rejection region made of levels 1..j (j = 0 for an empty
#   region). Splitting a merged gap inside level j + 1 could add part of it
#   to the region; merging the separating gap after level j would join level
#   j into the first excluded level and remove it from the region. Gaps
#   elsewhere cannot change which tables are rejected.
.bl_ambiguity_touches <- function(ord, j, kind = c("tail", "cut")) {
  kind <- match.arg(kind)
  amb <- ord$ambiguity
  if (is.null(amb)) return(FALSE)
  if (kind == "tail") {
    if (j <= 0L) return(FALSE)
    return(any((amb$merged & amb$level_low == j) | (!amb$merged & amb$level_low == j)))
  }
  any((amb$merged & amb$level_low == j + 1L) | (!amb$merged & amb$level_low == j))
}

# Ordering summary for one table (natural scale).
.bl_ordering_value <- function(ord, i) {
  list(value = exp(ord$lt[i]), complement = exp(ord$lq[i]), level = ord$level[i],
       n_levels = ord$n_levels)
}
