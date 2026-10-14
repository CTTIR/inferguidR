# Nuisance supremum with a rigorous bound -------------------------------------------
#
# For a fixed set A of tables, the null probability at the common nuisance
# probability pi is a Bernstein polynomial of degree N = n1 + n2:
#
#   f(pi) = P_pi(A) = sum_{k=0}^{N} c_k * dbinom(k, N, pi),
#
# where c_k = P_H0(A | X1 + X2 = k) lies in [0, 1] (Vandermonde's identity).
# Its second derivative is f''(pi) = N (N - 1) sum_{k=0}^{N-2} D2c_k
# dbinom(k, N - 2, pi) with D2c_k = c_{k+2} - 2 c_{k+1} + c_k. On an interval
# [a, b], dbinom(k, N - 2, pi) is largest at the point of [a, b] closest to
# its mode k / (N - 2), which gives the local bound
#
#   M(a, b) = N (N - 1) sum_k |D2c_k| max_{pi in [a, b]} dbinom(k, N - 2, pi).
#
# With h = b - a, s = (pi - a) / h and K = M h^2 / 2, f is bounded on [a, b]
# by the concave majorant f(a) + (f(b) - f(a)) s + K s (1 - s), whose maximum
# over s in [0, 1] is the interval bound. The global grid gives the starting
# intervals; intervals whose bound exceeds the best value found by more than
# the tolerance are bisected (branch and bound) until no interval can hide a
# larger value. The best computed value is a lower bound of the supremum and
# the largest remaining interval bound an upper bound, each up to a stated
# floating-point allowance for the evaluation of f.

# Bernstein basis: (N + 1) x length(pis) matrix of dbinom(k, N, pi).
.bl_basis <- function(N, pis) {
  matrix(stats::dbinom(rep.int(0:N, length(pis)), N, rep(pis, each = N + 1L)),
         nrow = N + 1L)
}

# (N - 1) x m matrix of max_{pi in [a_i, b_i]} dbinom(k, N - 2, pi), k = 0..N-2.
.bl_curvature_basis <- function(N, a, b) {
  M <- N - 2L
  if (M <= 0L) return(matrix(1, 1L, length(a)))
  k <- 0:M
  pis <- pmin(pmax(rep(k / M, length(a)), rep(a, each = M + 1L)), rep(b, each = M + 1L))
  matrix(stats::dbinom(rep.int(k, length(a)), M, pis), nrow = M + 1L)
}

# Grid matrices depend only on (N, grid size). They are large, so only the
# few most recent ones are kept (separately from the result cache).
.bl_grid_cache <- new.env(parent = emptyenv())
.bl_grid_cache_max <- 6L

.bl_grid <- function(N, G, control) {
  key <- paste(N, G, sep = "|")
  if (isTRUE(control$cache) && !is.null(.bl_grid_cache[[key]])) return(.bl_grid_cache[[key]])
  grid <- seq(0, 1, length.out = G)
  val <- list(grid = grid, B = .bl_basis(N, grid),
              C = .bl_curvature_basis(N, grid[-G], grid[-1L]))
  if (isTRUE(control$cache)) {
    keys <- ls(.bl_grid_cache, all.names = TRUE)
    if (length(keys) >= .bl_grid_cache_max) rm(list = keys, envir = .bl_grid_cache)
    assign(key, val, envir = .bl_grid_cache)
  }
  val
}

# Maximum of the concave quadratic majorant on each interval.
.bl_interval_bound <- function(fa, fb, m, h) {
  K <- m * h^2 / 2
  D <- fb - fa
  s <- ifelse(K > 0, pmin(1, pmax(0, (D + K) / (2 * K))), as.numeric(D > 0))
  fa + D * s + K * s * (1 - s)
}

# Floating-point allowance for evaluating f (a sum of N + 1 non-negative terms
# bounded by 1, each from dbinom with relative error far below 1e-13).
.bl_fp_allowance <- function(N) 1e-13 + 4 * (N + 1) * .Machine$double.eps

# Supremum over pi in [0, 1] of sum_k c_k dbinom(k, N, pi).
#
# Returns `value` (best computed value; a lower bound), `upper` (rigorous
# upper bound), `argmax`, convergence flags and work counters. With `alpha`,
# an interval [value, upper] that straddles alpha is refined further with
# tighter tolerances until the comparison with alpha is decided or the
# floating-point allowance is reached (`resolved`).
.bl_sup <- function(cvec, control, alpha = NULL) {
  out <- .bl_sup_once(cvec, control, control$nuisance_tol)
  if (!is.null(alpha)) {
    tol <- control$nuisance_tol
    # The bisection can only close the gap down to the evaluation allowance.
    floor_tol <- 2 * .bl_fp_allowance(length(cvec) - 1L)
    tries <- 0L
    while (out$value <= alpha && out$upper > alpha && tries < 8L && tol > floor_tol) {
      tol <- max(floor_tol, min(tol / 16, max((alpha - out$value) / 2, (out$upper - alpha) / 4)))
      out <- .bl_sup_once(cvec, control, tol)
      tries <- tries + 1L
    }
    out$resolved <- out$upper <= alpha || out$value > alpha
  }
  out
}

.bl_sup_once <- function(cvec, control, tol) {
  N <- length(cvec) - 1L
  fp <- .bl_fp_allowance(N)
  tol <- max(tol, 2 * fp)                       # the bound cannot close below the allowance
  base <- list(grid_points = control$grid_points, tolerance = tol, fp_allowance = fp)
  if (all(cvec == 0)) {
    return(c(list(value = 0, upper = 0, argmax = NA_real_, converged = TRUE, levels = 0L,
                  evaluations = 0L, curvature_bound = 0, method = "empty set"), base))
  }
  if (all(cvec == 1)) {
    return(c(list(value = 1, upper = 1, argmax = 0.5, converged = TRUE, levels = 0L,
                  evaluations = 0L, curvature_bound = 0, method = "all tables"), base))
  }
  d2 <- if (N >= 2L) abs(diff(cvec, differences = 2L)) else numeric()
  f_at <- function(pis) as.vector(crossprod(.bl_basis(N, pis), cvec))
  if (!length(d2) || max(d2) == 0) {
    # f is linear in pi: the supremum is attained at an endpoint.
    v <- max(cvec[1L], cvec[N + 1L])
    return(c(list(value = v, upper = min(1, v + fp),
                  argmax = if (cvec[1L] >= cvec[N + 1L]) 0 else 1, converged = TRUE,
                  levels = 0L, evaluations = 2L, curvature_bound = 0, method = "linear"), base))
  }
  G <- control$grid_points
  gr <- .bl_grid(N, G, control)
  fg <- as.vector(crossprod(gr$B, cvec))
  best <- max(fg)
  arg <- gr$grid[which.max(fg)]
  a <- gr$grid[-G]
  b <- gr$grid[-1L]
  fa <- fg[-G]
  fb <- fg[-1L]
  scale <- N * (N - 1)
  ub <- .bl_interval_bound(fa, fb, scale * as.vector(crossprod(gr$C, d2)), b - a) + fp
  u_discarded <- -Inf
  level <- 0L
  evaluations <- G
  repeat {
    keep <- ub > best + tol
    if (any(!keep)) u_discarded <- max(u_discarded, ub[!keep])
    a <- a[keep]; b <- b[keep]; fa <- fa[keep]; fb <- fb[keep]; ub <- ub[keep]
    # Stop at the level limit, or if the active set explodes (a plateau of
    # width far above the tolerance); the result is then flagged unconverged.
    if (!length(a) || level >= control$max_levels || length(a) > 2e5) break
    level <- level + 1L
    m <- (a + b) / 2
    fm <- f_at(m)
    evaluations <- evaluations + length(m)
    if (max(fm) > best) {
      best <- max(fm)
      arg <- m[which.max(fm)]
    }
    a2 <- c(a, m)
    b2 <- c(m, b)
    fa2 <- c(fa, fm)
    fb2 <- c(fm, fb)
    ub <- .bl_interval_bound(fa2, fb2,
                             scale * as.vector(crossprod(.bl_curvature_basis(N, a2, b2), d2)),
                             b2 - a2) + fp
    a <- a2; b <- b2; fa <- fa2; fb <- fb2
  }
  converged <- !length(a)
  upper <- max(best, u_discarded, if (length(ub)) max(ub) else -Inf)
  # Brent polish of the maximiser (improves the reported value only).
  h <- 2 / ((G - 1) * 2^level)
  lo <- max(0, arg - h)
  hi <- min(1, arg + h)
  if (hi > lo) {
    opt <- stats::optimize(function(p) sum(cvec * stats::dbinom(0:N, N, p)), c(lo, hi),
                           maximum = TRUE, tol = control$polish_tol)
    if (opt$objective > best) {
      best <- opt$objective
      arg <- opt$maximum
    }
  }
  upper <- max(upper, best)
  c(list(value = min(1, best), upper = min(1, upper), argmax = arg, converged = converged,
         levels = level, evaluations = as.integer(evaluations),
         curvature_bound = scale * max(d2), method = "grid + branch-and-bound + Brent"), base)
}

# Null supremum of an arbitrary set of tables of a design (n1, n2): `member`
# is a logical vector over the tables in the order of .bl_tables() (x1 varies
# fastest). This is the entry point for other unconditional tests that order
# the tables by a different statistic but share the equality null and its
# nuisance maximisation (specs/engine-decisions.md, G19).
.bl_set_sup <- function(n1, n2, member, control = boschloo_control(), alpha = NULL) {
  n1 <- as.integer(n1)
  n2 <- as.integer(n2)
  tab <- .bl_tables(n1, n2)
  if (!is.logical(member) || length(member) != length(tab$t) || anyNA(member)) {
    stop("`member` must be a logical vector with one element per table", call. = FALSE)
  }
  ord <- list(N = n1 + n2, t = tab$t, d = stats::dhyper(tab$x1, n1, n2, tab$t),
              stratum_size = tabulate(tab$t + 1L, n1 + n2 + 1L))
  .bl_sup(.bl_tail_c(ord, which(member)), control, alpha)
}

# Compact diagnostics record of a supremum (for results).
.bl_sup_record <- function(s) {
  list(value = s$value, upper = s$upper, gap = s$upper - s$value,
       argmax = s$argmax, converged = s$converged, bisection_levels = s$levels,
       evaluations = s$evaluations, grid_points = s$grid_points, tolerance = s$tolerance,
       fp_allowance = s$fp_allowance, global_curvature_bound = s$curvature_bound,
       method = s$method)
}

# Null probability curve of a tail set over a grid of pi (plot data).
.bl_null_curve <- function(cvec, pis) {
  as.vector(crossprod(.bl_basis(length(cvec) - 1L, pis), cvec))
}
