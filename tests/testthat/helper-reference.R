# Independent reference primitives (decision D18; specs/engine-decisions.md,
# Appendix A).
#
# These functions are written from the defining representations and share no
# numerical routine with the engines under test:
#
# * they never call pt()/qt()/pf()/pchisq() with a noncentrality parameter,
#   nor R's noncentral algorithms (nmath pnt, pnbeta, pnchisq);
# * they never call any inferguidR function (ref_reset_registry() at the end
#   is a test utility, not a reference primitive);
# * the noncentral t reference integrates the normal CDF over the chi-square
#   density (the engine's fallback integrates the chi-square CDF over the
#   normal density, and its main route is nmath's AS 243 series);
# * central t quantiles are found by root finding on the same integral, not
#   with qt();
# * noncentral chi-square uses the convolution (Z + sqrt(lambda))^2 +
#   chi2(k - 1), noncentral F integrates that convolution over the central
#   chi-square density of the denominator (the engine sums Poisson mixtures);
# * ref_pt() and the engine's integral route both use stats::integrate()
#   (QUADPACK), with different integrands; ref_pt_gl() is a second noncentral
#   t reference with its own Gauss-Legendre quadrature for that regime;
# * binomial probabilities are summed from lchoose() on the log scale.
#
# Only central primitives of base R are used (pnorm, dnorm, dchisq, pchisq
# without ncp, qchisq without ncp, lchoose). They are slow by design and are
# used on selected points in tests and on grids in validation/independent/t/.

# Integrate exp(logf(u)) over the real line, split around the mode of logf.
# Returns the log of the integral.
ref_log_integrate <- function(logf, lower = -Inf, upper = Inf, grid) {
  v <- logf(grid)
  if (!any(is.finite(v))) return(-Inf)
  u0 <- grid[which.max(v)]
  step <- if (length(grid) > 1L) diff(range(grid)) / (length(grid) - 1L) else 1
  opt <- stats::optimize(function(u) {
    w <- logf(u)
    if (is.finite(w)) -w else 1e300
  }, c(max(lower, u0 - 2 * step), min(upper, u0 + 2 * step)), tol = 1e-12)
  um <- if (logf(opt$minimum) > max(v)) opt$minimum else u0
  lm <- logf(um)
  # Width of the peak from the curvature of logf at the mode; the pieces are
  # cut at multiples of it so that a narrow peak is never missed.
  h <- 1e-4 * max(1, abs(um))
  curv <- (logf(um + h) - 2 * lm + logf(um - h)) / h^2
  width <- if (is.finite(curv) && curv < 0) 1 / sqrt(-curv) else 1
  g <- function(u) exp(logf(u) - lm)
  offsets <- width * c(1, 4, 16, 64, 256)
  cuts <- sort(unique(c(lower, um - offsets, um, um + offsets, upper)))
  cuts <- cuts[cuts >= lower & cuts <= upper]
  # The scaled integrand inherits the rounding error of logf, about |lm| eps,
  # so no quadrature can deliver a smaller relative error than that; asking
  # for 1e-12 where |lm| is ~1e6 only produces roundoff errors.
  rel_tol <- max(1e-12, 64 * .Machine$double.eps * abs(lm))
  total <- 0
  for (i in seq_len(length(cuts) - 1L)) {
    if (cuts[i] < cuts[i + 1L]) {
      total <- total + stats::integrate(g, cuts[i], cuts[i + 1L], rel.tol = rel_tol, abs.tol = 0,
                                        subdivisions = 5000L)$value
    }
  }
  log(total) + lm
}

# log(a log a - a - lgamma(a)) without cancellation for large a: by
# Stirling's formula it equals 0.5 log(a) - 0.5 log(2 pi) - e(a), where the
# Stirling error e(a) = lgamma(a) - ((a - 0.5) log a - a + 0.5 log(2 pi)) is
# taken from its asymptotic series for a >= 15.
ref_chisq_const <- function(a) {
  if (a < 15) return(a * log(a) - a - lgamma(a))
  stirling_error <- 1 / (12 * a) - 1 / (360 * a^3) + 1 / (1260 * a^5) - 1 / (1680 * a^7)
  0.5 * log(a) - 0.5 * log(2 * pi) - stirling_error
}

# Noncentral t tail probability, P(T <= t) (lower) or P(T > t) (upper):
#   P(T <= t) = int_0^Inf Phi(t sqrt(v / df) - ncp) f_chi2(v; df) dv.
# The integral is taken over x = log(v / df), where
#   f_chi2(v; df) dv = exp(a (x - expm1(x)) + C(a)) dx,  a = df / 2,
# which stays accurate for very large df and where exp() underflows.
ref_pt <- function(t, df, ncp, lower_tail = TRUE, log_p = FALSE) {
  if (is.infinite(t)) {
    one <- (t > 0) == lower_tail
    return(if (log_p) (if (one) 0 else -Inf) else as.double(one))
  }
  a <- df / 2
  const <- ref_chisq_const(a)
  logf <- function(x) {
    stats::pnorm(t * exp(x / 2) - ncp, lower.tail = lower_tail, log.p = TRUE) +
      a * (x - expm1(x)) + const
  }
  grid <- seq(log(stats::qchisq(1e-12, df) / df) - 30,
              log(stats::qchisq(1e-12, df, lower.tail = FALSE) / df) + 5, length.out = 400L)
  lp <- ref_log_integrate(logf, grid = grid)
  if (log_p) lp else exp(lp)
}

# Upper-tail quantile of the central t distribution from the integral CDF:
# the x with P(T > x) = p. f below is decreasing with f(0) = log(0.5 / p), so
# the root is positive for p < 1/2 and negative for p > 1/2 (a one-sided
# criterion alpha above one half).
ref_qt_upper <- function(p, df) {
  if (p == 0.5) return(0)
  f <- function(x) log(ref_pt(x, df, 0, lower_tail = FALSE)) - log(p)
  if (p < 0.5) {
    lo <- 0
    hi <- 1
    while (f(hi) > 0) hi <- hi * 2
  } else {
    hi <- 0
    lo <- -1
    while (f(lo) < 0) lo <- lo * 2
  }
  stats::uniroot(f, c(lo, hi), tol = 1e-14, maxiter = 1000L)$root
}

# Power of the two-sample t test (pooled variance) from its definition:
# df = n1 + n2 - 2, ncp = d sqrt(n1 n2 / (n1 + n2)); two-sided: alpha/2 per
# tail, both tails counted; one-sided: rejection in the direction of d.
ref_power_t2 <- function(d, n1, n2, alpha, tails) {
  df <- n1 + n2 - 2
  ncp <- d * sqrt(n1 * n2 / (n1 + n2))
  if (tails == "two") {
    tc <- ref_qt_upper(alpha / 2, df)
    power <- ref_pt(tc, df, ncp, lower_tail = FALSE) + ref_pt(-tc, df, ncp, lower_tail = TRUE)
    critical <- c(-tc, tc)
  } else {
    tc <- ref_qt_upper(alpha, df)
    if (d >= 0) {
      power <- ref_pt(tc, df, ncp, lower_tail = FALSE)
      critical <- tc
    } else {
      power <- ref_pt(-tc, df, ncp, lower_tail = TRUE)
      critical <- -tc
    }
  }
  list(power = power, critical = critical, df = df, ncp = ncp)
}

# Noncentral chi-square tail by convolution (k > 1 uses chi2(k - 1); k = 1
# is closed form). Only central pnorm/pchisq/dnorm are used.
ref_pchisq <- function(y, k, lambda, lower_tail = TRUE) {
  if (y <= 0) return(if (lower_tail) 0 else 1)
  s <- sqrt(y)
  mu <- sqrt(lambda)
  if (k == 1) {
    inside <- stats::pnorm(s - mu) - stats::pnorm(-s - mu)
    outside <- stats::pnorm(s - mu, lower.tail = FALSE) + stats::pnorm(-s - mu)
    return(if (lower_tail) inside else outside)
  }
  inner <- function(u) {
    stats::dnorm(u - mu) * stats::pchisq(pmax(y - u^2, 0), k - 1, lower.tail = lower_tail)
  }
  val <- stats::integrate(inner, -s, s, rel.tol = 1e-13, abs.tol = 0, subdivisions = 2000L)$value
  if (lower_tail) {
    val
  } else {
    val + stats::pnorm(s - mu, lower.tail = FALSE) + stats::pnorm(-s - mu)
  }
}

# Noncentral F tail by integrating the convolution chi-square over the
# central chi-square(df2) density of the denominator:
#   P(F > x) = int P(chi2(df1, lambda) > x df1 v / df2) f_chi2(v; df2) dv.
ref_pf <- function(x, df1, df2, lambda, lower_tail = TRUE) {
  outer <- function(v) {
    vapply(v, function(vi) ref_pchisq(x * df1 * vi / df2, df1, lambda, lower_tail), numeric(1)) *
      stats::dchisq(v, df2)
  }
  lo <- stats::qchisq(1e-14, df2)
  hi <- stats::qchisq(1e-14, df2, lower.tail = FALSE)
  md <- stats::qchisq(0.5, df2)
  stats::integrate(outer, lo, md, rel.tol = 1e-11, abs.tol = 0)$value +
    stats::integrate(outer, md, hi, rel.tol = 1e-11, abs.tol = 0)$value
}

# Binomial tails from lchoose() on the log scale: P(X <= q) or P(X > q).
ref_pbinom <- function(q, size, prob, lower_tail = TRUE) {
  q <- floor(q)
  idx <- if (lower_tail) seq.int(0, min(q, size)) else seq.int(max(q + 1, 0), size)
  if (!length(idx) || (lower_tail && q < 0) || (!lower_tail && q >= size)) return(0)
  lt <- lchoose(size, idx) + idx * log(prob) + (size - idx) * log1p(-prob)
  mx <- max(lt)
  exp(mx + log(sum(exp(lt - mx))))
}

# Second noncentral t reference without QUADPACK. The engine's integral route
# (R/num-dist.R) and ref_pt() above both use stats::integrate(); this one
# uses composite Gauss-Legendre quadrature with nodes computed here
# (Golub-Welsch, eigen()), over x = log S for the scaled chi variable
# S = sqrt(V / df):
#   P(T > t) = E[Phi(delta - t S)],  P(T <= t) = E[Phi(t S - delta)],
# where the density of x is proportional to exp(df x - (df / 2) expm1(2 x))
# and the expectation is normalised by the integral of that density computed
# with the same nodes (no lgamma or Stirling constant enters). It is used in
# the regime where the engine integrates (|ncp| > 37.62, df > 4e5, small
# tails), so that agreement there does not rest on a shared quadrature routine.
ref_gl_nodes <- local({
  n <- 24L
  k <- seq_len(n - 1L)
  b <- k / sqrt(4 * k^2 - 1)
  J <- matrix(0, n, n)
  J[cbind(k, k + 1L)] <- b
  J[cbind(k + 1L, k)] <- b
  e <- eigen(J, symmetric = TRUE)
  o <- order(e$values)
  list(x = e$values[o], w = 2 * e$vectors[1L, o]^2)
})

# log of the integral of exp(logdens(x) + logcond(x)) over x = log S.
ref_gl_log_integral <- function(logcond, df) {
  g <- function(x) df * x - (df / 2) * expm1(2 * x) + logcond(x)
  sdx <- 1 / sqrt(2 * df)
  grid <- sort(unique(c(seq(-745, 6, by = 0.01), sdx * seq(-80, 80, by = 0.1))))
  gv <- g(grid)
  gv[!is.finite(gv)] <- -Inf
  i <- which.max(gv)
  opt <- stats::optimize(function(x) {
    v <- g(x)
    if (is.finite(v)) -v else 1e300
  }, grid[c(max(1L, i - 1L), min(length(grid), i + 1L))], tol = 1e-13)
  xm <- if (-opt$objective > gv[i]) opt$minimum else grid[i]
  gm <- g(xm)
  h <- 1e-5 * max(sdx, 1e-3)
  curv <- (g(xm + h) - 2 * gm + g(xm - h)) / h^2
  w <- min(if (is.finite(curv) && curv < 0) 1 / sqrt(-curv) else sdx, 1)
  # Window: walk out until the log integrand is 110 below its maximum; keep
  # every scan point still within 110 of it.
  walk <- function(dir) {
    x <- xm
    step <- w / 2
    for (k in seq_len(5000L)) {
      x <- x + dir * step
      gx <- g(x)
      if (!is.finite(gx) || gx < gm - 110 || x < -745 || x > 8) break
      step <- step * 1.1
    }
    x
  }
  keep <- grid[gv > gm - 110]
  xl <- max(min(walk(-1), min(keep) - 0.01), -745)
  xr <- min(max(walk(1), max(keep) + 0.01), 8)
  npan <- min(60000L, max(50L, ceiling((xr - xl) / (w / 4))))
  edges <- seq(xl, xr, length.out = npan + 1L)
  half <- diff(edges) / 2
  X <- outer(half, ref_gl_nodes$x) + (edges[-1L] - half)
  W <- outer(half, ref_gl_nodes$w)
  lv <- g(X)
  top <- max(lv[is.finite(lv)])
  log(sum(W * exp(lv - top))) + top
}

ref_pt_gl <- function(t, df, ncp, lower_tail = TRUE, log_p = FALSE) {
  cond <- function(x) stats::pnorm(t * exp(x) - ncp, lower.tail = lower_tail, log.p = TRUE)
  lp <- ref_gl_log_integral(cond, df) - ref_gl_log_integral(function(x) 0 * x, df)
  if (log_p) lp else exp(lp)
}

# Comparison helper: absolute plus relative tolerance (VALIDATION_CONTRACT
# starting rule abs 1e-8 + rel 1e-6 is looser than these; the tighter values
# used in the tests are justified next to each use).
ref_close <- function(x, ref, abs_tol, rel_tol) {
  all(abs(x - ref) <= abs_tol + rel_tol * abs(ref))
}

# Re-create the registry state from the package's own providers so that a
# test file starts from the shipped definitions and solvers, whatever other
# test files registered or removed before.
ref_reset_registry <- function() {
  .ig_registry_init(ns = asNamespace("inferguidR"))
  invisible(TRUE)
}
