# Copyright (c) 2026 Raban Heller. MIT licence; see ../../../LICENSE.md.
# Independent log-density quadrature for variance tests.
g11r_gl_rule <- function(n) {
  x <- cos(pi * (seq_len(n) - 0.25) / (n + 0.5))
  legendre <- function(x) {
    p0 <- rep(1, length(x))
    p1 <- x
    for (k in 2:n) {
      p2 <- ((2 * k - 1) * x * p1 - (k - 1) * p0) / k
      p0 <- p1
      p1 <- p2
    }
    list(p = p1, dp = n * (x * p1 - p0) / (x^2 - 1))
  }
  for (it in 1:100) {
    L <- legendre(x)
    dx <- L$p / L$dp
    x <- x - dx
    if (max(abs(dx)) < 1e-16) break
  }
  L <- legendre(x)
  list(x = x, w = 2 / ((1 - x^2) * L$dp^2))
}

G11R_GL <- g11r_gl_rule(20)

## Integrals of exp(logf(x) - shift) over the panels [b_i, b_{i+1}].
g11r_panels <- function(logf, breaks, shift) {
  a <- breaks[-length(breaks)]
  b <- breaks[-1]
  half <- (b - a) / 2
  X <- outer(half, G11R_GL$x) + (a + b) / 2
  V <- exp(logf(X) - shift)
  V[!is.finite(V)] <- 0
  as.vector(V %*% G11R_GL$w) * half
}

## A distribution on [lo, hi] given by an unnormalised log density in some
## variable s, with panel breaks and optional analytic masses of the end pieces
## [lo, breaks[1]] and [breaks[n], hi] (already divided by exp(shift)).
g11r_make_dist <- function(logf, breaks, shift, end_lo = 0, end_hi = 0) {
  breaks <- sort(unique(breaks))
  pan <- g11r_panels(logf, breaks, shift)
  pre <- c(end_lo, end_lo + cumsum(pan))          # mass below breaks[i]
  suf <- c(rev(cumsum(rev(pan))) + end_hi, end_hi) # mass above breaks[i]
  total <- pre[length(pre)] + end_hi
  partial <- function(a, b) {
    half <- (b - a) / 2
    X <- outer(half, G11R_GL$x) + (a + b) / 2
    V <- exp(logf(X) - shift)
    V[!is.finite(V)] <- 0
    as.vector(V %*% G11R_GL$w) * half
  }
  # Mass below s (vectorised); s outside the break range uses the end pieces
  # (their partial masses are taken as 0 or total, adequate far in the tails).
  below <- function(s) {
    out <- numeric(length(s))
    nb <- length(breaks)
    i <- findInterval(s, breaks)
    lo_side <- i == 0
    hi_side <- i >= nb
    mid <- !lo_side & !hi_side
    out[lo_side] <- 0
    out[hi_side] <- total
    if (any(mid)) out[mid] <- pre[i[mid]] + partial(breaks[i[mid]], s[mid])
    out
  }
  above <- function(s) {
    out <- numeric(length(s))
    nb <- length(breaks)
    i <- findInterval(s, breaks)
    lo_side <- i == 0
    hi_side <- i >= nb
    mid <- !lo_side & !hi_side
    out[lo_side] <- total
    out[hi_side] <- 0
    if (any(mid)) out[mid] <- suf[i[mid] + 1L] + partial(s[mid], breaks[i[mid] + 1L])
    out
  }
  list(below = below, above = above, total = total, breaks = breaks)
}

## Chi-square(nu) in the variable u = sqrt(x). For nu > 1 the log density is
## taken relative to its mode m = sqrt(nu - 1), written without cancellation:
## (nu - 1) log(u / m) - (u - m)(u + m) / 2.
g11r_chisq_dist <- function(nu) {
  m <- sqrt(max(nu - 1, 0))
  logf <- if (m > 0) {
    function(u) (nu - 1) * log(u / m) - (u - m) * (u + m) / 2
  } else {
    function(u) (nu - 1) * log(u) - u^2 / 2
  }
  sd <- sqrt(0.5)
  lin <- m + sd * seq(-60, 60, by = 0.5)
  lin <- lin[lin > 0]
  first <- min(lin, if (m > 0) m else Inf)
  geo <- c(if (m > 0) m * 2^-(1:60), first * 2^-(1:100))
  breaks <- sort(unique(c(geo, lin)))
  eps <- breaks[1]
  # Analytic mass of [0, eps]: sum_j (-1/2)^j eps^(nu + 2j) / (j! (nu + 2j)),
  # relative to the mode for nu > 1 (density m^(nu - 1) exp(-m^2 / 2) there).
  lead <- nu * log(eps) - if (m > 0) (nu - 1) * log(m) - m^2 / 2 else 0
  end_lo <- exp(lead) * (1 / nu - eps^2 / (2 * (nu + 2)) + eps^4 / (8 * (nu + 4)))
  d <- g11r_make_dist(logf, breaks, 0, end_lo = end_lo, end_hi = 0)
  u <- function(x) sqrt(pmax(x, 0))
  list(lower = function(x) d$below(u(x)) / d$total,
       upper = function(x) d$above(u(x)) / d$total,
       mass = function(x1, x2) g11r_mass(d, u(x1), u(x2)),
       nu = nu)
}

## Mass between s1 < s2 of a distribution object, taken as a difference of
## the two smaller tails (no cancellation near 1).
g11r_mass <- function(d, s1, s2) {
  b2 <- d$below(s2)
  if (b2 <= 0.5 * d$total) (b2 - d$below(s1)) / d$total else (d$above(s1) - d$above(s2)) / d$total
}

## F(d1, d2) in the variable theta = atan(sqrt(d1 x / d2)). With an interior
## mode the log density is taken relative to it:
## a log(sin t / sin mode) + b log(cos t / cos mode), a = d1 - 1, b = d2 - 1.
g11r_f_dist <- function(d1, d2) {
  a <- d1 - 1
  b <- d2 - 1
  half_pi <- pi / 2
  interior <- a > 0 && b > 0
  if (interior) {
    mode <- atan(sqrt(a / b))
    curv <- a / sin(mode)^2 + b / cos(mode)^2
    sm <- sin(mode)
    cm <- cos(mode)
    logf <- function(t) a * log(sin(t) / sm) + b * log(cos(t) / cm)
  } else {
    mode <- if (a <= 0 && b > 0) 0 else if (a > 0 && b <= 0) half_pi else pi / 4
    curv <- max(abs(a) + abs(b), 1)
    logf0 <- function(t) a * log(sin(t)) + b * log(cos(t))
  }
  sd <- 1 / sqrt(curv)
  lin <- mode + sd * seq(-60, 60, by = 0.5)
  lin <- lin[lin > 0 & lin < half_pi]
  lo1 <- if (length(lin)) min(lin) else half_pi / 2
  hi1 <- if (length(lin)) max(lin) else half_pi / 2
  geo_lo <- lo1 * 2^-(1:60)
  geo_hi <- half_pi - (half_pi - hi1) * 2^-(1:40)
  if (interior) geo_lo <- c(geo_lo, mode * 2^-(1:40), mode)
  breaks <- sort(unique(c(geo_lo, lin, geo_hi)))
  e0 <- breaks[1]
  e1 <- half_pi - breaks[length(breaks)]
  if (interior) {
    ref <- a * log(sm) + b * log(cm)
    end_lo <- exp((a + 1) * log(e0) - ref) / (a + 1)
    end_hi <- exp((b + 1) * log(e1) - ref) / (b + 1)
    d <- g11r_make_dist(logf, breaks, 0, end_lo = end_lo, end_hi = end_hi)
  } else {
    shift <- max(logf0(breaks))
    end_lo <- exp((a + 1) * log(e0) - shift) / (a + 1)
    end_hi <- exp((b + 1) * log(e1) - shift) / (b + 1)
    d <- g11r_make_dist(logf0, breaks, shift, end_lo = end_lo, end_hi = end_hi)
  }
  th <- function(x) atan(sqrt(d1 * pmax(x, 0) / d2))
  list(lower = function(x) d$below(th(x)) / d$total,
       upper = function(x) d$above(th(x)) / d$total,
       mass = function(x1, x2) g11r_mass(d, th(x1), th(x2)),
       d1 = d1, d2 = d2)
}

## Quantile by bisection on log(x): P(X <= q) = p (lower) or P(X > q) = p.
g11r_quantile <- function(dist, p, lower, start = 1) {
  f <- if (lower) function(x) dist$lower(x) - p else function(x) p - dist$upper(x)
  lo <- log(start)
  hi <- log(start)
  while (f(exp(lo)) > 0) lo <- lo - 1
  while (f(exp(hi)) < 0) hi <- hi + 1
  for (i in 1:300) {
    mid <- (lo + hi) / 2
    if (f(exp(mid)) < 0) lo <- mid else hi <- mid
    if (hi - lo < 1e-15) break
  }
  exp((lo + hi) / 2)
}

## Normal distribution from chi-square(1).
G11R_CHI1 <- g11r_chisq_dist(1)
g11r_phi_upper <- function(a) {
  ifelse(a >= 0, 0.5 * G11R_CHI1$upper(a^2), 1 - 0.5 * G11R_CHI1$upper(a^2))
}
g11r_phi_lower <- function(a) g11r_phi_upper(-a)

## Noncentral chi-square(k, lambda) tails by convolution.
g11r_nchisq_tails <- function(x, k, lambda, inner = NULL) {
  r <- sqrt(lambda)
  s <- sqrt(x)
  out_lo <- g11r_phi_lower(s - r) - g11r_phi_lower(-s - r)   # P(|Z + r| <= s)
  out_hi <- g11r_phi_upper(s - r) + g11r_phi_lower(-s - r)   # P(|Z + r| > s)
  if (k == 1) return(c(lower = out_lo, upper = out_hi))
  if (is.null(inner)) inner <- g11r_chisq_dist(k - 1)
  lin <- r + seq(-12, 12, by = 0.25)
  lin <- lin[lin > -s & lin < s]
  gl <- if (length(lin)) min(lin) else 0
  gh <- if (length(lin)) max(lin) else 0
  br <- c(-s, s, lin, -s + (gl + s) * 2^-(0:50), s - (s - gh) * 2^-(0:50), -s + 2 * s * (1:64) / 65)
  br <- sort(unique(br[br >= -s & br <= s]))
  a <- br[-length(br)]
  b <- br[-1]
  half <- (b - a) / 2
  Z <- outer(half, G11R_GL$x) + (a + b) / 2
  z <- as.vector(Z)
  phi <- exp(-(z - r)^2 / 2) / sqrt(2 * pi)
  t <- pmax(x - z^2, 0)
  wv <- as.vector(outer(half, G11R_GL$w))
  lo <- sum(wv * phi * inner$lower(t))
  hi <- sum(wv * phi * inner$upper(t))
  c(lower = lo, upper = out_hi + hi)
}

## Power of the test with statistic ratio * X under H1 (X with the given
## distribution object); returns power, beta and the critical value(s).
g11r_scaled_power <- function(dist, ratio, alpha, tails, side, start = 1) {
  if (tails == "two") {
    lo <- g11r_quantile(dist, alpha / 2, TRUE, start)
    hi <- g11r_quantile(dist, alpha / 2, FALSE, start)
    if (ratio == 0) return(list(power = 1, beta = 0, critical = c(lo, hi)))
    pw <- dist$lower(lo / ratio) + dist$upper(hi / ratio)
    beta <- dist$mass(lo / ratio, hi / ratio)
    return(list(power = pw, beta = beta, critical = c(lo, hi)))
  }
  if (side == "upper") {
    c <- g11r_quantile(dist, alpha, FALSE, start)
    return(list(power = dist$upper(c / ratio), beta = dist$lower(c / ratio), critical = c))
  }
  c <- g11r_quantile(dist, alpha, TRUE, start)
  if (ratio == 0) return(list(power = 1, beta = 0, critical = c))
  list(power = dist$lower(c / ratio), beta = dist$upper(c / ratio), critical = c)
}

## F test of two variances: statistic s1^2 / s2^2 ~ ratio * F(n1 - 1, n2 - 1);
## one tail: upper for ratio >= 1, lower for ratio < 1 (G*Power 3.1.9.7).
g11r_power_f2 <- function(ratio, n1, n2, alpha, tails, dist = NULL) {
  if (is.null(dist)) dist <- g11r_f_dist(n1 - 1, n2 - 1)
  side <- if (ratio >= 1) "upper" else "lower"
  c(g11r_scaled_power(dist, ratio, alpha, tails, side), list(df = c(n1 - 1, n2 - 1)))
}

## Chi-square test of one variance: (N - 1) s^2 / var0 ~ ratio * chi2(N - 1);
## one tail: upper for ratio > 1, lower for ratio <= 1 (G*Power 3.1.9.7).
g11r_power_var1 <- function(ratio, N, alpha, tails, dist = NULL) {
  if (is.null(dist)) dist <- g11r_chisq_dist(N - 1)
  side <- if (ratio > 1) "upper" else "lower"
  c(g11r_scaled_power(dist, ratio, alpha, tails, side, start = N - 1), list(df = N - 1))
}

## Pearson chi-square test: chi2(df) under H0, noncentral chi2(df, N w^2) under H1.
g11r_power_gof <- function(w, N, df, alpha, central = NULL, inner = NULL) {
  if (is.null(central)) central <- g11r_chisq_dist(df)
  crit <- g11r_quantile(central, alpha, FALSE, start = df)
  lambda <- N * w^2
  if (lambda == 0) return(list(power = central$upper(crit), beta = central$lower(crit), critical = crit, ncp = 0))
  tl <- g11r_nchisq_tails(crit, df, lambda, inner)
  list(power = tl[["upper"]], beta = tl[["lower"]], critical = crit, ncp = lambda)
}

## Cohen's w from cell probabilities (definition, Cohen 1988 eq. 7.2.1).
g11r_cohen_w <- function(p0, p1) sqrt(sum((p1 - p0)^2 / p0))

## Plain bisection for an increasing (or decreasing) function g on [lo, hi]
## with g(lo) and g(hi) of opposite sign.
g11r_bisect <- function(g, lo, hi, iter = 200, tol = 1e-15) {
  glo <- g(lo)
  for (i in seq_len(iter)) {
    mid <- (lo + hi) / 2
    gm <- g(mid)
    if (sign(gm) == sign(glo)) {
      lo <- mid
      glo <- gm
    } else {
      hi <- mid
    }
    if (abs(hi - lo) <= tol * max(1, abs(lo))) break
  }
  (lo + hi) / 2
}

## Smallest integer n >= lo with f(n) >= target (f increasing in n).
g11r_first_integer <- function(f, lo, target, hi_limit = 2^31) {
  if (f(lo) >= target) return(lo)
  a <- lo
  b <- max(lo + 1, 2 * lo)
  while (f(b) < target) {
    a <- b
    b <- 2 * b
    if (b > hi_limit) return(NA_real_)
  }
  while (b - a > 1) {
    m <- floor((a + b) / 2)
    if (f(m) >= target) b <- m else a <- m
  }
  b
}

# Helpers --------------------------------------------------------------------------------
