# G11 engines against an independent reference (layer 2, decision D18):
# f_variances_two_groups, chi2_variance_one_sample and
# chi2_goodness_of_fit_contingency (R/methods-variances.R).
#
# The reference below is a compact copy of validation/independent/g11/ref_g11.R
# (the full check, 6790 comparisons, runs from validation/): probabilities are
# ratios of integrals of unnormalised densities computed with an own composite
# Gauss-Legendre rule, quantiles and inverse problems are solved by bisection,
# and the noncentral chi-square is the convolution (Z + sqrt(lambda))^2 +
# chi2(k - 1). It calls no stats distribution function, no root finder and no
# inferguidR function; the engine is reached through analysis_spec(),
# calculate_power() and effect_size().
#
# Tolerances (validation/independent/g11/summary.md, fixed before the first
# run): probabilities abs 1e-8 + rel 1e-6 (relative 1e-6 below 1e-8); critical
# values and noncentrality relative 1e-9; integers exact; solved quantities by
# the reference residual (1e-8; compromise relative 1e-6 on beta / (q alpha)).

ref_reset_registry()

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

g11_prob_ok <- function(ours, ref) {
  if (!is.finite(ours) || !is.finite(ref)) return(FALSE)
  if (ref == 0) return(ours < 1e-300)
  if (ref >= 1e-8) return(abs(ours - ref) <= 1e-8 + 1e-6 * abs(ref))
  ours > 0 && abs(log(ours) - log(ref)) <= 1e-6
}
g11_rel_ok <- function(ours, ref, tol = 1e-9) is.finite(ours) && is.finite(ref) && abs(ours / ref - 1) <= tol
g11_run <- function(method, mode, inputs, options = list()) {
  calculate_power(analysis_spec(method, mode, inputs = inputs, options = options))
}
g11_expect_forward <- function(res, ref, info) {
  expect_identical(res$status$code, "ok", info = info)
  expect_true(g11_prob_ok(res$achieved$power, ref$power),
              info = sprintf("%s power ours %.15g ref %.15g", info, res$achieved$power, ref$power))
  expect_true(g11_prob_ok(res$achieved$beta, ref$beta),
              info = sprintf("%s beta ours %.15g ref %.15g", info, res$achieved$beta, ref$beta))
  expect_identical(length(res$quantities$critical), length(ref$critical), info = info)
  for (i in seq_along(ref$critical)) {
    expect_true(g11_rel_ok(res$quantities$critical[i], ref$critical[i]),
                info = sprintf("%s critical[%d] ours %.15g ref %.15g", info, i, res$quantities$critical[i], ref$critical[i]))
  }
}

test_that("the reference reproduces closed forms (it does not call stats)", {
  f11 <- g11r_f_dist(1, 1)
  for (x in c(1e-6, 1, 161.4476)) expect_lt(abs(f11$lower(x) / (2 / pi * atan(sqrt(x))) - 1), 1e-12)
  c2 <- g11r_chisq_dist(2)
  for (x in c(0.5, 5, 50)) expect_lt(abs(c2$upper(x) / exp(-x / 2) - 1), 1e-12)
  f22 <- g11r_f_dist(2, 2)
  expect_lt(abs(f22$lower(3) / 0.75 - 1), 1e-12)
  expect_lt(abs(g11r_nchisq_tails(3, 2, 0)[["upper"]] / exp(-1.5) - 1), 1e-12)
})

# f_variances_two_groups -----------------------------------------------------------------

test_that("two-sample variance test: post hoc power, beta and critical values match the reference", {
  cases <- list(list(1.5, 50, 50, 0.05, "one"), list(1.5, 50, 50, 0.05, "two"), list(0.5, 2, 2, 0.05, "one"),
                list(0.5, 2, 2, 0.05, "two"), list(1, 20, 30, 0.05, "one"), list(0.25, 30, 60, 0.5, "two"),
                list(4, 5, 9, 1e-4, "one"), list(1.02, 1e5, 1e5, 1e-4, "two"), list(0, 20, 30, 0.05, "two"),
                list(1.02, 1647267, 164727, 0.001, "one"))
  for (cs in cases) {
    info <- paste(unlist(cs), collapse = " ")
    res <- g11_run("f_variances_two_groups", "post_hoc", list(ratio = cs[[1]], alpha = cs[[4]], n1 = cs[[2]], n2 = cs[[3]]),
                   list(tails = cs[[5]]))
    ref <- g11r_power_f2(cs[[1]], cs[[2]], cs[[3]], cs[[4]], cs[[5]])
    g11_expect_forward(res, ref, info)
    expect_identical(res$quantities$df, c(cs[[2]] - 1, cs[[3]] - 1), info = info)
  }
})

test_that("two-sample variance test: inverse modes satisfy the reference residual", {
  se <- g11_run("f_variances_two_groups", "sensitivity", list(alpha = 0.05, power = 0.8, n1 = 20, n2 = 60),
                list(tails = "one", effect_direction = "var1_le_var0"))
  expect_identical(se$status$code, "ok")
  expect_lt(se$solved$ratio, 1)
  expect_lt(abs(g11r_power_f2(se$solved$ratio, 20, 60, 0.05, "one")$power - 0.8), 1e-8)
  se2 <- g11_run("f_variances_two_groups", "sensitivity", list(alpha = 0.05, power = 0.8, n1 = 50, n2 = 50),
                 list(tails = "two", effect_direction = "var1_ge_var0"))
  expect_gt(se2$solved$ratio, 1)
  expect_lt(abs(g11r_power_f2(se2$solved$ratio, 50, 50, 0.05, "two")$power - 0.8), 1e-8)
  cr <- g11_run("f_variances_two_groups", "criterion", list(ratio = 1.5, power = 0.8, n1 = 50, n2 = 50), list(tails = "two"))
  expect_lt(abs(g11r_power_f2(1.5, 50, 50, cr$solved$alpha, "two")$power - 0.8), 1e-8)
  co <- g11_run("f_variances_two_groups", "compromise", list(ratio = 0.5, error_ratio = 1, n1 = 30, n2 = 60), list(tails = "one"))
  ref <- g11r_power_f2(0.5, 30, 60, co$solved$alpha, "one")
  expect_lt(abs(ref$beta / co$solved$alpha - 1), 1e-6)
  ap <- solve_design(analysis_spec("f_variances_two_groups", "a_priori", inputs = list(ratio = 1.5, alpha = 0.05, power = 0.8),
                                   options = list(tails = "one")))
  n <- g11r_first_integer(function(n) g11r_power_f2(1.5, n, n, 0.05, "one")$power, 2, 0.8)
  expect_identical(ap$solved$n_total, as.integer(2 * n))
  expect_true(g11_prob_ok(ap$achieved$power, g11r_power_f2(1.5, n, n, 0.05, "one")$power))
})

# chi2_variance_one_sample ---------------------------------------------------------------

test_that("one-sample variance test: post hoc power, beta and critical values match the reference", {
  cases <- list(list(1.5, 50, 0.05, "one"), list(1.5, 50, 0.05, "two"), list(0.5, 2, 0.05, "one"),
                list(0.5, 2, 0.05, "two"), list(1, 20, 0.05, "one"), list(1.02, 1e5, 1e-4, "two"),
                list(4, 3, 0.5, "two"), list(0, 20, 0.05, "one"), list(1e3, 10, 1e-4, "one"))
  for (cs in cases) {
    info <- paste(unlist(cs), collapse = " ")
    res <- g11_run("chi2_variance_one_sample", "post_hoc", list(ratio = cs[[1]], alpha = cs[[3]], n_total = cs[[2]]),
                   list(tails = cs[[4]]))
    g11_expect_forward(res, g11r_power_var1(cs[[1]], cs[[2]], cs[[3]], cs[[4]]), info)
    expect_identical(res$quantities$df, cs[[2]] - 1, info = info)
  }
})

test_that("one-sample variance test: inverse modes satisfy the reference residual", {
  ap <- solve_design(analysis_spec("chi2_variance_one_sample", "a_priori",
                                   inputs = list(ratio = 0.6666667, alpha = 0.05, power = 0.8), options = list(tails = "two")))
  N <- g11r_first_integer(function(N) g11r_power_var1(0.6666667, N, 0.05, "two")$power, 2, 0.8)
  expect_identical(ap$solved$n_total, as.integer(N))
  se <- g11_run("chi2_variance_one_sample", "sensitivity", list(alpha = 0.05, power = 0.8, n_total = 50), list(tails = "two"))
  expect_gt(se$solved$ratio, 1)
  expect_lt(abs(g11r_power_var1(se$solved$ratio, 50, 0.05, "two")$power - 0.8), 1e-8)
  # The root below 1 is reported as the other branch.
  other <- se$diagnostics$other_branch_roots
  expect_length(other, 1L)
  expect_lt(other, 1)
  expect_lt(abs(g11r_power_var1(other, 50, 0.05, "two")$power - 0.8), 1e-7)
  cr <- g11_run("chi2_variance_one_sample", "criterion", list(ratio = 0.5, power = 0.8, n_total = 30), list(tails = "one"))
  expect_lt(abs(g11r_power_var1(0.5, 30, cr$solved$alpha, "one")$power - 0.8), 1e-8)
  co <- g11_run("chi2_variance_one_sample", "compromise", list(ratio = 4, error_ratio = 4, n_total = 2), list(tails = "two"))
  expect_lt(abs(g11r_power_var1(4, 2, co$solved$alpha, "two")$beta / (4 * co$solved$alpha) - 1), 1e-6)
})

# chi2_goodness_of_fit_contingency -------------------------------------------------------

test_that("goodness of fit: post hoc power (noncentral chi-square by convolution) matches the reference", {
  cases <- list(list(0.3, 150, 5, 0.05), list(0.9, 2, 1, 0.05), list(0.01, 1e5, 200, 1e-4), list(0.3, 150, 2.5, 0.05),
                list(0, 150, 5, 0.05), list(0.3, 1, 5, 0.05), list(2.5, 150, 1, 0.5))
  for (cs in cases) {
    info <- paste(unlist(cs), collapse = " ")
    res <- g11_run("chi2_goodness_of_fit_contingency", "post_hoc", list(w = cs[[1]], alpha = cs[[4]], n_total = cs[[2]], df = cs[[3]]))
    ref <- g11r_power_gof(cs[[1]], cs[[2]], cs[[3]], cs[[4]])
    g11_expect_forward(res, ref, info)
    expect_equal(res$quantities$ncp, cs[[2]] * cs[[1]]^2, tolerance = 1e-14, info = info)
  }
})

test_that("goodness of fit: inverse modes satisfy the reference residual", {
  ap <- solve_design(analysis_spec("chi2_goodness_of_fit_contingency", "a_priori", inputs = list(w = 0.3, alpha = 0.05, power = 0.8, df = 5)))
  N <- g11r_first_integer(function(N) g11r_power_gof(0.3, N, 5, 0.05)$power, 1, 0.8)
  expect_identical(ap$solved$n_total, as.integer(N))
  se <- g11_run("chi2_goodness_of_fit_contingency", "sensitivity", list(alpha = 0.05, power = 0.8, n_total = 150, df = 5))
  expect_lt(abs(g11r_power_gof(se$solved$w, 150, 5, 0.05)$power - 0.8), 1e-8)
  cr <- g11_run("chi2_goodness_of_fit_contingency", "criterion", list(w = 0.9, power = 0.95, n_total = 2, df = 1))
  expect_lt(abs(g11r_power_gof(0.9, 2, 1, cr$solved$alpha)$power - 0.95), 1e-8)
  co <- g11_run("chi2_goodness_of_fit_contingency", "compromise", list(w = 0.3, error_ratio = 1, n_total = 150, df = 5))
  expect_lt(abs(g11r_power_gof(0.3, 150, 5, co$solved$alpha)$beta / co$solved$alpha - 1), 1e-6)
})

test_that("drawers agree with their definitions", {
  expect_equal(effect_size("f_variances_two_groups", "from_variances", var0 = 1.5, var1 = 1)$value, 1 / 1.5, tolerance = 1e-15)
  expect_equal(effect_size("chi2_variance_one_sample", "from_variances", var0 = 2, var1 = 3)$value, 1.5, tolerance = 1e-15)
  p0 <- c(0.25, 0.25, 0.25, 0.25)
  p1 <- c(0.4, 0.2, 0.2, 0.2)
  expect_equal(effect_size("chi2_goodness_of_fit_contingency", "cell_probabilities", p0 = p0, p1 = p1)$value,
               g11r_cohen_w(p0, p1), tolerance = 1e-14)
})
