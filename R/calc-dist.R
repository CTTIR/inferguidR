# Distribution calculator: distributions without a base R implementation -----------------
#
# Ownership: phase 10 utilities (group G12). Method notes:
# specs/method-notes/utilities.md, section "Distribution calculator".
#
# * Sample correlation coefficient r of N bivariate normal pairs with
#   population correlation rho (calculator corrcdf/corrpdf/corrinv).
#   Density (Hotelling 1953, eq. 23; Johnson, Kotz & Balakrishnan, vol. 2,
#   ch. 32, eq. 32.7):
#     f(r) = (N-2) Gamma(N-1) / (sqrt(2 pi) Gamma(N-1/2)) (1-rho^2)^((N-1)/2)
#            (1-r^2)^((N-4)/2) (1-rho r)^(3/2-N) 2F1(1/2, 1/2; N-1/2; (1+rho r)/2),
#   for rho = 0 the Student form (1-r^2)^((N-4)/2) / B(1/2, (N-2)/2). The CDF
#   integrates the density in Fisher's z = atanh(r) (stats::integrate,
#   pieces around the mode) and the quantile inverts the CDF.
# * Squared multiple correlation R^2 of a linear regression with k - 1
#   normally distributed (random) predictors, N cases and population rho^2
#   (calculator mr2cdf/mr2pdf/mr2inv; the k convention is the one of the
#   'G*Power' manual ch. 2 and was confirmed by the null-case capture
#   mr2pdf(0.5,0,2,10) = dbeta(0.5, 1/2, 4)). Fisher's (1928) mixture
#   representation (Johnson, Kotz & Balakrishnan, vol. 2, ch. 32, eq. 32.62):
#     R^2 ~ sum_j w_j Beta((k-1)/2 + j, (N-k)/2),
#     w_j = negative binomial(j; size (N-1)/2, prob 1 - rho^2).
# * Laplace (closed forms) and the continuous binomial quantile of the
#   'G*Power' calculator (binoinv inverts I_{1-pi}(N-x, x+1) in x; captured
#   live: binoinv(0.5,10,0.5) = 4.5).
#
# The independent checks of every routine here are in validation/independent/g12/
# (different representations and quadrature; D18).

# Gauss hypergeometric series sum_k (a)_k (b)_k / ((c)_k k!) x^k for |x| < 1,
# summed in vectorised blocks. The series is accepted when, beyond the index
# from which all factors are positive, the last term times 1 / (1 - ratio)
# (a bound of the remainder) is below 1e-17 of the sum; NaN when that is not
# reached within `max_terms` terms.
.ig_hyp2f1_series <- function(a, b, c, x, max_terms = 2e5) {
  if (x == 0) return(1)
  total <- 1
  term <- 1
  k0 <- 0
  block <- 16
  k_pos <- max(0, -a, -b, -c)
  repeat {
    k <- k0 + seq_len(block) - 1
    ratio <- (a + k) * (b + k) * x / ((c + k) * (k + 1))
    terms <- term * cumprod(ratio)
    total <- total + sum(terms)
    term <- terms[block]
    last <- abs(ratio[block])
    k0 <- k0 + block
    if (!is.finite(total)) return(NaN)
    if (k0 > k_pos && last < 1 && abs(term) <= abs(total) * 1e-17 * (1 - last)) break
    if (k0 >= max_terms) return(NaN)
    block <- min(2 * block, 65536)
  }
  total
}

# log 2F1(1/2, 1/2; c; x) for c > 1 and 0 <= x < 1, vectorised over x.
# The power series in x converges slowly when x is close to 1 and c is small
# (terms decay like k^(-c) x^k; for c >= 8 a few hundred terms suffice for
# every x). For c < 8 and x > 0.9 the connection formula in 1 - x is used
# (Abramowitz & Stegun 15.3.6 with a = b = 1/2, s = c - 1, y = 1 - x):
#   F(c; x) = Gamma(s+1) Gamma(s) / Gamma(s+1/2)^2 * 2F1(1/2, 1/2; 1-s; y)
#             - (y/x)^s / sin(pi s) * F(c; y),
# which needs s away from the integers (for the sample correlation with a
# whole number N, s = N - 3/2 is a half-integer). NaN when no representation
# converges.
.ig_log_hyp2f1_half <- function(c, x) {
  s <- c - 1
  frac_s <- s - floor(s)
  use_connection <- c < 8 && frac_s >= 0.25 && frac_s <= 0.75
  out <- numeric(length(x))
  for (i in seq_along(x)) {
    xi <- x[i]
    if (is.na(xi) || xi < 0 || xi >= 1) {
      out[i] <- NaN
    } else if (use_connection && xi > 0.9) {
      y <- 1 - xi
      lead <- exp(lgamma(s + 1) + lgamma(s) - 2 * lgamma(s + 0.5)) * .ig_hyp2f1_series(0.5, 0.5, 1 - s, y)
      corr <- (y / xi)^s / sin(pi * s) * .ig_hyp2f1_series(0.5, 0.5, c, y)
      out[i] <- log(lead - corr)
    } else {
      out[i] <- log(.ig_hyp2f1_series(0.5, 0.5, c, xi))
    }
  }
  out
}

# Log density of the sample correlation coefficient (N >= 3, |rho| < 1).
.ig_corr_logpdf <- function(r, rho, n) {
  out <- rep(-Inf, length(r))
  inside <- abs(r) < 1
  if (!any(inside)) return(out)
  ri <- r[inside]
  if (rho == 0) {
    out[inside] <- (n - 4) / 2 * log1p(-ri^2) - lbeta(0.5, (n - 2) / 2)
    return(out)
  }
  out[inside] <- log(n - 2) + lgamma(n - 1) - 0.5 * log(2 * pi) - lgamma(n - 0.5) +
    (n - 1) / 2 * log1p(-rho^2) + (n - 4) / 2 * log1p(-ri^2) + (1.5 - n) * log1p(-rho * ri) +
    .ig_log_hyp2f1_half(n - 0.5, (1 + rho * ri) / 2)
  out
}

.ig_corr_check <- function(rho, n) {
  if (!is.finite(rho) || abs(rho) >= 1) return("rho must lie in (-1, 1)")
  if (!is.finite(n) || n < 3) return("N must be at least 3")
  NULL
}

# Density of r.
.ig_corr_pdf <- function(r, rho, n) {
  bad <- .ig_corr_check(rho, n)
  if (!is.null(bad)) .ig_abort("invalid_input", bad)
  if (abs(r) > 1) return(0)
  if (abs(r) == 1) {
    # Limits at the boundary: 0 for N > 4, constant (rho = 0) or finite for N = 4,
    # unbounded for N < 4.
    if (n > 4) return(0)
    if (n < 4) return(Inf)
    return(exp(.ig_corr_logpdf(sign(r) * (1 - 1e-12), rho, n)))
  }
  exp(.ig_corr_logpdf(r, rho, n))
}

# Log density of z = atanh(r) (the density of r times dr/dz = sech(z)^2),
# written in z so that 1 - r^2 = sech(z)^2 keeps its relative accuracy in the
# tails: log(1 - r^2) = 2 (log 2 - |z| - log1p(exp(-2|z|))).
.ig_corr_logpdf_z <- function(z, rho, n) {
  az <- abs(z)
  log_sech2 <- 2 * (log(2) - az - log1p(exp(-2 * az)))
  r <- tanh(z)
  out <- if (rho == 0) {
    (n - 4) / 2 * log_sech2 - lbeta(0.5, (n - 2) / 2)
  } else {
    log(n - 2) + lgamma(n - 1) - 0.5 * log(2 * pi) - lgamma(n - 0.5) +
      (n - 1) / 2 * log1p(-rho^2) + (n - 4) / 2 * log_sech2 + (1.5 - n) * log1p(-rho * r) +
      .ig_log_hyp2f1_half(n - 0.5, (1 + rho * r) / 2)
  }
  out <- out + log_sech2
  out[!is.finite(z)] <- -Inf
  out
}

# The same log density relative to its value at a reference point zr (near
# the mode). For large N the absolute log density is a sum of terms of order
# N whose rounding (about N * eps) would limit the CDF to about 1e-12; near zr
# the differences are formed directly (d = z - zr):
#   log sech^2 z - log sech^2 zr = -2 log1p(2 sinh^2(d/2) + tanh(zr) sinh(d)),
#   log1p(-rho tanh z) - log1p(-rho tanh zr)
#     = log1p(-rho sinh(d) / (cosh z cosh zr (1 - rho tanh zr))).
# Far from zr (|d| >= 1) the absolute form is used; precision there is
# immaterial for the integral.
.ig_corr_logpdf_z_rel <- function(z, rho, n, zr) {
  out <- .ig_corr_logpdf_z(z, rho, n) - .ig_corr_logpdf_z(zr, rho, n)
  if (rho == 0) return(out)
  d <- z - zr
  near <- is.finite(z) & abs(d) < 1 & abs(z) < 300
  if (any(near)) {
    dn <- d[near]
    zn <- z[near]
    ls <- -2 * log1p(2 * sinh(dn / 2)^2 + tanh(zr) * sinh(dn))
    l1 <- log1p(-rho * sinh(dn) / (cosh(zn) * cosh(zr) * (1 - rho * tanh(zr))))
    h <- .ig_log_hyp2f1_half(n - 0.5, (1 + rho * tanh(zn)) / 2) - .ig_log_hyp2f1_half(n - 0.5, (1 + rho * tanh(zr)) / 2)
    out[near] <- (n - 2) / 2 * ls + (1.5 - n) * l1 + h
  }
  out
}

# Lower and upper tail of r at r0: integrals of the density of z = atanh(r)
# on the log scale (.ig_integrate_log_peak(), cut around the mode), so that
# tails far below 1e-16 keep their relative accuracy. The two tails are
# normalised by their sum, which must be 1 to 1e-8.
.ig_corr_tails <- function(r0, rho, n) {
  bad <- .ig_corr_check(rho, n)
  if (!is.null(bad)) .ig_abort("invalid_input", bad)
  if (r0 <= -1) return(c(lower = 0, upper = 1))
  if (r0 >= 1) return(c(lower = 1, upper = 0))
  if (rho == 0) {
    # Student t with N - 2 df: t = r sqrt(N-2) / sqrt(1 - r^2).
    t0 <- r0 * sqrt(n - 2) / sqrt((1 - r0) * (1 + r0))
    return(c(lower = stats::pt(t0, n - 2), upper = stats::pt(t0, n - 2, lower.tail = FALSE)))
  }
  m0 <- atanh(rho) + rho / (2 * (n - 1))
  s <- 1 / sqrt(max(n - 3, 1))
  opt <- stats::optimize(function(z) -.ig_corr_logpdf_z(z, rho, n), c(m0 - 20 * s, m0 + 20 * s), tol = 1e-12)
  zmode <- opt$minimum
  logf <- function(z) .ig_corr_logpdf_z_rel(z, rho, n, zmode)
  z0 <- atanh(r0)
  warn <- character()
  integ <- function(lo, hi, mode) {
    withCallingHandlers(.ig_integrate_log_peak(logf, lo, hi, mode),
                        warning = function(w) {
                          warn <<- c(warn, conditionMessage(w))
                          invokeRestart("muffleWarning")
                        })
  }
  ll <- integ(-Inf, z0, min(zmode, z0))
  lu <- integ(z0, Inf, max(zmode, z0))
  if (length(warn)) .ig_abort("numerical_failure", paste("correlation CDF:", warn[1]))
  # The integrals are relative to the density at the mode; their sum times
  # that density must be 1.
  total <- exp(.ig_logspace_add(ll, lu) + .ig_corr_logpdf_z(zmode, rho, n))
  if (!is.finite(total) || abs(total - 1) > 1e-8) {
    .ig_abort("numerical_failure", sprintf("the correlation density integrates to %s instead of 1",
                                           format(total, digits = 10)))
  }
  c(lower = stats::plogis(ll - lu), upper = stats::plogis(lu - ll))
}

.ig_corr_cdf <- function(r, rho, n) unname(.ig_corr_tails(r, rho, n)[["lower"]])

# Quantile of r: root of the CDF in z = atanh(r).
.ig_corr_inv <- function(p, rho, n) {
  bad <- .ig_corr_check(rho, n)
  if (!is.null(bad)) .ig_abort("invalid_input", bad)
  if (p <= 0) return(-1)
  if (p >= 1) return(1)
  if (rho == 0) {
    t <- stats::qt(p, n - 2)
    return(t / sqrt(n - 2 + t^2))
  }
  f <- function(z) {
    tl <- .ig_corr_tails(tanh(z), rho, n)
    if (p <= 0.5) tl[["lower"]] - p else (1 - p) - tl[["upper"]]
  }
  m <- atanh(rho)
  s <- 1 / sqrt(max(n - 3, 1))
  lo <- m - 2 * s
  hi <- m + 2 * s
  step <- 2 * s
  while (f(lo) > 0) {
    step <- step * 2
    lo <- lo - step
    if (lo < -40) break
  }
  step <- 2 * s
  while (f(hi) < 0) {
    step <- step * 2
    hi <- hi + step
    if (hi > 40) break
  }
  root <- stats::uniroot(f, c(lo, hi), tol = 1e-300, maxiter = 5000L)
  tanh(root$root)
}

# Squared multiple correlation (random predictors) ------------------------------------------

.ig_mr2_check <- function(rho2, k, n) {
  if (!is.finite(rho2) || rho2 < 0 || rho2 >= 1) return("rho^2 must lie in [0, 1)")
  if (!is.finite(k) || k < 2) return("k (number of predictors plus one) must be at least 2")
  if (!is.finite(n) || n <= k) return("N must be greater than k")
  NULL
}

# Negative binomial mixture weights j and log w_j covering all but 1e-18 of
# the mass.
.ig_mr2_weights <- function(rho2, n) {
  m <- (n - 1) / 2
  if (rho2 == 0) return(list(j = 0, logw = 0))
  q <- 1 - rho2
  lo <- stats::qnbinom(1e-18, size = m, prob = q)
  hi <- stats::qnbinom(1e-18, size = m, prob = q, lower.tail = FALSE)
  if (!is.finite(lo) || !is.finite(hi) || hi - max(0, lo) > .ig_calc_range$mixture_terms) {
    .ig_abort("unsupported", sprintf(paste("the mixture representation of the R-squared distribution needs more than %s terms",
                                           "for rho^2 = %s and N = %s (rho^2 too close to 1 for this N)"),
                                     format(.ig_calc_range$mixture_terms), format(rho2), format(n)))
  }
  j <- seq(max(0, lo), hi)
  list(j = j, logw = stats::dnbinom(j, size = m, prob = q, log = TRUE))
}

.ig_mr2_cdf <- function(x, rho2, k, n) {
  bad <- .ig_mr2_check(rho2, k, n)
  if (!is.null(bad)) .ig_abort("invalid_input", bad)
  if (x <= 0) return(0)
  if (x >= 1) return(1)
  a <- (k - 1) / 2
  b <- (n - k) / 2
  w <- .ig_mr2_weights(rho2, n)
  lt <- w$logw + stats::pbeta(x, a + w$j, b, log.p = TRUE)
  mx <- max(lt)
  if (!is.finite(mx)) return(0)
  min(1, exp(mx + log(sum(exp(lt - mx)))))
}

# Upper tail P(R^2 > x), summed directly (no cancellation near 1).
.ig_mr2_upper <- function(x, rho2, k, n) {
  if (x <= 0) return(1)
  if (x >= 1) return(0)
  a <- (k - 1) / 2
  b <- (n - k) / 2
  w <- .ig_mr2_weights(rho2, n)
  lt <- w$logw + stats::pbeta(x, a + w$j, b, lower.tail = FALSE, log.p = TRUE)
  mx <- max(lt)
  if (!is.finite(mx)) return(0)
  min(1, exp(mx + log(sum(exp(lt - mx)))))
}

.ig_mr2_pdf <- function(x, rho2, k, n) {
  bad <- .ig_mr2_check(rho2, k, n)
  if (!is.null(bad)) .ig_abort("invalid_input", bad)
  if (x < 0 || x > 1) return(0)
  a <- (k - 1) / 2
  b <- (n - k) / 2
  w <- .ig_mr2_weights(rho2, n)
  lt <- w$logw + stats::dbeta(x, a + w$j, b, log = TRUE)
  mx <- max(lt)
  if (is.infinite(mx) && mx > 0) return(Inf)
  if (!is.finite(mx)) return(0)
  exp(mx + log(sum(exp(lt - mx))))
}

.ig_mr2_inv <- function(p, rho2, k, n) {
  bad <- .ig_mr2_check(rho2, k, n)
  if (!is.null(bad)) .ig_abort("invalid_input", bad)
  if (p <= 0) return(0)
  if (p >= 1) return(1)
  if (rho2 == 0) return(stats::qbeta(p, (k - 1) / 2, (n - k) / 2))
  f <- if (p > 0.5) function(x) (1 - p) - .ig_mr2_upper(x, rho2, k, n) else function(x) .ig_mr2_cdf(x, rho2, k, n) - p
  stats::uniroot(f, c(0, 1), tol = 1e-300, maxiter = 5000L)$root
}

# Laplace (location m, scale s) ------------------------------------------------------------

.ig_lapl_cdf <- function(x, m, s) {
  if (x < m) 0.5 * exp((x - m) / s) else 1 - 0.5 * exp(-(x - m) / s)
}
.ig_lapl_pdf <- function(x, m, s) exp(-abs(x - m) / s) / (2 * s)
.ig_lapl_inv <- function(p, m, s) {
  if (p < 0.5) m + s * log(2 * p) else m - s * log(2 * (1 - p))
}

# Binomial: continuous extension in x ---------------------------------------------------

# P(X <= x) extended to real x in [0, N) as I_{1-pi}(N - x, x + 1)
# (= pbinom for integer x); 0 below 0 and 1 from N on (G*Power:
# binocdf(-0.5,10,0.3) = 0, captured).
.ig_binom_cdf_continuous <- function(x, size, prob) {
  if (x < 0) return(0)
  if (x >= size) return(1)
  stats::pbeta(1 - prob, size - x, x + 1)
}

# Inverse of the continuous extension on [0, N]: the x with
# I_{1-pi}(N - x, x + 1) = p, 0 when p does not exceed P(X <= 0) = (1 - pi)^N
# and N when p = 1 (G*Power: binoinv(0,10,0.5) = 0, binoinv(1,10,0.5) = 10,
# binoinv(1e-9,10,0.5) = 0, captured).
.ig_binom_inv_continuous <- function(p, size, prob) {
  if (p >= 1 || prob >= 1) return(size)
  if (p <= 0 || prob <= 0 || size == 0) return(0)
  f <- if (p > 0.5) {
    function(x) (1 - p) - stats::pbeta(1 - prob, size - x, x + 1, lower.tail = FALSE)
  } else {
    function(x) stats::pbeta(1 - prob, size - x, x + 1) - p
  }
  if (f(0) >= 0) return(0)
  hi <- size * (1 - 1e-15)
  if (f(hi) <= 0) return(size)
  stats::uniroot(f, c(0, hi), tol = 1e-300, maxiter = 5000L)$root
}

# Noncentral t -----------------------------------------------------------------------------
#
# The calculator displays 15 significant digits. nmath's pnt (AS 243) has an
# absolute error target of 1e-12 and its dnt is a difference of two pnt
# values, which is visible in the displayed digits (for example
# pt(2, 10, 1.5) = 0.659154072442026, the integral below and G*Power both give
# 0.65915407244219). The calculator therefore evaluates noncentral t
# probabilities with the integral of R/num-dist.R (.ig_p_t_integral(), the
# chi-square CDF integrated over the normal numerator on the log scale) and
# the density with the analogous integral over the chi-square variable:
#   f(t) = int_0^Inf phi(t sqrt(v/df) - nc) sqrt(v/df) dchisq(v, df) dv.

.ig_calc_capture_warnings <- function(expr, what) {
  withCallingHandlers(expr, warning = function(w) {
    .ig_abort("numerical_failure", sprintf("%s: %s", what, conditionMessage(w)))
  })
}

.ig_calc_nct_cdf <- function(x, df, nc) {
  if (nc == 0) return(stats::pt(x, df))
  if (is.infinite(x)) return(as.double(x > 0))
  if (x == 0) return(stats::pnorm(-nc))
  .ig_calc_capture_warnings(exp(.ig_p_t_integral(x, df, nc, TRUE)), "noncentral t CDF")
}

.ig_calc_nct_pdf <- function(x, df, nc) {
  if (nc == 0) return(stats::dt(x, df))
  if (is.infinite(x)) return(0)
  # Integrate in the positive scale S=sqrt(V/df), rescaled to its mode.
  # The original V integral loses its narrow near-zero peak for low df and large x.
  scale <- max(abs(x),sqrt(df))
  h <- scale*sqrt((x/scale)^2+(sqrt(df)/scale)^2)
  b <- (x/h)*nc
  disc <- sqrt(b*b+4*df)
  u <- if(b>=0) (b+disc)/2 else 2*df/(disc-b)
  mode <- u/h
  if(!is.finite(mode)||mode<=0) .ig_abort("numerical_failure","noncentral t density: unresolved integration scale")
  constant <- log(2)+(df/2)*log(df/2)-lgamma(df/2)
  logf <- function(z) {
    out <- constant+df*(log(mode)+log(z)) + log(mode) +
      stats::dnorm((x*mode)*z-nc,log=TRUE) - (sqrt(df)*mode*z)^2/2
    out[z<=0] <- -Inf
    out
  }
  .ig_calc_capture_warnings(exp(.ig_integrate_log_peak(logf,0,Inf,1)),"noncentral t density")
}

# Quantile by inversion of a continuous, increasing lower-tail CDF. For
# p > 0.5 the upper tail is inverted instead when `upper` is given, so that
# quantiles near p = 1 are not limited by the absolute rounding of the lower
# tail near 1. The bracket is grown from `start`; uniroot() with tolerance 0
# stops at about 2 ulp of the root.
.ig_calc_invert <- function(cdf, p, start, support = c(-Inf, Inf), upper = NULL) {
  f <- if (!is.null(upper) && p > 0.5) function(x) (1 - p) - upper(x) else function(x) cdf(x) - p
  step <- max(1e-3, 1e-3 * abs(start))
  lo <- max(start - step, support[1])
  hi <- min(start + step, support[2])
  flo <- f(lo)
  fhi <- f(hi)
  for (i in seq_len(400L)) {
    if (flo <= 0 && fhi >= 0) break
    step <- step * 2
    if (flo > 0) {
      lo <- if (is.finite(support[1])) support[1] + (lo - support[1]) / 2 else lo - step
      flo <- f(lo)
    }
    if (fhi < 0) {
      hi <- hi + step
      fhi <- f(hi)
    }
  }
  if (!(flo <= 0 && fhi >= 0)) .ig_abort("numerical_failure", "the quantile could not be bracketed")
  if (flo == 0) return(lo)
  if (fhi == 0) return(hi)
  stats::uniroot(f, c(lo, hi), f.lower = flo, f.upper = fhi, tol = 1e-300, maxiter = 5000L)$root
}

.ig_calc_nct_upper <- function(x, df, nc) {
  if (nc == 0) return(stats::pt(x, df, lower.tail = FALSE))
  if (is.infinite(x)) return(as.double(x < 0))
  if (x == 0) return(stats::pnorm(nc))
  .ig_calc_capture_warnings(exp(.ig_p_t_integral(x, df, nc, FALSE)), "noncentral t CDF")
}

# Noncentral chi-square and F: Poisson mixtures on the log scale
# (.ig_poisson_mixture_log() of R/num-dist.R widens its window until the terms
# are negligible), here with a cap on the number of terms so that a far tail
# of a distribution with a very large noncentrality parameter ends in a typed
# failure instead of a long computation.
.ig_calc_mixture <- function(m, log_term) {
  lp <- .ig_poisson_mixture_log(m, log_term, max_terms = .ig_calc_range$mixture_terms)
  if (is.nan(lp)) {
    .ig_abort("numerical_failure", sprintf("the Poisson mixture needs more than %s terms (a far tail with a very large noncentrality parameter)",
                                           format(.ig_calc_range$mixture_terms)))
  }
  lp
}

# Chernoff bounds of the noncentral chi-square tails, on the log scale:
#   log P(X >= x) <= min_{0 < t < 1/2} K(t) - t x,   log P(X <= x) <= min_{t < 0} K(t) - t x,
# K(t) = nc t / (1 - 2t) - (df / 2) log(1 - 2t) (cumulant generating function).
# Used only to recognise tails below the double range (exp(-745)), where the
# mixture would need a window far from the Poisson mode.
.ig_calc_ncx2_log_bound <- function(x, df, nc, upper) {
  obj <- function(u) {
    t <- if (upper) 0.5 * u else -u / (1 - u)
    nc * t / (1 - 2 * t) - df / 2 * log1p(-2 * t) - t * x
  }
  if (upper && x <= df + nc || !upper && x >= df + nc) return(0)
  stats::optimize(obj, c(0, 1), tol = 1e-10)$objective
}

.ig_calc_ncx2_cdf <- function(x, df, nc, lower_tail = TRUE) {
  if (nc == 0) return(stats::pchisq(x, df, lower.tail = lower_tail))
  if (x <= 0 || is.infinite(x)) return(as.double((x > 0) == lower_tail))
  # A tail below the double range is 0 (and its complement 1 to the last bit).
  if (.ig_calc_ncx2_log_bound(x, df, nc, upper = FALSE) < -750) return(as.double(!lower_tail))
  if (.ig_calc_ncx2_log_bound(x, df, nc, upper = TRUE) < -750) return(as.double(lower_tail))
  exp(.ig_calc_mixture(nc / 2, function(j) stats::pchisq(x, df + 2 * j, lower.tail = lower_tail, log.p = TRUE)))
}

# Noncentral chi-square density. nmath's dnchisq is 2.6% off at
# dchisq(1500, 10, ncp = 1000) (1.33482e-13 against 1.37114e-13 from the
# 40-digit reference), so the mixture of central densities is summed here.
# Beyond the mode the density is decreasing, so f(x) <= P(x - 1 < X <= x);
# below the mode f(x) <= P(x < X <= x + 1): where the Chernoff bound of that
# probability is below the double range the density is 0.
.ig_calc_ncx2_pdf <- function(x, df, nc) {
  if (nc == 0) return(stats::dchisq(x, df))
  if (x < 0 || is.infinite(x)) return(0)
  if (x > df + nc + 1 && .ig_calc_ncx2_log_bound(x - 1, df, nc, upper = TRUE) < -750) return(0)
  if (df >= 2 && x < df + nc - 1 && .ig_calc_ncx2_log_bound(x + 1, df, nc, upper = FALSE) < -750) return(0)
  exp(.ig_calc_mixture(nc / 2, function(j) stats::dchisq(x, df + 2 * j, log = TRUE)))
}

.ig_calc_ncf_cdf <- function(x, df1, df2, nc, lower_tail = TRUE) {
  if (nc == 0) return(stats::pf(x, df1, df2, lower.tail = lower_tail))
  if (x <= 0 || is.infinite(x)) return(as.double((x > 0) == lower_tail))
  # Beta argument and its complement, both formed without cancellation.
  y <- df1 * x / (df1 * x + df2)
  y_c <- df2 / (df1 * x + df2)
  exp(.ig_calc_mixture(nc / 2, function(j) {
    if (lower_tail) stats::pbeta(y, df1 / 2 + j, df2 / 2, log.p = TRUE)
    else stats::pbeta(y_c, df2 / 2, df1 / 2 + j, log.p = TRUE)
  }))
}

# Log-normal density; the log form is used where the direct formula of nmath
# is 0/0 (an SD so small that x * s underflows).
.ig_logn_pdf <- function(x, m, s) {
  if (x <= 0 || is.infinite(x)) return(0)
  v <- stats::dlnorm(x, m, s)
  if (!is.nan(v)) return(v)
  y <- (log(x) - m) / s
  exp(-(0.5 * log(2 * pi) + 0.5 * y^2 + log(x) + log(s)))
}

.ig_calc_nct_inv <- function(p, df, nc) {
  if (nc == 0) return(if (p == 0.5) 0 else stats::qt(p, df))
  start <- suppressWarnings(stats::qt(p, df, nc))
  if (!is.finite(start)) start <- nc
  .ig_calc_invert(function(x) .ig_calc_nct_cdf(x, df, nc), p, start,
                  upper = function(x) .ig_calc_nct_upper(x, df, nc))
}

.ig_calc_ncx2_inv <- function(p, df, nc) {
  if (nc == 0) return(stats::qchisq(p, df))
  if (p == 0) return(0)
  start <- suppressWarnings(stats::qchisq(p, df, nc))
  if (!is.finite(start) || start <= 0) start <- df + nc
  .ig_calc_invert(function(x) .ig_calc_ncx2_cdf(x, df, nc), p, start, support = c(0, Inf),
                  upper = function(x) .ig_calc_ncx2_cdf(x, df, nc, lower_tail = FALSE))
}

.ig_calc_ncf_inv <- function(p, df1, df2, nc) {
  if (nc == 0) return(stats::qf(p, df1, df2))
  if (p == 0) return(0)
  start <- suppressWarnings(stats::qf(p, df1, df2, nc))
  if (!is.finite(start) || start <= 0) start <- (df1 + nc) / df1
  .ig_calc_invert(function(x) .ig_calc_ncf_cdf(x, df1, df2, nc), p, start, support = c(0, Inf),
                  upper = function(x) .ig_calc_ncf_cdf(x, df1, df2, nc, lower_tail = FALSE))
}
