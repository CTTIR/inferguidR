# Copyright (c) 2026 Raban Heller. MIT licence; see ../../../LICENSE.md.
# Pure mathematical references: own Gauss-Legendre integrals and information matrices.
ck_gl_nodes <- function(n) {
  i <- seq_len(n - 1L)
  b <- i / sqrt(4 * i^2 - 1)
  J <- matrix(0, n, n)
  J[cbind(i, i + 1L)] <- b
  J[cbind(i + 1L, i)] <- b
  e <- eigen(J, symmetric = TRUE)
  list(x = e$values, w = 2 * e$vectors[1, ]^2)
}
CK_GL <- ck_gl_nodes(64L)

## Composite Gauss-Legendre integral of a vectorised f over [a, b].
ck_int <- function(f, a, b, panels = 160L) {
  edges <- seq(a, b, length.out = panels + 1L)
  h <- diff(edges) / 2
  mid <- (edges[-1L] + edges[-length(edges)]) / 2
  x <- outer(CK_GL$x, h) + rep(mid, each = length(CK_GL$x))
  fx <- matrix(f(as.vector(x)), nrow = length(CK_GL$x))
  sum(colSums(fx * CK_GL$w) * h)
}

ck_phi <- function(z) exp(-z^2 / 2) / sqrt(2 * pi)

## Lower normal tail Phi(x), by integrating the density (relative accuracy
## near double precision for x <= 0; upper tails by symmetry).
ck_Phi <- function(x) {
  if (!is.finite(x)) return(if (x > 0) 1 else 0)
  if (x > 0) return(1 - ck_Phi(-x))
  if (x < -38) return(0)
  ck_int(ck_phi, max(x - 12, -40), x, panels = 24L)
}
ck_Q <- function(x) ck_Phi(-x) # upper tail

## z with upper tail probability a (bisection).
ck_zq <- function(a) {
  lo <- -40
  hi <- 40
  for (i in 1:200) {
    m <- (lo + hi) / 2
    if (ck_Q(m) > a) lo <- m else hi <- m
    if (hi - lo < 1e-15 * max(1, abs(m))) break
  }
  (lo + hi) / 2
}

ck_pois_pmf <- function(lambda, kmax) {
  if (!is.finite(kmax) || kmax < 0 || kmax > 1e6) stop("reference Poisson grid exceeds 1,000,001 terms")
  p <- numeric(kmax + 1L)
  p[1] <- exp(-lambda)
  for (k in seq_len(kmax)) p[k + 1L] <- p[k] * lambda / k
  p
}

## Expectations E[g(X)] for the covariate distributions (centred variable
## and wide fixed ranges; `tilt` shifts the range for exponential weights).
ck_expect <- function(xd, par, g, tilt = 0) {
  switch(xd,
    binomial = (1 - par$pi) * g(0) + par$pi * g(1),
    poisson = {
      kmax <- ceiling(par$lambda * exp(max(tilt, 0)) * 3 + 80)
      p <- ck_pois_pmf(par$lambda, kmax)
      k <- 0:kmax
      v <- p * g(k)
      sum(v[is.finite(v)])
    },
    normal = {
      s <- tilt * par$sigma
      ck_int(function(z) { v <- ck_phi(z) * g(par$mu + par$sigma * z); v[is.nan(v)] <- 0; v },
             min(-40, s - 40), max(40, s + 40), panels = 200L)
    },
    lognormal = ck_int(function(z) { v <- ck_phi(z) * g(exp(par$mu + par$sigma * z)); v[is.nan(v)] <- 0; v },
                       -40, 40 + 2 * par$sigma, panels = 300L),
    exponential = {
      r <- 1 - tilt / par$rate
      ck_int(function(e) { v <- exp(-e) * g(e / par$rate); v[is.nan(v)] <- 0; v }, 0, 800 / r, panels = 500L)
    },
    uniform = ck_int(function(x) g(x) / (par$b - par$a), par$a, par$b, panels = 40L)
  )
}

ck_mean <- function(xd, par) switch(xd, binomial = par$pi, poisson = par$lambda, normal = par$mu,
                                    lognormal = exp(par$mu + par$sigma^2 / 2), exponential = 1 / par$rate,
                                    uniform = (par$a + par$b) / 2)
ck_var <- function(xd, par) switch(xd, binomial = par$pi * (1 - par$pi), poisson = par$lambda,
                                   normal = par$sigma^2, lognormal = (exp(par$sigma^2) - 1) * exp(2 * par$mu + par$sigma^2),
                                   exponential = 1 / par$rate^2, uniform = (par$b - par$a)^2 / 12)

## v = [I^{-1}]_{11} for the information matrix of (1, X - c) with weight w.
ck_v <- function(xd, par, w, tilt = 0) {
  c0 <- ck_mean(xd, par)
  a00 <- ck_expect(xd, par, function(x) w(x), tilt)
  a01 <- ck_expect(xd, par, function(x) (x - c0) * w(x), tilt)
  a11 <- ck_expect(xd, par, function(x) (x - c0)^2 * w(x), tilt)
  a00 / (a00 * a11 - a01^2)
}

## ---------------------------------------------------------------------------
## reference power models
## ---------------------------------------------------------------------------

ck_zpower <- function(m1, s1, alpha, tails) {
  m <- abs(m1)
  if (tails == "one") {
    z <- ck_zq(alpha)
    list(power = ck_Phi((m - z) / s1), beta = ck_Phi((z - m) / s1), crit = z)
  } else {
    z <- ck_zq(alpha / 2)
    list(power = ck_Phi((m - z) / s1) + ck_Phi((-m - z) / s1),
         beta = ck_Phi((z - m) / s1) - ck_Phi((-z - m) / s1), crit = z)
  }
}

ck_logit_v <- function(xd, par, b0, b1) {
  w <- function(x) { e <- exp(-abs(b0 + b1 * x)); e / (1 + e)^2 }
  ck_v(xd, par, w)
}

ck_logistic_power <- function(cs) {
  b0 <- log(cs$p1 / (1 - cs$p1))
  b1 <- log(cs$or)
  p2 <- 1 / (1 + exp(-(b0 + b1)))
  N <- cs$n * (1 - cs$r2)
  if (cs$proc == "hsieh" && (cs$xd == "binomial" || (cs$xd == "normal" && cs$par$mu == 0 && cs$par$sigma == 1))) {
    z <- if (cs$tails == "one") ck_zq(cs$alpha) else ck_zq(cs$alpha / 2)
    if (cs$xd == "binomial") {
      B <- cs$par$pi
      pb <- (1 - B) * cs$p1 + B * p2
      arg <- (sqrt(N * (1 - B)) * abs(p2 - cs$p1) - z * sqrt(pb * (1 - pb) / B)) /
        sqrt(cs$p1 * (1 - cs$p1) + p2 * (1 - p2) * (1 - B) / B)
    } else {
      arg <- sqrt(N * cs$p1 * (1 - cs$p1)) * abs(b1) - z
    }
    return(list(power = ck_Phi(arg), beta = ck_Phi(-arg), crit = z))
  }
  vc <- if (cs$proc == "hsieh") "on" else cs$vc
  if (cs$xd == "manual") {
    v1 <- cs$par$v1
    v0 <- cs$par$v0
    a <- if (vc == "on") 1 else 0
  } else if (cs$xd == "binomial") {
    B <- cs$par$pi
    v1 <- 1 / ((1 - B) * cs$p1 * (1 - cs$p1)) + 1 / (B * p2 * (1 - p2))
    mu <- (1 - B) * cs$p1 + B * p2
    v0 <- 1 / (mu * (1 - mu) * B * (1 - B))
    a <- if (vc == "on") 0.85 else 0
  } else {
    v1 <- ck_logit_v(cs$xd, cs$par, b0, b1)
    mu <- ck_expect(cs$xd, cs$par, function(x) 1 / (1 + exp(-(b0 + b1 * x))))
    v0 <- 1 / (mu * (1 - mu) * ck_var(cs$xd, cs$par))
    a <- if (vc == "on") (if (cs$xd == "lognormal") 0.75 else 1) else 0
  }
  m1 <- sqrt(N / v1) * b1
  s1 <- sqrt((a * v0 + (1 - a) * v1) / v1)
  ck_zpower(m1, s1, cs$alpha, cs$tails)
}

## Poisson regression variances: closed forms where they exist, Gauss-Legendre
## otherwise (lognormal); `route = "quadrature"` forces the quadrature.
ck_pois_v <- function(xd, par, b0, b1, route = "closed") {
  # The antiderivative difference loses its leading terms near beta1=0.
  # Use the independent quadrature already provided for this distribution.
  if (xd == "uniform" && abs(b1) * (par$b - par$a) < 0.1) route <- "quadrature"
  if (route == "closed" && xd != "lognormal") {
    v <- switch(xd,
      binomial = 1 / (1 - par$pi) + 1 / (par$pi * exp(b1)),
      normal = exp(-(b1 * par$mu + (b1 * par$sigma)^2 / 2)) / par$sigma^2,
      exponential = (par$rate - b1)^3 / par$rate,
      poisson = exp(-b1 + par$lambda - exp(b1) * par$lambda) / par$lambda,
      uniform = {
        if (b1 == 0) 12 / (par$b - par$a)^2 else {
          A <- par$a; B <- par$b; D <- B - A
          F0 <- function(x) exp(b1 * x) / b1
          F1 <- function(x) exp(b1 * x) * (x / b1 - 1 / b1^2)
          F2 <- function(x) exp(b1 * x) * (x^2 / b1 - 2 * x / b1^2 + 2 / b1^3)
          i0 <- (F0(B) - F0(A)) / D; i1 <- (F1(B) - F1(A)) / D; i2 <- (F2(B) - F2(A)) / D
          i0 / (i0 * i2 - i1^2)
        }
      })
    if (xd == "binomial" && b1 == 0) v <- 1 / (par$pi * (1 - par$pi))
    return(v / exp(b0))
  }
  ck_v(xd, par, function(x) exp(b0 + b1 * x), tilt = b1)
}

ck_poisson_power <- function(cs, route = "closed") {
  b0 <- log(cs$base)
  b1 <- log(cs$eb1)
  N <- cs$n * (1 - cs$r2)
  t <- cs$t
  if (cs$xd == "manual") {
    v1 <- cs$par$v1; v0 <- cs$par$v0; v0s <- cs$par$v0
  } else {
    v1 <- ck_pois_v(cs$xd, cs$par, b0, b1, route)
    v0 <- ck_pois_v(cs$xd, cs$par, b0, 0, route)
    mu <- if (route == "closed" && cs$xd != "lognormal") {
      switch(cs$xd, binomial = exp(b0) * (1 - cs$par$pi + cs$par$pi * exp(b1)),
             normal = exp(b0 + b1 * cs$par$mu + (b1 * cs$par$sigma)^2 / 2),
             exponential = exp(b0) * cs$par$rate / (cs$par$rate - b1),
             poisson = exp(b0 + cs$par$lambda * (exp(b1) - 1)),
             uniform = {
               z <- b1 * (cs$par$b - cs$par$a)
               exp(b0 + b1 * cs$par$a) * (if (z == 0) 1 else expm1(z) / z)
             })
    } else ck_expect(cs$xd, cs$par, function(x) exp(b0 + b1 * x), tilt = b1)
    v0s <- 1 / (mu * ck_var(cs$xd, cs$par))
  }
  if (cs$proc == "signorini" && cs$vc == "off" && !(cs$xd %in% c("binomial", "manual"))) {
    # box unchecked: the null variance is replaced by v1 (method notes); no
    # effect for binomial and Manual X
    m1 <- b1 * sqrt(N * t / v1); s1 <- 1
  } else if (cs$proc == "signorini") {
    m1 <- b1 * sqrt(N * t / v0); s1 <- sqrt(v1 / v0)
  } else {
    a <- if (cs$vc == "on") (if (cs$xd == "lognormal") 0.75 else 1) else 0
    m1 <- b1 * sqrt(N * t / v1); s1 <- sqrt((a * v0s + (1 - a) * v1) / v1)
  }
  ck_zpower(m1, s1, cs$alpha, cs$tails)
}


# Enumeration references are capped and used only for small synthetic grids.
ck_grid <- function(xd, par, ng, family = "logistic") {
  q <- (seq_len(ng) - 3 / 8) / (ng + 1 / 4)
  upper <- family == "logistic"
  switch(xd,
    binomial = if (upper) ifelse(q >= 1 - par$pi, 1, 0) else ifelse(q > 1 - par$pi, 1, 0),
    poisson = {
      kmax <- ceiling(par$lambda * 4 + 80)
      Fk <- cumsum(ck_pois_pmf(par$lambda, kmax))
      vapply(q, function(qq) (if (upper) which(Fk > qq)[1] else which(Fk >= qq)[1]) - 1, numeric(1))
    },
    exponential = -log(1 - q) / par$rate,
    lognormal = exp(par$mu + par$sigma * stats::qnorm(q)),
    normal = par$mu + par$sigma * stats::qnorm(q),
    uniform = par$a + (par$b - par$a) * q)
}

## Poisson regression: exemplary data set written out row by row. Count range
## L..K: K is the first count at which the cumulative probability of the
## largest expected count reaches 1 - 1e-5, L the first count at which the
## cumulative probability of the smallest expected count reaches 1e-5.
## Returns the rows (x, y, w) aggregated over equal covariate values.
ck_pois_rows <- function(x, m, tail = 1e-5) {
  cum_first <- function(mean, level) {
    if (!is.finite(mean) || mean > 700) stop("reference count recursion cannot start reliably above mean 700")
    p <- exp(-mean); cum <- p; k <- 0
    while (cum < level) { k <- k + 1; p <- p * mean / k; cum <- cum + p; if (k > 1e6) stop("count limit not found") }
    k
  }
  K <- cum_first(max(m), 1 - tail)
  L <- cum_first(min(m), tail)
  ux <- sort(unique(x))
  if (length(ux) * (K - L + 1) > 2e6) stop("reference Poisson enumeration exceeds 2,000,000 rows")
  rows <- do.call(rbind, lapply(ux, function(u) {
    mi <- m[match(u, x)]
    p <- ck_pois_pmf(mi, K)
    data.frame(x = u, y = L:K, w = sum(x == u) * p[(L:K) + 1L])
  }))
  attr(rows, "range") <- c(L, K)
  rows
}

## Weighted Poisson regression of y on x fitted to rows by iteratively
## reweighted least squares; returns the coefficients, the information matrix
## at the fit and the log likelihood (without the factorial term).
ck_pois_irls <- function(rows, intercept_only = FALSE) {
  X <- if (intercept_only) matrix(1, nrow(rows), 1) else cbind(1, rows$x - mean(rows$x))
  y <- rows$y; w <- rows$w
  eta <- rep(log(sum(w * y) / sum(w)), nrow(rows))
  for (it in 1:200) {
    mu <- exp(eta)
    z <- eta + (y - mu) / mu
    ww <- sqrt(w * mu)
    b <- base::qr.solve(X * ww, z * ww, tol = 1e-13)
    eta_new <- drop(X %*% b)
    done <- max(abs(eta_new - eta)) < 1e-13
    eta <- eta_new
    if (done) break
  }
  mu <- exp(eta)
  list(b = b, info = crossprod(X * sqrt(w * mu)), loglik = sum(w * (y * eta - mu)), sum_wy = sum(w * y))
}

ck_enum_lambda <- function(cs, family) {
  ng <- floor(cs$n * (1 - cs$r2) + 0.5)
  if (ng < 2) return(0)
  x <- ck_grid(cs$xd, cs$par, ng, family)
  # A constant covariate contributes no information about its slope.
  if (length(unique(x)) < 2) return(0)
  if (family != "logistic") {
    b1 <- log(cs$eb1)
    if (length(unique(x)) < 2) return(0)
    rows <- ck_pois_rows(x, exp(log(cs$base) + b1 * x) * cs$t)
    if (cs$stat == "wald") return(b1^2 / base::solve(ck_pois_irls(rows)$info)[2, 2])
    # likelihood ratio: log likelihood of the rows at the hypothesised means
    # minus that of the fitted intercept-only model
    null <- ck_pois_irls(rows, intercept_only = TRUE)
    m_row <- exp(log(cs$base) + b1 * rows$x) * cs$t
    return(2 * (sum(rows$w * (rows$y * log(m_row) - m_row)) - null$loglik))
  }
  {
    b0 <- log(cs$p1 / (1 - cs$p1)); b1 <- log(cs$or)
    mu <- 1 / (1 + exp(-(b0 + b1 * x)))
    if (cs$stat == "wald") {
      w <- mu * (1 - mu)
    } else {
      nu <- 1 / (1 + exp(b0 + b1 * x))
      entropy_term <- function(p) ifelse(p > 0, p * log(p), 0)
      ll1 <- sum(entropy_term(mu) + entropy_term(nu))
      ll0 <- ng * (entropy_term(mean(mu)) + entropy_term(mean(nu)))
      return(max(2 * (ll1 - ll0), 0))
    }
  }
  xc <- x - mean(x)
  I <- matrix(c(sum(w), sum(w * xc), sum(w * xc), sum(w * xc^2)), 2)
  V <- base::solve(I)
  b1^2 / V[2, 2]
}

ck_enum_power <- function(lam, alpha, tails) {
  a <- if (tails == "one") 2 * alpha else alpha
  if (a >= 1) return(list(power = 1, beta = 0, crit = 0))
  z <- ck_zq(a / 2)
  r <- sqrt(lam)
  list(power = ck_Phi(r - z) + ck_Phi(-r - z), beta = ck_Phi(z - r) - ck_Phi(-z - r), crit = z^2, lambda = lam)
}
