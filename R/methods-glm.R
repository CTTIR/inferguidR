# Method definitions: logistic and Poisson regression (z tests) -----------------------
#
# Ownership: group G6 (blueprint phases 06-10 for z_logistic_regression and
# z_poisson_regression). Effect-size drawer: R/effect-glm.R. Method notes:
# specs/method-notes/z_logistic_regression.md and
# specs/method-notes/z_poisson_regression.md (the mathematical contract,
# G*Power's conventions as found in the emulator captures, and every manual
# or reference defect met).
#
# Both procedures test H0: beta1 = 0 for one covariate X1 of a generalized
# linear model (logit link with a binary response, log link with a count
# response and exposure t). The covariate X has one of six predefined
# distributions, or its variances under H0 and H1 are given directly
# ("Manual"). G*Power offers two kinds of procedure:
#
# * Large-sample approximations of the Wald test (manual 29.6.2, 30.5.2):
#   the statistic is N(0, 1) under H0 and N(m1, s1) under H1.
#     logistic, Demidenko (2007):  m1 = sqrt(N (1 - R2) / v1) beta1,
#                                  s1 = sqrt((a v0 + (1 - a) v1) / v1)
#     logistic, Hsieh et al. (1998): their sample-size equations (1) (standard
#                                  normal X) and (2) (binary X) solved for
#                                  power; any other X switches to Demidenko
#     Poisson, Signorini (1991):   m1 = beta1 sqrt(N (1 - R2) t / v0),
#                                  s1 = sqrt(v1 / v0)
#     Poisson, Demidenko:          m1 = beta1 sqrt(N (1 - R2) t / v1),
#                                  s1 = sqrt((a v0* + (1 - a) v1) / v1)
#   v = I00 / (I00 I11 - I01^2) is the per-observation variance of beta1-hat
#   from the information integrals I_jk = E[X^(j+k) w(X)], w = p (1 - p)
#   (logistic) or exp(beta0 + beta1 X) (Poisson). v1 is taken at the H1
#   parameters; v0 at beta1 = 0 with the intercept matched to the marginal
#   mean (Demidenko) or kept at beta0 (Signorini). Variance correction:
#   a = 0.75 (lognormal X), 0.85 (binomial X, logistic only), 1 otherwise;
#   without correction a = 0 (s1 = 1).
# * The enumeration procedure of Lyles et al. (2007): an exemplary data set
#   of Ng = N (1 - R2) covariate values at the quantiles of the X
#   distribution, p_i = (i - 3/8) / (Ng + 1/4) (Blom positions; identified
#   from the emulator captures, which it reproduces to all 8 displayed
#   digits for logistic regression), carries the expected responses. The
#   Wald noncentrality is Ng beta1^2 / v_grid and the likelihood-ratio
#   noncentrality the expected deviance difference between the H1 model and
#   the fitted intercept-only model; the test statistic is chi-square(1).
#   G*Power reports "Critical chi^2" at level 2 alpha for one tail and alpha
#   for two tails, and the power counts both tails of the noncentral
#   chi-square distribution.
#
# Integrals: stats::integrate (QUADPACK) on the standardised covariate with
# the integration range split around the peak of the integrand; discrete X
# by summation. The independent checker (validation/independent/g6/) uses
# closed-form moment generating functions and its own Gauss-Legendre rule.

.ig_define_methods_z_g6 <- function() {
  list(.ig_method_z_logistic_regression(), .ig_method_z_poisson_regression())
}

# Shared parameters and options ------------------------------------------------------------

.ig_glm_x_codes <- c("binomial", "exponential", "lognormal", "normal", "poisson", "uniform", "manual")

.ig_glm_x_params <- function() {
  when <- function(codes) list(option = "x_distribution", "in" = codes)
  list(
    ig_param("x_pi", role = "nuisance",
             domain = ig_domain_real(0, 1, lower_inclusive = FALSE, upper_inclusive = FALSE),
             label = list(en = "X parm \u03c0"),
             description = "Binomial (Bernoulli) covariate: Pr(X = 1).",
             applies_when = when("binomial")),
    ig_param("x_lambda", role = "nuisance",
             domain = ig_domain_real(0, Inf, lower_inclusive = FALSE, upper_inclusive = FALSE),
             label = list(en = "X parm \u03bb"),
             description = paste("Exponential covariate: rate (mean 1 / lambda); Poisson covariate: mean.",
                                 "The rate reading of the exponential parameter is confirmed by the captures."),
             applies_when = when(c("exponential", "poisson"))),
    ig_param("x_mu", role = "nuisance", domain = ig_domain_real(), default = 0,
             label = list(en = "X parm \u03bc"),
             description = paste("Normal covariate: mean; lognormal covariate: mean of log X. Default 0 is",
                                 "the G*Power field default (live capture 2026-10-05)."),
             applies_when = when(c("lognormal", "normal"))),
    ig_param("x_sigma", role = "nuisance",
             domain = ig_domain_real(0, Inf, lower_inclusive = FALSE, upper_inclusive = FALSE), default = 1,
             label = list(en = "X parm \u03c3"),
             description = paste("Normal covariate: SD; lognormal covariate: SD of log X. Default 1 is the",
                                 "G*Power field default (live capture 2026-10-05)."),
             applies_when = when(c("lognormal", "normal"))),
    ig_param("x_left", role = "nuisance", domain = ig_domain_real(),
             label = list(en = "X left border"), description = "Uniform covariate: lower border a.",
             applies_when = when("uniform")),
    ig_param("x_right", role = "nuisance", domain = ig_domain_real(),
             label = list(en = "X right border"), description = "Uniform covariate: upper border b > a.",
             applies_when = when("uniform")),
    ig_param("x_var_h0", role = "nuisance",
             domain = ig_domain_real(0, Inf, lower_inclusive = FALSE, upper_inclusive = FALSE),
             label = list(en = "X Var(\u03b2|H0)"),
             description = "Manual mode: per-observation variance v0 of beta1-hat under H0.",
             modes = c("a_priori", "post_hoc", "criterion", "compromise"),
             applies_when = when("manual")),
    ig_param("x_var_h1", role = "nuisance",
             domain = ig_domain_real(0, Inf, lower_inclusive = FALSE, upper_inclusive = FALSE),
             label = list(en = "X Var(\u03b2|H1)"),
             description = paste("Manual mode: per-observation variance v1 of beta1-hat under H1. Not",
                                 "available in sensitivity mode (v1 depends on the effect being solved;",
                                 "manual 30.5.2)."),
             modes = c("a_priori", "post_hoc", "criterion", "compromise"),
             applies_when = when("manual"))
  )
}

.ig_glm_common_params <- function() {
  list(
    ig_param_alpha(),
    ig_param_power(),
    ig_param_error_ratio(),
    ig_param("n_total", role = "design", domain = ig_domain_integer(lower = 2), unit = "subjects",
             label = list(en = "Total sample size", de = "Gesamtstichprobenumfang"),
             description = "Total sample size N (G*Power: 'Sample size must be greater than or equal to 2')."),
    ig_param("r2_other_x", role = "nuisance", domain = ig_domain_real(0, 1), default = 0,
             label = list(en = "R\u00b2 other X"),
             description = paste("Squared multiple correlation of X1 with the other covariates; the",
                                 "effective sample size is N (1 - R2) (Hsieh et al. 1998). Default 0 is the",
                                 "G*Power field default (live capture 2026-10-05).")),
    ig_param("enumeration_threshold", role = "control", domain = ig_domain_integer(lower = 1),
             default = 10000L, unit = "subjects",
             label = list(en = "Use enumeration procedure if N <"),
             description = paste("Options dialog: the enumeration procedure is used for N below this value;",
                                 "for larger N the selected (default) large-sample approximation is used",
                                 "(manual 29.2). Default 10000 is the G*Power field default (live capture)."),
             applies_when = list(option = "computation", "in" = "enumeration"))
  )
}

.ig_glm_common_options <- function(procedures, procedure_labels, vc_procedures,
                                   direction_codes, direction_labels, direction_default) {
  list(
    ig_option("tails", values = c("one", "two"), default = "one", labels = c(one = "One", two = "Two"),
              label = list(en = "Tail(s)", de = "Seiten")),
    ig_option("computation", values = c("enumeration", "large_sample"), default = "large_sample",
              labels = c(enumeration = "Use enumeration procedure (Lyles et al., 2007) if N < x",
                         large_sample = "Use large sample approximation"),
              label = list(en = "Computation")),
    ig_option("test_statistic", values = c("wald", "likelihood_ratio"), default = "wald",
              labels = c(wald = "Wald", likelihood_ratio = "Likelihood ratio"),
              applies_when = list(option = "computation", "in" = "enumeration"),
              label = list(en = "Test statistic")),
    ig_option("test_procedure", values = procedures, default = "demidenko", labels = procedure_labels,
              applies_when = list(option = "computation", "in" = "large_sample"),
              label = list(en = "Test procedure")),
    ig_option("variance_correction", values = c("on", "off"), default = "on",
              labels = c(on = "checked", off = "unchecked"),
              applies_when = list(option = "test_procedure", "in" = vc_procedures),
              label = list(en = "With variance correction")),
    ig_option("x_distribution", values = .ig_glm_x_codes, default = "normal",
              labels = c(binomial = "Binomial", exponential = "Exponential", lognormal = "Lognormal",
                         normal = "Normal", poisson = "Poisson", uniform = "Uniform", manual = "Manual"),
              label = list(en = "X distribution")),
    ig_option("effect_direction", values = direction_codes, default = direction_default,
              labels = direction_labels, applies_when = list(mode = "sensitivity"),
              label = list(en = "Effect direction (sensitivity)"))
  )
}

# Cross-parameter checks shared by both methods.
.ig_glm_constraints <- function(p) {
  out <- character()
  if (!is.null(p$x_left) && !is.null(p$x_right) && !(p$x_left < p$x_right)) {
    out <- c(out, "the uniform covariate needs X left border < X right border")
  }
  if (identical(p$computation, "enumeration") && identical(p$x_distribution, "manual")) {
    out <- c(out, paste("the Manual X distribution is available only with the large-sample procedures",
                        "(manual 29.1, 30.1; G*Power does not offer it with the enumeration procedure)"))
  }
  if (!is.null(p$alpha) && !is.null(p$power) && p$power < p$alpha) {
    out <- c(out, "the target power must not be below alpha (G*Power: 'Power must be greater than or equal to alpha.')")
  }
  out
}

# Covariate distributions --------------------------------------------------------------------

.ig_glm_xspec <- function(p) {
  switch(p$x_distribution,
         binomial = list(type = "binomial", pi = p$x_pi),
         exponential = list(type = "exponential", rate = p$x_lambda),
         lognormal = list(type = "lognormal", mu = p$x_mu, sigma = p$x_sigma),
         normal = list(type = "normal", mu = p$x_mu, sigma = p$x_sigma),
         poisson = list(type = "poisson", lambda = p$x_lambda),
         uniform = list(type = "uniform", a = p$x_left, b = p$x_right),
         manual = {
           if (is.null(p$x_var_h0) || is.null(p$x_var_h1)) {
             .ig_abort("invalid_input", paste(
               "the Manual X distribution is not available in sensitivity mode: the variance of beta1-hat under H1",
               "depends on the effect size being solved (manual 30.5.2; G*Power does not offer it either)"))
           }
           list(type = "manual", v0 = p$x_var_h0, v1 = p$x_var_h1)
         },
         stop("unknown X distribution", call. = FALSE))
}

.ig_glm_x_mean <- function(xs) {
  switch(xs$type,
         binomial = xs$pi,
         exponential = 1 / xs$rate,
         lognormal = exp(xs$mu + xs$sigma^2 / 2),
         normal = xs$mu,
         poisson = xs$lambda,
         uniform = (xs$a + xs$b) / 2)
}

.ig_glm_x_var <- function(xs) {
  switch(xs$type,
         binomial = xs$pi * (1 - xs$pi),
         exponential = 1 / xs$rate^2,
         lognormal = expm1(xs$sigma^2) * exp(2 * xs$mu + xs$sigma^2),
         normal = xs$sigma^2,
         poisson = xs$lambda,
         uniform = (xs$b - xs$a)^2 / 12)
}

# Quantile function used for the enumeration grid. A Blom position can fall
# exactly on a jump of the distribution function of a discrete covariate (in
# practice: a binomial covariate with pi = 0.5 and an odd grid, whose middle
# position is exactly 0.5). G*Power's two modules resolve the tie differently,
# as two dedicated captures of the smallest odd design show:
#   logistic regression: the upper value, x = min{k : F(k) > q}
#     (`tie = "upper"`; SUP-G6-Z-LOGIT-POS-E02: N = 3 has grid 0, 1, 1; the
#     campaign case GPW-Z-LOGIT-APR-B001-O01 has 308 ones among 615 values);
#   Poisson regression: the lower value, x = min{k : F(k) >= q}
#     (`tie = "lower"`; SUP-G6-Z-POIS-POS-Q25: N = 3 has grid 0, 0, 1;
#     GPW-Z-POIS-APR-B001-O01 has 327 ones among 655 values).
.ig_glm_x_quantile <- function(xs, q, tie = c("upper", "lower")) {
  tie <- match.arg(tie)
  upper <- identical(tie, "upper")
  switch(xs$type,
         binomial = as.double(if (upper) q >= 1 - xs$pi else q > 1 - xs$pi),
         exponential = stats::qexp(q, xs$rate),
         lognormal = stats::qlnorm(q, xs$mu, xs$sigma),
         normal = stats::qnorm(q, xs$mu, xs$sigma),
         poisson = {
           k <- stats::qpois(q, xs$lambda)
           # qpois returns the smallest k with F(k) >= q up to a small fuzz;
           # make the rule exact, then move to the strict rule if requested.
           lo <- pmax(k - 1, 0)
           k <- ifelse(k > 0 & stats::ppois(lo, xs$lambda) >= q, lo, k)
           k <- k + (stats::ppois(k, xs$lambda) < q)
           if (upper) k + (stats::ppois(k, xs$lambda) <= q) else k
         },
         uniform = stats::qunif(q, xs$a, xs$b))
}

# One-dimensional integral of a vectorised integrand over a range split at
# the given points (each piece by QUADPACK). Returns NA when QUADPACK fails.
.ig_glm_integrate <- function(f, breaks) {
  breaks <- sort(unique(breaks[is.finite(breaks)]))
  total <- 0
  for (k in seq_len(length(breaks) - 1L)) {
    a <- breaks[k]
    b <- breaks[k + 1L]
    if (b <= a) next
    r <- tryCatch(stats::integrate(f, a, b, rel.tol = 1e-12, abs.tol = 0, subdivisions = 2000L,
                                   stop.on.error = FALSE),
                  error = function(e) NULL)
    if (is.null(r) || !is.finite(r$value) ||
        !(r$message %in% c("OK", "roundoff error was detected", "roundoff error is detected in the extrapolation table"))) {
      r <- tryCatch(stats::integrate(f, a, b, rel.tol = 1e-10, abs.tol = 0, subdivisions = 5000L),
                    error = function(e) NULL)
      if (is.null(r) || !is.finite(r$value)) return(NA_real_)
    }
    total <- total + r$value
  }
  total
}

# E[g(X) w(X)] for the covariate distribution. `log_w(x)` is the log of the
# weight function and `g` multiplies it (for example (x - c)^k). The
# integrand is formed on the log scale (log density + log weight) so that an
# exponential weight cannot overflow before it is damped by the density.
# `tilt` is the coefficient of an exponential weight exp(tilt * x) (0 for
# bounded weights) and `peak` an optional list(x, scale) locating a narrow
# feature of the weight (the logistic weight p (1 - p) peaks where the linear
# predictor is 0, with width 1 / |beta1|); both only place the breakpoints of
# the integration range.
.ig_glm_expect <- function(xs, log_w, g, tilt = 0, peak = NULL) {
  val <- function(x, logd) {
    v <- g(x) * exp(logd + log_w(x))
    v[is.nan(v)] <- 0
    v
  }
  # Breakpoints around a narrow peak, mapped to the integration variable.
  peak_breaks <- function(to_u) {
    if (is.null(peak) || !is.finite(peak$x) || !is.finite(peak$scale)) return(numeric())
    xb <- peak$x + peak$scale * c(-20, -5, -1, 0, 1, 5, 20)
    u <- suppressWarnings(to_u(xb))
    u[is.finite(u)]
  }
  switch(xs$type,
    binomial = sum(val(c(0, 1), log(c(1 - xs$pi, xs$pi)))),
    poisson = {
      lam_t <- xs$lambda * exp(max(tilt, 0))
      kmax <- ceiling(lam_t + 40 * sqrt(lam_t) + 60)
      # Bound the explicitly enumerated integration grid before allocation.
      # A missing integral becomes the existing typed numerical failure;
      # do not truncate the grid and return an inaccurate finite value.
      if (!is.finite(kmax) || kmax > 1e6) return(NA_real_)
      k <- 0:kmax
      sum(val(k, stats::dpois(k, xs$lambda, log = TRUE)))
    },
    normal = {
      s <- tilt * xs$sigma
      h <- function(u) val(xs$mu + xs$sigma * u, stats::dnorm(u, log = TRUE))
      pb <- peak_breaks(function(x) (x - xs$mu) / xs$sigma)
      lo <- min(-39, s - 39)
      hi <- max(39, s + 39)
      .ig_glm_integrate(h, c(lo, s - 10, s - 4, s, s + 4, s + 10, hi, pb[pb > lo & pb < hi]))
    },
    lognormal = {
      sg <- xs$sigma
      h <- function(u) val(exp(xs$mu + sg * u), stats::dnorm(u, log = TRUE))
      pb <- peak_breaks(function(x) (log(x) - xs$mu) / sg)
      hi <- 2 * sg + 39
      .ig_glm_integrate(h, c(-39, -10, -4, 0, sg, 2 * sg, 2 * sg + 4, 2 * sg + 10, hi, pb[pb > -39 & pb < hi]))
    },
    exponential = {
      r <- 1 - tilt / xs$rate # rate of the weighted density in units of the covariate's rate
      if (!(r > 0)) return(Inf)
      h <- function(u) val(u / xs$rate, -u)
      pb <- peak_breaks(function(x) x * xs$rate)
      hi <- 800 / r
      .ig_glm_integrate(h, c(0, 1 / r, 4 / r, 15 / r, 60 / r, 200 / r, hi, pb[pb > 0 & pb < hi]))
    },
    uniform = {
      h <- function(u) val(u, rep(-log(xs$b - xs$a), length(u)))
      pb <- peak_breaks(function(x) x)
      .ig_glm_integrate(h, c(seq(xs$a, xs$b, length.out = 9L), pb[pb > xs$a & pb < xs$b]))
    },
    stop("unknown X distribution", call. = FALSE))
}

# Cache of the N- and alpha-independent part of the large-sample models
# (variances under H0 and H1), keyed by the inputs that determine it. The
# a priori, criterion and compromise searches evaluate the same effect many
# times; the cached values are those a fresh evaluation returns.
.ig_glm_cache <- new.env(parent = emptyenv())

.ig_glm_cached <- function(key, compute) {
  hit <- .ig_glm_cache[[key]]
  if (!is.null(hit)) return(hit)
  val <- compute()
  if (length(ls(.ig_glm_cache, all.names = TRUE)) >= 500L) rm(list = ls(.ig_glm_cache, all.names = TRUE), envir = .ig_glm_cache)
  assign(key, val, envir = .ig_glm_cache)
  val
}

.ig_glm_key <- function(...) {
  x <- list(...)
  paste(vapply(x, function(v) if (is.numeric(v)) paste(sprintf("%.17g", v), collapse = ",") else
    paste(as.character(v), collapse = ","), ""), collapse = "|")
}

# Per-observation variance of beta1-hat, v = I00 / (I00 I11 - I01^2), from the
# weighted moments, computed as 1 / (I00 * weighted variance of X) with a
# two-pass centring so that a covariate with a large mean relative to its
# spread loses no precision. Returns list(v, i00) or NA values.
.ig_glm_v <- function(xs, log_w, tilt = 0, peak = NULL) {
  one <- function(x) rep(1, length(x))
  i00 <- .ig_glm_expect(xs, log_w, one, tilt, peak)
  if (!is.finite(i00) || i00 <= 0) return(list(v = NA_real_, i00 = i00))
  c0 <- .ig_glm_x_mean(xs)
  j1 <- .ig_glm_expect(xs, log_w, function(x) x - c0, tilt, peak)
  mw <- c0 + j1 / i00
  j2 <- .ig_glm_expect(xs, log_w, function(x) (x - mw)^2, tilt, peak)
  j1b <- .ig_glm_expect(xs, log_w, function(x) x - mw, tilt, peak)
  varw <- j2 / i00 - (j1b / i00)^2
  if (!is.finite(varw) || varw <= 0) return(list(v = NA_real_, i00 = i00))
  list(v = 1 / (i00 * varw), i00 = i00)
}

# Logistic regression: v1 (H1), the marginal event probability mu and v0
# (H0 with beta0* = logit(mu)) for a covariate distribution.
.ig_logistic_variances <- function(xs, b0, b1) {
  .ig_glm_cached(.ig_glm_key("logit", xs$type, unlist(xs[-1]), b0, b1), function() {
    lw1 <- function(x) {
      eta <- b0 + b1 * x
      stats::plogis(eta, log.p = TRUE) + stats::plogis(-eta, log.p = TRUE)
    }
    peak <- if (b1 != 0) list(x = -b0 / b1, scale = 1 / abs(b1)) else NULL
    v1 <- .ig_glm_v(xs, lw1, peak = peak)$v
    mu <- .ig_glm_expect(xs, function(x) stats::plogis(b0 + b1 * x, log.p = TRUE),
                         function(x) rep(1, length(x)), peak = peak)
    list(v1 = v1, mu = mu, v0 = 1 / (mu * (1 - mu) * .ig_glm_x_var(xs)))
  })
}

# Poisson regression: v1 (H1), v0 (beta1 = 0, beta0; Signorini), the marginal
# rate mu = E[exp(beta0 + beta1 X)] and v0* (beta1 = 0, beta0* = log mu;
# Demidenko).
.ig_poisson_variances <- function(xs, b0, b1) {
  .ig_glm_cached(.ig_glm_key("pois", xs$type, unlist(xs[-1]), b0, b1), function() {
    r1 <- .ig_glm_v(xs, function(x) b0 + b1 * x, tilt = b1)
    vx <- .ig_glm_x_var(xs)
    list(v1 = r1$v, mu = r1$i00, v0 = 1 / (exp(b0) * vx), v0_star = 1 / (r1$i00 * vx))
  })
}

# Normal-statistic power: H0 N(0, 1), H1 N(m1, s1). One tail rejects on the
# side of m1 (the upper side for m1 = 0); two tails count both rejection
# regions. Returns power, beta (computed directly) and critical values.
.ig_glm_z_power <- function(m1, s1, alpha, tails) {
  m <- abs(m1)
  side <- if (m1 < 0) -1 else 1
  if (identical(tails, "one")) {
    z <- stats::qnorm(alpha, lower.tail = FALSE)
    power <- stats::pnorm((m - z) / s1)
    beta <- stats::pnorm((z - m) / s1)
    critical <- side * z
  } else {
    z <- stats::qnorm(alpha / 2, lower.tail = FALSE)
    power <- stats::pnorm((m - z) / s1) + stats::pnorm((-m - z) / s1)
    beta <- stats::pnorm((z - m) / s1) - stats::pnorm((-z - m) / s1)
    critical <- c(-z, z)
  }
  list(power = power, beta = beta, critical = critical, side = side, z = z)
}

# Chi-square(1) power of the enumeration procedure: critical value at level
# 2 alpha (one tail) or alpha (two tails); both tails of the noncentral
# chi-square count.
.ig_glm_chisq_power <- function(lambda, alpha, tails) {
  a_eff <- if (identical(tails, "one")) 2 * alpha else alpha
  if (a_eff >= 1) {
    return(list(power = 1, beta = 0, critical = 0, a_eff = a_eff))
  }
  z <- stats::qnorm(a_eff / 2, lower.tail = FALSE)
  r <- sqrt(max(lambda, 0))
  power <- stats::pnorm(r - z) + stats::pnorm(-r - z)
  beta <- stats::pnorm(z - r) - stats::pnorm(-z - r)
  list(power = power, beta = beta, critical = z^2, a_eff = a_eff)
}

# Size of the enumeration grid: N (1 - R2) rounded half up, evaluated as
# floor(N * (1 - R2) + 0.5) in double precision. The captures need both
# details: N = 30925, R2 = 0.5 gives 15463 points (GPW-Z-LOGIT-APR-B066-D01:
# lambda 29.3408463; R's round() would give 15462), and N = 705, R2 = 0.9
# gives 705 * 0.09999999999999998 = 70.49999999999999, i.e. 70 points
# (GPW-Z-LOGIT-APR-B026-E01: lambda 6.11312); 20 * (1 - 0.9) gives 2.
.ig_glm_grid_size <- function(n, r2) {
  floor(as.double(n) * (1 - r2) + 0.5)
}

.ig_glm_grid <- function(xs, ng, tie = "upper") {
  q <- (seq_len(ng) - 0.375) / (ng + 0.25)
  .ig_glm_x_quantile(xs, q, tie)
}

.ig_glm_side_note <- function(tails) {
  if (identical(tails, "two")) {
    "two tails: alpha/2 in each tail; the power counts both rejection regions"
  } else {
    "one tail: the rejection region lies on the side of the sign of beta1 (the upper side for beta1 = 0)"
  }
}

# Effect direction restriction of the sensitivity mode: effects on the other
# side of the null are outside the requested analysis (power undefined there,
# reported as NA with a reason so that the effect-size scan skips them).
.ig_glm_outside_direction <- function(p, beta1) {
  d <- p$effect_direction
  if (is.null(d)) return(FALSE)
  le <- d %in% c("p2_le_p1", "beta1_le_0")
  (le && beta1 > 0) || (!le && beta1 < 0)
}

# Sensitivity search (G*Power convention) -----------------------------------------------------
#
# G*Power solves the sensitivity mode of these procedures within the effect
# direction chosen in the main window ("p2 <= p1" / "p2 >= p1",
# "beta1 <= 0" / "beta1 >= 0"). The power need not be monotone in the
# effect (binary or skewed covariates: the power rises and falls again away
# from the null, manual 29.3 and Fig. 36), so the convention matters:
#
# * the reported effect is the crossing of the target power closest to the
#   null in the selected direction (the captures, e.g.
#   GPW-Z-LOGIT-SEN-B012-O01, binary X, odds ratio 0.4832796, where a second
#   crossing far from the null exists);
# * when the target is not reached in that direction, the effect that
#   maximises the power there is returned with its actual power (manual 29.3:
#   "From version 3.1.8 on G*Power returns in these cases the effect size
#   which maximizes the power in the selected direction"). The result then
#   carries a warning, diagnostics$target_reached = FALSE and the maximum
#   power; it is not a solution of power = target.
#
# The search runs on a fixed grid of 201 points from the null outward
# (logit spacing of the distance to the boundary of the effect domain: p2 in
# (0, p1) or (p1, 1) for logistic regression; |beta1| in (0, 50), or below the
# rate of an exponential covariate, for Poisson regression), refines the first
# sign change with uniroot() and a maximum with optimize(). It returns the
# shape required of an ig_method sensitivity hook:
# list(status = "ok", value, diagnostics, warnings) or list(status, reason,
# diagnostics). It is registered as `sensitivity_search` when the core
# ig_method() accepts that argument (change request CR-G6-01); until then
# the generic scan of R/num-modes.R is used, restricted to the selected
# direction by the forward models (effects on the other side give NA power).
.ig_glm_sensitivity_search <- function(power_model) {
  force(power_model)
  function(p, target_name) {
    target <- p$power
    le <- p$effect_direction %in% c("p2_le_p1", "beta1_le_0")
    d <- if (le) -1 else 1
    logistic <- identical(target_name, "odds_ratio") || identical(target_name, "p2")
    if (logistic) {
      b0 <- stats::qlogis(p$p1)
      bound <- if (le) 0 else 1
      # Effect value at distance t in (0, 1) from the null towards the bound.
      at_t <- function(t) {
        p2 <- p$p1 + t * (bound - p$p1)
        if (identical(target_name, "p2")) p2 else exp(stats::qlogis(p2) - b0)
      }
      domain_text <- sprintf("p2 in %s", if (le) sprintf("(0, %s)", .ig_format_double(p$p1)) else
        sprintf("(%s, 1)", .ig_format_double(p$p1)))
    } else {
      umax <- 50
      if (identical(p$x_distribution, "exponential") && !le && !identical(p$computation, "enumeration")) {
        umax <- min(umax, p$x_lambda)
      }
      if (identical(p$x_distribution, "lognormal") && !le && !identical(p$computation, "enumeration")) {
        return(list(status = "invalid_input", reason = paste(
          "for a lognormal covariate every Exp(\u03b21) > 1 has infinite information integrals (no moment",
          "generating function for positive arguments), so the sensitivity analysis in the direction",
          "beta1 >= 0 has no defined solution"), diagnostics = list(search = "glm_sensitivity_first_crossing")))
      }
      at_t <- function(t) exp(d * t * umax)
      domain_text <- sprintf("beta1 in %s", if (le) sprintf("(-%s, 0)", .ig_format_double(umax)) else
        sprintf("(0, %s)", .ig_format_double(umax)))
    }
    pw <- function(t) {
      q <- p
      q[[target_name]] <- at_t(t)
      v <- power_model(q)$power
      if (is.numeric(v) && length(v) == 1L && is.finite(v)) v else NA_real_
    }
    s_grid <- seq(-25, 25, by = 0.25)
    t_grid <- stats::plogis(s_grid)
    vals <- vapply(t_grid, pw, numeric(1))
    ok <- is.finite(vals)
    diag <- list(search = "glm_sensitivity_first_crossing", domain = domain_text,
                 direction = p$effect_direction, grid_points = length(t_grid),
                 failed_evaluations = sum(!ok), target_power = target)
    if (!any(ok)) {
      return(list(status = "numerical_failure", reason = paste(
        "the power could not be evaluated anywhere in the selected effect direction"), diagnostics = diag))
    }
    f <- function(s) pw(stats::plogis(s)) - target
    above <- which(ok & vals >= target)
    if (length(above)) {
      k <- above[1]
      prev <- if (k > 1L) max(which(ok[seq_len(k - 1L)]), -Inf) else -Inf
      if (!is.finite(prev)) {
        # Reached at the first evaluable point next to the null.
        return(list(status = "ok", value = at_t(t_grid[k]), diagnostics = c(diag, list(
          target_reached = TRUE, root_kind = "grid_point", residual = vals[k] - target)), warnings = character()))
      }
      r <- stats::uniroot(f, c(s_grid[prev], s_grid[k]), f.lower = vals[prev] - target, f.upper = vals[k] - target,
                          tol = 1e-13, maxiter = 1000L)
      value <- at_t(stats::plogis(r$root))
      other <- which(ok[-length(ok)] & ok[-1] & (vals[-length(vals)] - target) * (vals[-1] - target) < 0)
      diag <- c(diag, list(target_reached = TRUE, root_kind = "crossing", iterations = r$iter,
                           bracket = c(at_t(t_grid[prev]), at_t(t_grid[k])),
                           sign_changes_on_grid = length(other),
                           rule = "the crossing of the target power closest to the null in the selected direction"))
      return(list(status = "ok", value = value, diagnostics = diag, warnings = if (length(other) > 1L)
        sprintf("%d crossings of the target power in the selected direction; the one closest to the null is reported",
                length(other)) else character()))
    }
    # Target not reached: the power-maximising effect in the selected direction.
    k <- which.max(replace(vals, !ok, -Inf))
    lo <- s_grid[max(1L, k - 1L)]
    hi <- s_grid[min(length(s_grid), k + 1L)]
    opt <- stats::optimize(function(s) { v <- pw(stats::plogis(s)); if (is.na(v)) -Inf else v },
                           c(lo, hi), maximum = TRUE, tol = 1e-10)
    s_best <- if (is.finite(opt$objective) && opt$objective >= vals[k]) opt$maximum else s_grid[k]
    value <- at_t(stats::plogis(s_best))
    best <- max(opt$objective, vals[k])
    diag <- c(diag, list(target_reached = FALSE, max_power = best,
                         rule = paste("target not reached in the selected direction: the effect maximising the power",
                                      "there is reported (G*Power 3.1.8+ convention, manual 29.3)")))
    list(status = "ok", value = value, diagnostics = diag, warnings = sprintf(paste(
      "the target power %s is not reached in the selected effect direction; the reported effect maximises the",
      "power there (actual power %s), following G*Power (manual 29.3)"), format(target),
      format(signif(best, 7))))
  }
}

# ig_method() arguments for the sensitivity hook, used when the core accepts
# it (forward compatible with change request CR-G6-01).
.ig_glm_sensitivity_hook <- function(power_model) {
  if ("sensitivity_search" %in% names(formals(ig_method))) {
    list(sensitivity_search = .ig_glm_sensitivity_search(power_model))
  } else {
    list()
  }
}

.ig_glm_na_result <- function(reason, n_total) {
  list(power = NA_real_, beta = NA_real_, critical = NA_real_, n_total = as.double(n_total),
       notes = list(power = reason, beta = reason, critical = reason))
}

# Logistic regression ---------------------------------------------------------------------

.ig_method_z_logistic_regression <- function() {
  args <- list(
    id = "z_logistic_regression",
    family = "z",
    label = list(en = "Logistic regression", de = "Logistische Regression"),
    params = ig_params(
      ig_param("odds_ratio", role = "effect",
               domain = ig_domain_real(0, Inf, lower_inclusive = FALSE, upper_inclusive = FALSE),
               label = list(en = "Odds ratio", de = "Odds Ratio"),
               description = paste("OR = [p2 / (1 - p2)] / [p1 / (1 - p1)] = exp(beta1): the odds ratio for a",
                                   "unit increase of X (manual 29.1)."),
               applies_when = list(option = "effect_size_input_mode", "in" = "odds_ratio")),
      ig_param("p2", role = "effect",
               domain = ig_domain_real(0, 1, lower_inclusive = FALSE, upper_inclusive = FALSE),
               label = list(en = "Pr(Y=1|X=1) H1"),
               description = "Event probability at X = 1 under H1, p2 = plogis(beta0 + beta1) (manual 29.1).",
               applies_when = list(option = "effect_size_input_mode", "in" = "two_probabilities")),
      ig_param("p1", role = "nuisance",
               domain = ig_domain_real(0, 1, lower_inclusive = FALSE, upper_inclusive = FALSE),
               label = list(en = "Pr(Y=1|X=1) H0"),
               description = paste("Event probability under H0, p1 = plogis(beta0), i.e. the event probability",
                                   "at X = 0 of the H1 model (manual 29.1).")),
      .ig_glm_common_params(),
      .ig_glm_x_params()
    ),
    options = ig_options(
      .ig_glm_common_options(
        procedures = c("demidenko", "hsieh"),
        procedure_labels = c(demidenko = "Demidenko (2007) - recommended",
                             hsieh = "Hsieh et al. (1998) - for binary and N(0,1) X-Distribution"),
        vc_procedures = "demidenko",
        direction_codes = c("p2_le_p1", "p2_ge_p1"),
        direction_labels = c(p2_le_p1 = "p2 <= p1", p2_ge_p1 = "p2 >= p1"),
        direction_default = "p2_le_p1"),
      ig_option("effect_size_input_mode", values = c("odds_ratio", "two_probabilities"), default = "odds_ratio",
                labels = c(odds_ratio = "Odds ratio", two_probabilities = "Two probabilities"),
                label = list(en = "Input effect size as ..."))
    ),
    modes = c("a_priori", "post_hoc", "sensitivity", "criterion", "compromise"),
    targets = list(sensitivity = c("odds_ratio", "p2")),
    power = .ig_power_z_logistic_regression,
    design = NULL,
    discrete = FALSE,
    effect_direction = "unknown",
    engine = list(package = "stats",
                  fun = paste("integrate (QUADPACK) for the information integrals; pnorm/qnorm;",
                              "Blom-quantile enumeration grid (Lyles et al. 2007)"),
                  version = NA_character_),
    checker = "validation/independent/g6/check_g6.R",
    reference = list(
      gpower_path = "Tests > Correlation and regression > Logistic regression",
      gpower_label = "z tests - Logistic regression",
      manual_pages = "ch. 29 (PDF pp. 70-74)",
      method_notes = "specs/method-notes/z_logistic_regression.md"
    ),
    status = "independently_checked",
    constraints = .ig_glm_constraints
  )
  do.call(ig_method, c(args, .ig_glm_sensitivity_hook(.ig_power_z_logistic_regression)))
}

# Logistic forward model. `p` holds inputs and options (and the solver's
# n_total / alpha / effect).
.ig_power_z_logistic_regression <- function(p) {
  n <- as.double(p$n_total)
  r2 <- p$r2_other_x
  b0 <- stats::qlogis(p$p1)
  b1 <- if (!is.null(p$odds_ratio)) log(p$odds_ratio) else stats::qlogis(p$p2) - b0
  or <- exp(b1)
  p2 <- stats::plogis(b0 + b1)
  if (.ig_glm_outside_direction(p, b1)) {
    return(.ig_glm_na_result(sprintf("the effect (odds ratio %s) lies outside the selected effect direction %s",
                                     .ig_format_double(or), p$effect_direction), n))
  }
  xs <- .ig_glm_xspec(p)
  tails <- p$tails
  outputs <- list(odds_ratio = or, p2 = p2, beta0 = b0, beta1 = b1)
  assumptions <- c(
    "logit P(Y = 1 | X) = beta0 + beta1 X; H0: beta1 = 0 (other coefficients unrestricted)",
    sprintf("beta0 = logit(p1) = %s, beta1 = log(odds ratio) = %s", .ig_format_double(b0), .ig_format_double(b1)),
    sprintf("effective sample size N (1 - R2) for R\u00b2 other X = %s (Hsieh et al. 1998)", .ig_format_double(r2))
  )

  enumeration <- identical(p$computation, "enumeration") && n < p$enumeration_threshold
  if (enumeration) {
    ng <- .ig_glm_grid_size(n, r2)
    if (ng < 2) {
      lambda <- 0
      grid_note <- "the exemplary data set has fewer than 2 covariate values, so beta1 is not identified (lambda = 0)"
    } else {
      x <- .ig_glm_grid(xs, ng)
      pr <- stats::plogis(b0 + b1 * x)
      outputs$marginal_event_probability <- mean(pr)
      if (identical(p$test_statistic, "wald")) {
        w <- pr * (1 - pr)
        s0 <- mean(w)
        mw <- sum(w * x) / sum(w)
        varw <- sum(w * (x - mw)^2) / sum(w)
        lambda <- ng * b1^2 * s0 * varw
      } else {
        # Expected deviance difference per observation between the H1 model
        # and the fitted intercept-only model (pbar), with 0 log 0 = 0.
        pbar <- mean(pr)
        eta <- b0 + b1 * x
        qr <- stats::plogis(-eta)
        qbar <- mean(qr)
        lp <- stats::plogis(eta, log.p = TRUE)
        lq <- stats::plogis(-eta, log.p = TRUE)
        # Keep both event and non-event means: either complement can round
        # to zero even while its log probability is finite.
        lpbar <- if (qbar < 0.5) log1p(-qbar) else log(pbar)
        lqbar <- if (pbar < 0.5) log1p(-pbar) else log(qbar)
        term <- ifelse(pr > 0, pr * (lp - lpbar), 0) + ifelse(qr > 0, qr * (lq - lqbar), 0)
        lambda <- max(2 * ng * mean(term), 0)
      }
      grid_note <- NULL
    }
    cp <- .ig_glm_chisq_power(lambda, p$alpha, tails)
    outputs$grid_size <- ng
    outputs$procedure <- paste0("enumeration_", p$test_statistic)
    return(list(
      power = cp$power, beta = cp$beta, alpha_attained = p$alpha,
      critical = cp$critical, df = 1, ncp = lambda, n_total = n,
      statistic_distribution = list(h0 = list(family = "chisq", df = 1, ncp = 0),
                                    h1 = list(family = "chisq", df = 1, ncp = lambda),
                                    rejection_side = "upper"),
      outputs = outputs,
      assumptions = c(assumptions, grid_note, sprintf(paste(
        "enumeration procedure (Lyles et al. 2007): exemplary data set of %s covariate values at the",
        "quantiles (i - 3/8) / (Ng + 1/4) of the X distribution; %s noncentrality; chi-square(1) test with",
        "the critical value at level %s (G*Power: 2 alpha for one tail); the power counts both tails"),
        .ig_format_double(ng), if (identical(p$test_statistic, "wald")) "Wald" else "likelihood-ratio",
        .ig_format_double(cp$a_eff)))
    ))
  }

  # Large-sample approximations.
  procedure <- if (identical(p$computation, "enumeration")) "demidenko" else p$test_procedure
  vc <- if (identical(p$computation, "enumeration")) "on" else (p$variance_correction %||% "on")
  if (identical(p$computation, "enumeration")) {
    assumptions <- c(assumptions, sprintf(paste(
      "N = %s is not below the enumeration threshold %s: the large-sample approximation of Demidenko",
      "(2007) with variance correction is used (manual 29.2)"), .ig_format_double(n),
      .ig_format_double(p$enumeration_threshold)))
  }
  if (identical(procedure, "hsieh")) {
    std_normal <- identical(xs$type, "normal") && xs$mu == 0 && xs$sigma == 1
    if (identical(xs$type, "binomial") || std_normal) {
      return(.ig_logistic_hsieh(n, r2, p$p1, p2, b1, xs, p$alpha, tails, outputs, assumptions))
    }
    procedure <- "demidenko"
    vc <- "on"
    assumptions <- c(assumptions, paste(
      "Hsieh et al. (1998) covers only binary and standard normal covariates; as in G*Power the",
      "Demidenko (2007) procedure with variance correction is used instead (manual 29.2)"))
  }

  if (identical(xs$type, "manual")) {
    v1 <- xs$v1
    v0 <- xs$v0
    a <- if (identical(vc, "on")) 1 else 0
  } else {
    vv <- .ig_logistic_variances(xs, b0, b1)
    v1 <- vv$v1
    v0 <- vv$v0
    a <- if (identical(vc, "on")) switch(xs$type, lognormal = 0.75, binomial = 0.85, 1) else 0
    if (is.finite(vv$mu)) outputs$marginal_event_probability <- vv$mu
  }
  if (!is.finite(v1) || (a > 0 && !is.finite(v0))) {
    return(.ig_glm_na_result("the information integrals could not be evaluated (numerical integration failed)", n))
  }
  neff <- n * (1 - r2)
  m1 <- sqrt(neff / v1) * b1
  s1 <- if (a > 0) sqrt((a * v0 + (1 - a) * v1) / v1) else 1
  zp <- .ig_glm_z_power(m1, s1, p$alpha, tails)
  outputs <- c(outputs, list(v1 = v1, m1 = m1, s1 = s1, variance_correction_a = a, procedure = "demidenko"))
  if (a > 0 && is.finite(v0)) outputs$v0 <- v0
  list(
    power = zp$power, beta = zp$beta, alpha_attained = p$alpha, critical = zp$critical,
    ncp = m1, n_total = n,
    statistic_distribution = list(h0 = list(family = "normal", mean = 0, sd = 1),
                                  h1 = list(family = "normal", mean = m1, sd = s1),
                                  rejection_side = if (identical(tails, "two")) "both" else if (zp$side < 0) "lower" else "upper"),
    outputs = outputs,
    assumptions = c(assumptions, sprintf(paste(
      "Demidenko (2007) large-sample Wald test: H1 statistic N(m1, s1) with m1 = sqrt(N (1 - R2) / v1) beta1",
      "and s1 = sqrt((a v0 + (1 - a) v1) / v1), a = %s (%s)"), .ig_format_double(a),
      if (a > 0) "with variance correction" else "without variance correction"),
      .ig_glm_side_note(tails))
  )
}

# Hsieh et al. (1998) equations (1) and (2) solved for power. As in the
# published sample-size equations only the rejection region on the side of
# the effect is counted (also for two tails, with alpha/2).
.ig_logistic_hsieh <- function(n, r2, p1, p2, b1, xs, alpha, tails, outputs, assumptions) {
  z <- if (identical(tails, "one")) stats::qnorm(alpha, lower.tail = FALSE) else stats::qnorm(alpha / 2, lower.tail = FALSE)
  side <- if (b1 < 0) -1 else 1
  neff <- n * (1 - r2)
  if (identical(xs$type, "binomial")) {
    bb <- xs$pi
    pbar <- (1 - bb) * p1 + bb * p2
    s0 <- sqrt(pbar * (1 - pbar) / bb)
    s1 <- sqrt(p1 * (1 - p1) + p2 * (1 - p2) * (1 - bb) / bb)
    m <- sqrt(neff * (1 - bb)) * abs(p2 - p1)
    arg <- (m - z * s0) / s1
    form <- "Hsieh et al. (1998) eq. (2), binary X with Pr(X = 1) = B"
    outputs <- c(outputs, list(hsieh_sd0 = s0, hsieh_sd1 = s1, hsieh_mean = m))
  } else {
    arg <- sqrt(neff * p1 * (1 - p1)) * abs(b1) - z
    form <- "Hsieh et al. (1998) eq. (1), standard normal X"
  }
  power <- stats::pnorm(arg)
  beta <- stats::pnorm(-arg)
  crit <- if (identical(tails, "one")) side * z else c(-z, z)
  outputs$procedure <- "hsieh"
  outputs$hsieh_argument <- arg
  list(
    power = power, beta = beta, alpha_attained = alpha, critical = crit, n_total = n,
    statistic_distribution = list(h0 = list(family = "normal", mean = 0, sd = 1),
                                  h1 = list(family = "normal", mean = arg + z, sd = 1),
                                  rejection_side = if (side < 0) "lower" else "upper"),
    outputs = outputs,
    assumptions = c(assumptions, paste0(form, " solved for power; only the rejection region on the side",
                                        " of the effect is counted (as in the published sample-size equations)"))
  )
}

# Poisson regression ---------------------------------------------------------------------

.ig_method_z_poisson_regression <- function() {
  args <- list(
    id = "z_poisson_regression",
    family = "z",
    label = list(en = "Poisson regression", de = "Poisson-Regression"),
    params = ig_params(
      ig_param("exp_beta1", role = "effect",
               domain = ig_domain_real(0, Inf, lower_inclusive = FALSE, upper_inclusive = FALSE),
               label = list(en = "Exp(\u03b21)"),
               description = paste("Rate ratio exp(beta1) for a unit increase of X (manual 30.1); 1 means no",
                                   "effect.")),
      ig_param("base_rate", role = "nuisance",
               domain = ig_domain_real(0, Inf, lower_inclusive = FALSE, upper_inclusive = FALSE),
               label = list(en = "Base rate exp(\u03b20)"),
               description = "Mean event rate exp(beta0) under H0 (manual 30.1)."),
      ig_param("mean_exposure", role = "nuisance",
               domain = ig_domain_real(0, Inf, lower_inclusive = FALSE, upper_inclusive = FALSE), default = 1,
               label = list(en = "Mean exposure"),
               description = paste("Exposure time t during which events are counted (manual 30.1). Default 1",
                                   "is the G*Power field default (live capture 2026-10-05).")),
      .ig_glm_common_params(),
      .ig_glm_x_params()
    ),
    options = ig_options(
      .ig_glm_common_options(
        procedures = c("demidenko", "signorini"),
        procedure_labels = c(demidenko = "Demidenko (2007) - recommended",
                             signorini = "Signorini (1991) - not recommended"),
        vc_procedures = c("demidenko", "signorini"),
        direction_codes = c("beta1_le_0", "beta1_ge_0"),
        direction_labels = c(beta1_le_0 = "\u03b21 <= 0", beta1_ge_0 = "\u03b21 >= 0"),
        direction_default = "beta1_ge_0")
    ),
    modes = c("a_priori", "post_hoc", "sensitivity", "criterion", "compromise"),
    power = .ig_power_z_poisson_regression,
    design = NULL,
    discrete = FALSE,
    effect_direction = "unknown",
    engine = list(package = "stats",
                  fun = paste("integrate (QUADPACK) for the information integrals; pnorm/qnorm;",
                              "Blom-quantile enumeration grid (Lyles et al. 2007)"),
                  version = NA_character_),
    checker = "validation/independent/g6/check_g6.R",
    reference = list(
      gpower_path = "Tests > Correlation and regression > Poisson regression",
      gpower_label = "z tests - Poisson regression",
      manual_pages = "ch. 30 (PDF pp. 75-79)",
      method_notes = "specs/method-notes/z_poisson_regression.md"
    ),
    status = "independently_checked",
    constraints = function(p) {
      out <- .ig_glm_constraints(p)
      # The large-sample procedures need E[exp(beta1 X)]; the enumeration
      # procedure works on a finite exemplary design and stays defined.
      if (!is.null(p$exp_beta1) && !identical(p$computation, "enumeration")) {
        b1 <- log(p$exp_beta1)
        div <- .ig_poisson_divergence(p, b1)
        if (!is.null(div)) out <- c(out, div)
      }
      out
    }
  )
  do.call(ig_method, c(args, .ig_glm_sensitivity_hook(.ig_power_z_poisson_regression)))
}

# The Poisson information integrals E[X^k exp(beta1 X)] are infinite for a
# lognormal covariate with beta1 > 0 (no moment generating function) and for
# an exponential covariate with beta1 >= its rate. Returns the reason or NULL.
.ig_poisson_divergence <- function(p, b1) {
  xd <- p$x_distribution
  if (identical(xd, "lognormal") && b1 > 0) {
    return(paste("for a lognormal covariate and Exp(\u03b21) > 1 the expected count E[exp(beta0 + beta1 X)] and",
                 "the information integrals are infinite (the lognormal distribution has no moment generating",
                 "function for positive arguments), so the large-sample variance of beta1-hat is undefined"))
  }
  if (identical(xd, "exponential") && !is.null(p$x_lambda) && b1 >= p$x_lambda) {
    return(sprintf(paste("for an exponential covariate with rate %s the information integrals are infinite when",
                         "beta1 = log(Exp(\u03b21)) >= the rate (beta1 = %s)"),
                   .ig_format_double(p$x_lambda), .ig_format_double(b1)))
  }
  NULL
}

# Enumeration procedure for counts (Lyles et al. 2007, as G*Power computes it).
#
# The exemplary data set holds, for every covariate value x_i of the grid, a
# range of counts y = L..K with weights P(Y = y | m_i), m_i = t exp(beta0 +
# beta1 x_i). The count range must be cut off somewhere (Lyles et al.: where
# the remaining probability is negligible). G*Power's rule was identified
# with dedicated captures (SUP-G6-Z-POIS-POS-P01, Q01 ... Q13 and the campaign
# cases with mean exposure 10; method notes): one range for all covariate
# values, from L = the 1e-5 quantile of the smallest expected count to K = the
# 1 - 1e-5 quantile of the largest expected count; the weights are not
# renormalised. The Poisson regression is then fitted to these weighted data
# and
#   Wald:             lambda = beta1^2 / Var(beta1-hat), the variance from the
#                     information matrix at the fitted means (beta1 is the
#                     hypothesised value, not the fitted one);
#   likelihood ratio: lambda = 2 [l(hypothesised parameters) - l(fitted
#                     intercept-only model)]; the full model is not refitted
#                     to the cut-off data (33 likelihood-ratio captures agree
#                     with this form to half a displayed digit; with a
#                     refitted full model four of them do not, and a grid of
#                     equal covariate values shows G*Power's small negative
#                     value, -1.256486e-09 in GPW-Z-POIS-POS-B004-E01, which
#                     only this form produces; reported as 0 here).
# Without the cut-off the fit returns the true parameters and both quantities
# have closed forms (`lambda_exact`); the cut-off lowers them by about 1e-5
# relative. Only the sufficient statistics per covariate value are needed:
# S_i = P(L <= Y <= K | m_i) (total weight) and T_i = m_i P(L - 1 <= Y <= K - 1
# | m_i) (weighted sum of the counts). The fit is a Newton iteration from the
# true parameters on the centred covariate.
#
# Returns list(lambda, lambda_exact, count_range, expected_events, note) or
# list(reason) when the quantity is undefined.
.ig_glm_count_tail <- 1e-5

.ig_poisson_enumeration <- function(x, m, b1, statistic) {
  wald <- identical(statistic, "wald")
  if (!all(is.finite(m)) || any(m <= 0)) {
    return(list(reason = paste("the expected counts of the exemplary data set are not all positive and finite",
                               "(overflow or underflow of exp(beta0 + beta1 x))")))
  }
  sm <- sum(m)
  cx <- sum(m * x) / sm
  xc <- x - cx
  lambda_exact <- if (wald) b1^2 * sum(m * xc^2) else {
    mbar <- sm / length(m)
    2 * sum(m * log(m / mbar) - (m - mbar))
  }
  k <- stats::qpois(1 - .ig_glm_count_tail, max(m))
  l <- stats::qpois(.ig_glm_count_tail, min(m))
  s <- stats::ppois(k, m) - stats::ppois(l - 1, m)
  tt <- m * (stats::ppois(k - 1, m) - stats::ppois(l - 2, m))
  out <- list(lambda_exact = lambda_exact, count_range = c(l, k), expected_events = sum(tt), note = NULL)
  if (!(sum(tt) > 0)) {
    return(list(reason = sprintf(paste(
      "the largest expected count (%s) is so small that the count range of the enumeration procedure, cut off",
      "at the 1 - 1e-5 quantile, holds only the count 0: the exemplary data set contains no events and no",
      "model can be fitted (the expected-data noncentrality without the cut-off would be %s)"),
      .ig_format_double(max(m)), .ig_format_double(lambda_exact))))
  }
  if (max(x) == min(x)) {
    out$lambda <- 0
    out$note <- "all covariate values of the exemplary data set are equal, so beta1 is not identified (lambda = 0)"
    return(out)
  }
  if (!wald) {
    # Log likelihood at the hypothesised parameters minus the maximised log
    # likelihood of the intercept-only model (its fitted mean is sum(T) / sum(S)).
    mu0 <- sum(tt) / sum(s)
    out$lambda <- max(2 * sum(tt * log(m / mu0) - s * (m - mu0)), 0)
    return(out)
  }
  g1 <- b1
  g0 <- mean(log(m) - b1 * xc)
  converged <- FALSE
  for (it in seq_len(50L)) {
    mu <- exp(g0 + g1 * xc)
    w <- s * mu
    r <- tt - w
    i00 <- sum(w); i01 <- sum(w * xc); i11 <- sum(w * xc^2)
    det <- i00 * i11 - i01^2
    if (!is.finite(det) || det <= 0) break
    u0 <- sum(r); u1 <- sum(r * xc)
    d0 <- (i11 * u0 - i01 * u1) / det
    d1 <- (i00 * u1 - i01 * u0) / det
    g0 <- g0 + d0
    g1 <- g1 + d1
    if (max(abs(d0), abs(d1)) < 1e-13) {
      converged <- TRUE
      break
    }
  }
  if (!converged) {
    return(list(reason = "the Poisson regression fit to the exemplary data set of the enumeration procedure did not converge"))
  }
  w <- s * exp(g0 + g1 * xc)
  xw <- sum(w * xc) / sum(w)
  out$lambda <- b1^2 * sum(w * (xc - xw)^2)
  out
}

.ig_power_z_poisson_regression <- function(p) {
  n <- as.double(p$n_total)
  r2 <- p$r2_other_x
  t <- p$mean_exposure
  b0 <- log(p$base_rate)
  b1 <- log(p$exp_beta1)
  if (.ig_glm_outside_direction(p, b1)) {
    return(.ig_glm_na_result(sprintf("the effect (Exp(beta1) %s) lies outside the selected effect direction %s",
                                     .ig_format_double(p$exp_beta1), p$effect_direction), n))
  }
  div <- .ig_poisson_divergence(p, b1)
  xs <- .ig_glm_xspec(p)
  tails <- p$tails
  outputs <- list(exp_beta1 = p$exp_beta1, beta0 = b0, beta1 = b1)
  assumptions <- c(
    "log E[Y | X] = beta0 + beta1 X + log(t) (Poisson response, exposure t); H0: beta1 = 0",
    sprintf("beta0 = log(base rate) = %s, beta1 = log(Exp(beta1)) = %s, mean exposure t = %s",
            .ig_format_double(b0), .ig_format_double(b1), .ig_format_double(t)),
    sprintf("effective sample size N (1 - R2) for R\u00b2 other X = %s (Hsieh et al. 1998)", .ig_format_double(r2))
  )

  enumeration <- identical(p$computation, "enumeration") && n < p$enumeration_threshold
  if (enumeration) {
    ng <- .ig_glm_grid_size(n, r2)
    if (ng < 2) {
      lambda <- 0
      grid_note <- "the exemplary data set has fewer than 2 covariate values, so beta1 is not identified (lambda = 0)"
    } else {
      x <- .ig_glm_grid(xs, ng, tie = "lower")
      en <- .ig_poisson_enumeration(x, t * exp(b0 + b1 * x), b1, p$test_statistic)
      if (!is.null(en$reason)) return(.ig_glm_na_result(en$reason, n))
      lambda <- en$lambda
      outputs$ncp_untruncated <- en$lambda_exact
      outputs$enumeration_count_range <- en$count_range
      outputs$enumeration_expected_events <- en$expected_events
      grid_note <- en$note
      if (!is.null(div)) {
        grid_note <- c(grid_note, paste(
          "the covariate distribution has no finite E[exp(beta1 X)] for this effect, so the large-sample",
          "procedures are undefined; the enumeration result describes the finite exemplary design only and is",
          "dominated by its largest covariate values"))
      }
    }
    cp <- .ig_glm_chisq_power(lambda, p$alpha, tails)
    outputs$grid_size <- ng
    outputs$procedure <- paste0("enumeration_", p$test_statistic)
    return(list(
      power = cp$power, beta = cp$beta, alpha_attained = p$alpha,
      critical = cp$critical, df = 1, ncp = lambda, n_total = n,
      statistic_distribution = list(h0 = list(family = "chisq", df = 1, ncp = 0),
                                    h1 = list(family = "chisq", df = 1, ncp = lambda),
                                    rejection_side = "upper"),
      outputs = outputs,
      assumptions = c(assumptions, grid_note, sprintf(paste(
        "enumeration procedure (Lyles et al. 2007): exemplary data set of %s covariate values at the",
        "quantiles (i - 3/8) / (Ng + 1/4) of the X distribution, each with the counts L..K weighted by their",
        "Poisson probabilities, L the 1e-5 quantile of the smallest and K the 1 - 1e-5 quantile of the largest",
        "expected count (G*Power's cut-off); %s noncentrality from these data",
        "(outputs$ncp_untruncated: the value without the cut-off); chi-square(1) test with the",
        "critical value at level %s (G*Power: 2 alpha for one tail); the power counts both tails"),
        .ig_format_double(ng),
        if (identical(p$test_statistic, "wald")) "Wald" else "likelihood-ratio", .ig_format_double(cp$a_eff)))
    ))
  }

  if (!is.null(div) && !identical(p$x_distribution, "manual")) return(.ig_glm_na_result(div, n))
  procedure <- if (identical(p$computation, "enumeration")) "demidenko" else p$test_procedure
  vc <- if (identical(p$computation, "enumeration")) "on" else (p$variance_correction %||% "on")
  if (identical(p$computation, "enumeration")) {
    assumptions <- c(assumptions, sprintf(paste(
      "N = %s is not below the enumeration threshold %s: the large-sample approximation of Demidenko",
      "(2007) with variance correction is used (manual 30.2)"), .ig_format_double(n),
      .ig_format_double(p$enumeration_threshold)))
  }
  if (identical(xs$type, "manual")) {
    v1 <- xs$v1
    v0 <- xs$v0
    v0s <- xs$v0
  } else {
    vv <- .ig_poisson_variances(xs, b0, b1)
    v1 <- vv$v1
    v0 <- vv$v0
    v0s <- vv$v0_star
    if (is.finite(vv$mu)) outputs$marginal_rate <- vv$mu
  }
  if (!is.finite(v1) || !is.finite(v0) || !is.finite(v0s)) {
    return(.ig_glm_na_result("the information integrals could not be evaluated (numerical integration failed)", n))
  }
  neff <- n * (1 - r2)
  if (identical(procedure, "signorini") && identical(vc, "off") && !(xs$type %in% c("binomial", "manual"))) {
    # "With variance correction" unchecked under Signorini: G*Power replaces
    # the null variance v0 by v1, so the statistic is N(m1, 1) under H1 and
    # the result equals the Demidenko procedure without correction. G*Power
    # does this for exponential, lognormal, normal, Poisson and uniform X; for
    # binomial and Manual X the box has no effect there (method notes,
    # "Signorini and the variance-correction box"). The option is a G*Power
    # control without a counterpart in Signorini (1991); its scope is taken
    # from the reference and stated in the result.
    m1 <- b1 * sqrt(neff * t / v1)
    s1 <- 1
    a <- NA_real_
    form <- paste("Signorini (1991) with the variance-correction box unchecked: as in G*Power the null variance v0 is",
                  "replaced by v1, i.e. m1 = beta1 sqrt(N (1 - R2) t / v1), s1 = 1 (equal to the Demidenko procedure",
                  "without correction)")
  } else if (identical(procedure, "signorini")) {
    m1 <- b1 * sqrt(neff * t / v0)
    s1 <- sqrt(v1 / v0)
    a <- NA_real_
    form <- "Signorini (1991): m1 = beta1 sqrt(N (1 - R2) t / v0), s1 = sqrt(v1 / v0), v0 at beta1 = 0 and beta0"
    if (identical(vc, "off")) {
      form <- paste0(form, "; the unchecked variance-correction box has no effect for a ",
                     if (identical(xs$type, "manual")) "Manual" else "binomial",
                     " covariate (as in G*Power, which changes the result only for the other X distributions)")
    }
  } else {
    a <- if (identical(vc, "on")) (if (identical(xs$type, "lognormal")) 0.75 else 1) else 0
    m1 <- b1 * sqrt(neff * t / v1)
    s1 <- sqrt((a * v0s + (1 - a) * v1) / v1)
    form <- sprintf(paste("Demidenko (2007): m1 = beta1 sqrt(N (1 - R2) t / v1), s1 = sqrt((a v0* + (1 - a) v1) / v1),",
                          "v0* at beta1 = 0 and beta0* = log E[exp(beta0 + beta1 X)], a = %s (%s)"),
                    .ig_format_double(a), if (a > 0) "with variance correction" else "without variance correction")
  }
  zp <- .ig_glm_z_power(m1, s1, p$alpha, tails)
  outputs <- c(outputs, list(v1 = v1, v0 = v0, v0_star = v0s, m1 = m1, s1 = s1, procedure = procedure))
  if (!is.na(a)) outputs$variance_correction_a <- a
  list(
    power = zp$power, beta = zp$beta, alpha_attained = p$alpha, critical = zp$critical,
    ncp = m1, n_total = n,
    statistic_distribution = list(h0 = list(family = "normal", mean = 0, sd = 1),
                                  h1 = list(family = "normal", mean = m1, sd = s1),
                                  rejection_side = if (identical(tails, "two")) "both" else if (zp$side < 0) "lower" else "upper"),
    outputs = outputs,
    assumptions = c(assumptions, paste("large-sample Wald test;", form), .ig_glm_side_note(tails))
  )
}
