# Distribution calculator: the allowlisted function table --------------------------------
#
# Ownership: phase 10 utilities (group G12). Method notes:
# specs/method-notes/utilities.md, section "Distribution calculator".
#
# Every function of the 'G*Power' 3.1.9.7 calculator (manual ch. 2, the
# Windows function-name table and the live capture of 2026-10-05,
# specs/gpower-utilities.json util_calculator u006-u079) with its argument
# order, argument domain and implementation. The calculator evaluates nothing
# else: a name that is not in this table is an unknown function.
#
# Conventions (decided from the live captures; method notes):
#
# * Arguments are given in the G*Power order, e.g. hygecdf(x, N, ns, nt) with
#   sample size N, ns successes in a population of nt; atan2(y, x) is the
#   angle of the point (x, y) (atan2(1, 2) = atan(1/2)).
# * A parameter outside its mathematical domain (an SD <= 0, df <= 0, a
#   probability outside [0, 1], ...) is a math error: the calculator refuses
#   with status invalid_input and no number, as G*Power shows 'Math error in
#   formula ...'. The value of a CDF or density outside the support of the
#   distribution is defined (0 or 1) and returned.
# * Quantile functions of distributions with unbounded support refuse p = 0
#   and p = 1 (the quantile is infinite; G*Power: 'Math error').
# * Elementary functions follow IEEE arithmetic as G*Power does: log(-1) and
#   sqrt(-1) are NaN, log(0) is -Inf, 1/0 is Inf (no message in G*Power).
# * Where G*Power 3.1.9.7 computes a different function than its name says
#   (chi2pdf and ncx2pdf evaluate the quantile functions chi2inv and ncx2inv;
#   unipdf evaluates the CDF), the calculator computes the function the name
#   and the manual define and documents the reference defect (`gpower_note`).
#
# Each entry: name, args (argument names), category, impl (function of the
# numeric argument vector), domain (function returning NULL or the reason of a
# math error), description, gpower_note.

.ig_calc_fn <- function(name, args, category, impl, description, domain = NULL, gpower_note = NA_character_,
                        range = NULL) {
  list(name = name, args = args, category = category, impl = impl, domain = domain, range = range,
       description = description, gpower_note = gpower_note)
}

# Largest parameter values for which the calculator's own numerical routines
# (Poisson and negative binomial mixtures, quadrature, discrete searches) are
# used. Larger values are valid mathematically; they are refused with status
# unsupported instead of a slow or unchecked computation. Functions that go
# straight to R's nmath have no such limit.
.ig_calc_range <- list(nc = 1e8, nct_nc = 1e6, df = 1e9, count = 1e9, search = 1e7, mixture_terms = 5e5, corr_n = 1e7)
.ig_calc_at_most <- function(v, limit, what, fun) {
  list(is.finite(v) && abs(v) <= limit,
       sprintf("%s: %s beyond %s is not computed", fun, what, format(limit, scientific = TRUE)))
}
.ig_calc_range_nc <- function(fun, nc, ..., limit = .ig_calc_range$nc) {
  dfs <- list(...)
  function(a) {
    checks <- c(list(.ig_calc_at_most(a[nc], limit, "a noncentrality parameter", fun)),
                lapply(dfs, function(i) .ig_calc_at_most(a[i], .ig_calc_range$df, "degrees of freedom", fun)))
    do.call(.ig_calc_need, checks)
  }
}
.ig_calc_range_corr <- function(fun) {
  function(a) .ig_calc_need(.ig_calc_at_most(a[3], .ig_calc_range$corr_n, "a sample size N", fun))
}
.ig_calc_range_hyge <- function(fun, search = FALSE) {
  function(a) {
    .ig_calc_need(.ig_calc_at_most(a[4], .ig_calc_range$count, "a population size nt", fun),
                  if (search) .ig_calc_at_most(min(.ig_calc_round_count(a[2]), a[3]), .ig_calc_range$search,
                                               "a quantile search over min(N, ns)", fun))
  }
}

# Domain helpers (return the reason of a math error or NULL). Each argument
# of .ig_calc_need() is a pair list(condition, reason), or the result of
# another domain helper: NULL (nothing to report) or the reason itself.
.ig_calc_need <- function(...) {
  for (cn in list(...)) {
    if (is.null(cn)) next
    if (is.character(cn)) return(cn[1])
    if (!isTRUE(cn[[1]])) return(cn[[2]])
  }
  NULL
}
.ig_calc_prob_open <- function(p, fun) {
  list(is.finite(p) && p > 0 && p < 1,
       sprintf("%s: p must lie strictly between 0 and 1 (the quantile at 0 or 1 is infinite)", fun))
}
.ig_calc_prob_closed <- function(p, fun) list(is.finite(p) && p >= 0 && p <= 1, sprintf("%s: p must lie in [0, 1]", fun))
.ig_calc_pos <- function(v, what, fun) list(is.finite(v) && v > 0, sprintf("%s: %s must be greater than 0", fun, what))
.ig_calc_nonneg <- function(v, what, fun) list(is.finite(v) && v >= 0, sprintf("%s: %s must be non-negative", fun, what))
.ig_calc_finite <- function(v, what, fun) list(is.finite(v), sprintf("%s: %s must be finite", fun, what))

# Round half up, as G*Power does for the count argument of the binomial mass
# and the hypergeometric functions (captured: binopdf(2.5, ...) = binopdf(3, ...)).
.ig_calc_round_count <- function(x) floor(x + 0.5)

.ig_calc_function_table <- function() {
  fns <- list(
    # Elementary functions ----------------------------------------------------------------
    .ig_calc_fn("abs", "x", "function", function(a) abs(a[1]), "absolute value"),
    .ig_calc_fn("sin", "x", "function", function(a) .ig_calc_trig(a[1], "sin"), "sine (radians)"),
    .ig_calc_fn("asin", "x", "function", function(a) asin(a[1]), "arc sine (radians); NaN outside [-1, 1]"),
    .ig_calc_fn("cos", "x", "function", function(a) .ig_calc_trig(a[1], "cos"), "cosine (radians)"),
    .ig_calc_fn("acos", "x", "function", function(a) acos(a[1]), "arc cosine (radians); NaN outside [-1, 1]"),
    .ig_calc_fn("tan", "x", "function", function(a) .ig_calc_trig(a[1], "tan"), "tangent (radians)"),
    .ig_calc_fn("atan", "x", "function", function(a) atan(a[1]), "arc tangent (radians)"),
    .ig_calc_fn("atan2", c("y", "x"), "function", function(a) atan2(a[1], a[2]),
                "angle of the point (x, y) in (-pi, pi]; atan2(y, x) = atan(y/x) for x > 0",
                gpower_note = "manual ch. 2: 'Arc tangent of y/x'; captured live: atan2(1,2) = atan(1/2), atan2(0,-1) = pi"),
    .ig_calc_fn("exp", "x", "function", function(a) exp(a[1]), "exponential function"),
    .ig_calc_fn("log", "x", "function", function(a) log(a[1]), "natural logarithm; -Inf at 0, NaN below 0"),
    .ig_calc_fn("sqrt", "x", "function", function(a) sqrt(a[1]), "square root; NaN below 0"),
    .ig_calc_fn("sqr", "x", "function", function(a) a[1] * a[1], "square x * x"),
    .ig_calc_fn("sign", "x", "function", function(a) sign(a[1]), "sign: -1, 0 or 1"),
    .ig_calc_fn("lngamma", "x", "function", .ig_calc_lngamma,
                "natural logarithm of the gamma function; Inf at the poles 0, -1, -2, ...; NaN where Gamma(x) < 0",
                gpower_note = "G*Power 3.1.9.7 returns NaN for every negative x (captured: lngamma(-1.5) = -nan(ind), although Gamma(-1.5) = 2.363 > 0)"),
    .ig_calc_fn("frac", "x", "function", function(a) a[1] - trunc(a[1]), "fractional part x - int(x) (sign of x)"),
    .ig_calc_fn("int", "x", "function", function(a) trunc(a[1]), "integer part (truncation toward 0)"),
    .ig_calc_fn("min", c("x", "y"), "function", function(a) min(a[1], a[2]), "minimum of two values"),
    .ig_calc_fn("max", c("x", "y"), "function", function(a) max(a[1], a[2]), "maximum of two values"),
    .ig_calc_fn("uround", c("x", "m"), "function", function(a) {
      m <- trunc(a[2])
      ceiling(a[1] / m) * m
    }, paste("round x up to a multiple of int(m): ceiling(x / int(m)) * int(m); NaN when int(m) = 0",
             "(the multiple is truncated to an integer, as in G*Power)"),
    gpower_note = paste("captured live: uround(2.3,2.5) = 4 and uround(5,2.5) = 6 (multiple truncated to 2),",
                        "uround(2.3,0.5) and uround(2.3,0) = -nan(ind), uround(-2.3,1) = -2, uround(2.3,-1) = 2")),

    # Normal -------------------------------------------------------------------------------
    .ig_calc_fn("zcdf", "x", "distribution", function(a) stats::pnorm(a[1]), "standard normal CDF"),
    .ig_calc_fn("zpdf", "x", "distribution", function(a) stats::dnorm(a[1]), "standard normal density"),
    .ig_calc_fn("zinv", "p", "distribution", function(a) stats::qnorm(a[1]), "standard normal quantile",
                domain = function(a) .ig_calc_need(.ig_calc_prob_open(a[1], "zinv"))),
    .ig_calc_fn("normcdf", c("x", "m", "s"), "distribution", function(a) stats::pnorm(a[1], a[2], a[3]),
                "normal CDF with mean m and SD s",
                domain = function(a) .ig_calc_need(.ig_calc_finite(a[2], "m", "normcdf"), .ig_calc_pos(a[3], "s", "normcdf"))),
    .ig_calc_fn("normpdf", c("x", "m", "s"), "distribution", function(a) stats::dnorm(a[1], a[2], a[3]),
                "normal density with mean m and SD s",
                domain = function(a) .ig_calc_need(.ig_calc_finite(a[2], "m", "normpdf"), .ig_calc_pos(a[3], "s", "normpdf"))),
    .ig_calc_fn("norminv", c("p", "m", "s"), "distribution", function(a) stats::qnorm(a[1], a[2], a[3]),
                "normal quantile with mean m and SD s",
                domain = function(a) .ig_calc_need(.ig_calc_prob_open(a[1], "norminv"), .ig_calc_finite(a[2], "m", "norminv"), .ig_calc_pos(a[3], "s", "norminv"))),

    # Chi-square ---------------------------------------------------------------------------
    .ig_calc_fn("chi2cdf", c("x", "df"), "distribution", function(a) stats::pchisq(a[1], a[2]), "chi-square CDF",
                domain = function(a) .ig_calc_need(.ig_calc_pos(a[2], "df", "chi2cdf"))),
    .ig_calc_fn("chi2pdf", c("x", "df"), "distribution", function(a) stats::dchisq(a[1], a[2]), "chi-square density",
                domain = function(a) .ig_calc_need(.ig_calc_pos(a[2], "df", "chi2pdf")),
                gpower_note = paste("reference defect: G*Power 3.1.9.7 evaluates the quantile chi2inv(x, df) under this name",
                                    "(captured: chi2pdf(0.5,2) = 1.38629436111989 = chi2inv(0.5,2); 'Math error' for x > 1)")),
    .ig_calc_fn("chi2inv", c("p", "df"), "distribution", function(a) stats::qchisq(a[1], a[2]), "chi-square quantile",
                domain = function(a) .ig_calc_need(.ig_calc_pos(a[2], "df", "chi2inv"),
                                                   list(is.finite(a[1]) && a[1] >= 0 && a[1] < 1,
                                                        "chi2inv: p must lie in [0, 1) (the quantile at 1 is infinite)"))),
    .ig_calc_fn("ncx2cdf", c("x", "df", "nc"), "distribution", function(a) .ig_calc_ncx2_cdf(a[1], a[2], a[3]),
                "noncentral chi-square CDF (noncentrality nc)",
                domain = function(a) .ig_calc_need(.ig_calc_pos(a[2], "df", "ncx2cdf"), .ig_calc_nonneg(a[3], "nc", "ncx2cdf")),
                range = .ig_calc_range_nc("ncx2cdf", 3, 2)),
    .ig_calc_fn("ncx2pdf", c("x", "df", "nc"), "distribution", function(a) .ig_calc_ncx2_pdf(a[1], a[2], a[3]),
                "noncentral chi-square density (noncentrality nc)",
                domain = function(a) .ig_calc_need(.ig_calc_pos(a[2], "df", "ncx2pdf"), .ig_calc_nonneg(a[3], "nc", "ncx2pdf")),
                range = .ig_calc_range_nc("ncx2pdf", 3, 2),
                gpower_note = paste("reference defect: G*Power 3.1.9.7 evaluates the quantile ncx2inv(x, df, nc) under this name",
                                    "(captured: ncx2pdf(0.5,3,2) = 4.13751512339913 = ncx2inv(0.5,3,2); 'Math error' for x > 1)")),
    .ig_calc_fn("ncx2inv", c("p", "df", "nc"), "distribution", function(a) .ig_calc_ncx2_inv(a[1], a[2], a[3]),
                "noncentral chi-square quantile (noncentrality nc)",
                domain = function(a) .ig_calc_need(.ig_calc_pos(a[2], "df", "ncx2inv"), .ig_calc_nonneg(a[3], "nc", "ncx2inv"),
                                                   list(is.finite(a[1]) && a[1] >= 0 && a[1] < 1,
                                                        "ncx2inv: p must lie in [0, 1) (the quantile at 1 is infinite)")),
                range = .ig_calc_range_nc("ncx2inv", 3, 2)),

    # F ----------------------------------------------------------------------------------------
    .ig_calc_fn("fcdf", c("x", "df1", "df2"), "distribution", function(a) stats::pf(a[1], a[2], a[3]), "F CDF",
                domain = function(a) .ig_calc_need(.ig_calc_pos(a[2], "df1", "fcdf"), .ig_calc_pos(a[3], "df2", "fcdf"))),
    .ig_calc_fn("fpdf", c("x", "df1", "df2"), "distribution", function(a) stats::df(a[1], a[2], a[3]), "F density",
                domain = function(a) .ig_calc_need(.ig_calc_pos(a[2], "df1", "fpdf"), .ig_calc_pos(a[3], "df2", "fpdf")),
                gpower_note = paste("G*Power 3.1.9.7 returns 0 at x = 0 (captured: fpdf(0,2,20) = 0, where the density is 1,",
                                    "and fpdf(0,1,5) = 0, where it is unbounded)")),
    .ig_calc_fn("finv", c("p", "df1", "df2"), "distribution", function(a) stats::qf(a[1], a[2], a[3]), "F quantile",
                domain = function(a) .ig_calc_need(.ig_calc_pos(a[2], "df1", "finv"), .ig_calc_pos(a[3], "df2", "finv"),
                                                   list(is.finite(a[1]) && a[1] >= 0 && a[1] < 1,
                                                        "finv: p must lie in [0, 1) (the quantile at 1 is infinite)"))),
    .ig_calc_fn("ncfcdf", c("x", "df1", "df2", "nc"), "distribution", function(a) .ig_calc_ncf_cdf(a[1], a[2], a[3], a[4]),
                "noncentral F CDF (noncentrality nc)",
                domain = function(a) .ig_calc_need(.ig_calc_pos(a[2], "df1", "ncfcdf"), .ig_calc_pos(a[3], "df2", "ncfcdf"),
                                                   .ig_calc_nonneg(a[4], "nc", "ncfcdf")),
                range = .ig_calc_range_nc("ncfcdf", 4, 2, 3)),
    .ig_calc_fn("ncfpdf", c("x", "df1", "df2", "nc"), "distribution", function(a) {
      if (a[1] < 0 || is.infinite(a[1])) 0 else stats::df(a[1], a[2], a[3], ncp = a[4])
    },
                "noncentral F density (noncentrality nc)",
                domain = function(a) .ig_calc_need(.ig_calc_pos(a[2], "df1", "ncfpdf"), .ig_calc_pos(a[3], "df2", "ncfpdf"),
                                                   .ig_calc_nonneg(a[4], "nc", "ncfpdf")),
                range = .ig_calc_range_nc("ncfpdf", 4, 2, 3)),
    .ig_calc_fn("ncfinv", c("p", "df1", "df2", "nc"), "distribution", function(a) .ig_calc_ncf_inv(a[1], a[2], a[3], a[4]),
                "noncentral F quantile (noncentrality nc)",
                domain = function(a) .ig_calc_need(.ig_calc_pos(a[2], "df1", "ncfinv"), .ig_calc_pos(a[3], "df2", "ncfinv"),
                                                   .ig_calc_nonneg(a[4], "nc", "ncfinv"),
                                                   list(is.finite(a[1]) && a[1] >= 0 && a[1] < 1,
                                                        "ncfinv: p must lie in [0, 1) (the quantile at 1 is infinite)")),
                range = .ig_calc_range_nc("ncfinv", 4, 2, 3)),

    # t ----------------------------------------------------------------------------------------
    .ig_calc_fn("tcdf", c("x", "df"), "distribution", function(a) stats::pt(a[1], a[2]), "Student t CDF",
                domain = function(a) .ig_calc_need(.ig_calc_pos(a[2], "df", "tcdf"))),
    .ig_calc_fn("tpdf", c("x", "df"), "distribution", function(a) stats::dt(a[1], a[2]), "Student t density",
                domain = function(a) .ig_calc_need(.ig_calc_pos(a[2], "df", "tpdf"))),
    .ig_calc_fn("tinv", c("p", "df"), "distribution", function(a) if (a[1] == 0.5) 0 else stats::qt(a[1], a[2]),
                "Student t quantile (0 at p = 0.5 by symmetry)",
                domain = function(a) .ig_calc_need(.ig_calc_prob_open(a[1], "tinv"), .ig_calc_pos(a[2], "df", "tinv"))),
    .ig_calc_fn("nctcdf", c("x", "df", "nc"), "distribution", function(a) .ig_calc_nct_cdf(a[1], a[2], a[3]),
                "noncentral t CDF (noncentrality nc)",
                domain = function(a) .ig_calc_need(.ig_calc_pos(a[2], "df", "nctcdf"),
                                                   list(is.finite(a[3]), "nctcdf: nc must be finite")),
                range = .ig_calc_range_nc("nctcdf", 3, 2, limit = .ig_calc_range$nct_nc)),
    .ig_calc_fn("nctpdf", c("x", "df", "nc"), "distribution", function(a) .ig_calc_nct_pdf(a[1], a[2], a[3]),
                "noncentral t density (noncentrality nc)",
                domain = function(a) .ig_calc_need(.ig_calc_pos(a[2], "df", "nctpdf"),
                                                   list(is.finite(a[3]), "nctpdf: nc must be finite")),
                range = .ig_calc_range_nc("nctpdf", 3, 2, limit = .ig_calc_range$nct_nc)),
    .ig_calc_fn("nctinv", c("p", "df", "nc"), "distribution", function(a) .ig_calc_nct_inv(a[1], a[2], a[3]),
                "noncentral t quantile (noncentrality nc)",
                domain = function(a) .ig_calc_need(.ig_calc_prob_open(a[1], "nctinv"), .ig_calc_pos(a[2], "df", "nctinv"),
                                                   list(is.finite(a[3]), "nctinv: nc must be finite")),
                range = .ig_calc_range_nc("nctinv", 3, 2, limit = .ig_calc_range$nct_nc)),

    # Beta -------------------------------------------------------------------------------------
    .ig_calc_fn("betacdf", c("x", "a", "b"), "distribution", function(a) stats::pbeta(a[1], a[2], a[3]), "beta CDF",
                domain = function(a) .ig_calc_need(.ig_calc_pos(a[2], "a", "betacdf"), .ig_calc_pos(a[3], "b", "betacdf"))),
    .ig_calc_fn("betapdf", c("x", "a", "b"), "distribution", function(a) stats::dbeta(a[1], a[2], a[3]), "beta density",
                domain = function(a) .ig_calc_need(.ig_calc_pos(a[2], "a", "betapdf"), .ig_calc_pos(a[3], "b", "betapdf"))),
    .ig_calc_fn("betainv", c("p", "a", "b"), "distribution", function(a) stats::qbeta(a[1], a[2], a[3]), "beta quantile",
                domain = function(a) .ig_calc_need(.ig_calc_prob_closed(a[1], "betainv"), .ig_calc_pos(a[2], "a", "betainv"),
                                                   .ig_calc_pos(a[3], "b", "betainv"))),

    # Poisson ----------------------------------------------------------------------------------
    .ig_calc_fn("poisscdf", c("x", "lambda"), "distribution", function(a) stats::ppois(floor(a[1]), a[2]),
                "Poisson CDF P(X <= x) with mean lambda (x rounded down, as in G*Power)",
                domain = function(a) .ig_calc_need(.ig_calc_nonneg(a[2], "lambda", "poisscdf")),
                gpower_note = "captured live: poisscdf(2.5,2) = poisscdf(2,2)"),
    .ig_calc_fn("poisspdf", c("x", "lambda"), "distribution", function(a) {
      if (a[1] != floor(a[1])) 0 else stats::dpois(a[1], a[2])
    }, "Poisson probability P(X = x) with mean lambda (0 for a non-integer x, as in G*Power)",
    domain = function(a) .ig_calc_need(.ig_calc_nonneg(a[2], "lambda", "poisspdf")),
    gpower_note = "captured live: poisspdf(2.5,2) = 0"),
    .ig_calc_fn("poissinv", c("p", "lambda"), "distribution", function(a) stats::qpois(a[1], a[2]),
                "Poisson quantile: the smallest integer x with P(X <= x) >= p (Inf at p = 1)",
                domain = function(a) .ig_calc_need(.ig_calc_prob_closed(a[1], "poissinv"),
                                                   .ig_calc_nonneg(a[2], "lambda", "poissinv")),
                gpower_note = "G*Power 3.1.9.7 returns the largest double 1.79769313486232e+308 at p = 1 (captured), where the quantile is infinite"),
    .ig_calc_fn("poissmean", c("x", "lambda"), "distribution", function(a) a[2],
                "mean lambda of the Poisson distribution (x is not used)",
                domain = function(a) .ig_calc_need(.ig_calc_nonneg(a[2], "lambda", "poissmean")),
                gpower_note = "captured live: poissmean(x, lambda) returns lambda for every x (poissmean(3,2.5) = 2.5, poissmean(0.5,7) = 7)"),

    # Binomial ---------------------------------------------------------------------------------
    .ig_calc_fn("binocdf", c("x", "N", "pi"), "distribution", function(a) {
      .ig_binom_cdf_continuous(a[1], .ig_calc_round_count(a[2]), a[3])
    }, paste("binomial CDF P(X <= x) for N trials with success probability pi, continuous in x:",
             "I_{1-pi}(N - x, x + 1) for 0 <= x < N, 0 below 0, 1 from N on (N rounded half up)"),
    domain = function(a) .ig_calc_need(.ig_calc_prob_closed(a[3], "binocdf"), .ig_calc_nonneg(a[2], "N", "binocdf")),
    gpower_note = paste("captured live: binocdf(2.5,10,0.3) = 0.519251402167188 = I_0.7(7.5, 3.5) (continuous in x),",
                        "binocdf(-0.5,10,0.3) = 0, binocdf(3,10.5,0.3) = binocdf(3,11,0.3)")),
    .ig_calc_fn("binopdf", c("x", "N", "pi"), "distribution", function(a) {
      stats::dbinom(.ig_calc_round_count(a[1]), .ig_calc_round_count(a[2]), a[3])
    }, "binomial probability P(X = x) (x and N rounded half up, as in G*Power)",
    domain = function(a) .ig_calc_need(.ig_calc_prob_closed(a[3], "binopdf"), .ig_calc_nonneg(a[2], "N", "binopdf")),
    gpower_note = "captured live: binopdf(2.5,10,0.3) = binopdf(3,10,0.3), binopdf(2.2,10,0.3) = binopdf(2,10,0.3)"),
    .ig_calc_fn("binoinv", c("p", "N", "pi"), "distribution", function(a) {
      .ig_binom_inv_continuous(a[1], .ig_calc_round_count(a[2]), a[3])
    },
                paste("continuous binomial quantile: the x in [0, N] with I_{1-pi}(N - x, x + 1) = p (inverse of binocdf;",
                      "0 when p is at most the probability (1 - pi)^N of x = 0); not an integer in general"),
                domain = function(a) .ig_calc_need(.ig_calc_prob_closed(a[1], "binoinv"), .ig_calc_prob_closed(a[3], "binoinv"),
                                                   .ig_calc_nonneg(a[2], "N", "binoinv")),
                gpower_note = "captured live: binoinv(0.5,10,0.5) = 4.5, binoinv(0,10,0.5) = 0"),

    # Hypergeometric -----------------------------------------------------------------------------
    .ig_calc_fn("hygecdf", c("x", "N", "ns", "nt"), "distribution", function(a) {
      stats::phyper(.ig_calc_round_count(a[1]), a[3], a[4] - a[3], .ig_calc_round_count(a[2]))
    }, paste("hypergeometric CDF P(X <= x): sample of N from a population of nt with ns successes",
             "(x and N rounded half up)"),
    domain = function(a) .ig_calc_hyge_domain(a, "hygecdf"),
    range = .ig_calc_range_hyge("hygecdf"),
    gpower_note = "captured live: hygecdf(2.5,10,5,50) = hygecdf(3,10,5,50), hygecdf(2.2,...) = hygecdf(2,...), hygecdf(2,10.5,5,50) = hygecdf(2,11,5,50)"),
    .ig_calc_fn("hygepdf", c("x", "N", "ns", "nt"), "distribution", function(a) {
      stats::dhyper(.ig_calc_round_count(a[1]), a[3], a[4] - a[3], .ig_calc_round_count(a[2]))
    }, "hypergeometric probability P(X = x) (x and N rounded half up)",
    domain = function(a) .ig_calc_hyge_domain(a, "hygepdf"),
    range = .ig_calc_range_hyge("hygepdf")),
    .ig_calc_fn("hygeinv", c("p", "N", "ns", "nt"), "distribution", function(a) {
      stats::qhyper(a[1], a[3], a[4] - a[3], .ig_calc_round_count(a[2]))
    }, "hypergeometric quantile: the smallest integer x with P(X <= x) >= p (N rounded half up)",
    domain = function(a) .ig_calc_need(.ig_calc_prob_closed(a[1], "hygeinv"), .ig_calc_hyge_domain(a, "hygeinv")),
    range = .ig_calc_range_hyge("hygeinv", search = TRUE)),

    # Correlation and squared multiple correlation -----------------------------------------------
    .ig_calc_fn("corrcdf", c("r", "rho", "N"), "distribution", function(a) {
      if (abs(a[2]) == 1) return(as.double(a[1] >= a[2]))
      .ig_corr_cdf(a[1], a[2], a[3])
    }, paste("CDF of the sample correlation coefficient of N bivariate normal pairs with correlation rho",
             "(N >= 3; for rho = -1 or 1 the distribution is a point mass at rho)"),
    domain = function(a) .ig_corr_check(if (abs(a[2]) == 1) 0 else a[2], a[3]),
    range = .ig_calc_range_corr("corrcdf"),
    gpower_note = paste("G*Power 3.1.9.7 evaluates N < 4 as N = 4 when rho != 0 (captured: corrcdf(0.3,0.3,2) and",
                        "corrcdf(0.3,0.3,3.5) both equal corrcdf(0.3,0.3,4) = 0.45881379627718); here N is used as given",
                        "and N < 3 is refused")),
    .ig_calc_fn("corrpdf", c("r", "rho", "N"), "distribution", function(a) .ig_corr_pdf(a[1], a[2], a[3]),
                "density of the sample correlation coefficient",
                domain = function(a) .ig_corr_check(a[2], a[3]),
                range = .ig_calc_range_corr("corrpdf")),
    .ig_calc_fn("corrinv", c("p", "rho", "N"), "distribution", function(a) .ig_corr_inv(a[1], a[2], a[3]),
                "quantile of the sample correlation coefficient",
                gpower_note = "G*Power 3.1.9.7 returns 0.984375 for corrinv(1,0.3,20) (captured), where the quantile is 1",
                domain = function(a) .ig_calc_need(.ig_calc_prob_closed(a[1], "corrinv"), .ig_corr_check(a[2], a[3])),
                range = .ig_calc_range_corr("corrinv")),
    .ig_calc_fn("mr2cdf", c("R2", "rho2", "k", "N"), "distribution", function(a) .ig_mr2_cdf(a[1], a[2], a[3], a[4]),
                "CDF of the squared multiple correlation R^2 with k - 1 random (normal) predictors, N cases, population rho^2",
                domain = function(a) .ig_mr2_check(a[2], a[3], a[4])),
    .ig_calc_fn("mr2pdf", c("R2", "rho2", "k", "N"), "distribution", function(a) .ig_mr2_pdf(a[1], a[2], a[3], a[4]),
                "density of the squared multiple correlation R^2 (k - 1 random predictors)",
                domain = function(a) .ig_mr2_check(a[2], a[3], a[4])),
    .ig_calc_fn("mr2inv", c("p", "rho2", "k", "N"), "distribution", function(a) .ig_mr2_inv(a[1], a[2], a[3], a[4]),
                "quantile of the squared multiple correlation R^2 (k - 1 random predictors)",
                domain = function(a) .ig_calc_need(.ig_calc_prob_closed(a[1], "mr2inv"), .ig_mr2_check(a[2], a[3], a[4]))),

    # Log-normal, Laplace, exponential, uniform ----------------------------------------------------
    .ig_calc_fn("logncdf", c("x", "m", "s"), "distribution", function(a) stats::plnorm(a[1], a[2], a[3]),
                "log-normal CDF; m and s are the mean and SD of the associated normal (0 for x <= 0)",
                domain = function(a) .ig_calc_need(.ig_calc_finite(a[2], "m", "logncdf"), .ig_calc_pos(a[3], "s", "logncdf")),
                gpower_note = "G*Power 3.1.9.7 refuses x < 0 with 'Math error' (captured: logncdf(-1,0,1)), where the CDF is 0"),
    .ig_calc_fn("lognpdf", c("x", "m", "s"), "distribution", function(a) .ig_logn_pdf(a[1], a[2], a[3]),
                "log-normal density (0 for x <= 0)",
                domain = function(a) .ig_calc_need(.ig_calc_finite(a[2], "m", "lognpdf"), .ig_calc_pos(a[3], "s", "lognpdf"))),
    .ig_calc_fn("logninv", c("p", "m", "s"), "distribution", function(a) stats::qlnorm(a[1], a[2], a[3]),
                "log-normal quantile",
                domain = function(a) .ig_calc_need(.ig_calc_finite(a[2], "m", "logninv"), .ig_calc_pos(a[3], "s", "logninv"),
                                                   list(is.finite(a[1]) && a[1] >= 0 && a[1] < 1,
                                                        "logninv: p must lie in [0, 1) (the quantile at 1 is infinite)"))),
    .ig_calc_fn("laplcdf", c("x", "m", "s"), "distribution", function(a) .ig_lapl_cdf(a[1], a[2], a[3]),
                "Laplace CDF with location m and scale s",
                domain = function(a) .ig_calc_need(.ig_calc_finite(a[2], "m", "laplcdf"), .ig_calc_pos(a[3], "s", "laplcdf"))),
    .ig_calc_fn("laplpdf", c("x", "m", "s"), "distribution", function(a) .ig_lapl_pdf(a[1], a[2], a[3]),
                "Laplace density with location m and scale s",
                domain = function(a) .ig_calc_need(.ig_calc_finite(a[2], "m", "laplpdf"), .ig_calc_pos(a[3], "s", "laplpdf"))),
    .ig_calc_fn("laplinv", c("p", "m", "s"), "distribution", function(a) .ig_lapl_inv(a[1], a[2], a[3]),
                "Laplace quantile with location m and scale s",
                domain = function(a) .ig_calc_need(.ig_calc_prob_open(a[1], "laplinv"), .ig_calc_finite(a[2], "m", "laplinv"), .ig_calc_pos(a[3], "s", "laplinv"))),
    .ig_calc_fn("expcdf", c("x", "lambda"), "distribution", function(a) stats::pexp(a[1], a[2]),
                "exponential CDF with rate lambda (0 for x < 0)",
                domain = function(a) .ig_calc_need(.ig_calc_pos(a[2], "lambda", "expcdf")),
                gpower_note = "G*Power 3.1.9.7 computes 1 - exp(-lambda x) with cancellation (captured: expcdf(1e-10,1e-5) = 9.99200722162641e-16 instead of 1e-15)"),
    .ig_calc_fn("exppdf", c("x", "lambda"), "distribution", function(a) stats::dexp(a[1], a[2]),
                "exponential density with rate lambda",
                domain = function(a) .ig_calc_need(.ig_calc_pos(a[2], "lambda", "exppdf"))),
    .ig_calc_fn("expinv", c("p", "lambda"), "distribution", function(a) stats::qexp(a[1], a[2]),
                "exponential quantile with rate lambda",
                domain = function(a) .ig_calc_need(.ig_calc_pos(a[2], "lambda", "expinv"),
                                                   list(is.finite(a[1]) && a[1] >= 0 && a[1] < 1,
                                                        "expinv: p must lie in [0, 1) (the quantile at 1 is infinite)"))),
    .ig_calc_fn("unicdf", c("x", "a", "b"), "distribution", function(a) stats::punif(a[1], a[2], a[3]),
                "uniform CDF on [a, b]",
                domain = function(a) .ig_calc_uni_domain(a, "unicdf")),
    .ig_calc_fn("unipdf", c("x", "a", "b"), "distribution", function(a) stats::dunif(a[1], a[2], a[3]),
                "uniform density on [a, b]: 1/(b - a) inside, 0 outside",
                domain = function(a) .ig_calc_uni_domain(a, "unipdf"),
                gpower_note = paste("reference defect: G*Power 3.1.9.7 evaluates the CDF unicdf under this name",
                                    "(captured: unipdf(0.3,0,1) = 0.3, unipdf(2,0,1) = 1, unipdf(0.5,0,2) = 0.25)")),
    .ig_calc_fn("uniinv", c("p", "a", "b"), "distribution", function(a) stats::qunif(a[1], a[2], a[3]),
                "uniform quantile on [a, b]",
                domain = function(a) .ig_calc_need(.ig_calc_prob_closed(a[1], "uniinv"), .ig_calc_uni_domain(a, "uniinv")))
  )
  stats::setNames(fns, vapply(fns, `[[`, "", "name"))
}

# N is rounded half up (captured); ns and nt must be whole numbers (their
# rounding was not captured, so a fractional value is refused rather than
# guessed).
.ig_calc_hyge_domain <- function(a, fun) {
  n <- .ig_calc_round_count(a[2])
  ns <- a[3]
  nt <- a[4]
  ok <- all(is.finite(c(n, ns, nt))) && n >= 0 && ns >= 0 && nt >= 0 &&
    ns == round(ns) && nt == round(nt) && n <= nt && ns <= nt
  if (ok) return(NULL)
  sprintf("%s: ns and nt must be whole numbers, with N <= nt and ns <= nt", fun)
}

.ig_calc_uni_domain <- function(a, fun) {
  if (is.finite(a[2]) && is.finite(a[3]) && a[2] < a[3]) return(NULL)
  sprintf("%s: the interval [a, b] needs a < b", fun)
}

.ig_calc_lngamma <- function(a) {
  x <- a[1]
  if (is.na(x)) return(x)
  if (x > 0 || is.infinite(x)) return(lgamma(x))
  if (x == floor(x)) return(Inf)
  # Gamma(x) for negative non-integer x is positive when floor(x) is even.
  if (floor(x) %% 2 == 0) lgamma(x) else NaN
}

# Cached table (built once per session; the closures are package code).
.ig_calc_table_cache <- new.env(parent = emptyenv())
.ig_calc_table <- function() {
  if (is.null(.ig_calc_table_cache$table)) .ig_calc_table_cache$table <- .ig_calc_function_table()
  .ig_calc_table_cache$table
}

#' Functions of the distribution calculator
#'
#' Lists the functions that [distribution_calculator()] accepts, with their
#' arguments in the order of the 'G\*Power' 3.1.9.7 calculator, a description
#' and notes on reference behaviour where 'G\*Power' differs.
#'
#' @return A data frame with columns `name`, `arguments`, `n_args`,
#'   `category` (`function` or `distribution`), `description` and
#'   `gpower_note` (`NA` when there is nothing to note).
#' @examples
#' fns <- calculator_functions()
#' fns[fns$category == "distribution", c("name", "arguments")][1:6, ]
#' @seealso [distribution_calculator()]
#' @export
calculator_functions <- function() {
  tab <- .ig_calc_table()
  data.frame(
    name = names(tab),
    arguments = vapply(tab, function(e) paste(e$args, collapse = ", "), ""),
    n_args = vapply(tab, function(e) length(e$args), integer(1)),
    category = vapply(tab, `[[`, "", "category"),
    description = vapply(tab, `[[`, "", "description"),
    gpower_note = vapply(tab, function(e) e$gpower_note, ""),
    row.names = NULL, stringsAsFactors = FALSE
  )
}
