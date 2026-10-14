# Distribution calculator: the function table (R/calc-functions.R,
# R/calc-dist.R). Group G12; method notes specs/method-notes/utilities.md.
#
# Evidence kinds kept apart:
# * closed forms and identities computed here with independent arithmetic;
# * reference values of the independent mpmath checker
#   (validation/independent/g12/calc_reference.py, 30 significant digits;
#   shares no routine with the engine, D18), cited as "mpmath";
# * G*Power conventions captured live (record ids cited); G*Power numbers
#   themselves are compared in test-calc-fixtures.R.

cv <- function(expr) {
  r <- distribution_calculator(expr)
  expect_identical(r$status$code, "ok", info = paste(expr, r$status$reason))
  r$value
}
refuses <- function(expr) {
  r <- distribution_calculator(expr)
  expect_identical(r$status$code, "invalid_input", info = expr)
  expect_identical(r$status$details$kind, "math", info = expr)
  expect_true(is.na(r$value), info = expr)
}
close_rel <- function(x, ref, rel = 1e-13) expect_lte(abs(x - ref), rel * max(abs(ref), 1e-300))

test_that("the table holds every function of the G*Power 3.1.9.7 calculator with its arity", {
  # specs/gpower-utilities.json util_calculator u006-u079 (certified live 2026-10-05).
  expected <- c(abs = 1, sin = 1, asin = 1, cos = 1, acos = 1, tan = 1, atan = 1, atan2 = 2, exp = 1, log = 1,
                sqrt = 1, sqr = 1, sign = 1, lngamma = 1, frac = 1, int = 1, min = 2, max = 2, uround = 2,
                zcdf = 1, zpdf = 1, zinv = 1, normcdf = 3, normpdf = 3, norminv = 3, chi2cdf = 2, chi2pdf = 2,
                chi2inv = 2, fcdf = 3, fpdf = 3, finv = 3, tcdf = 2, tpdf = 2, tinv = 2, ncx2cdf = 3, ncx2pdf = 3,
                ncx2inv = 3, ncfcdf = 4, ncfpdf = 4, ncfinv = 4, nctcdf = 3, nctpdf = 3, nctinv = 3, betacdf = 3,
                betapdf = 3, betainv = 3, poisscdf = 2, poisspdf = 2, poissinv = 2, poissmean = 2, binocdf = 3,
                binopdf = 3, binoinv = 3, hygecdf = 4, hygepdf = 4, hygeinv = 4, corrcdf = 3, corrpdf = 3,
                corrinv = 3, mr2cdf = 4, mr2pdf = 4, mr2inv = 4, logncdf = 3, lognpdf = 3, logninv = 3,
                laplcdf = 3, laplpdf = 3, laplinv = 3, expcdf = 2, exppdf = 2, expinv = 2, unicdf = 3,
                unipdf = 3, uniinv = 3)
  fns <- calculator_functions()
  expect_setequal(fns$name, names(expected))
  expect_identical(stats::setNames(fns$n_args, fns$name)[names(expected)], stats::setNames(as.integer(expected), names(expected)))
  expect_true(all(nzchar(fns$description)))
  for (f in c("chi2pdf", "ncx2pdf", "unipdf", "uround", "binoinv", "poissmean")) {
    expect_false(is.na(fns$gpower_note[fns$name == f]), info = f)
  }
})

test_that("elementary functions and the captured G*Power conventions", {
  expect_identical(cv("frac(1.56)"), 1.56 - 1)                # DOC-024
  expect_identical(cv("int(1.56)"), 1)                        # DOC-025
  expect_identical(cv("frac(-1.56)"), -1.56 + 1)              # U020-UE01
  expect_identical(cv("int(-1.56)"), -1)                      # U021-UE01
  expect_identical(cv("atan2(1,2)"), atan(1 / 2))             # U013-UO01: atan2(y, x)
  expect_identical(cv("atan2(0,-1)"), pi)                     # U013-UD01
  expect_identical(cv("sqr(-3)"), 9)
  expect_identical(cv("sign(-1e-300)"), -1)
  # uround(x, m) = ceiling(x / int(m)) * int(m) (SUP-G12-UT-CALC-U01..U07, DOC-026/027):
  expect_identical(cv("uround(2.3, 1)"), 3)
  expect_identical(cv("uround(2.3, 2)"), 4)
  expect_identical(cv("uround(4,2)"), 4)
  expect_identical(cv("uround(-2.3,1)"), -2)
  expect_identical(cv("uround(2.3,2.5)"), 4)
  expect_identical(cv("uround(5,2.5)"), 6)
  expect_identical(cv("uround(2.3,-1)"), 2)
  expect_true(is.nan(cv("uround(2.3,0.5)")))
  expect_true(is.nan(cv("uround(2.3,0)")))
  # lngamma: log Gamma(x); Inf at the poles; NaN where Gamma(x) < 0.
  expect_identical(cv("lngamma(5)"), log(24))
  expect_identical(cv("lngamma(0)"), Inf)                     # SUP-G12-UT-CALC-P30
  close_rel(cv("lngamma(-1.5)"), log(4 * sqrt(pi) / 3))       # Gamma(-1.5) = 4 sqrt(pi)/3 > 0
  expect_true(is.nan(cv("lngamma(-0.5)")))                    # Gamma(-0.5) < 0
  expect_identical(cv("poissmean(3,2.5)"), 2.5)               # U055-UO01: returns lambda
  expect_identical(cv("poissmean(0.5,7)"), 7)                 # SUP-G12-UT-CALC-D31
})

test_that("closed forms of the central distributions", {
  close_rel(cv("chi2cdf(3,2)"), 1 - exp(-1.5))
  close_rel(cv("chi2pdf(3,2)"), exp(-1.5) / 2)
  close_rel(cv("chi2inv(0.5,2)"), 2 * log(2))
  close_rel(cv("tcdf(2,1)"), 0.5 + atan(2) / pi)
  close_rel(cv("tpdf(0,1)"), 1 / pi)
  close_rel(cv("tinv(0.75,1)"), 1)
  close_rel(cv("fcdf(1,2,2)"), 0.5)
  close_rel(cv("fpdf(0,2,20)"), 1)                            # density at 0 for df1 = 2
  close_rel(cv("finv(0.5,1,1)"), 1, 1e-15)                    # the median of F(1, 1) is 1
  close_rel(cv("betacdf(0.3,2,5)"), 1 - 0.7^6 - 6 * 0.3 * 0.7^5)
  close_rel(cv("betapdf(0.3,2,5)"), 30 * 0.3 * 0.7^4)
  close_rel(cv("normcdf(5,5,2)"), 0.5)
  close_rel(cv("norminv(0.5,5,2)"), 5)
  close_rel(cv("expcdf(1,2)"), 1 - exp(-2))
  close_rel(cv("expcdf(1e-10,1e-5)"), -expm1(-1e-15))         # no cancellation (G*Power: U074-UD01)
  close_rel(cv("expinv(0.5,2)"), log(2) / 2)
  close_rel(cv("laplcdf(1,0,1)"), 1 - exp(-1) / 2)
  close_rel(cv("laplcdf(-700,0,1)"), exp(-700) / 2)
  close_rel(cv("laplinv(1e-300,0,1)"), log(2e-300))
  close_rel(cv("laplpdf(1,0,1)"), exp(-1) / 2)
  close_rel(cv("logncdf(1,0,1)"), 0.5)
  close_rel(cv("lognpdf(1,0,1)"), 1 / sqrt(2 * pi))
  close_rel(cv("logninv(0.975,1,0.5)"), exp(1 + 0.5 * qnorm(0.975)))
  expect_identical(cv("unicdf(0.3,0,2)"), 0.15)
  expect_identical(cv("unipdf(0.3,0,2)"), 0.5)                # the density (G*Power returns the CDF: D26, U078)
  expect_identical(cv("unipdf(2,0,1)"), 0)
  expect_identical(cv("uniinv(0.3,0,2)"), 0.6)
})

test_that("discrete distributions and their captured argument conventions", {
  ref_pmf <- function(k, n, p) choose(n, k) * p^k * (1 - p)^(n - k)
  close_rel(cv("binocdf(3,10,0.3)"), sum(ref_pmf(0:3, 10, 0.3)))
  close_rel(cv("binopdf(3,10,0.3)"), ref_pmf(3, 10, 0.3))
  # binocdf is continuous in x: I_{1-pi}(N - x, x + 1) (SUP-G12-UT-CALC-D01) ...
  close_rel(cv("binocdf(2.5,10,0.3)"), pbeta(0.7, 7.5, 3.5))
  expect_identical(cv("binocdf(-0.5,10,0.3)"), 0)             # SUP-G12-UT-CALC-R09
  # ... while N is rounded half up (SUP-G12-UT-CALC-B54) and the mass rounds x (R01, R02, D02):
  close_rel(cv("binocdf(3,10.5,0.3)"), sum(ref_pmf(0:3, 11, 0.3)))
  close_rel(cv("binopdf(2.2,10,0.3)"), ref_pmf(2, 10, 0.3))
  close_rel(cv("binopdf(2.5,10,0.3)"), ref_pmf(3, 10, 0.3))
  # binoinv inverts the continuous CDF on [0, N] (U058-UO01, SUP-G12-UT-CALC-D15, B31).
  x <- cv("binoinv(0.5,10,0.5)")
  close_rel(x, 4.5, 1e-12)
  x <- cv("binoinv(0.3,20,0.4)")
  close_rel(pbeta(0.6, 20 - x, x + 1), 0.3, 1e-14)
  expect_identical(cv("binoinv(0,10,0.5)"), 0)
  expect_identical(cv("binoinv(1e-9,10,0.5)"), 0)
  expect_identical(cv("binoinv(1,10,0.5)"), 10)
  # Poisson: CDF rounds x down (D03, R06); the mass is 0 for a non-integer x (D04).
  close_rel(cv("poisscdf(3,2.5)"), sum(exp(-2.5) * 2.5^(0:3) / factorial(0:3)))
  close_rel(cv("poisscdf(2.7,2)"), sum(exp(-2) * 2^(0:2) / factorial(0:2)))
  expect_identical(cv("poisspdf(2.5,2)"), 0)
  expect_identical(cv("poissinv(0.5,2.5)"), 2)
  expect_identical(cv("poissinv(1,2.5)"), Inf)                # G*Power: 1.79769313486232e+308 (D16)
  # Hypergeometric hygecdf(x, N, ns, nt); x and N rounded half up (D05, D06, R03-R05, B53).
  ref_h <- function(k, n, ns, nt) choose(ns, k) * choose(nt - ns, n - k) / choose(nt, n)
  close_rel(cv("hygecdf(2,10,5,50)"), sum(ref_h(0:2, 10, 5, 50)))
  close_rel(cv("hygecdf(2.5,10,5,50)"), sum(ref_h(0:3, 10, 5, 50)))
  close_rel(cv("hygecdf(2.2,10,5,50)"), sum(ref_h(0:2, 10, 5, 50)))
  close_rel(cv("hygecdf(2,10.5,5,50)"), sum(ref_h(0:2, 11, 5, 50)))
  close_rel(cv("hygepdf(2,10,5,50)"), ref_h(2, 10, 5, 50))
  expect_identical(cv("hygeinv(0.5,10,5,50)"), 1)
})

test_that("sample correlation distribution", {
  # rho = 0: Student t with N - 2 df.
  r <- 0.3
  t <- r * sqrt(48) / sqrt(1 - r^2)
  close_rel(cv("corrcdf(0.3,0,50)"), pt(t, 48))
  close_rel(cv("corrcdf(0.3,0,3)"), 0.5 + asin(0.3) / pi)     # N = 3: arcsine law
  close_rel(cv("corrpdf(0,0,4)"), 0.5)
  # rho != 0: mpmath reference values (30 digits).
  close_rel(cv("corrcdf(-0.2,0.3,20)"), 0.0163603305820813122594377738518, 1e-12)
  close_rel(cv("corrcdf(0.6,0.6,1e5)"), 0.499621526166630435848090326164, 1e-13)
  close_rel(cv("corrcdf(-0.9,0.9,50)"), 1.91055700451306965323645093267e-49, 1e-9)
  close_rel(cv("corrcdf(0.3,0.3,3)"), 0.445069676907227278845261382102, 1e-12)
  close_rel(cv("corrinv(0.5,0.5,10)"), 0.523673585072965421156871485088, 1e-12)
  # symmetry: F(r; rho) = 1 - F(-r; -rho); density integrates the CDF.
  close_rel(cv("corrcdf(0.25,0.4,30)") + cv("corrcdf(-0.25,-0.4,30)"), 1, 1e-12)
  h <- 1e-5
  close_rel((cv("corrcdf(0.30001,0.4,30)") - cv("corrcdf(0.29999,0.4,30)")) / (2 * h), cv("corrpdf(0.3,0.4,30)"), 1e-7)
  expect_identical(cv("corrcdf(0.3,1,20)"), 0)                # point mass at rho = 1 (G*Power: SUP-G12-UT-CALC-B21)
  expect_identical(cv("corrinv(1,0.3,20)"), 1)                # G*Power: 0.984375 (SUP-G12-UT-CALC-B09)
  refuses("corrcdf(0.3,0.3,2)")                                # N >= 3 (G*Power uses N = 4: B22)
  refuses("corrcdf(0.3,1.2,20)")
})

test_that("squared multiple correlation distribution (k - 1 random predictors)", {
  # rho^2 = 0: Beta((k-1)/2, (N-k)/2) (captured null case U066-UE01: mr2pdf(0.5,0,2,10) = dbeta(0.5, 1/2, 4)).
  close_rel(cv("mr2pdf(0.5,0,2,10)"), dbeta(0.5, 0.5, 4))
  close_rel(cv("mr2cdf(0.3,0,6,40)"), pbeta(0.3, 2.5, 17))
  close_rel(cv("mr2inv(0.95,0,6,100)"), qbeta(0.95, 2.5, 47))
  # rho^2 > 0: mpmath reference values (Fisher's 2F1 density, tanh-sinh quadrature).
  close_rel(cv("mr2cdf(0.2,0.1,4,50)"), 0.735455057097931921471059064661, 1e-13)
  close_rel(cv("mr2pdf(0.2,0.1,4,50)"), 3.28706002558171716972019004534, 1e-13)
  close_rel(cv("mr2inv(0.8,0.2,4,30)"), 0.384976187040376850486674502227, 1e-12)
  refuses("mr2cdf(0.3,0.2,1,30)")                              # k >= 2 (G*Power: Math error, B23)
  refuses("mr2cdf(0.3,0.2,30,30)")                             # N > k (G*Power shows 0: B24)
})

test_that("noncentral distributions agree with the independent mpmath reference", {
  close_rel(cv("nctcdf(2,10,1.5)"), 0.659154072442190815401477895068, 1e-13)
  close_rel(cv("nctpdf(1,10,1.5)"), 0.341074931926083828672953199296, 1e-13)
  close_rel(cv("nctpdf(-5,3,10)"), 2.75514744207997702215242084416e-28, 1e-9)
  close_rel(cv("finv(0.95,3,20)"), 3.0983912121407791534055844585, 1e-14)
  expect_identical(cv("nctcdf(2,10,0)"), pt(2, 10))
  x <- cv("nctinv(0.8,10,1.5)")
  close_rel(cv(sprintf("nctcdf(%.17g,10,1.5)", x)), 0.8, 1e-14)
  x <- cv("ncx2inv(0.8,3,5)")
  close_rel(cv(sprintf("ncx2cdf(%.17g,3,5)", x)), 0.8, 1e-14)
  x <- cv("ncfinv(0.8,3,20,5)")
  close_rel(cv(sprintf("ncfcdf(%.17g,3,20,5)", x)), 0.8, 1e-14)
  # The densities are the documented ones (G*Power evaluates chi2inv / ncx2inv under these names).
  close_rel(cv("chi2pdf(0.5,2)"), exp(-0.25) / 2)              # G*Power: 1.38629436111989 (SUP-G12-UT-CALC-D19)
  close_rel(cv("ncx2pdf(5,3,0)"), dchisq(5, 3))                # G*Power: Math error (U041-UE01)
})

test_that("parameters outside their domain are math errors; values outside the support are defined", {
  refuses("normcdf(1,0,-1)")                                   # SUP-G12-UT-CALC-D09
  refuses("normpdf(0,0,0)")                                    # G*Power returns 1 (B29)
  refuses("tcdf(1,0)")                                         # D10
  refuses("tcdf(1,-1)")                                        # B13
  refuses("zinv(0)")                                           # D12
  refuses("zinv(1)")                                           # D13
  refuses("tinv(1,10)")                                        # D14
  refuses("finv(1,3,10)")                                      # B02
  refuses("norminv(1.5,0,1)")                                  # B08
  refuses("uniinv(1.5,0,1)")                                   # B07
  refuses("binocdf(3,10,1.2)")                                 # B12
  refuses("chi2cdf(1,0)")                                      # B39
  refuses("ncfcdf(2,3,10,-1)")                                 # B34
  refuses("ncx2cdf(1,2,-1)")                                   # G*Power treats nc < 0 as 0 (B18)
  refuses("hygecdf(2,10,5,8)")                                 # B25
  refuses("unicdf(0.5,1,0)")                                   # B26
  refuses("expcdf(1,-1)")                                      # B27
  refuses("laplcdf(0,0,0)")                                    # B28
  refuses("laplinv(0,0,1)")                                    # G*Power: -inf (B06)
  refuses("expinv(1,2)")                                       # G*Power: nan (B04)
  expect_identical(cv("chi2cdf(-1,2)"), 0)                     # G*Power: Math error (B14)
  expect_identical(cv("fcdf(-1,2,3)"), 0)                      # G*Power: Math error (B15)
  expect_identical(cv("logncdf(-1,0,1)"), 0)                   # G*Power: Math error (D29)
  expect_identical(cv("lognpdf(0,0,1)"), 0)                    # G*Power: inf (B36)
  expect_identical(cv("betacdf(1.5,2,3)"), 1)                  # G*Power: 3.24e+282 (B16)
  expect_identical(cv("expcdf(-1,2)"), 0)                      # D30
  expect_identical(cv("binocdf(11,10,0.3)"), 1)                # D08
  expect_identical(cv("chi2inv(0,3)"), 0)                      # B01
  expect_identical(cv("betainv(1,2,3)"), 1)                    # B03
  expect_identical(cv("betainv(0,2,3)"), 0)                    # G*Power: 5e-51 (B41)
  expect_identical(cv("logninv(0,0,1)"), 0)                    # G*Power: 1 (B05)
  expect_identical(cv("fpdf(0,1,5)"), Inf)                     # G*Power: 0 (D33)
  expect_identical(cv("tpdf(0,0.5)"), dt(0, 0.5))              # df < 1 allowed (B51)
})

test_that("combined domain rules give math errors, not R errors", {
  # A probability check followed by another domain helper (regression: the
  # second helper's reason was read as a condition pair).
  refuses("uniinv(0.5,1,0)")                                   # interval needs a < b (unicdf: B26)
  refuses("hygeinv(0.5,10,5,8)")                               # N > nt (hygecdf: B25)
  refuses("corrinv(0.5,1.2,20)")
  refuses("mr2inv(0.5,0.2,1,30)")
  refuses("normcdf(1,1/0,1)")                                  # an infinite location has no distribution
  refuses("laplcdf(1/0,1/0,1)")
  refuses("zcdf(0/0)")                                         # NaN argument of a distribution function
  expect_true(is.nan(cv("max(1,0/0)")))                        # elementary functions propagate NaN
  r <- distribution_calculator("uniinv(0.5,1,0)")
  expect_match(r$status$reason, "^math error: uniinv: ")
  r <- distribution_calculator("mr2cdf(0.3,0.2,1,30)")
  expect_match(r$status$reason, "^math error: mr2cdf: ")       # the function is always named
})

test_that("parameters beyond the computable range are typed as unsupported", {
  unsupported <- function(expr) {
    r <- distribution_calculator(expr)
    expect_identical(r$status$code, "unsupported", info = expr)
    expect_true(is.na(r$value), info = expr)
    expect_match(r$status$reason, "^outside the computable range: ", info = expr)
  }
  rg <- inferguidR:::.ig_calc_range
  expect_identical(rg$nc, 1e8)
  unsupported("ncx2cdf(1,2,1e9)")
  unsupported("ncx2pdf(1,2,1e300)")
  unsupported("ncx2inv(0.5,1e10,3)")
  unsupported("ncfcdf(1,2,3,1e9)")
  unsupported("ncfinv(0.5,2,1e10,3)")
  unsupported("nctcdf(1,5,1e7)")
  unsupported("nctpdf(1,5,-1e7)")
  unsupported("corrcdf(0.3,0.3,1e8)")
  unsupported("corrinv(0.5,0.3,1e300)")
  unsupported("hygecdf(1,10,5,1e10)")
  unsupported("hygeinv(0.5,1e8,1e8,1e9)")                      # a linear search over 1e8 values
  unsupported("mr2cdf(0.5,0.999999,5,100)")                    # mixture with more than 5e5 terms
  # Inside the range the same functions compute.
  expect_identical(distribution_calculator("ncx2cdf(1e8,5,1e8)")$status$code, "ok")
  expect_identical(distribution_calculator("corrcdf(0.3,0.3,1e7)")$status$code, "ok")
  expect_identical(cv("hygeinv(0.5,1e7,1e8,1e9)"), 1e6)
})

test_that("tails below the double range are 0 or 1 without a long mixture", {
  expect_identical(cv("ncx2cdf(1,5,1e8)"), 0)
  expect_identical(cv("ncx2cdf(3e8,5,1e8)"), 1)
  expect_identical(cv("ncx2pdf(1e300,0.999999,1.5)"), 0)
  expect_identical(cv("ncx2pdf(1,5,1e8)"), 0)
  expect_identical(cv("ncfpdf(1/0,2,3,1)"), 0)
  expect_identical(cv("lognpdf(1e-300,0.999999,1e-300)"), 0)   # nmath's direct formula is 0/0 here
  # The Chernoff bound that decides this is a bound: where it does not apply the mixture is summed.
  b <- inferguidR:::.ig_calc_ncx2_log_bound
  for (x in c(30, 60, 200)) {
    expect_gte(b(x, 4, 20, upper = TRUE), pchisq(x, 4, 20, lower.tail = FALSE, log.p = TRUE) - 1e-9)
  }
  for (x in c(0.5, 3, 10)) {
    expect_gte(b(x, 4, 20, upper = FALSE), pchisq(x, 4, 20, log.p = TRUE) - 1e-9)
  }
  expect_identical(b(24, 4, 20, upper = TRUE), 0)             # at the mean the bound is trivial
})

test_that("the hypergeometric function of the correlation density has its closed forms", {
  h <- function(c, x) exp(inferguidR:::.ig_log_hyp2f1_half(c, x))
  # 2F1(1/2, 1/2; 3/2; x) = asin(sqrt(x)) / sqrt(x);
  # 2F1(1/2, 1/2; 5/2; x) = 3/(4x) ((2x - 1) asin(sqrt(x)) / sqrt(x) + sqrt(1 - x))  (N = 3).
  f15 <- function(x) asin(sqrt(x)) / sqrt(x)
  f25 <- function(x) 3 / (4 * x) * ((2 * x - 1) * asin(sqrt(x)) / sqrt(x) + sqrt(1 - x))
  for (x in c(0.01, 0.3, 0.75, 0.9, 0.9000001, 0.95, 0.999, 0.9999995)) {
    close_rel(h(1.5, x), f15(x), 1e-14)
    close_rel(h(2.5, x), f25(x), 1e-14)
  }
  expect_identical(h(2.5, 0), 1)
  # The series and the connection formula in 1 - x agree where both converge.
  s <- function(c, x) inferguidR:::.ig_hyp2f1_series(0.5, 0.5, c, x)
  for (c in c(2.5, 3.5, 4.4, 7.5)) close_rel(h(c, 0.93), s(c, 0.93), 1e-14)
  # A series that cannot converge within its term limit gives NaN, never a wrong number.
  expect_true(is.nan(inferguidR:::.ig_hyp2f1_series(0.5, 0.5, 3, 1 - 1e-9, max_terms = 1e4)))
})

test_that("every function returns a typed outcome for edge-case arguments", {
  # A deterministic sample of the checks in validation/independent/g12/fuzz_calculator.R.
  vals <- c("0", "1", "-1", "0.5", "2", "3.5", "1e5", "1e300", "-1e300", "1e-300", "(0/0)", "(1/0)", "(-1/0)", "0.999999", "-0")
  fns <- calculator_functions()
  set.seed(20261006)
  for (i in seq_len(nrow(fns))) {
    for (j in seq_len(6L)) {
      e <- sprintf("%s(%s)", fns$name[i], paste(sample(vals, fns$n_args[i], replace = TRUE), collapse = ","))
      r <- distribution_calculator(e)
      expect_true(r$status$code %in% c("ok", "invalid_input", "unsupported", "numerical_failure"), info = e)
      if (identical(r$status$code, "ok") && fns$category[i] == "distribution") {
        expect_false(is.nan(r$value), info = e)
        if (grepl("cdf$", fns$name[i])) expect_true(r$value >= 0 && r$value <= 1, info = e)
        if (grepl("pdf$", fns$name[i])) expect_true(r$value >= 0, info = e)
      } else if (!identical(r$status$code, "ok")) {
        expect_true(is.na(r$value), info = e)
      }
    }
  }
})
