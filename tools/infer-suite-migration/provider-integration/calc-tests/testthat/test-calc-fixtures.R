# Distribution calculator against G*Power 3.1.9.7 (emulator_capture, Wine):
# a curated, fast subset of the fixture comparison of
# validation/compare/utilities/util_calculator (full run: summary.md there).
# Every G*Power string below is the displayed Result of the cited fixture
# record (validation/fixtures/emulator_capture/ or
# emulator_capture_supplement/g12/). Rule MR-DISPLAY: closed rounding
# interval of the last displayed digit plus 8 eps max(1, |x|)
# (validation/tolerances.md section 4); text displays are compared exactly.

g12_display_interval <- function(displayed, eps_slack = 8) {
  s <- trimws(displayed)
  if (grepl("e", s, fixed = TRUE)) {
    mant <- sub("e.*$", "", s)
    ex <- as.integer(sub("^.*e", "", s))
    d <- if (grepl(".", mant, fixed = TRUE)) nchar(sub("^-?[0-9]*\\.", "", mant)) else 0L
    half <- 0.5 * 10^(ex - d)
  } else {
    d <- if (grepl(".", s, fixed = TRUE)) nchar(sub("^-?[0-9]*\\.", "", s)) else 0L
    half <- 0.5 * 10^(-d)
  }
  x <- as.numeric(s)
  a <- eps_slack * .Machine$double.eps * max(1, abs(x))
  c(x - half - a, x + half + a)
}
expect_display <- function(expr, displayed, record, eps_slack = 8) {
  r <- distribution_calculator(expr)
  expect_identical(r$status$code, "ok", info = record)
  if (displayed %in% c("inf", "-inf", "-nan(ind)")) {
    expect_identical(r$display, displayed, info = record)
  } else {
    iv <- g12_display_interval(displayed, eps_slack)
    expect_true(r$value >= iv[1] && r$value <= iv[2],
                info = sprintf("%s: %s = %s, G*Power %s", record, expr, format(r$value, digits = 17), displayed))
  }
}

test_that("operators, elementary functions and scripts agree with the G*Power display", {
  cases <- list(
    c("(-0.0667*7.5/4)^2", "0.01564062890625", "GPW-UT-CALC-UTL-U001-UO02--r01"),
    c("2/3", "0.666666666666667", "GPW-UT-CALC-UTL-U003-UO02--r01"),
    c("1/7e-310", "inf", "GPW-UT-CALC-UTL-U003-UD01--r01"),
    c("sin(1e10)", "-0.487506025087511", "GPW-UT-CALC-UTL-U007-UD01--r01"),
    c("acos(0.9999999999)", "1.41421362089116e-05", "GPW-UT-CALC-UTL-U010-UD01--r01"),
    c("tan(1.5707963)", "37320539.6343548", "GPW-UT-CALC-UTL-U011-UD01--r01"),
    c("log(1e-300)", "-690.775527898214", "GPW-UT-CALC-UTL-U015-UD01--r01"),
    c("lngamma(1e10)", "220258509288.811", "GPW-UT-CALC-UTL-U019-UD01--r01"),
    c("frac(1e15+0.5)", "0.5", "GPW-UT-CALC-UTL-U020-UD01--r01"),
    c("uround(7.1,5)", "10", "GPW-UT-CALC-UTL-U024-UO03--r01"),
    c("uround(2.3,0.1)", "-nan(ind)", "GPW-UT-CALC-UTL-U024-UD01--r01"),
    c("z = (0.509-0.397)/0.0788; 1-normcdf(z,0,1)", "0.0776119053969239", "GPW-UT-CALC-UTL-U028-UO01--r01"),
    c("alpha1 = 0.05; tinv(1-alpha1,10)", "1.81246112281168", "GPW-UT-CALC-UTL-U082-UD01--r01"),
    c("x = 0.25\ny = x*4\ny", "1", "GPW-UT-CALC-UTL-U085-UO01--r01"),
    c("sqrt(-1)", "-nan(ind)", "GPW-UT-CALC-UTL-U086-I06--r01")
  )
  for (k in cases) expect_display(k[1], k[2], k[3])
})

test_that("distribution functions agree with the G*Power display", {
  cases <- list(
    c("zcdf(-37)", "5.72557122252458e-300", "GPW-UT-CALC-UTL-U025-UD01--r01"),
    c("zpdf(38)", "1.09722105199497e-314", "GPW-UT-CALC-UTL-U026-UD01--r01"),
    c("zinv(1e-300)", "-37.0470962993612", "GPW-UT-CALC-UTL-U027-UD01--r01"),
    c("norminv(1e-12,0,1)", "-7.03448382530113", "GPW-UT-CALC-UTL-U030-UD01--r01"),
    c("chi2cdf(100000,100000)", "0.500594708104793", "GPW-UT-CALC-UTL-U031-UD01--r01"),
    c("chi2inv(1e-10,1000)", "741.268071719353", "GPW-UT-CALC-UTL-U033-UD01--r01"),
    # df = 1e6 makes this evaluation ill-conditioned: at 50 digits the value is
    # 0.69137445042890253; G*Power and R on Linux are 1.5e-14 above it, R on
    # macOS 2.5e-14. This record alone gets 64 eps instead of MR-DISPLAY's 8.
    c("fcdf(1.001,1e6,1e6)", "0.691374450428913", "GPW-UT-CALC-UTL-U034-UD01--r01", "64"),
    c("finv(0.999999,1,1)", "405284734545.376", "GPW-UT-CALC-UTL-U036-UD01--r01"),
    c("tcdf(-50,3)", "8.80857602063598e-06", "GPW-UT-CALC-UTL-U037-UD01--r01"),
    c("tinv(0.9999999,1)", "3183098.86351325", "GPW-UT-CALC-UTL-U039-UD01--r01"),
    c("ncx2cdf(10,3,5)", "0.706648647777452", "GPW-UT-CALC-UTL-U040-UO01--r01"),
    c("ncfcdf(3,2,30,10)", "0.191313461482448", "GPW-UT-CALC-UTL-U043-UO01--r01"),
    c("ncfpdf(2,3,20,5)", "0.240123085349335", "GPW-UT-CALC-UTL-U044-UO01--r01"),
    c("nctcdf(2,10,1.5)", "0.65915407244219", "GPW-UT-CALC-UTL-U046-UO01--r01"),
    c("nctpdf(1,10,1.5)", "0.341074931926085", "GPW-UT-CALC-UTL-U047-UO01--r01"),
    c("nctinv(0.8,10,1.5)", "2.52752212662531", "GPW-UT-CALC-UTL-U048-UO01--r01"),
    c("betacdf(0.3,2,5)", "0.579825", "GPW-UT-CALC-UTL-U049-UO01--r01"),
    c("betapdf(0.001,0.5,0.5)", "10.0708791199471", "GPW-UT-CALC-UTL-U050-UD01--r01"),
    c("poisspdf(500,1000)", "4.16015547553561e-69", "GPW-UT-CALC-UTL-U053-UD01--r01"),
    c("poissinv(0.999999,1e5)", "101507", "GPW-UT-CALC-UTL-U054-UD01--r01"),
    c("binocdf(5000,10000,0.5)", "0.503989323069691", "GPW-UT-CALC-UTL-U056-UD01--r01"),
    c("binoinv(1e-9,10000,0.5)", "4699.69708710705", "GPW-UT-CALC-UTL-U058-UD01--r01"),
    c("hygecdf(500,1000,5000,10000)", "0.513294714980647", "GPW-UT-CALC-UTL-U059-UD01--r01"),
    c("hygeinv(0.999999,1000,5000,10000)", "571", "GPW-UT-CALC-UTL-U061-UD01--r01"),
    c("corrcdf(0.3,0,50)", "0.982856909983535", "GPW-UT-CALC-UTL-U062-UO01--r01"),
    c("corrpdf(0.99,0.95,500)", "5.3078224384403e-63", "GPW-UT-CALC-UTL-U063-UD01--r01"),
    c("corrcdf(-0.2,0.3,20)", "0.0163603305820813", "SUP-G12-UT-CALC-D21--r01"),
    c("corrpdf(-0.2,0.3,20)", "0.167203987735204", "SUP-G12-UT-CALC-D37--r01"),
    c("mr2cdf(0.2,0.1,4,50)", "0.735455057097933", "GPW-UT-CALC-UTL-U065-UO02--r01"),
    c("mr2cdf(0.3,0.2,4,30)", "0.598278054073064", "SUP-G12-UT-CALC-D22--r01"),
    c("mr2pdf(0.5,0,2,10)", "0.193349510480697", "GPW-UT-CALC-UTL-U066-UE01--r01"),
    c("logninv(1e-12,0,1)", "0.00088097278344686", "GPW-UT-CALC-UTL-U070-UD01--r01"),
    c("laplinv(1e-300,0,1)", "-690.082380717654", "GPW-UT-CALC-UTL-U073-UD01--r01"),
    c("expinv(0.999999999,0.001)", "20723.2658652283", "GPW-UT-CALC-UTL-U076-UD01--r01"),
    c("binocdf(2.5,10,0.3)", "0.519251402167188", "SUP-G12-UT-CALC-D01--r01"),
    c("hygecdf(2,10.5,5,50)", "0.936001717985992", "SUP-G12-UT-CALC-B53--r01"),
    c("tcdf(2,10.5)", "0.963983099921516", "SUP-G12-UT-CALC-D11--r01")
  )
  for (k in cases) expect_display(k[1], k[2], k[3], if (length(k) > 3L) as.numeric(k[4]) else 8)
})

test_that("documented G*Power defects give the documented function, not the G*Power value", {
  # Each pair: the G*Power display of the record and the property the engine
  # satisfies instead (triage: validation/compare/utilities/util_calculator/triage.R).
  r <- distribution_calculator("chi2pdf(0.5,2)")               # SUP-G12-UT-CALC-D19: 1.38629436111989 = chi2inv(0.5,2)
  expect_false(abs(r$value - 1.38629436111989) < 1e-6)
  expect_equal(r$value, exp(-0.25) / 2, tolerance = 1e-15)
  r <- distribution_calculator("unipdf(0.3,0,1)")              # GPW-UT-CALC-UTL-U078-UO01: 0.3 (the CDF)
  expect_identical(r$value, 1)
  r <- distribution_calculator("finv(0.5,1,1)")                # GPW-UT-CALC-UTL-U036-UE01: 0.999999999999997
  expect_lt(abs(r$value - 1), 1e-15)
  r <- distribution_calculator("betainv(1e-10,0.1,10)")        # GPW-UT-CALC-UTL-U051-UD01: 5e-51 (inversion failure)
  expect_equal(r$value, 6.35634284678601199211251961922e-102, tolerance = 1e-12)  # mpmath
  r <- distribution_calculator("1e+3")                         # SUP-G12-UT-CALC-F05: 1003
  expect_identical(r$value, 1000)
})
