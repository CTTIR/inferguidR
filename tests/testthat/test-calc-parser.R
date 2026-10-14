# Distribution calculator: syntax, evaluation and safety (R/calc-parse.R,
# R/calc-eval.R). Group G12; method notes specs/method-notes/utilities.md.
#
# The syntax conventions were captured from G*Power 3.1.9.7 under Wine
# (emulator_capture); the record ids are cited. The expected values here are
# exact arithmetic facts (64, 4, 0.5, ...), not G*Power numbers.

calc <- function(expr) distribution_calculator(expr)
calc_value <- function(expr) {
  r <- calc(expr)
  expect_identical(r$status$code, "ok", info = expr)
  r$value
}
calc_refuses <- function(expr, kind) {
  r <- calc(expr)
  expect_identical(r$status$code, "invalid_input", info = expr)
  expect_identical(r$status$details$kind, kind, info = expr)
  expect_true(is.na(r$value), info = expr)
  expect_true(is.na(r$display), info = expr)
  invisible(r)
}

test_that("operators follow the captured G*Power precedence and associativity", {
  expect_identical(calc_value("2^3"), 8)                      # GPW-UT-CALC-UTL-U001-UO01
  expect_identical(calc_value("2^3^2"), 64)                   # U001-UO03: '^' is left-associative
  expect_identical(calc_value("2*3^2"), 18)                   # SUP-G12-UT-CALC-P01
  expect_identical(calc_value("-2^2"), 4)                     # inventory u001 note: unary minus binds tighter
  expect_identical(calc_value("-3^2"), 9)                     # SUP-G12-UT-CALC-P03
  expect_identical(calc_value("2*-3^2"), 18)                  # SUP-G12-UT-CALC-P02
  expect_identical(calc_value("2^-1"), 0.5)                   # U001-UE01
  expect_identical(calc_value("2^-1^2"), 0.25)                # SUP-G12-UT-CALC-P04
  expect_identical(calc_value("1-2-3"), -4)                   # SUP-G12-UT-CALC-P05
  expect_identical(calc_value("8/2/2"), 2)                    # SUP-G12-UT-CALC-P06
  expect_identical(calc_value("2--3"), 5)                     # SUP-G12-UT-CALC-P08
  expect_identical(calc_value("2+-3"), -1)                    # SUP-G12-UT-CALC-P09
  expect_identical(calc_value("-3*-2"), 6)                    # U002-UO02
  expect_identical(calc_value("-(-3)"), 3)                    # U005-UO02
  expect_identical(calc_value("(2+3)*4"), 20)
})

test_that("number literals are read as in G*Power", {
  expect_identical(calc_value(".5+5."), 5.5)                  # SUP-G12-UT-CALC-P10
  expect_identical(calc_value("1e3+1E3+2.5e-3"), 2000.0025)
  expect_identical(calc_value("1e-3"), 0.001)                 # SUP-G12-UT-CALC-F06
  # An explicit '+' in the exponent is read as a decimal exponent here; G*Power
  # 3.1.9.7 reads the digits twice (1e+3 = 1003, SUP-G12-UT-CALC-F05), which
  # the result notes as a reference defect.
  r <- calc("1e+3")
  expect_identical(r$value, 1000)
  expect_true(any(grepl("1e\\+3 is 1003", r$notes)))
  expect_identical(calc_value("1e400"), Inf)                  # SUP-G12-UT-CALC-P28
})

test_that("statements, assignments and comments", {
  expect_identical(calc_value("1\n2\n3"), 3)                  # U080-UE01
  expect_identical(calc_value("1;2;3"), 3)                    # U081-UE01
  expect_identical(calc_value("a = 2\nb = a*3\nb+1"), 7)      # U080-UO02
  expect_identical(calc_value("x = 1; x = x + 1; x"), 2)      # U082-UE01
  expect_identical(calc_value("a=2"), 2)                      # SUP-G12-UT-CALC-P15: an assignment is the last value
  expect_identical(calc_value("1+1;"), 2)                     # SUP-G12-UT-CALC-P16: trailing separator
  expect_identical(calc_value("2+2 # trailing comment"), 4)   # U083-UE01
  expect_identical(calc_value("# comment\n2+2"), 4)           # U083-UO01
  expect_identical(calc_value("x1=2; y2=3; x1*y2"), 6)
  r <- calc("a = 1; b = a + 1; b * 10")
  expect_identical(r$statements$kind, c("assign", "assign", "expr"))
  expect_identical(r$variables, list(a = 1, b = 2))
  expect_identical(r$statements$value, c(1, 2, 20))
})

test_that("the G*Power matched-pairs script of the manual evaluates", {
  # DOC-018 / GPW-UT-CALC-UTL-U080-UO01. Value of the independent mpmath
  # reference (validation/independent/g12, 30 digits: 0.815566411124118980...).
  script <- paste("# Power matched pairs",
                  "m1 = 140; m2=0; s1 = 210; s2 = 210; rho=0.5; N = 16; alpha=0.05",
                  "dz = abs(m1-m2)/sqrt(s1^2 + s2^2 - 2*rho*s1*s2)",
                  "d = dz*sqrt(N)", "tc = tinv(1-alpha,N-1)", "power = 1-nctcdf(tc,N-1,d)", sep = "\n")
  r <- calc(script)
  expect_identical(r$status$code, "ok")
  expect_lt(abs(r$value - 0.815566411124118980), 1e-14)
  expect_setequal(r$functions_used, c("abs", "sqrt", "tinv", "nctcdf"))
})

test_that("syntax errors are typed refusals without a number", {
  calc_refuses("1+", "syntax")                                # U086-I01: 'Syntax error: Missing operand at:'
  calc_refuses("(2+3", "syntax")                              # U086-I02: unbalanced parentheses
  calc_refuses("2+3)", "syntax")                              # SUP-G12-UT-CALC-P24
  calc_refuses("foo(1)", "syntax")                            # U086-I03: unknown function
  calc_refuses("max(1:2)", "syntax")                          # U084-UE01: ':' is not an argument separator
  calc_refuses("normcdf(1)", "syntax")                        # SUP-G12-UT-CALC-P18: too few arguments
  calc_refuses("sin(1,2)", "syntax")                          # SUP-G12-UT-CALC-P19: too many arguments
  calc_refuses("y+1", "syntax")                               # SUP-G12-UT-CALC-P13: unknown variable
  calc_refuses("a=1; A", "syntax")                            # SUP-G12-UT-CALC-P12: names are case-sensitive
  calc_refuses("zCDF(1)", "syntax")                           # inventory u084: function names are case-sensitive
  calc_refuses("pi", "syntax")                                # SUP-G12-UT-CALC-P21: no constants
  calc_refuses("e+1", "syntax")                               # SUP-G12-UT-CALC-P22
  calc_refuses("+3", "syntax")                                # SUP-G12-UT-CALC-P07: no unary plus
  calc_refuses("2(3)", "syntax")                              # SUP-G12-UT-CALC-P23: no implicit multiplication
  calc_refuses("sin=2; sin", "syntax")                        # SUP-G12-UT-CALC-P20: function names are reserved
  calc_refuses("my_var=3", "syntax")                          # SUP-G12-UT-CALC-P29: '_' is not part of a name
  r <- calc_refuses("1+", "syntax")
  expect_match(r$status$reason, "missing operand")
  # A script without an expression has no result (G*Power shows 0, SUP-G12-UT-CALC-P17).
  r <- calc("# only a comment")
  expect_identical(r$status$code, "invalid_input")
  expect_true(is.na(r$value))
})

test_that("input text is never evaluated as R code", {
  for (txt in c("system('echo hi')", "q()", "Sys.time()", "eval(1)", "library(stats)", "`+`(1,2)", "x <- 1",
                "1; quit()", "get('x')", "\"text\"", "1 | 2", "!1", "c(1,2)", "1[1]", "x$y", "{1}")) {
    r <- calc(txt)
    expect_identical(r$status$code, "invalid_input", info = txt)
    expect_true(is.na(r$value), info = txt)
  }
  expect_error(distribution_calculator(c("1", "2")), class = "inferguidR_invalid_input")
  expect_error(distribution_calculator(NA_character_), class = "inferguidR_invalid_input")
  expect_error(distribution_calculator(1), class = "inferguidR_invalid_input")
  long <- paste(rep("1+", 20000), collapse = "")
  expect_identical(calc(paste0(long, "1"))$status$code, "invalid_input")
  deep <- paste0(strrep("(", 500), "1", strrep(")", 500))
  expect_identical(calc(deep)$status$code, "invalid_input")
})

test_that("IEEE results are values with notes, as in G*Power", {
  r <- calc("log(-1)")                                        # U086-I04: -nan(ind), no message
  expect_identical(r$status$code, "ok")
  expect_true(is.nan(r$value))
  expect_identical(r$display, "-nan(ind)")
  expect_true(any(grepl("non-finite result: log\\(-1\\)", r$notes)))
  expect_identical(calc("1/0")$display, "inf")                # U086-I05
  expect_identical(calc("1/-0")$display, "-inf")              # SUP-G12-UT-CALC-B49
  expect_identical(calc("log(0)")$display, "-inf")            # SUP-G12-UT-CALC-P26
  expect_identical(calc("0/0")$display, "-nan(ind)")          # SUP-G12-UT-CALC-P25
  expect_identical(calc("10^308*10")$display, "inf")          # U001-UD01
  expect_identical(calc("(-8)^(1/3)")$display, "-nan(ind)")   # SUP-G12-UT-CALC-P33
})

test_that("the display follows the 15 significant digits of the G*Power result field", {
  expect_identical(calc("2/3")$display, "0.666666666666667")                  # U003-UO02
  expect_identical(calc("1e16+1")$display, "1e+16")                           # U004-UD01
  expect_identical(calc("0.1+0.2")$display, "0.3")                            # U004-UO02
  expect_identical(calc("123456789012345678")$display, "1.23456789012346e+17") # SUP-G12-UT-CALC-F02
  expect_identical(calc("0.000012345")$display, "1.2345e-05")                 # SUP-G12-UT-CALC-F03
  expect_identical(calc("-0.5*0")$display, "-0")                              # SUP-G12-UT-CALC-F04
  expect_identical(calc("int(-0.5)")$display, "-0")                           # SUP-G12-UT-CALC-B30
  expect_identical(calc("exp(709)")$display, "8.21840746155497e+307")         # U014-UD01
})

test_that("the result object is complete and reproducible", {
  r <- calc("a = 2; tcdf(a, 10)")
  expect_s3_class(r, "ig_calc_result")
  expect_identical(r$object_type, "calculator_result")
  expect_identical(r$r_call, "inferguidR::distribution_calculator(\"a = 2; tcdf(a, 10)\")")
  expect_identical(r$functions_used, "tcdf")
  expect_true(any(grepl("result:", format.ig_calc_result(r))))
  r2 <- calc("x = 1\ny = \"q\"")
  expect_identical(r2$status$code, "invalid_input")
  expect_identical(calc("a = 1\nb = 2")$r_call, "inferguidR::distribution_calculator(\"a = 1\\nb = 2\")")
})
