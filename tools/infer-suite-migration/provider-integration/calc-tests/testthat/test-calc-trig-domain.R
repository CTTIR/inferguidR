test_that("calculator reduction covers medium arguments and represented pole neighbours", {
 # Independently evaluated at 500 digits from exact binary64 input rationals.
 d<-read.csv(test_path("fixtures","calculator-trig-domain-independent.csv"),colClasses="character")
 d[]<-lapply(d,.ig_parse_decimal)
 for(i in seq_len(nrow(d))) for(kind in c("sin","cos","tan")) {
  expression<-sprintf("%s(%.17g)",kind,d$x[i])
  actual<-distribution_calculator(expression)
  expect_identical(actual$status$code,"ok",info=expression)
  expected<-d[[kind]][i]
  if(expected==0) expect_identical(actual$value,0,info=expression) else
   expect_equal(actual$value/expected,1,tolerance=2e-15,info=expression)
 }
})
