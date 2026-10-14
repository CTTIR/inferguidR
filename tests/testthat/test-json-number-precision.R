test_that("canonical decimal tokens use the same conversion when written and read", {
  x <- c(as.numeric("0.04537547795318898"), .1 + .2, 1 / 3, -1 / 7,
    .Machine$double.xmin, .Machine$double.xmax, 5e-324, -5e-324,
    2^53 + 2, 1e-300, 1e300, 1e-30)
  expect_identical(unlist(.ig_json_parse(ig_to_json(I(x))), use.names = FALSE), x)
  for (v in x) expect_identical(as.double(.ig_json_parse(ig_to_json(v))), v)
  expect_identical(.ig_json_parse('{"a":1,"b":1.5,"c":null,"d":true}'),
    list(a = 1L, b = 1.5, c = NULL, d = TRUE))
  expect_identical(.ig_json_parse("[]"), list())
})

test_that("number recovery never consumes digits inside escaped strings or object keys", {
  txt <- '{"123":[1.25,"0.04537547795318898",null,"quote: \\\" 42",true],"\\u0031":"\\u0032", "n":-6.25e-7}'
  got <- .ig_json_parse(txt)
  expect_identical(got[["123"]][[1]], 1.25)
  expect_identical(got[["123"]][[2]], "0.04537547795318898")
  expect_null(got[["123"]][[3]])
  expect_identical(got[["123"]][[4]], 'quote: " 42')
  expect_true(got[["123"]][[5]])
  expect_identical(got[["1"]], "2")
  expect_identical(got$n, -6.25e-7)
  strings <- c('123 -4.5e6', 'a"12\\34', '\\"7', 'x\n8\t9', '\u03b1 0.25')
  expect_identical(unlist(.ig_json_parse(ig_to_json(I(strings)))), strings)
  long <- paste(rep('a"12\\34 ', 10000), collapse = "")
  expect_identical(.ig_json_parse(ig_to_json(long)), long)
})

test_that("a fixed binary floating-point corpus survives canonical round trips", {
  withr::local_seed(1739)
  bytes <- as.raw(sample.int(256L, 8L * 1024L, replace = TRUE) - 1L)
  values <- readBin(bytes, "double", n = 1024L, size = 8L, endian = "little")
  values <- values[is.finite(values)]
  expect_gt(length(values), 1000L)
  expect_identical(unlist(.ig_json_parse(ig_to_json(I(values))), use.names = FALSE), values)
  nested <- list(z = list(values[1:20], NULL), a = values[21:40])
  expect_identical(as.character(ig_to_json(.ig_json_parse(ig_to_json(nested)))), as.character(ig_to_json(nested)))
})
