# Exercise the foundation contracts independently of later method providers.
# Restore the live registry and solvers when this file finishes.
th_registry_only()

test_that("canonical JSON sorts keys and is compact", {
  expect_identical(as.character(ig_to_json(list(b = 1, a = list(d = TRUE, c = "x")))),
                   '{"a":{"c":"x","d":true},"b":1}')
})

test_that("numbers use the shortest representation that round-trips exactly", {
  x <- c(0.1 + 0.2, 0.3, 1 / 3, 1e-300, 2^53 + 2, 5e-324, 1e21, 123456789.123456789)
  txt <- as.character(ig_to_json(I(x)))
  back <- unlist(jsonlite::parse_json(txt))
  expect_identical(back, x)
  expect_identical(as.character(ig_to_json(0.1 + 0.2)), "0.30000000000000004")
  expect_identical(as.character(ig_to_json(0.3)), "0.3")
  expect_identical(as.character(ig_to_json(-0)), "0")
  expect_identical(as.character(ig_to_json(128L)), as.character(ig_to_json(128)))
})

test_that("missing and non-finite values are encoded explicitly", {
  txt <- as.character(ig_to_json(list(a = NA, b = c(1, Inf, -Inf, NaN, NA))))
  expect_identical(txt, paste0(
    '{"a":null,"b":[1,{"nonfinite":"positive_infinity","value":null},',
    '{"nonfinite":"negative_infinity","value":null},',
    '{"nonfinite":"not_a_number","value":null},null]}'))
  back <- .ig_json_parse(txt)
  expect_identical(back$b[[2]], Inf)
  expect_identical(back$b[[3]], -Inf)
  expect_true(is.nan(back$b[[4]]))
  expect_null(back$b[[5]])
})

test_that("boxing rules: I() arrays, named vectors as objects, empty containers", {
  expect_identical(as.character(ig_to_json(list(a = I("x"), b = "x"))), '{"a":["x"],"b":"x"}')
  expect_identical(as.character(ig_to_json(c(n2 = 3L, n1 = 4L))), '{"n1":4,"n2":3}')
  expect_identical(as.character(ig_to_json(list(e = .ig_empty_object(), l = list(), n = NULL))),
                   '{"e":{},"l":[],"n":null}')
  expect_identical(as.character(ig_to_json(character())), "[]")
})

test_that("strings are escaped and kept as UTF-8", {
  txt <- as.character(ig_to_json(list(s = "quote \" backslash \\ newline \n \u03b1")))
  expect_identical(jsonlite::parse_json(txt)$s, "quote \" backslash \\ newline \n \u03b1")
})

test_that("values that are not data cannot be serialised", {
  expect_error(ig_to_json(list(f = function(x) x)), "cannot serialise")
  expect_error(ig_to_json(list(e = new.env())), "cannot serialise")
  expect_error(ig_to_json(stats::setNames(list(1, 2), c("a", "a"))), "duplicate key")
  expect_error(ig_to_json(stats::setNames(list(1, 2), c("a", ""))), "partially named")
})

test_that("hashes are SHA-256 and independent of key order", {
  h1 <- ig_hash(list(a = 1, b = list(c = "x", d = 2L)))
  h2 <- ig_hash(list(b = list(d = 2, c = "x"), a = 1))
  expect_identical(h1, h2)
  expect_match(h1, "^[0-9a-f]{64}$")
  expect_false(identical(h1, ig_hash(list(a = 1 + 1e-15, b = list(c = "x", d = 2L)))))
})

test_that("parsing rejects duplicate keys and invalid JSON without evaluating", {
  expect_identical(th_condition_code(.ig_json_parse('{"a":1,"a":2}')), "invalid_input")
  expect_identical(th_condition_code(.ig_json_parse('{"a":{"b":1,"b":1}}')), "invalid_input")
  expect_identical(th_condition_code(.ig_json_parse("not json")), "invalid_input")
  # A URL is data, never fetched: it is simply not valid JSON.
  expect_identical(th_condition_code(.ig_json_parse("https://example.org/spec.json")),
                   "invalid_input")
  expect_identical(th_condition_code(.ig_json_parse('"system(\'echo hi\')"')), NA_character_)
  expect_identical(.ig_json_parse('"system(\'echo hi\')"'), "system('echo hi')")
})

test_that("ig_from_json refuses unknown types, versions and mismatched types", {
  spec <- th_minimal_spec("t_means_one_sample", "post_hoc")
  lst <- .ig_spec_to_list(spec)
  expect_identical(th_condition_code(ig_from_json('{"object_type":"power_magic"}')), "invalid_input")
  expect_identical(th_condition_code(ig_from_json("[1,2]")), "invalid_input")
  future <- lst
  future$schema_version <- "inferguidR.power_spec/9.0.0"
  cond <- tryCatch(ig_from_json(ig_to_json(future)), inferguidR_invalid_input = function(e) e)
  expect_match(cond$status$reason, "unsupported schema_version")
  expect_identical(th_condition_code(ig_from_json(ig_to_json(spec), type = "power_result")),
                   "invalid_input")
  expect_identical(th_condition_code(ig_from_json(json = "{}", file = "x.json")), "invalid_input")
  expect_identical(th_condition_code(ig_from_json(file = "https://example.org/x.json")),
                   "invalid_input")
})

test_that("ig_from_json reads local files", {
  spec <- th_minimal_spec("t_means_one_sample", "post_hoc")
  path <- withr::local_tempfile(fileext = ".json")
  writeLines(ig_to_json(spec, pretty = TRUE), path, useBytes = TRUE)
  expect_identical(ig_from_json(file = path), spec)
  expect_identical(th_condition_code(ig_from_json(file = paste0(path, ".missing"))), "invalid_input")
})

test_that("pretty output contains the same tokens as the canonical form", {
  spec <- th_minimal_spec("t_means_one_sample", "post_hoc")
  pretty <- ig_to_json(spec, pretty = TRUE)
  expect_identical(as.character(jsonlite::minify(pretty)), as.character(ig_to_json(spec)))
})
