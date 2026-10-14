rr_rehash <- function(x) {
  x <- unclass(x); x$snapshot_hash <- NULL
  x$snapshot_hash <- .ig_hash_value(x)
  .ig_canonical_json(x)
}

test_that("the current teaching catalog identifies the shipped lesson bytes", {
  record <- .ig_report_resource_catalog()$entries[[1]]
  for (entry in record$files) {
    path <- system.file("lessons", "boschloo", entry$name, package = "inferguidR")
    expect_identical(digest::digest(file = path, algo = "sha256"), entry$sha256)
  }
})

test_that("verified lessons freeze exact topic metadata in bilingual reports", {
  result <- boschloo_test(1, 2, 0, 2, variant = "minlike")
  e <- analysis_envelope(result$spec, result)
  request <- list(list(analysis = 1, resource_id = "inferguidR_boschloo_minlike"))
  for (lang in c("en", "de")) {
    x <- create_report_snapshot(e, locale = lang, resources = request)
    expect_identical(x$schema_version, "inferguidR.report_document/0.3.0")
    expect_identical(x$resources$status, "attached")
    expect_identical(x$resources$entries[[1]]$locale, lang)
    expect_identical(ig_to_json(ig_from_json(ig_to_json(x))), ig_to_json(x))
    expect_match(paste(.ig_report_summary_lines(x), collapse = "\n"), "10.1111/j.1467-9574.1970.tb00104.x", fixed = TRUE)
    for (field in c("title", "scope", "source", "access", "match", "analysis_variant_id")) {
      bad <- unclass(x); bad$resources$entries[[1]][[field]] <- "fabricated"
      expect_error(ig_from_json(rr_rehash(bad)), class = "inferguidR_invalid_input")
    }
    bad <- unclass(x); bad$resources$entries[[1]]$references[[1]]$url <- "https://example.invalid/"
    expect_error(ig_from_json(rr_rehash(bad)), class = "inferguidR_invalid_input")
    bad <- unclass(x); bad$resources$entries[[1]]$files[[1]]$sha256 <- strrep("0", 64)
    expect_error(ig_from_json(rr_rehash(bad)), class = "inferguidR_invalid_input")
    bad <- unclass(x); bad$resources$catalog_version <- "2099.1"
    expect_error(ig_from_json(rr_rehash(bad)), class = "inferguidR_invalid_input")
    for (version in c("0.1.0", "0.2.0")) {
      bad <- unclass(x); bad$schema_version <- paste0("inferguidR.report_document/", version)
      expect_error(ig_from_json(rr_rehash(bad)), class = "inferguidR_invalid_input")
    }
    path <- tempfile(fileext = ".txt"); on.exit(unlink(path), add = TRUE)
    export_report(x, path, "txt")
    expect_match(paste(readLines(path), collapse = "\n"), x$resources$entries[[1]]$title, fixed = TRUE)
    if (requireNamespace("jsonvalidate", quietly = TRUE)) {
      schema <- system.file("schemas", "report_document.schema.json", package = "inferguidR")
      expect_true(jsonvalidate::json_validate(ig_to_json(x), schema, engine = "ajv"))
      expect_false(jsonvalidate::json_validate(rr_rehash(bad), schema, engine = "ajv"))
    }
  }
  old <- create_report_snapshot(e)
  expect_identical(old$schema_version, "inferguidR.report_document/0.1.0")
  expect_identical(ig_to_json(ig_from_json(ig_to_json(old))), ig_to_json(old))
  expect_error(create_report_snapshot(e, resources = c(request, request)), "Duplicate")
  for (index in list(0, 2, 1.1, NA_real_, Inf, "1", NULL)) {
    bad <- request; bad[[1]]$analysis <- index
    expect_error(create_report_snapshot(e, resources = bad), class = "inferguidR_invalid_input")
  }
  bad <- request; bad[[1]]$resource_id <- "not_verified"
  expect_error(create_report_snapshot(e, resources = bad), "Unknown verified")
  central <- boschloo_test(1, 2, 0, 2, variant = "central")
  expect_error(create_report_snapshot(analysis_envelope(central$spec, central), resources = request), "selected analysis variant")
})

test_that("power resources and sweeps retain independent associations", {
  s <- analysis_spec("boschloo_proportions_two_independent", "post_hoc",
    list(n1 = 2, n2 = 2, p1 = .8, p2 = .2, alpha = .05), list(variant = "minlike"))
  e <- analysis_envelope(s, calculate_power(s))
  request <- list(list(analysis = 2, resource_id = "inferguidR_boschloo_minlike"))
  x <- create_report_snapshot(list(e, e), report_kind = "comparison", resources = request)
  expect_identical(x$resources$entries[[1]]$analysis, 2L)
  expect_identical(ig_to_json(ig_from_json(ig_to_json(x))), ig_to_json(x))
  sw <- sweep_design(s, "p1", c(.7, .8))
  request[[1]]$analysis <- 1
  x <- create_report_snapshot(e, resources = request, sweeps = list(list(analysis = 1, sweep = sw)))
  expect_length(x$sweeps, 1)
  expect_identical(ig_to_json(ig_from_json(ig_to_json(x))), ig_to_json(x))
})

test_that("app freezes eligible references and keeps frozen selection after changes", {
  skip_if_not_installed("shiny")
  r <- boschloo_test(1, 2, 0, 2, variant = "minlike")
  e <- analysis_envelope(r$spec, r)
  expect_length(.ig_report_envelope_resources(e), 1)
  r2 <- boschloo_test(1, 2, 0, 2, variant = "central")
  expect_length(.ig_report_envelope_resources(analysis_envelope(r2$spec, r2)), 0)
  shiny::testServer(.ig_report_server, args = list(envelope = shiny::reactive(e), language = shiny::reactive("en")), {
    session$flushReact()
    session$setInputs(include_resources = TRUE, freeze = 1)
    expect_identical(snapshot()$resources$status, "attached")
    frozen <- snapshot()$snapshot_hash
    session$setInputs(include_resources = FALSE)
    expect_identical(snapshot()$snapshot_hash, frozen)
    session$setInputs(freeze = 2)
    expect_identical(snapshot()$resources$status, "not_attached")
    expect_identical(snapshot()$schema_version, "inferguidR.report_document/0.1.0")
  })
})
