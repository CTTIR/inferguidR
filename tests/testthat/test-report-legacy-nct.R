legacy_nct_json <- function(x) {
  x <- unclass(x); x$snapshot_hash <- NULL
  x$snapshot_hash <- .ig_hash_value(x)
  .ig_canonical_json(x)
}

test_that("schema 0.1 reports retain the known previous successful noncentral-t renderer", {
  raw <- paste(readLines(test_path("fixtures", "report-legacy-nct.json")), collapse = "\n")
  saved <- .ig_json_parse(raw)
  doc <- ig_from_json(raw)
  expect_s3_class(doc, "report_document")
  expect_identical(doc$visualisations, saved$visualisations)
  expect_identical(as.character(ig_to_json(doc)), raw)
  expect_identical(as.character(ig_to_json(ig_from_json(ig_to_json(doc)))), raw)
  result <- .ig_result_from_list(doc$results[[1]])
  modern <- power_plot_data(result)
  expect_identical(modern$status, "ok")
  expect_false(identical(modern$data$x, vapply(doc$visualisations[[1]]$data, `[[`, 0, "x")))
  for (field in c("x", "y", "reject", "hypothesis")) {
    bad <- unclass(doc)
    row <- bad$visualisations[[1]]$data[[10]]
    row[[field]] <- switch(field, x = row$x + .01, y = row$y + .01,
      reject = !row$reject, hypothesis = "altered")
    bad$visualisations[[1]]$data[[10]] <- row
    expect_error(ig_from_json(legacy_nct_json(bad)), class = "inferguidR_invalid_input")
  }
  bad <- unclass(doc); bad$visualisations[[1]]$data <- head(bad$visualisations[[1]]$data, -1L)
  expect_error(ig_from_json(legacy_nct_json(bad)), class = "inferguidR_invalid_input")
  bad <- unclass(doc); bad$schema_version <- "inferguidR.report_document/0.2.0"
  expect_error(ig_from_json(legacy_nct_json(bad)), class = "inferguidR_invalid_input")
  new <- create_report_snapshot(analysis_envelope(result$spec, result))
  expect_identical(ig_to_json(ig_from_json(ig_to_json(new))), ig_to_json(new))
})
