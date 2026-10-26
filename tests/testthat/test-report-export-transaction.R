.report_export_clone <- function(bindings) {
  fun <- export_report
  environment(fun) <- list2env(bindings, parent = environment(fun))
  fun
}

.report_export_fixture <- function() {
  spec <- analysis_spec("t_means_one_sample", "post_hoc",
    list(d = .5, alpha = .05, n_total = 30), list(tails = "two"))
  create_report_snapshot(analysis_envelope(spec, calculate_power(spec)))
}

test_that("report replacement commits a complete sibling file", {
  snapshot <- .report_export_fixture()
  directory <- tempfile(); dir.create(directory)
  on.exit(unlink(directory, recursive = TRUE), add = TRUE)
  formats <- c("json", "txt", if (report_pdf_capabilities()$available) "pdf")
  for (format in formats) {
    file <- file.path(directory, paste0("report.", format))
    first <- export_report(snapshot, file, format)
    expect_identical(first$sha256, digest::digest(file = file, algo = "sha256"))
    writeLines("old destination", file)
    rename <- function(from, to) {
      expect_identical(normalizePath(dirname(from)), normalizePath(dirname(to)))
      expect_identical(readLines(to), "old destination")
      expect_gt(file.size(from), 0)
      base::file.rename(from, to)
    }
    receipt <- .report_export_clone(list(file.rename = rename))(snapshot, file, format, overwrite = TRUE)
    if (format != "pdf") expect_identical(receipt$sha256, first$sha256)
    expect_identical(receipt$sha256, digest::digest(file = file, algo = "sha256"))
    expect_identical(list.files(directory, all.files = TRUE, no.. = TRUE), basename(file))
    unlink(file)
  }
})

test_that("failed staging rendering and commit preserve an existing report", {
  snapshot <- .report_export_fixture()
  directory <- tempfile(); dir.create(directory)
  on.exit(unlink(directory, recursive = TRUE), add = TRUE)
  file <- file.path(directory, "report.txt")
  writeLines("old destination", file)
  original <- readBin(file, "raw", n = file.size(file))
  faults <- list(
    stage = list(file.create = function(...) FALSE),
    permissions = list(Sys.chmod = function(...) FALSE),
    partial_render = list(writeLines = function(text, con, ...) {
      base::writeLines("partial staged output", con)
      stop("render cancelled")
    }),
    empty_render = list(writeLines = function(...) invisible(NULL)),
    commit = list(file.rename = function(...) FALSE),
    commit_error = list(file.rename = function(...) stop("rename cancelled"))
  )
  for (fault in faults) {
    expect_error(.report_export_clone(fault)(snapshot, file, "txt", overwrite = TRUE))
    expect_identical(readBin(file, "raw", n = file.size(file)), original)
    expect_identical(list.files(directory, all.files = TRUE, no.. = TRUE), "report.txt")
  }
  cancel <- function(...) stop(structure(list(message = "user interrupt"), class = c("interrupt", "condition")))
  interrupted <- tryCatch({
    .report_export_clone(list(file.rename = cancel))(snapshot, file, "txt", overwrite = TRUE)
    FALSE
  }, interrupt = function(e) TRUE)
  expect_true(interrupted)
  expect_identical(readBin(file, "raw", n = file.size(file)), original)
  expect_identical(list.files(directory, all.files = TRUE, no.. = TRUE), "report.txt")
})

test_that("a destination created during rendering is retained without overwrite", {
  snapshot <- .report_export_fixture()
  directory <- tempfile(); dir.create(directory)
  on.exit(unlink(directory, recursive = TRUE), add = TRUE)
  file <- file.path(directory, "new.txt")
  render <- function(text, con, ...) {
    base::writeLines(text, con, ...)
    base::writeLines("another writer", file)
  }
  fun <- .report_export_clone(list(writeLines = render))
  expect_error(fun(snapshot, file, "txt"), class = "inferguidR_invalid_input")
  expect_identical(readLines(file), "another writer")
  expect_identical(list.files(directory, all.files = TRUE, no.. = TRUE), "new.txt")
})

test_that("dangling destination links require explicit replacement", {
  snapshot <- .report_export_fixture()
  directory <- tempfile(); dir.create(directory)
  on.exit(unlink(directory, recursive = TRUE), add = TRUE)
  file <- file.path(directory, "report.txt")
  target <- file.path(directory, "missing.txt")
  skip_if_not(suppressWarnings(file.symlink(target, file)), "symbolic links unavailable")
  expect_error(export_report(snapshot, file, "txt"), class = "inferguidR_invalid_input")
  expect_identical(Sys.readlink(file), target)
  expect_false(file.exists(target))
  receipt <- export_report(snapshot, file, "txt", overwrite = TRUE)
  expect_identical(Sys.readlink(file), "")
  expect_false(file.exists(target))
  expect_identical(receipt$sha256, digest::digest(file = file, algo = "sha256"))
  unlink(file)
  render <- function(text, con, ...) {
    base::writeLines(text, con, ...)
    stopifnot(file.symlink(target, file))
  }
  expect_error(.report_export_clone(list(writeLines = render))(snapshot, file, "txt"),
    class = "inferguidR_invalid_input")
  expect_identical(Sys.readlink(file), target)
  expect_false(file.exists(target))
  expect_identical(list.files(directory, all.files = TRUE, no.. = TRUE), "report.txt")
  unlink(file)
  writeLines("retained linked target", target)
  stopifnot(file.symlink(target, file))
  export_report(snapshot, file, "txt", overwrite = TRUE)
  expect_identical(Sys.readlink(file), "")
  expect_identical(readLines(target), "retained linked target")
})
