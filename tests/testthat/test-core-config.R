test_that("configuration options have documented defaults", {
  withr::local_options(inferguidR.reference_archive_dir = NULL,
                       inferguidR.validate_json_schema = NULL)
  cfg <- ig_config()
  expect_identical(names(cfg), c("inferguidR.reference_archive_dir",
                                 "inferguidR.validate_json_schema"))
  expect_null(cfg$inferguidR.reference_archive_dir)
  expect_true(cfg$inferguidR.validate_json_schema)
  expect_null(.ig_reference_archive_dir())
  expect_error(.ig_get_config("unknown_option"), "unknown inferguidR option")
})

test_that("a configured reference archive must exist", {
  dir <- withr::local_tempdir()
  withr::local_options(inferguidR.reference_archive_dir = dir)
  expect_identical(.ig_reference_archive_dir(), normalizePath(dir, winslash = "/"))
  withr::local_options(inferguidR.reference_archive_dir = file.path(dir, "missing"))
  expect_identical(th_condition_code(.ig_reference_archive_dir()), "invalid_input")
  withr::local_options(inferguidR.reference_archive_dir = 42)
  expect_identical(th_condition_code(.ig_reference_archive_dir()), "invalid_input")
})

test_that("config initialisation never overrides user settings", {
  withr::local_options(inferguidR.validate_json_schema = FALSE)
  .ig_config_init()
  expect_false(getOption("inferguidR.validate_json_schema"))
})

test_that("launch_app locates resources in the installed package and fails typed without an app", {
  app <- .ig_app_dir()
  if (is.null(app)) {
    cond <- tryCatch(launch_app(launch.browser = FALSE), inferguidR_unsupported = function(e) e)
    expect_s3_class(cond, "inferguidR_unsupported")
    expect_identical(cond$status$code, "unsupported")
    expect_match(conditionMessage(cond), "system.file", fixed = TRUE)
  } else {
    expect_true(file.exists(file.path(app, "app.R")) || file.exists(file.path(app, "ui.R")))
  }
})

test_that("bundled resources are found through system.file", {
  expect_true(nzchar(system.file("registry", package = "inferguidR")))
  expect_true(nzchar(system.file("schemas", "power_spec.schema.json", package = "inferguidR")))
})
