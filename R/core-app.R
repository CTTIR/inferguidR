# Application launcher ------------------------------------------------------------------

#' Launch the inferguidR application
#'
#' Starts the bundled 'shiny' application. The application directory is
#' located inside the installed package (`system.file("app", package =
#' "inferguidR")`), never through a development path.
#'
#' While the application is not part of the installed build, the function
#' stops with an error of class `inferguidR_unsupported` (its `$status` is an
#' [ig_status()] with code `"unsupported"`).
#'
#' @param ... Further arguments passed to `shiny::runApp()`.
#' @param launch.browser Passed to `shiny::runApp()`.
#' @return Called for its side effect; returns the value of
#'   `shiny::runApp()` when the application exists.
#' @examples
#' app_bundled <- nzchar(system.file("app", "app.R", package = "inferguidR"))
#' if (!app_bundled) {
#'   # Without a bundled application the call fails with a typed condition.
#'   tryCatch(launch_app(), inferguidR_unsupported = function(e) e$status)
#' }
#' if (interactive() && app_bundled) launch_app()
#' @export
launch_app <- function(..., launch.browser = interactive()) {
  app_dir <- .ig_app_dir()
  if (is.null(app_dir)) {
    .ig_abort("unsupported", paste(
      "the inferguidR application is not included in this build:",
      "no 'app' directory with app.R (or ui.R and server.R) was found in the installed package",
      "(system.file(\"app\", package = \"inferguidR\"))"))
  }
  if (!requireNamespace("shiny", quietly = TRUE)) {
    .ig_abort("unsupported", "launch_app() needs the 'shiny' package; install it with install.packages(\"shiny\")")
  }
  shiny::runApp(app_dir, launch.browser = launch.browser, ...)
}

# Directory of the bundled application, or NULL when absent.
.ig_app_dir <- function() {
  dir <- system.file("app", package = "inferguidR")
  if (!nzchar(dir)) return(NULL)
  has_app <- file.exists(file.path(dir, "app.R")) ||
    (file.exists(file.path(dir, "ui.R")) && file.exists(file.path(dir, "server.R")))
  if (has_app) dir else NULL
}
