.onLoad <- function(libname, pkgname) {
  .ig_config_init()
  .ig_registry_init(ns = asNamespace(pkgname))
  problems <- c(if (!is.null(.ig_state$load_error)) .ig_state$load_error, .ig_state$problems)
  if (length(problems)) {
    warning("inferguidR registry problems:\n", paste0("* ", problems, collapse = "\n"),
            call. = FALSE)
  }
  invisible(NULL)
}
