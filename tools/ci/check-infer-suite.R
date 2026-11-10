# Build and check exact local providers before consumers. Run with Rscript.
args <- commandArgs(TRUE)
if (length(args) != 4L) stop("Usage: Rscript check-infer-suite.R CALC_SOURCE GUID_SOURCE REPORT_SOURCE OUTPUT")
sources <- normalizePath(args[1:3], mustWork = TRUE)
names(sources) <- c("infercalcR", "inferguidR", "infereportR")
output <- normalizePath(args[4], mustWork = FALSE)
dir.create(output, recursive = TRUE, showWarnings = FALSE)
output <- normalizePath(output)
lib <- file.path(output, "library")
minimal <- file.path(output, "minimal-library")
empty <- file.path(output, "empty-site")
for (p in c(lib, minimal, empty)) dir.create(p, showWarnings = FALSE)
original_libs <- .libPaths()
R <- file.path(R.home("bin"), "R")
Rscript <- file.path(R.home("bin"), "Rscript")
results <- data.frame(check_id = character(), package = character(), status = character(),
  exit_code = integer(), log = character(), stringsAsFactors = FALSE)
run <- function(id, package, command, arguments, env = character()) {
  log <- file.path(output, paste0(id, "-", package, ".log"))
  status <- system2(command, vapply(arguments, shQuote, ""), stdout = log, stderr = log, env = env)
  results[nrow(results) + 1L, ] <<- list(id, package, if (status == 0L) "PASS" else "FAIL", status, log)
  utils::write.csv(results, file.path(output, "results.csv"), row.names = FALSE)
  if (status != 0L) stop("Failed ", id, " for ", package, "; inspect ", log)
}
writeLines(c(capture.output(sessionInfo()), "", ".libPaths:", original_libs), file.path(output, "environment.txt"))
old <- getwd(); setwd(output)
archives <- character()
for (pkg in names(sources)) {
  version <- read.dcf(file.path(sources[[pkg]], "DESCRIPTION"))[1, "Version"]
  archive <- file.path(output, paste0(pkg, "_", version, ".tar.gz"))
  # Build-vignette installation replaces its user library; R_LIBS must retain
  # the exact provider and the declared integration dependencies in child R.
  libraries <- shQuote(paste(c(lib, original_libs), collapse = .Platform$path.sep))
  full_env <- c(paste0("R_LIBS=", libraries), paste0("R_LIBS_USER=", libraries),
                "_R_CHECK_FORCE_SUGGESTS_=false", "INFER_SUITE_FULL=true")
  run("build", pkg, R, c("CMD", "build", sources[[pkg]]), full_env)
  stopifnot(file.exists(archive))
  archives[pkg] <- archive
  run("install", pkg, R, c("CMD", "INSTALL", "--install-tests", paste0("--library=", lib), archive), full_env)
  # Check built resources, examples and installed tests; record every warning.
  run("check", pkg, R, c("CMD", "check", "--no-manual", "--no-vignettes", paste0("--library=", lib), archive), full_env)
  checklog <- file.path(output, paste0(pkg, ".Rcheck"), "00check.log")
  lines <- readLines(checklog, warn = FALSE)
  writeLines(lines[grepl("Status:|ERROR|WARNING|NOTE", lines)], file.path(output, paste0(pkg, "-check-disposition.txt")))
  if (any(grepl(" ERROR| WARNING", lines))) stop("Package check requires repair or a recorded disposition: ", pkg)
}
minimal_env <- c(paste0("R_LIBS=", shQuote(minimal)), paste0("R_LIBS_USER=", shQuote(minimal)), paste0("R_LIBS_SITE=", shQuote(empty)),
                 "_R_CHECK_FORCE_SUGGESTS_=false", "INFER_SUITE_FULL=false")
for (pkg in c("infercalcR", "infereportR")) {
  run("minimal-install", pkg, R, c("CMD", "INSTALL", "--install-tests", paste0("--library=", minimal), archives[[pkg]]), minimal_env)
  dir.create(file.path(output, "minimal-check"), showWarnings = FALSE)
  run("minimal-check", pkg, R, c("CMD", "check", "--no-manual", "--no-vignettes",
    paste0("--output=", file.path(output, "minimal-check")),
    paste0("--library=", minimal), archives[[pkg]]), minimal_env)
  checklog <- file.path(output, "minimal-check", paste0(pkg, ".Rcheck"), "00check.log")
  lines <- readLines(checklog, warn = FALSE)
  writeLines(lines[grepl("Status:|ERROR|WARNING|NOTE", lines)],
    file.path(output, paste0(pkg, "-minimal-disposition.txt")))
  if (any(grepl(" ERROR| WARNING", lines))) stop("Minimal package check needs disposition: ", pkg)
}
probe <- file.path(output, "minimal-probe.R")
writeLines(c("stopifnot(!any(vapply(c('digest','jsonlite','ggplot2','testthat','DescrTab2'), requireNamespace, logical(1), quietly=TRUE)))",
  "library(infercalcR); library(infereportR)",
  "x <- describe_data(data.frame(value=c(1,2,NA)))",
  "stopifnot(is.data.frame(report_table(x)), is.character(report_text(x)))",
  "print(.libPaths()); print(sessionInfo())"), probe)
run("minimal-probe", "suite", Rscript, c("--vanilla", probe), minimal_env)
setwd(old)
cat("Installed suite lanes completed; review every recorded check disposition.\n")
