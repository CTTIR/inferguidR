#!/usr/bin/env Rscript

# Install one verified CI source bundle without using a cached package library.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L || any(!nzchar(args))) {
  stop("Usage: Rscript tools/ci/install-candidate.R CANDIDATE_DIRECTORY INSTALL_LOG", call. = FALSE)
}
bundle <- normalizePath(args[[1L]], winslash = "/", mustWork = TRUE)
log <- args[[2L]]
runner_temp <- Sys.getenv("RUNNER_TEMP")
environment_file <- Sys.getenv("GITHUB_ENV")
source_commit <- Sys.getenv("GITHUB_SHA")
if (!nzchar(runner_temp) || !dir.exists(runner_temp) || !nzchar(environment_file) || !nzchar(source_commit)) {
  stop("RUNNER_TEMP, GITHUB_ENV and GITHUB_SHA are required", call. = FALSE)
}
archives <- list.files(bundle, pattern = "^inferguidR_.*[.]tar[.]gz$", full.names = TRUE)
if (length(archives) != 1L) stop("The source bundle must contain exactly one candidate archive", call. = FALSE)
manifest <- readLines(file.path(bundle, "SOURCE-SHA256.txt"), warn = FALSE)
if (length(manifest) != 1L) stop("Expected one source checksum record", call. = FALSE)
fields <- strsplit(trimws(manifest), "[[:space:]]+")[[1L]]
if (length(fields) != 2L || !grepl("^[a-f0-9]{64}$", fields[[1L]]) ||
    !identical(fields[[2L]], basename(archives))) {
  stop("Source checksum manifest does not identify the candidate archive", call. = FALSE)
}
expected_sha <- fields[[1L]]
if (!identical(digest::digest(file = archives, algo = "sha256"), expected_sha)) {
  stop("Candidate archive SHA-256 mismatch", call. = FALSE)
}
if (!identical(readLines(file.path(bundle, "SOURCE-COMMIT.txt"), warn = FALSE), source_commit)) {
  stop("Candidate archive belongs to a different source commit", call. = FALSE)
}

candidate_library <- tempfile("inferguidR-candidate-", tmpdir = runner_temp)
if (!dir.create(candidate_library)) stop("Could not create a fresh candidate library", call. = FALSE)
candidate_library <- normalizePath(candidate_library, winslash = "/", mustWork = TRUE)
dir.create(dirname(log), recursive = TRUE, showWarnings = FALSE)
r <- file.path(R.home("bin"), if (.Platform$OS.type == "windows") "R.exe" else "R")
status <- system2(r, c("CMD", "INSTALL", shQuote(paste0("--library=", candidate_library)),
                       shQuote(archives)), stdout = log, stderr = log)
installed_description <- file.path(candidate_library, "inferguidR", "DESCRIPTION")
installed <- status == 0L && file.exists(installed_description)
if (installed) {
  installed <- identical(unname(read.dcf(installed_description, fields = "Package")[1L, 1L]), "inferguidR")
}
jsonlite::write_json(list(
  source_sha256 = expected_sha,
  source_commit = source_commit,
  archive = archives,
  candidate_library = candidate_library,
  r_version = R.version.string,
  install_exit_code = status,
  success = installed
), paste0(log, ".json"), pretty = TRUE, auto_unbox = TRUE)
if (!installed) {
  cat(paste(tail(readLines(log, warn = FALSE), 30L), collapse = "\n"), "\n")
  stop(sprintf("Candidate installation failed (exit %s); no candidate environment exported", status), call. = FALSE)
}

old_libraries <- Sys.getenv("R_LIBS")
libraries <- paste(c(candidate_library, if (nzchar(old_libraries)) old_libraries), collapse = .Platform$path.sep)
con <- file(environment_file, open = "a")
writeLines(c(paste0("R_LIBS=", libraries),
             paste0("INFERGUIDR_EXPECTED_LIBRARY=", candidate_library),
             paste0("INFERGUIDR_SOURCE_SHA256=", expected_sha)), con)
close(con)
cat("Candidate installed into fresh library:", candidate_library, "\n")
