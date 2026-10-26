#!/usr/bin/env Rscript

# Install only mandatory dependencies and hide dependency caches in later steps.
stopifnot(getRversion() == "4.1.0")
runner_temp <- Sys.getenv("RUNNER_TEMP")
environment_file <- Sys.getenv("GITHUB_ENV")
if (!nzchar(runner_temp) || !dir.exists(runner_temp) || !nzchar(environment_file))
  stop("RUNNER_TEMP and GITHUB_ENV are required", call. = FALSE)
dependency_library <- tempfile("inferguidR-core-dependencies-", tmpdir = runner_temp)
if (!dir.create(dependency_library)) stop("Could not create a fresh dependency library", call. = FALSE)
dependency_library <- normalizePath(dependency_library, winslash = "/", mustWork = TRUE)
required <- c("digest", "jsonlite")
install.packages(required, lib = dependency_library, repos = "https://cloud.r-project.org",
                 dependencies = NA, Ncpus = 2L)
installed <- rownames(installed.packages(lib.loc = dependency_library))
if (!setequal(installed, required))
  stop("Isolated core dependency installation did not provide exactly digest and jsonlite", call. = FALSE)

# Nonexistent explicit paths prevent R from restoring its usual site/user paths.
hidden <- tempfile("inferguidR-hidden-libraries-", tmpdir = runner_temp)
con <- file(environment_file, open = "a")
writeLines(c(paste0("R_LIBS=", dependency_library),
             paste0("R_LIBS_USER=", file.path(hidden, "user")),
             paste0("R_LIBS_SITE=", file.path(hidden, "site"))), con)
close(con)
cat("Isolated mandatory dependencies installed into:", dependency_library, "\n")
