library(testthat)
library(inferguidR)

local({
  # Respect a one-worker request and cap the test pool at two workers.
  requested <- suppressWarnings(as.integer(getOption("Ncpus", Sys.getenv("TESTTHAT_CPUS", "2"))))
  workers <- if (length(requested) == 1L && !is.na(requested) && requested >= 1L) min(2L, requested) else 1L
  old <- options(Ncpus = workers)
  on.exit(options(old), add = TRUE)
  test_check("inferguidR")
})
