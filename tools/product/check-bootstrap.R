#!/usr/bin/env Rscript
# Copyright (c) 2026 Raban Heller. MIT licence; see ../../LICENSE.md.
args <- commandArgs(TRUE)
if (!length(args) || length(args) > 2L || (length(args) == 2L && args[2] != "--negative-control"))
  stop("usage: Rscript tools/product/check-bootstrap.R OUTPUT_DIRECTORY [--negative-control]")
out <- path.expand(args[1])
if (dir.exists(out) && length(list.files(out, all.files = TRUE, no.. = TRUE))) stop("output directory must be new or empty")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
if (!dir.exists(out)) stop("cannot create output directory")
script <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)))
checks <- list(); failure <- NULL; coverage <- NULL
record <- function(case, quantity, observed, reference, tolerance = 1e-12) {
  if (length(args) == 2L && case == "class" && quantity == "tpr.lower") observed <- observed + .01
  same <- if (is.numeric(reference)) isTRUE(is.finite(observed) && abs(observed-reference) <= tolerance) else identical(observed, reference)
  checks[[length(checks)+1L]] <<- data.frame(case, quantity, observed = as.character(observed),
    reference = as.character(reference), tolerance, passed = same)
}
tryCatch({
  library(inferguidR)
  installed <- find.package("inferguidR")
  if (!file.exists(file.path(installed, "Meta", "package.rds"))) stop("installed package required")
  expected <- Sys.getenv("INFERGUIDR_EXPECTED_LIBRARY")
  if (nzchar(expected) && normalizePath(file.path(expected,"inferguidR")) != normalizePath(installed)) stop("unexpected package library")
  source(file.path(dirname(script), "reference-bootstrap.R"), local = TRUE)
  d <- data.frame(y = rep(c(1,0), 6), p = c(.9,.8,.7,.6,.5,.5,.4,.3,.2,.1,.75,.25),
    w = c(1,2,3,1,2,1,1,2,1,3,2,1), cluster = rep(letters[1:6], each = 2))
  for (design in c("class", "none", "cluster")) {
    refs <- reference_bootstrap(d, design, 250, 742)
    unc <- classification_uncertainty(proportion_methods = character(), auc_method = "none",
      bootstrap_replicates = 250, seed = 742, bootstrap_strata = if (design == "cluster") "class" else design,
      bootstrap_metrics = vapply(refs, `[[`, "", "metric_id"))
    sp <- classification_spec(d, truth = "y", positive_label = 1, scores = "p", score_scale = "probability",
      cutoff = .5, weights = "w", cluster = if (design == "cluster") "cluster" else NULL, uncertainty = unc)
    before <- .Random.seed
    result <- evaluate_cutoff(sp)
    record(design, "rng_restored", identical(before, .Random.seed), TRUE)
    observed <- classification_table(result, "intervals")
    for (ref in refs) {
      row <- observed[observed$metric_id == ref$metric_id & observed$method == "bootstrap_percentile", ]
      if (nrow(row) != 1L) stop("missing or duplicate bootstrap interval")
      record(design, paste0(ref$metric_id,".status"), row$status, ref$status)
      if (ref$status == "ok") for (bound in c("lower","upper"))
        record(design, paste0(ref$metric_id,".",bound), row[[bound]], ref[[bound]])
    }
  }
  coverage <- expand.grid(n_event = c(5,10,20), true_sensitivity = c(.01,.1,.5,.9,.99))
  coverage$nominal <- .95
  coverage$coverage_infinite_replicates <- mapply(reference_binomial_coverage, coverage$n_event, coverage$true_sensitivity)
  utils::write.csv(coverage, file.path(out,"coverage-diagnostic.csv"), row.names = FALSE)
}, error = function(e) failure <<- conditionMessage(e))
table <- if (length(checks)) do.call(rbind, checks) else data.frame()
utils::write.csv(table, file.path(out,"comparisons.csv"), row.names = FALSE)
hash <- function(path) tryCatch(digest::digest(file = path, algo = "sha256"), error = function(e) NA_character_)
summary <- list(status = if (is.null(failure) && nrow(table) > 0L && all(table$passed)) "passed" else "failed",
  comparisons = nrow(table), failed = if (nrow(table)) sum(!table$passed) else 0L, error = failure,
  negative_control = length(args) == 2L, package_path = tryCatch(find.package("inferguidR"), error=function(e) NULL),
  package_version = tryCatch(as.character(utils::packageVersion("inferguidR")), error=function(e) NULL),
  source_sha256 = Sys.getenv("INFERGUIDR_SOURCE_SHA256"), r_version = R.version.string,
  checker_sha256 = hash(script), reference_sha256 = hash(file.path(dirname(script),"reference-bootstrap.R")),
  scope = "Synthetic weighted rows: six metrics, class/unstratified/cluster resampling, 250 draws, one seed; independent arithmetic and type-7 interpolation share R multinomial RNG.",
  coverage_scope = "Exact binomial repeated-sampling coverage for limiting infinite-replicate percentile sensitivity intervals only; not finite-B coverage or general classifier certification.",
  validation_status = "development; no nominal-coverage certification")
jsonlite::write_json(summary, file.path(out,"summary.json"), auto_unbox = TRUE, pretty = TRUE, null = "null", na = "null")
if (summary$status != "passed") quit(status = 1L)
