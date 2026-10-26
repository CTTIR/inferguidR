rp_envelope <- function() {
  s <- classification_spec(data.frame(y=c(0,1,0,1,1,0),p=c(.1,.2,.4,.6,.8,.9)),
    truth="y",positive_label=1,scores="p",score_scale="probability",cutoff=.5)
  analysis_envelope(s,evaluate_cutoff(s,threshold_probability=c(.2,.5,.8),calibration_bins=3))
}

test_that("probability report aggregates require explicit opt-in and preserve old schemas", {
  e <- rp_envelope()
  default <- create_report_snapshot(e)
  expect_identical(default$schema_version,"inferguidR.report_document/0.1.0")
  expect_null(default$probability_summaries)
  expect_identical(ig_from_json(ig_to_json(default)),default)
  for (lang in c("en","de")) {
    doc <- create_report_snapshot(e,locale=lang,include_probability_summaries=TRUE)
    expect_identical(doc$schema_version,"inferguidR.report_document/0.4.0")
    expect_true(doc$privacy$probability_aggregates_included)
    expect_false(doc$privacy$raw_rows_included)
    expect_false(doc$privacy$score_curve_tables_included)
    expect_length(doc$probability_summaries,1)
    rows <- doc$probability_summaries[[1]]
    expect_equal(vapply(rows$calibration,`[[`,0,"mean_predicted"),c(.15,.5,.85))
    expect_equal(vapply(rows$calibration,`[[`,0,"events"),rep(1,3))
    expect_equal(vapply(rows$net_benefit,`[[`,0,"net_benefit"),c(3/6-2/6/4,2/6-1/6,1/6-4/6))
    expect_identical(ig_from_json(ig_to_json(doc)),doc)
    serialized <- as.character(ig_to_json(doc))
    for (private in c('"min_predicted"','"max_predicted"','"raw_rows"','"boundary_score"','"ig_cls_data"'))
      expect_false(grepl(private,serialized,fixed=TRUE))
    files <- tempfile(fileext=c(".txt",".pdf"))
    on.exit(unlink(files),add=TRUE)
    export_report(doc,files[1],"txt")
    expect_true(any(grepl("mean_predicted",readLines(files[1]),fixed=TRUE)))
    if (report_pdf_capabilities()$available) {
      export_report(doc,files[2],"pdf")
      expect_gt(file.info(files[2])$size,1000)
    }
  }
  expect_error(create_report_snapshot(e,include_probability_summaries="yes"),class="inferguidR_invalid_input")
})

test_that("resealed probability reports reject invalid arithmetic, thresholds, duplicated plots and private fields", {
  doc <- unclass(create_report_snapshot(rp_envelope(),include_probability_summaries=TRUE))
  changed <- function(mutator) {
    d <- mutator(doc); d$snapshot_hash <- NULL; d$snapshot_hash <- ig_hash(d)
    ig_from_json(ig_to_json(d),validate=FALSE)
  }
  for (mutator in list(
    function(d) {d$probability_summaries[[1]]$calibration[[1]]$events <- 0;d},
    function(d) {d$probability_summaries[[1]]$calibration[[1]]$mean_predicted <- .9;d},
    function(d) {d$probability_summaries[[1]]$net_benefit[[1]]$net_benefit <- 99;d},
    function(d) {d$probability_summaries[[1]]$net_benefit[[1]]$threshold_probability <- .3;d},
    function(d) {d$probability_summaries[[1]]$net_benefit <- list();d},
    function(d) {d$probability_summaries <- list();d},
    function(d) {
      b <- d$probability_summaries[[1]]$calibration
      b[[1]]$mean_predicted <- .5; b[[2]]$mean_predicted <- .15
      d$probability_summaries[[1]]$calibration <- b;d
    },
    function(d) {
      b <- d$probability_summaries[[1]]$net_benefit
      b[[1]]$tp <- 1; b[[1]]$net_benefit <- (1-b[[1]]$fp*b[[1]]$harm_weight)/b[[1]]$n_total
      d$probability_summaries[[1]]$net_benefit <- b;d
    },
    function(d) {d$probability_summaries[[1]]$calibration[[1]]$min_predicted <- .1;d},
    function(d) {d$visualisations <- c(d$visualisations,tail(d$visualisations,1));d},
    function(d) {d$privacy$probability_aggregates_included <- FALSE;d},
    function(d) {d$schema_version <- "inferguidR.report_document/0.3.0";d}
  )) expect_error(changed(mutator),class="inferguidR_invalid_input")
})

test_that("report opt-in freezes only current probability settings", {
  shiny::testServer(.ig_classifier_workspace_server,args=list(language=shiny::reactive("en")),{
    session$setInputs(example="scored_synthetic",selection="manual",cutoff="0.3",scenario=FALSE,
      threshold_probability="0.2; 0.5",calibration_bins="3",`report-include_probability_summaries`=TRUE)
    session$setInputs(evaluate=1)
    session$setInputs(`report-freeze`=1)
    old <- report_module$snapshot()
    expect_identical(old$schema_version,"inferguidR.report_document/0.4.0")
    expect_true(length(old$probability_summaries)>0)
    session$setInputs(threshold_probability="0.2; 0.6")
    session$setInputs(`report-freeze`=2)
    expect_identical(report_module$snapshot(),old)
    session$setInputs(evaluate=2)
    session$setInputs(`report-freeze`=3)
    fresh <- report_module$snapshot()
    expect_false(identical(fresh$snapshot_hash,old$snapshot_hash))
    expect_equal(vapply(fresh$probability_summaries[[1]]$net_benefit,`[[`,0,"threshold_probability"),c(.2,.6))
    expect_identical(ig_from_json(ig_to_json(old)),old)
  })
})
