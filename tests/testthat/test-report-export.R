.report_test_envelope <- function() {
 s <- analysis_spec("t_means_one_sample","post_hoc",list(d=.5,alpha=.05,n_total=30),list(tails="two"))
 analysis_envelope(s,calculate_power(s))
}

test_that("portable reports round-trip without evaluating literal content", {
 e <- .report_test_envelope()
 x <- create_report_snapshot(e,title="Literal <script>stop('bad')</script> \\input{secret}")
 expect_s3_class(x,"report_snapshot")
 expect_identical(ig_to_json(ig_from_json(ig_to_json(x))),ig_to_json(x))
 old <- x$snapshot_hash
 e <- unclass(e); bad <- unclass(e$spec); bad$inputs$d <- .8
 e$spec <- structure(bad,class=class(.report_test_envelope()$spec)); e <- structure(e,class="analysis_envelope")
 expect_identical(x$snapshot_hash,old)
 expect_error(create_report_snapshot(e),class="inferguidR_invalid_input")
 y <- unclass(x);y$title <- "modified";y <- structure(y,class=class(x))
 expect_error(export_report(y,tempfile(),"txt"),class="inferguidR_invalid_input")
 for(fmt in c("json","txt")) {
  file <- tempfile();on.exit(unlink(file),add=TRUE)
  receipt <- export_report(x,file,fmt)
  expect_identical(receipt$snapshot_id,x$snapshot_id)
  txt <- paste(readLines(file),collapse="\n")
  expect_match(txt,x$snapshot_id,fixed=TRUE)
  expect_match(txt,"stop('bad')",fixed=TRUE)
  expect_error(export_report(x,file,fmt),class="inferguidR_invalid_input")
 }
 expect_error(create_report_snapshot(.report_test_envelope(),display_timezone="made_up"),class="inferguidR_invalid_input")
})

test_that("classifier reports omit row-level attachments and score-derived curves", {
 s <- classification_examples("scored_synthetic")$spec
 r <- evaluate_cutoff(s)
 x <- create_report_snapshot(analysis_envelope(s,r),report_kind="model_cutoff")
 raw <- .ig_report_document_to_list(x)
 expect_false(x$privacy$raw_rows_included)
 expect_null(raw$results[[1]]$models[[1]]$curve)
 expect_null(raw$results[[1]]$models[[1]]$operating_point$boundary_score)
 expect_null(raw$results[[1]]$models[[1]]$probability)
 expect_true(raw$reproducibility[[1]]$private_data_required)
 expect_null(attr(raw$specifications[[1]]$specification,"ig_cls_data"))
 expect_identical(ig_to_json(x),ig_to_json(ig_from_json(ig_to_json(x))))
 expect_equal(raw$results[[1]]$models[[1]]$operating_point$tp,r$models[[1]]$operating_point$tp)
})

test_that("incomplete reports are explicitly labelled", {
 e <- .report_test_envelope();e <- analysis_envelope(e$spec)
 expect_error(create_report_snapshot(e),class="inferguidR_invalid_input")
 x <- create_report_snapshot(e,report_kind="incomplete_settings")
 expect_identical(x$completeness,"incomplete")
 expect_identical(x$results[[1]]$status$code,"not_computed")
 expect_error(ig_from_json(sub("report_document/0.1.0","report_document/9.0.0",ig_to_json(x),fixed=TRUE)),class="inferguidR_invalid_input")
})

test_that("PDF renderer creates a real document and preserves active devices", {
 skip_if_not(report_pdf_capabilities()$available)
 x <- create_report_snapshot(.report_test_envelope())
 file <- tempfile(fileext=".pdf");on.exit(unlink(file),add=TRUE)
 prior <- grDevices::dev.cur()
 expect_no_error(export_report(x,file,"pdf"))
 expect_identical(grDevices::dev.cur(),prior)
 expect_identical(rawToChar(readBin(file,"raw",5L)),"%PDF-")
 expect_gt(file.size(file),1000)
})

test_that("report sessions retain frozen snapshots across analysis changes", {
 skip_if_not_installed("shiny")
 source <- shiny::reactiveVal(.report_test_envelope())
 eligible <- shiny::reactiveVal(TRUE)
 shiny::testServer(.ig_report_server,args=list(envelope=source,language=shiny::reactive("en"),eligible=eligible), {
  session$flushReact();session$setInputs(freeze=1)
  old <- snapshot()$snapshot_hash
  e <- source(); s <- update(e$spec,d=.8);source(analysis_envelope(s,calculate_power(s)))
  session$flushReact()
  expect_identical(snapshot()$snapshot_hash,old)
  eligible(FALSE);session$setInputs(freeze=2)
  expect_identical(snapshot()$snapshot_hash,old)
  expect_match(message(),"current settings")
  eligible(TRUE);session$setInputs(freeze=3)
  expect_false(identical(snapshot()$snapshot_hash,old))
 })
})

test_that("private score and identifier sentinels never enter report exports", {
 data <- data.frame(outcome=c(0,1,0,1),score=c(.123456789123456,.856789123456789,.23,.91),
                    person=c("private-person-A","private-person-B","private-person-C","private-person-D"))
 s <- classification_spec(data=data,truth="outcome",positive_label=1,scores="score",id="person",cutoff=.5)
 x <- create_report_snapshot(analysis_envelope(s,evaluate_cutoff(s)))
 text <- as.character(ig_to_json(x))
 expect_false(grepl("private-person",text,fixed=TRUE))
 expect_false(grepl("0.123456789123456",text,fixed=TRUE))
 expect_false(grepl("0.856789123456789",text,fixed=TRUE))
 forged <- unclass(x);forged$results[[1]]$scores <- list(.123)
 forged$snapshot_hash <- NULL;forged$snapshot_hash <- .ig_hash_value(forged)
 expect_error(ig_from_json(.ig_canonical_json(forged),validate=FALSE),class="inferguidR_invalid_input")
})
