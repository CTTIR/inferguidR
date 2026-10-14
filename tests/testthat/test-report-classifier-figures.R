test_that("classifier report mosaics and prevalence curves retain aggregate mathematics", {
 s<-classification_spec(counts=c(tp=40,fn=10,fp=20,tn=130),cutoff=.5,prevalence=.1)
 x<-create_report_snapshot(analysis_envelope(s,evaluate_cutoff(s)))
 mosaics<-Filter(function(v) v$plot_type=="classifier_mosaic",x$visualisations)
 expect_length(mosaics,2L)
 areas<-vapply(mosaics[[1]]$data$rectangles,function(r)(r$x1-r$x0)*(r$y1-r$y0),0)
 expect_equal(areas,c(.2,.05,.1,.65),tolerance=1e-14)
 expect_identical(mosaics[[2]]$data$count_kind,"expected")
 curve<-Filter(function(v)v$plot_type=="classifier_prevalence",x$visualisations)[[1]]
 at<-function(pi) curve$data$prevalence_grid[[which(vapply(curve$data$prevalence_grid,`[[`,0,"prevalence")==pi)]]
 expect_equal(at(.1)$ppv,.4,tolerance=1e-14)
 expect_equal(at(.1)$npv,.975,tolerance=1e-14)
 expect_equal(at(.25)$ppv,2/3,tolerance=1e-14)
 expect_equal(at(.25)$npv,13/14,tolerance=1e-14)
 expect_identical(ig_to_json(ig_from_json(ig_to_json(x))),ig_to_json(x))
 expect_false(x$privacy$raw_rows_included)
 expect_false(x$privacy$score_curve_tables_included)
})

test_that("classifier report curves preserve undefined denominators and absent classes", {
 for(counts in list(c(tp=0,fn=10,fp=0,tn=20),c(tp=10,fn=0,fp=20,tn=0),c(tp=0,fn=0,fp=3,tn=7))) {
  s<-classification_spec(counts=counts,cutoff=.5)
  x<-create_report_snapshot(analysis_envelope(s,evaluate_cutoff(s)))
  curve<-Filter(function(v)v$plot_type=="classifier_prevalence",x$visualisations)[[1]]
  if(sum(counts[c("tp","fn")])==0) {
   expect_identical(curve$status,"not_applicable")
   expect_length(curve$data$prevalence_grid,0)
  } else {
   metric<-if(counts[["tp"]]==0) "ppv" else "npv"
   expect_true(all(vapply(curve$data$prevalence_grid,function(r)is.null(r[[metric]]),FALSE)))
  }
  expect_identical(ig_to_json(ig_from_json(ig_to_json(x))),ig_to_json(x))
 }
})

test_that("classifier aggregate figures reject altered geometry or scenario values", {
 s<-classification_spec(counts=c(tp=40,fn=10,fp=20,tn=130),cutoff=.5,prevalence=.1)
 x<-create_report_snapshot(analysis_envelope(s,evaluate_cutoff(s)))
 for(plot in c("classifier_mosaic","classifier_prevalence")) {
  y<-unclass(x);j<-which(vapply(y$visualisations,function(v)v$plot_type==plot,FALSE))[1]
  if(plot=="classifier_mosaic") y$visualisations[[j]]$data$rectangles[[1]]$x1<-.9 else y$visualisations[[j]]$data$prevalence_grid[[11]]$ppv<-.9
  y$snapshot_hash<-NULL;y$snapshot_hash<-.ig_hash_value(y)
  expect_error(ig_from_json(.ig_canonical_json(y)),class="inferguidR_invalid_input")
 }
})

test_that("long model subtitles retain basis without entering the title line", {
  for (basis in c("Original data", "Ausgangsdaten", "Scenario", "Szenario",
                  "Theoretical expected values", "Theoretische Erwartungswerte")) {
    label <- .ig_report_classifier_subtitle(strrep("LongModelLabel", 30), basis)
    expect_false(grepl("\\n", label))
    expect_true(endsWith(label, paste("[...] |", basis)))
  }
  expect_identical(.ig_report_classifier_subtitle("model_a", "Original data"), "model_a | Original data")
})
