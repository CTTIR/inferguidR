.binormal_template_env<-function() {
 e<-new.env(parent=baseenv());invisible(capture.output(sys.source(system.file("advisor","templates","tpl.classifier_binormal_theory.v1.R",package="inferguidR"),e)));e
}
test_that("binormal template retains independently computed hypothetical rates and counts", {
 e<-.binormal_template_env();refs<-jsonlite::fromJSON(test_path("fixtures","binormal-template-independent.json"),simplifyVector=FALSE)
 for(ref in refs) {
  r<-do.call(e$analyse,ref[c("auc","sigma","cutoff","prevalence","n_total")])
  expect_identical(r$status,"ok");expect_identical(r$count_kind,"expected")
  expect_equal(r$model$d,as.numeric(ref$d),tolerance=2e-12)
  for(key in c("tpr","fpr")) expect_equal(r$analytic_rates[[key]],as.numeric(ref[[key]]),tolerance=2e-12)
  for(key in c("tp","fn","fp","tn")) expect_equal(r$operating_point[[key]],as.numeric(ref$counts[[key]]),tolerance=2e-12)
  expect_equal(r$summaries$value[match("roc_auc",r$summaries$metric_id)],ref$auc)
  expect_identical(r$interval_status,"not_applicable");expect_null(r$confidence_interval)
  expect_length(r$result$models[[1]]$intervals$metric_id,0L)
 }
})
test_that("binormal scenarios distinguish scaling, prevalence endpoints and limits", {
 e<-.binormal_template_env();r<-e$analyse(.8,1.5,1,.1,1000);larger<-e$analyse(.8,1.5,1,.1,2000)
 expect_identical(r$analytic_rates,larger$analytic_rates)
 for(key in c("tp","fn","fp","tn")) expect_equal(larger$operating_point[[key]],2*r$operating_point[[key]])
 for(prev in c(0,1)) {
  z<-e$analyse(.8,1.5,1,prev,1000)
  expect_identical(z$analytic_rates,r$analytic_rates);expect_identical(z$count_kind,"expected")
  expect_equal(z$operating_point$n_events,prev*1000)
  expect_identical(z$metrics$status[match(if(prev==0) "tpr" else "fpr",z$metrics$metric_id)],"undefined")
 }
 for(cut in c(-Inf,Inf)) expect_error(e$analyse(.8,1.5,cut,.1,1000))
 for(auc in list(0,1,NA,Inf,c(.5,.8))) expect_error(e$analyse(auc,1,1,.1,1000))
 expect_error(e$analyse(.8,0,1,.1,1000));expect_error(e$analyse(.8,1,NA,.1,1000))
 expect_error(e$analyse(.8,1,1,-.1,1000));expect_error(e$analyse(.8,1,1,.1,0))
 expect_error(e$analyse(.8,1,c(0,1),.1,1000))
 expect_equal(e$analyse(.5,1,0,.1,1000)$analytic_rates$tpr,.5)
})
test_that("theoretical classifier exports cannot imply empirical inference", {
 for(lang in c("en","de")) {
  s<-advisor_replay(c(q.aim="performance",q.stage="analysis",q.perf_inputs="theoretical_model"),language=lang)
  expect_identical(advisor_recommendation(s)$preferred$analysis_variant_id,"classifier_binormal_theory")
  e<-new.env();invisible(capture.output(eval(parse(text=advisor_r_code(s)),e)))
  expect_identical(e$example_result$status,"ok");expect_identical(e$example_result$count_kind,"expected")
  expect_false(exists("alternative",e,inherits=FALSE));expect_false(exists("confidence_level",e,inherits=FALSE))
  expect_error(advisor_r_code(s,alternative="two.sided"),class="inferguidR_invalid_input")
 }
})
