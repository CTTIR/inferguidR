.counts_template_env <- function() {
 e<-new.env(parent=baseenv());e$confidence_level<-.95
 invisible(capture.output(sys.source(system.file("advisor","templates","tpl.diagnostic_2x2_fixed_cutoff.v1.R",package="inferguidR"),e)));e
}
.counts_template_call <- function(e,counts,...) do.call(e$analyse,c(counts,list(sampling="cross_sectional",cutoff_origin="predeclared_rule",evaluation_split="same_data"),list(...)))
test_that("diagnostic counts preserve independent fractions and exact binomial intervals", {
 e<-.counts_template_env();refs<-jsonlite::fromJSON(test_path("fixtures","classifier-counts-independent.json"),simplifyVector=FALSE)
 for(ref in refs) {
  r<-.counts_template_call(e,ref$counts);expect_identical(r$status,"ok")
  for(key in names(ref$metrics)) {
   j<-match(key,r$metrics$metric_id);value<-ref$metrics[[key]]
   if(is.null(value)) {expect_identical(r$metrics$status[j],"undefined");expect_true(is.na(r$metrics$value[j]))}
   else if(value=="infinity") {expect_identical(r$metrics$status[j],"positive_infinity");expect_true(is.na(r$metrics$value[j]))}
   else {expect_identical(r$metrics$status[j],"ok");expect_equal(r$metrics$value[j],as.numeric(value),tolerance=2e-12)}
  }
  for(key in names(ref$intervals)) {
   j<-match(key,r$intervals$metric_id);value<-ref$intervals[[key]]
   if(is.null(value)) expect_identical(r$intervals$status[j],"undefined") else {
    expect_identical(r$intervals$status[j],"ok")
    expect_equal(c(r$intervals$lower[j],r$intervals$upper[j]),as.numeric(unlist(value)),tolerance=2e-12)
   }
  }
  expect_true(all(r$summaries$status=="not_applicable"))
 }
})
test_that("counts templates require fixed-rule provenance and distinguish prevalence scenarios", {
 e<-.counts_template_env();cells<-list(tp=40,fn=10,fp=20,tn=130)
 for(v in list(-1,NA,Inf,.5,c(1,2),"2")) {a<-cells;a$tp<-v;expect_error(.counts_template_call(e,a))}
 expect_error(e$analyse(1,2,3,4))
 expect_error(e$analyse(1,2,3,4,"cross_sectional","same_data","same_data"))
 expect_error(e$analyse(1,2,3,4,"cross_sectional","held_out_tuning","same_data"))
 expect_error(.counts_template_call(e,cells,prevalence=.1))
 expect_error(.counts_template_call(e,cells,prevalence_source="external"))
 expect_error(.counts_template_call(e,cells,conf.level=1))
 r<-.counts_template_call(e,cells,prevalence=.1,prevalence_source="assumed")
 expect_identical(r$specification$prevalence$source,"assumed");expect_false(is.null(r$scenario))
 expect_equal(unname(r$counts),unname(unlist(cells)))
 expect_identical(r$result$models[[1]]$count_kind,"observed")
 for(sampling in c("case_control","unknown")) {
  r<-e$analyse(40,10,20,130,sampling,"held_out_tuning","external")
  expect_identical(r$specification$prevalence$sampling,sampling)
  expect_identical(r$specification$operating_rule$selection_data,"held_out_tuning")
 }
 expect_identical(e$analyse(40,10,20,130,"cross_sectional","held_out_tuning","held_out_test")$status,"ok")
 r<-.counts_template_call(e,cells,conf.level=1-.Machine$double.eps/2)
 selected<-r$intervals$status=="ok";expect_true(all(is.finite(r$intervals$lower[selected])));expect_true(all(is.finite(r$intervals$upper[selected])))
})
test_that("counts-only advisor exports explicit fixed-rule evaluation in both languages", {
 for(lang in c("en","de")) {
  s<-advisor_replay(c(q.aim="performance",q.stage="analysis",q.perf_inputs="counts_2x2",q.missingness="none"),language=lang)
  expect_identical(advisor_recommendation(s)$preferred$analysis_variant_id,"diagnostic_2x2_fixed_cutoff")
  code<-advisor_r_code(s);e<-new.env();invisible(capture.output(eval(parse(text=code),e)))
  expect_identical(e$example_result$status,"ok");expect_equal(e$confidence_level,.95)
  expect_false(exists("alternative",e,inherits=FALSE))
  expect_error(advisor_r_code(s,alternative="two.sided"),class="inferguidR_invalid_input")
 }
})
