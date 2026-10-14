.ancova_template <- function() {
 e<-new.env();e$alternative<-"two.sided";e$confidence_level<-.95;e$example_group_count<-3L
 invisible(capture.output(sys.source(system.file("advisor","templates","tpl.ancova_baseline_adjusted.v1.R",package="inferguidR"),e)));e
}
.ancova_case <- function(c) {
 x<-do.call(rbind,lapply(c$x,unlist));colnames(x)<-paste0("baseline",seq_len(ncol(x)))
 data.frame(id=seq_len(nrow(x)),y=unlist(c$y),group=paste0("g",unlist(c$groups)),x)
}
test_that("additive ANCOVA matches independent normal-equation contrasts and tails", {
 e<-.ancova_template()
 cases<-jsonlite::fromJSON(test_path("fixtures","ancova-template-cases.json"),simplifyVector=FALSE)
 refs<-jsonlite::fromJSON(test_path("fixtures","ancova-template-independent.json"),simplifyVector=FALSE)
 for(i in seq_along(cases)) {
  c<-cases[[i]];d<-.ancova_case(c);covariates<-paste0("baseline",seq_len(c$nc));groups<-paste0("g",0:(c$k-1))
  r<-e$analyse(d,covariates,groups,TRUE);ref<-refs[[i]]
  for(key in names(ref$group_test)) expect_equal(r$group_test[[key]],as.numeric(ref$group_test[[key]]),tolerance=2e-11)
  for(j in seq_along(ref$adjusted_means)) for(key in names(ref$adjusted_means[[j]]))
   expect_equal(r$adjusted_means[[key]][j],as.numeric(ref$adjusted_means[[j]][[key]]),tolerance=2e-11)
  for(j in seq_along(ref$differences)) for(key in names(ref$differences[[j]]))
   expect_equal(r$differences[[key]][j],as.numeric(ref$differences[[j]][[key]]),tolerance=2e-11)
  expect_equal(unname(r$reference_covariates),as.numeric(unlist(ref$reference)),tolerance=2e-13)
  expect_equal(r$group_test$df_numerator,c$k-1)
  expect_equal(r$group_test$df_denominator,nrow(d)-c$k-c$nc)
  expect_identical(r$n_excluded,0L)
 }
})
test_that("ANCOVA design coding and complete-case policy preserve the declared estimand", {
 e<-.ancova_template();d<-e$example_data;levels<-c("control","intervention","comparison")
 r<-e$analyse(d,"baseline",levels,TRUE)
 reordered<-e$analyse(d[nrow(d):1,],"baseline",levels,TRUE)
 expect_equal(r$group_test,reordered$group_test,tolerance=1e-12)
 expect_equal(r$adjusted_means,reordered$adjusted_means,tolerance=1e-12)
 reverse<-e$analyse(d,"baseline",rev(levels),TRUE)
 expect_equal(r$group_test,reverse$group_test,tolerance=1e-12)
 expect_equal(r$adjusted_means$estimate,rev(reverse$adjusted_means$estimate),tolerance=1e-12)
 changed<-d;changed$baseline<-100+3*d$baseline
 transformed<-e$analyse(changed,"baseline",levels,TRUE)
 expect_equal(r$differences,transformed$differences,tolerance=1e-11)
 expect_equal(r$adjusted_means,transformed$adjusted_means,tolerance=1e-11)
 missing<-d;missing$baseline[2]<-NA
 expect_error(e$analyse(missing,"baseline",levels,TRUE),"complete_cases")
 excluded<-e$analyse(missing,"baseline",levels,TRUE,missing="complete_cases")
 direct<-e$analyse(d[-2,],"baseline",levels,TRUE)
 expect_equal(excluded$group_test,direct$group_test)
 expect_equal(excluded$adjusted_means,direct$adjusted_means)
 expect_identical(excluded$n_excluded,1L)
 expect_error(e$analyse(d,"baseline",levels,FALSE))
 expect_error(e$analyse(d,"baseline",levels,TRUE,alternative="greater"))
 expect_error(e$analyse(d,"baseline",levels,TRUE,conf.level=1))
 duplicate<-d;duplicate$id[2]<-duplicate$id[1]
 expect_error(e$analyse(duplicate,"baseline",levels,TRUE))
 aliased<-d;aliased$copy<-aliased$baseline
 expect_error(e$analyse(aliased,c("baseline","copy"),levels,TRUE))
 exact<-d;exact$y<-2+exact$baseline
 expect_error(e$analyse(exact,"baseline",levels,TRUE),"exact fit")
 absent<-d[d$group!="comparison",]
 expect_error(e$analyse(absent,"baseline",levels,TRUE))
 hostile<-d;names(hostile)[names(hostile)=="baseline"]<-"stop('not code')"
 expect_equal(e$analyse(hostile,"stop('not code')",levels,TRUE)$group_test,r$group_test)
 expect_match(r$multiplicity,"unadjusted")
 expect_match(r$assumptions,"not individual prediction intervals")
})
test_that("public ANCOVA examples execute for ENDE without directional substitution", {
 answers<-c(q.aim="compare",q.stage="analysis",q.outcome_type="continuous",q.unit="individuals",q.structure="independent",
  q.independence_check="no",q.n_factors="one",q.n_groups="two",q.covariates="yes",q.missingness="none")
 for(groups in c("two","more_than_two")) for(lang in c("en","de")) {
  answers[["q.n_groups"]]<-groups
  s<-advisor_replay(answers,language=lang)
  expect_identical(advisor_recommendation(s)$preferred$analysis_variant_id,"ancova_baseline_adjusted")
  e<-new.env();invisible(capture.output(eval(parse(text=advisor_r_code(s)),e)))
  expect_identical(e$example_result$status,"ok")
  expect_length(e$example_result$group_levels,if(groups=="two") 2L else 3L)
  expect_error(advisor_r_code(s,alternative="greater"),class="inferguidR_invalid_input")
 }
})
