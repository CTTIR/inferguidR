.paired_nominal_env <- function(method) {
 e<-new.env(parent=baseenv());e$alternative<-"two.sided";e$confidence_level<-.95
 invisible(capture.output(sys.source(system.file("advisor","templates",paste0("tpl.",method,".v1.R"),package="inferguidR"),e)));e
}
test_that("paired nominal examples match independent transition covariance calculations", {
 cases<-jsonlite::fromJSON(test_path("fixtures","paired-nominal-cases.json"),simplifyVector=FALSE)
 refs<-read.csv(test_path("fixtures","paired-nominal-independent.csv"))
 for(method in unique(refs$method)) {
  e<-.paired_nominal_env(method)
  for(i in which(refs$method==method)) {
   x<-do.call(rbind,lapply(cases[[refs$id[i]]],unlist));dimnames(x)<-rep(list(paste0("category",seq_len(nrow(x)))),2)
   r<-e$analyse(x);expect_identical(r$status,"ok");expect_equal(r$df,refs$df[i])
   expect_lte(abs(r$statistic-refs$statistic[i]),1e-12);expect_lte(abs(r$p_value-refs$p_value[i]),1e-12)
  }
 }
})
test_that("paired nominal examples preserve category identity and singularity boundaries", {
 for(method in c("stuart_maxwell_marginal_homogeneity","bowker_symmetry")) {
  e<-.paired_nominal_env(method);x<-e$example_counts
  expect_error(e$analyse(unname(x)))
  a<-x;colnames(a)<-rev(colnames(a));expect_error(e$analyse(a))
  a<-x;dimnames(a)<-rep(list(c("A","A","C")),2);expect_error(e$analyse(a))
  for(v in c(-1,.5,Inf,NA)) {a<-x;a[1]<-v;expect_error(e$analyse(a))}
  expect_error(e$analyse(x[1:2,1:2]));expect_error(e$analyse(x*0))
  expect_error(e$analyse(x,alternative="less"));expect_error(e$analyse(x,conf.level=1))
  for(v in c(-1,.5,NA,Inf)) expect_error(e$analyse(x,n_excluded=v))
  r<-e$analyse(x,n_excluded=3L);expect_identical(r$n_excluded,3L);expect_null(r$confidence_interval)
  expect_equal(e$analyse(t(x))$p_value,r$p_value)
  expect_equal(e$analyse(x[c(3,1,2),c(3,1,2)])$p_value,r$p_value)
  a<-x;a[]<-0;diag(a)<-10;r<-e$analyse(a)
  expect_identical(r$status,"unavailable");expect_null(r$p_value)
 }
 # Equal margins do not imply symmetry: a directed cycle has homogeneous margins.
 x<-matrix(c(10,9,1,1,10,9,9,1,10),3,byrow=TRUE,dimnames=rep(list(c("A","B","C")),2))
 expect_equal(.paired_nominal_env("stuart_maxwell_marginal_homogeneity")$analyse(x)$p_value,1)
 expect_lt(.paired_nominal_env("bowker_symmetry")$analyse(x)$p_value,.001)
})
test_that("paired nominal primary and override export in both languages", {
 a<-c(q.aim="compare",q.stage="analysis",q.outcome_type="nominal",q.unit="individuals",
   q.structure="paired",q.pairing_id="yes",q.missingness="none")
 for(lang in c("en","de")) for(method in c("stuart_maxwell_marginal_homogeneity","bowker_symmetry")) {
  s<-advisor_replay(a,language=lang);if(method=="bowker_symmetry") s<-advisor_override(s,method)
  expect_identical(advisor_recommendation(s)$preferred$analysis_variant_id,method)
  e<-new.env();invisible(capture.output(eval(parse(text=advisor_r_code(s)),e)))
  expect_identical(e$example_result$status,"ok");expect_null(e$example_result$confidence_interval)
  expect_error(advisor_r_code(s,alternative="greater"),class="inferguidR_invalid_input")
 }
})
