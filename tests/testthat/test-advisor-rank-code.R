.rank_code_env <- function(id,alt) {
 e<-new.env();e$alternative<-alt;e$confidence_level<-.95
 file<-system.file("advisor","templates",paste0("tpl.",id,".v1.R"),package="inferguidR")
 expect_true(nzchar(file));invisible(capture.output(sys.source(file,e)));e
}
.rank_expected_p <- function(z,alt) switch(alt,less=pnorm(z),greater=pnorm(z,lower.tail=FALSE),two.sided=2*pnorm(-abs(z)))

test_that("signed-rank examples agree with sign enumeration and Walsh location", {
 for(id in c("signed_rank_wilcoxon_one_sample","signed_rank_wilcoxon")) for(alt in c("less","greater","two.sided")) {
  e<-.rank_code_env(id,alt);d<-c(-1,2,-3,4,5,-6)
  signs<-expand.grid(rep(list(c(0,1)),6));support<-drop(as.matrix(signs)%*%(1:6))
  observed<-sum(rank(abs(d))[d>0]);lo<-mean(support<=observed);hi<-mean(support>=observed)
  expected<-switch(alt,less=lo,greater=hi,two.sided=min(1,2*min(lo,hi)))
  r<-e$signed_analysis(d,0,alt,.8,TRUE)
  walsh<-outer(d,d,"+")/2;walsh<-walsh[upper.tri(walsh,diag=TRUE)]
  expect_equal(r$p_value,expected,tolerance=1e-12)
  expect_equal(r$pseudomedian,median(walsh),tolerance=1e-10)
  expect_identical(r$exact,TRUE);expect_equal(r$requested_confidence_level,.8)
  expect_true(r$confidence_interval[1]<=r$confidence_interval[2])
  expect_error(e$signed_analysis(c(1,1,2),0,alt,.95,TRUE))
  expect_error(e$signed_analysis(c(0,1,2),0,alt,.95,TRUE))
  expect_error(e$signed_analysis(c(NA,1),0,alt,.95,FALSE))
  expect_error(e$signed_analysis(matrix(1:4,2),0,alt,.95,FALSE))
  expect_identical(e$signed_analysis(c(0,0),0,alt,.95,FALSE)$status,"not_applicable")
  d<-c(-1,1,-2,3,4,0);a<-d[d!=0];n<-length(a);ranks<-rank(abs(a))
  delta<-sum(ranks[a>0])-n*(n+1)/4;t<-as.numeric(table(abs(a)))
  se<-sqrt(n*(n+1)*(2*n+1)/24-sum(t^3-t)/48)
  correction<-switch(alt,less=-.5,greater=.5,two.sided=sign(delta)*.5)
  expected<-.rank_expected_p((delta-correction)/se,alt)
  r<-e$signed_analysis(d,0,alt,.8,FALSE)
  expect_equal(r$p_value,expected,tolerance=1e-12)
  expect_equal(r$n_zero,1);expect_true(r$tied_absolute_ranks)
  if(id=="signed_rank_wilcoxon") {
   expect_error(e$analyse(data.frame(id=c(1,1),x=c(1,2),y=0)))
   expect_error(e$analyse(data.frame(id=1:2,x=c(1,Inf),y=0)))
  }
 }
})

test_that("rank-sum examples distinguish descriptive superiority from shift intervals", {
 for(alt in c("less","greater","two.sided")) {
  e<-.rank_code_env("rank_sum_mann_whitney",alt)
  x<-c(1,4,6);y<-c(2,3,5,7);nx<-length(x);ny<-length(y)
  support<-colSums(combn(seq_len(nx+ny),nx))-nx*(nx+1)/2
  u<-sum(outer(x,y,">"));lo<-mean(support<=u);hi<-mean(support>=u)
  expected<-switch(alt,less=lo,greater=hi,two.sided=min(1,2*min(lo,hi)))
  r<-e$analyse(x,y,alt)
  expect_equal(r$p_value,expected,tolerance=1e-12)
  expect_equal(r$probability_superiority,mean(outer(x,y,">")),tolerance=1e-12)
  expect_null(r$confidence_interval);expect_null(r$shift_estimate)
  expect_identical(r$superiority_interval_status,"not_computed")
  r<-e$analyse(x,y,alt,location_shift=TRUE)
  expect_equal(r$shift_estimate,median(outer(x,y,"-")),tolerance=1e-10)
  expect_identical(r$interval_status,"shift_model_only")
  expect_error(e$analyse(c(1,1),c(2,3),alt))
  expect_error(e$analyse(c(1,NA),c(2,3),alt,exact=FALSE))
  expect_identical(e$analyse(c(1,1),c(1,1),alt,exact=FALSE)$status,"not_applicable")
  x<-c(1,2,2);y<-c(2,3,4,5);n<-length(x)+length(y);nx<-length(x);ny<-length(y)
  u<-sum(outer(x,y,">"))+.5*sum(outer(x,y,"=="));delta<-u-nx*ny/2;t<-as.numeric(table(c(x,y)))
  se<-sqrt(nx*ny/12*(n+1-sum(t^3-t)/(n*(n-1))))
  correction<-switch(alt,less=-.5,greater=.5,two.sided=sign(delta)*.5)
  r<-e$analyse(x,y,alt,exact=FALSE)
  expect_equal(r$p_value,.rank_expected_p((delta-correction)/se,alt),tolerance=1e-12)
  expect_equal(r$probability_superiority,u/(nx*ny),tolerance=1e-12)
 }
})

test_that("rank examples export through permitted alternatives without changing the estimand", {
 one<-c(q.aim="compare",q.stage="analysis",q.outcome_type="continuous",q.unit="individuals",q.structure="one_sample",q.reference_value="yes",q.direction="two_sided",q.missingness="none")
 paired<-c(q.aim="compare",q.stage="analysis",q.outcome_type="continuous",q.unit="individuals",q.structure="paired",q.pairing_id="yes",q.estimand_continuous="mean",q.direction="two_sided",q.missingness="none")
 independent<-c(q.aim="compare",q.stage="analysis",q.outcome_type="continuous",q.unit="individuals",q.structure="independent",q.independence_check="no",q.n_factors="one",q.n_groups="two",q.covariates="no",q.estimand_continuous="mean",q.variance="not_justified",q.direction="two_sided",q.missingness="none")
 cases<-list(signed_rank_wilcoxon_one_sample=one,signed_rank_wilcoxon=paired,rank_sum_mann_whitney=independent)
 for(id in names(cases)) for(direction in c("two_sided","one_sided_prespecified")) {
  ans<-cases[[id]];ans["q.direction"]<-direction
  s<-advisor_override(advisor_replay(ans),id);r<-advisor_recommendation(s)
  expect_identical(r$preferred$analysis_variant_id,id)
  expect_identical(r$r_implementation$template_status,"available")
  for(alt in if(direction=="two_sided") "two.sided" else c("greater","less")) {
   code<-advisor_r_code(s,alternative=alt);expect_match(code,r$recommendation_hash,fixed=TRUE)
   e<-new.env();invisible(capture.output(eval(parse(text=code),e)))
   expect_identical(e$example_result$alternative,alt);expect_true(is.finite(e$example_result$p_value))
  }
 }
})

test_that("exact rank intervals agree with independent permutation inversion", {
 for(alt in c("two.sided","less","greater")) {
  e<-.rank_code_env("signed_rank_wilcoxon_one_sample",alt)
  expected<-switch(alt,two.sided=c(-.5,7.5),less=c(-Inf,7),greater=c(.5,Inf))
  expect_equal(e$example_result$confidence_interval,expected,tolerance=1e-12)
  e<-.rank_code_env("rank_sum_mann_whitney",alt)
  expected<-switch(alt,two.sided=c(-5,7),less=c(-Inf,5),greater=c(-3,Inf))
  expect_equal(e$analyse(e$example_x,e$example_y,alt,location_shift=TRUE)$confidence_interval,expected,tolerance=1e-12)
 }
})
