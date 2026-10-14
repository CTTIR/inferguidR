.omnibus_code_env <- function(id) {
 e<-new.env();e$alternative<-"two.sided";e$confidence_level<-.95
 invisible(capture.output(sys.source(system.file("advisor","templates",paste0("tpl.",id,".v1.R"),package="inferguidR"),e)));e
}

test_that("one-way examples match independent closed-form balanced references", {
 x<-1:9;g<-rep(c("A","B","C"),each=3)
 refs<-list(anova_one_way=c(27,.001),anova_one_way_welch=c(162/7,(7/88)^2),kruskal_wallis=c(7.2,exp(-3.6)))
 for(id in names(refs)) {
  e<-.omnibus_code_env(id);r<-e$analyse(x,g)
  expect_equal(r$statistic,refs[[id]][1]);expect_equal(r$p_value,refs[[id]][2])
  expect_null(r$confidence_interval);expect_identical(r$interval_status,"not_computed")
  expect_equal(e$analyse(3*x+10,g)$statistic,r$statistic)
  expect_equal(e$analyse(x,factor(g,levels=c("C","B","A","unused")))$p_value,r$p_value)
  expect_error(e$analyse(x,g,"greater"))
  for(bad in list(c(NA,x[-1]),c(Inf,x[-1]),matrix(x,3))) expect_error(e$analyse(bad,g))
  for(bad in list(rep("A",9),1:9,g[-1],c(g[-1],NA),c(g[-1],""))) expect_error(e$analyse(x,bad))
  for(conf in list(0,1,NA,Inf,c(.9,.95))) expect_error(e$analyse(x,g,conf.level=conf))
  deg<-e$analyse(rep(2,9),g);expect_identical(deg$status,"not_applicable");expect_null(deg$p_value)
 }
})

test_that("rank omnibus applies tie correction and does not imply median inference", {
 e<-.omnibus_code_env("kruskal_wallis")
 # Average ranks 1.5,1.5,4,4,4,6; group sums 7 and 14.
 r<-e$analyse(c(1,1,2,2,2,3),rep(c("A","B"),each=3))
 h<-(12/(6*7)*(49/3+196/3)-21)/(1-(6+24)/(216-6))
 expect_equal(r$statistic,h);expect_equal(r$p_value,2*pnorm(-sqrt(h)))
 expect_match(r$assumptions,"identical distributions",fixed=TRUE)
})

test_that("variance degeneracy is explicit for pooled and Welch mean comparisons", {
 x<-c(1,1,1,2,3,4);g<-rep(c("A","B"),each=3)
 pooled<-.omnibus_code_env("anova_one_way")$analyse(x,g)
 welch<-.omnibus_code_env("anova_one_way_welch")$analyse(x,g)
 expect_identical(pooled$status,"ok");expect_true(is.finite(pooled$p_value))
 expect_identical(welch$status,"not_applicable");expect_null(welch$p_value)
})

test_that("omnibus methods export actual routes without inventing directional contrasts", {
 ids<-c("anova_one_way","anova_one_way_welch","kruskal_wallis")
 for(id in ids) {
  a<-c(q.aim="compare",q.stage="analysis",q.outcome_type="continuous",q.unit="individuals",
   q.structure="independent",q.independence_check="no",q.n_factors="one",q.n_groups="more_than_two",
   q.covariates="no",q.estimand_continuous=if(id=="kruskal_wallis") "rank" else "mean")
  if(id!="kruskal_wallis") a<-c(a,q.variance=if(id=="anova_one_way") "justified" else "not_justified")
  a<-c(a,q.missingness="none");s<-advisor_replay(a)
  expect_identical(advisor_recommendation(s)$preferred$analysis_variant_id,id)
  code<-advisor_r_code(s);e<-new.env();invisible(capture.output(eval(parse(text=code),e)))
  expect_true(is.finite(e$example_result$p_value))
  expect_error(advisor_r_code(s,alternative="less"),"omnibus")
 }
})
