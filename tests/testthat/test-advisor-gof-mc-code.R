.gof_mc_env<-function() {
 e<-new.env(parent=baseenv());e$confidence_level<-.95
 invisible(capture.output(sys.source(system.file("advisor","templates","tpl.chisq_goodness_of_fit_mc.v1.R",package="inferguidR"),e)));e
}
test_that("Monte Carlo goodness of fit agrees with independently enumerated tails", {
 e<-.gof_mc_env();refs<-jsonlite::fromJSON(test_path("fixtures","gof-mc-independent.json"),simplifyVector=FALSE)
 for(ref in refs) {
  x<-unlist(ref$observed);p<-unlist(ref$probabilities);names(x)<-names(p)<-LETTERS[seq_along(x)]
  r<-e$analyse(x,p,seed=104729L,B=19999L)
  expect_equal(r$statistic,ref$statistic,tolerance=1e-12)
  bound<-6*sqrt(ref$exact_p*(1-ref$exact_p)/r$B)+2/(r$B+1)
  expect_lte(abs(r$p_value-ref$exact_p),bound)
  expect_equal(r$p_value,(r$simulated_exceedances+1)/(r$B+1))
  expect_gte(r$p_value,r$minimum_p_value);expect_lte(r$p_value,1)
  expect_true(all(r$monte_carlo_interval>=0&r$monte_carlo_interval<=1))
  expect_identical(r$interval_status,"not_computed");expect_null(r$confidence_interval)
 }
})
test_that("Monte Carlo examples restore RNG state and validate fixed null probabilities", {
 e<-.gof_mc_env();x<-c(A=1,B=3,C=6);p<-c(A=.2,B=.3,C=.5)
 local({
  kind<-RNGkind();had<-exists(".Random.seed",.GlobalEnv,inherits=FALSE);if(had) saved<-get(".Random.seed",.GlobalEnv)
  on.exit({do.call(RNGkind,as.list(kind));if(had) assign(".Random.seed",saved,.GlobalEnv) else if(exists(".Random.seed",.GlobalEnv,inherits=FALSE)) rm(".Random.seed",envir=.GlobalEnv)})
  RNGkind("L'Ecuyer-CMRG");set.seed(123);before<-get(".Random.seed",.GlobalEnv);rng<-RNGkind()
  a<-e$analyse(x,p,seed=5,B=999);expect_identical(get(".Random.seed",.GlobalEnv),before);expect_identical(RNGkind(),rng)
  expect_identical(e$analyse(x,p,seed=5,B=999),a)
  rm(".Random.seed",envir=.GlobalEnv);invisible(e$analyse(x,p,seed=5,B=999));expect_false(exists(".Random.seed",.GlobalEnv,inherits=FALSE))
 })
 for(B in c(1,1000001,999.5,NA,Inf)) expect_error(e$analyse(x,p,seed=1,B=B))
 for(seed in list(-1,.5,NA,Inf,"1")) expect_error(e$analyse(x,p,seed=seed))
 expect_error(e$analyse(unname(x),p,seed=1));expect_error(e$analyse(x,p[c(2,1,3)],seed=1))
 expect_error(e$analyse(x,p*.99,seed=1));expect_error(e$analyse(c(A=0,B=0,C=0),p,seed=1))
 expect_error(e$analyse(x,c(A=0,B=.5,C=.5),seed=1));expect_error(e$analyse(c(A=1.5,B=3,C=6),p,seed=1))
 expect_error(e$analyse(x,p,seed=1,mc.level=1));expect_error(e$analyse(x,p,seed=1,conf.level=1))
 r<-e$analyse(c(A=1,B=0),c(A=.5,B=.5),seed=7,B=999,mc.level=1-.Machine$double.eps/2)
 expect_equal(r$p_value,1);expect_equal(r$monte_carlo_standard_error,0)
 expect_lt(r$monte_carlo_interval[1],1);expect_equal(r$monte_carlo_interval[2],1)
 expect_true(all(is.finite(r$monte_carlo_interval)))
})
test_that("Monte Carlo GOF remains an explicit omnibus advisor alternative", {
 for(lang in c("en","de")) {
  a<-c(q.aim="compare",q.stage="analysis",q.outcome_type="nominal",q.unit="individuals",q.structure="one_sample",q.reference_value="yes",q.direction="two_sided",q.missingness="none")
  s<-advisor_override(advisor_replay(a,language=lang),"chisq_goodness_of_fit_mc")
  expect_identical(advisor_recommendation(s)$preferred$analysis_variant_id,"chisq_goodness_of_fit_mc")
  e<-new.env();invisible(capture.output(eval(parse(text=advisor_r_code(s)),e)));expect_identical(e$example_result$status,"ok")
  expect_error(advisor_r_code(s,alternative="greater"),class="inferguidR_invalid_input")
  expect_equal(e$example_result$seed,104729L)
 }
})
test_that("GOF Monte Carlo is not a transitive replacement for an ordinal sign test", {
 s<-advisor_replay(c(q.aim="compare",q.stage="analysis",q.outcome_type="ordinal",q.unit="individuals",q.structure="one_sample",q.reference_value="yes",q.direction="two_sided",q.missingness="none"))
 expect_error(advisor_override(s,"chisq_goodness_of_fit_mc"),class="inferguidR_override_refused")
})
