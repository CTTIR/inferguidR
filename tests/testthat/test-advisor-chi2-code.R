.chi2_code_env <- function(id, alternative="two.sided") {
 e<-new.env();e$alternative<-alternative;e$confidence_level<-.95
 file<-system.file("advisor","templates",paste0("tpl.",id,".v1.R"),package="inferguidR")
 expect_true(nzchar(file));invisible(capture.output(sys.source(file,e)));e
}

test_that("Pearson examples preserve omnibus statistics, effects and sparse-cell warnings", {
 e<-.chi2_code_env("chisq_goodness_of_fit")
 r<-e$analyse(c(10,20,30),c(.2,.3,.5));s<-4/12+4/18
 expect_equal(r$statistic,s);expect_equal(r$p_value,exp(-s/2));expect_equal(r$cohen_w,sqrt(s/60))
 expect_null(r$confidence_interval);expect_equal(r$df,2)
 expect_true(length(e$analyse(c(1,2,3),c(.2,.3,.5))$warnings)>0)
 for(x in list(c(-1,2,3),c(NA,2,3),c(Inf,2,3),c(1.5,2,3),c(0,0,0),matrix(1:3,1)))
  expect_error(e$analyse(x,c(.2,.3,.5)))
 for(p in list(c(0,.5,.5),c(.2,.2,.2),c(NA,.3,.5),matrix(c(.2,.3,.5),1)))
  expect_error(e$analyse(c(10,20,30),p))
 expect_error(e$analyse(c(10,20,30),c(.2,.3,.5),"greater"))
 e<-.chi2_code_env("chisq_independence_rxc")
 x<-rbind(c(10,20,30),c(25,20,15));ex<-rbind(c(17.5,20,22.5),c(17.5,20,22.5))
 r<-e$analyse(x);s<-2*(7.5^2/17.5+7.5^2/22.5)
 expect_equal(r$expected,ex);expect_equal(r$statistic,s);expect_equal(r$p_value,exp(-s/2))
 expect_equal(r$cramer_v,sqrt(s/120));expect_false(r$continuity_correction)
 for(x in list(matrix(0,2,2),matrix(1:4,1),rbind(c(1,0),c(2,0)),rbind(c(1,2),c(3,NA)),rbind(c(1,2),c(3,Inf)),matrix(.5,2,2)))
  expect_error(e$analyse(x))
})

test_that("asymptotic McNemar variants preserve discordance, correction and undefined cases", {
 for(id in c("mcnemar_asymptotic_uncorrected","mcnemar_asymptotic_cc")) {
  e<-.chi2_code_env(id);correct<-id=="mcnemar_asymptotic_cc"
  for(pair in list(c(20,10),c(10,10),c(1,0),c(0,12))) {
   b<-pair[1];cc<-pair[2];r<-e$analyse(rbind(c(40,b),c(cc,30)))
   s<-(if(correct) max(0,abs(b-cc)-1) else abs(b-cc))^2/(b+cc)
   expect_equal(r$statistic,s);expect_equal(r$p_value,2*pnorm(-sqrt(s)))
   expect_equal(r$marginal_difference,(b-cc)/(70+b+cc));expect_identical(r$continuity_correction,correct)
  }
  r<-e$analyse(diag(c(10,20)));expect_identical(r$status,"not_applicable");expect_null(r$p_value)
  for(x in list(matrix(0,2,2),matrix(1,2,3),matrix(.5,2,2),matrix(Inf,2,2),matrix(NA_real_,2,2)))
   expect_error(e$analyse(x))
  expect_error(e$analyse(matrix(1,2,2),"greater"))
 }
})

test_that("independent proportion example separates score test from difference interval", {
 for(alt in c("two.sided","less","greater")) {
  e<-.chi2_code_env("prop_two_independent_pearson_chisq",alt)
  r<-e$analyse(rbind(c(30,70),c(45,55)),alt,.9)
  delta<--.15;z<-delta/sqrt(.375*.625*.02);se<-sqrt((.3*.7+.45*.55)/100)
  q<-qnorm(if(alt=="two.sided") .95 else .9)
  ci<-switch(alt,two.sided=delta+c(-1,1)*q*se,less=c(-1,delta+q*se),greater=c(delta-q*se,1))
  prob<-switch(alt,two.sided=2*pnorm(-abs(z)),less=pnorm(z),greater=pnorm(z,lower.tail=FALSE))
  expect_equal(r$p_value,prob);expect_equal(r$confidence_interval,ci);expect_equal(r$proportion_difference,delta)
  expect_identical(r$alternative,alt);expect_false(r$continuity_correction)
  for(x in list(matrix(1,2,3),matrix(.5,2,2),rbind(c(0,10),c(0,20)),matrix(NA_real_,2,2)))
   expect_error(e$analyse(x))
 }
})

test_that("count templates export through actual advisor recommendations and overrides", {
 routes<-list(chisq_goodness_of_fit = c(q.aim = "compare", q.stage = "planning", 
q.outcome_type = "ordinal", q.unit = "individuals", q.structure = "one_sample", 
q.reference_value = "yes", q.direction = "two_sided", q.missingness = "some"
), prop_two_independent_pearson_chisq = c(q.aim = "compare", 
q.stage = "planning", q.outcome_type = "binary", q.unit = "individuals", 
q.structure = "independent", q.independence_check = "no", q.n_factors = "one", 
q.n_groups = "two", q.covariates = "no", q.binary_margins = "group_sizes_fixed", 
q.direction = "two_sided", q.missingness = "some"), mcnemar_asymptotic_uncorrected = c(q.aim = "compare", 
q.stage = "planning", q.outcome_type = "binary", q.unit = "individuals", 
q.structure = "paired", q.pairing_id = "yes", q.direction = "two_sided", 
q.missingness = "some"), mcnemar_asymptotic_cc = c(q.aim = "compare", 
q.stage = "planning", q.outcome_type = "binary", q.unit = "individuals", 
q.structure = "paired", q.pairing_id = "yes", q.direction = "two_sided", 
q.missingness = "some"))
 routes$chisq_independence_rxc<-c(q.aim="compare",q.stage="analysis",q.outcome_type="binary",q.unit="individuals",q.structure="independent",q.independence_check="no",q.n_factors="one",q.n_groups="more_than_two",q.covariates="no",q.small_expected="no",q.missingness="none")
 for(id in names(routes)) {
  chosen<-advisor_override(advisor_replay(routes[[id]]),id)
  expect_false(is.null(chosen),info=id)
  if(is.null(chosen)) next
  r<-advisor_recommendation(chosen);expect_identical(r$preferred$analysis_variant_id,id)
  code<-advisor_r_code(chosen);expect_match(code,r$recommendation_hash,fixed=TRUE)
  e<-new.env();invisible(capture.output(eval(parse(text=code),e)))
  expect_true(is.finite(e$example_result$p_value))
  if(id!="prop_two_independent_pearson_chisq") {
   expect_error(advisor_r_code(chosen,alternative="greater"),class="inferguidR_invalid_input")
   if("q.direction" %in% names(routes[[id]])) {
    a<-routes[[id]];a["q.direction"]<-"one_sided_prespecified"
    directional<-advisor_override(advisor_replay(a),id)
    expect_error(advisor_r_code(directional,alternative="greater"),class="inferguidR_invalid_input")
   }
  } else {
   a<-routes[[id]];a["q.direction"]<-"one_sided_prespecified"
   directional<-advisor_override(advisor_replay(a),id)
   for(alt in c("less","greater")) {
    code<-advisor_r_code(directional,alternative=alt);env<-new.env()
    invisible(capture.output(eval(parse(text=code),env)))
    expect_identical(env$example_result$alternative,alt)
   }
  }
 }
})
