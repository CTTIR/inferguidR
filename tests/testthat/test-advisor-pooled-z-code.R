e<-new.env(parent=baseenv());e$alternative<-'two.sided';e$confidence_level<-.95
invisible(capture.output(sys.source(system.file("advisor","templates","tpl.unconditional_exact_z_pooled.v1.R",package="inferguidR"),e)))
test_that('pooled score probabilities enclose independent rational polynomial maxima',{
 refs<-jsonlite::read_json(test_path("fixtures","pooled-z-polynomial.json"))
 for(x in refs) {
  tab<-matrix(c(x$x1,x$n1-x$x1,x$x2,x$n2-x$x2),2,byrow=TRUE)
  z<-e$analyse(tab,x$alternative,independent_binomial_groups=TRUE);truth<-as.numeric(x$p_value)
  expect_identical(z$status,'ok');expect_lte(z$p_lower,truth+1e-15)
  expect_gte(z$p_upper,truth-1e-15);expect_lte(z$p_upper-z$p_lower,1e-8)
  expect_equal(z$p_value,z$p_upper);expect_null(z$confidence_interval)
  pooled<-(x$x1+x$x2)/(x$n1+x$n2)
  expected_z<-if(pooled %in% c(0,1)) 0 else
   (x$x1/x$n1-x$x2/x$n2)/sqrt(pooled*(1-pooled)*(1/x$n1+1/x$n2))
  expect_equal(z$statistic,expected_z,tolerance=1e-14)
  opposite<-switch(x$alternative,greater='less',less='greater',two.sided='two.sided')
  swapped<-e$analyse(tab[2:1,,drop=FALSE],opposite,independent_binomial_groups=TRUE)
  complemented<-e$analyse(tab[,2:1,drop=FALSE],opposite,independent_binomial_groups=TRUE)
  expect_equal(swapped$p_value,z$p_value,tolerance=1e-8)
  expect_equal(complemented$p_value,z$p_value,tolerance=1e-8)
 }
})
test_that('pooled score boundaries include nuisance endpoints and resource refusals',{
 run<-function(tab=matrix(c(4,1,1,4),2,byrow=TRUE),...)e$analyse(tab,independent_binomial_groups=TRUE,...)
 for(alt in c('two.sided','greater','less')) {
  for(tab in list(matrix(c(0,5,0,5),2,byrow=TRUE),matrix(c(5,0,5,0),2,byrow=TRUE))) {
   z<-run(tab,alternative=alt);expect_equal(z$p_value,1);expect_identical(z$decision,'do_not_reject')
  }
 }
 fail<-run(matrix(c(0,5,1,4),2,byrow=TRUE),tolerance=1e-10,max_subdivisions=1)
 expect_identical(fail$status,'numerical_limit');expect_true(is.na(fail$p_value))
 expect_identical(fail$decision,'numerical_limit');expect_gt(fail$p_upper,fail$p_lower)
 z<-run();alpha<-(z$p_lower+z$p_upper)/2
 expect_identical(run(conf.level=1-alpha)$decision,'numerically_indeterminate')
 expect_error(e$analyse(matrix(1,2,2)));expect_error(run(matrix(0,2,2)))
 expect_error(run(matrix(8,2,2)));expect_error(run(matrix(c(1,1,0,0),2,byrow=TRUE)))
 expect_error(run(matrix(c(1,1,NA,1),2)));expect_error(run(matrix(c(1,1,Inf,1),2)))
 expect_error(run(matrix(.5,2,2)));expect_error(run(matrix(-1,2,2)))
 expect_error(run(alternative='central'));expect_error(run(conf.level=0))
 expect_error(run(conf.level=1));expect_error(run(tolerance=1e-11))
 expect_error(run(max_subdivisions=4097));expect_error(run(max_subdivisions=1.5))
 expect_error(run(matrix(1,2,3)))
})
test_that('small-N binomial multiplicities are exact integers on this platform',{
 ref<-read.csv(test_path("fixtures","pooled-z-binomial-integers.csv"))
 expect_identical(round(choose(ref$n,ref$k)),as.numeric(ref$value))
})

test_that('pooled-Z overrides preserve the exact event and direction contract in EN and DE',{
 base<-c(q.aim='compare',q.stage='analysis',q.outcome_type='binary',q.unit='individuals',
  q.structure='independent',q.independence_check='no',q.n_factors='one',q.n_groups='two',
  q.covariates='no',q.binary_margins='group_sizes_fixed',q.direction='two_sided',q.missingness='none')
 for(lang in c('en','de')) for(direction in c('two_sided','one_sided_prespecified')) {
  ans<-base;ans['q.direction']<-direction;s<-advisor_replay(ans,language=lang)
  s<-advisor_override(s,'unconditional_exact_z_pooled');r<-advisor_recommendation(s)
  expect_identical(r$preferred$analysis_variant_id,'unconditional_exact_z_pooled')
  expect_identical(r$r_implementation$template_status,'available')
  for(alt in if(direction=='two_sided') 'two.sided' else c('less','greater')) {
   code<-advisor_r_code(s,alternative=alt);expect_match(code,r$recommendation_hash,fixed=TRUE)
   env<-new.env(parent=baseenv());invisible(capture.output(eval(parse(text=code),env)))
   expect_identical(env$example_result$method,'Unconditional exact pooled-Z test')
   expect_identical(env$example_result$alternative,alt);expect_identical(env$example_result$status,'ok')
  }
  expect_error(advisor_r_code(s,alternative=if(direction=='two_sided') 'greater' else 'two.sided'))
  if(direction=='one_sided_prespecified')expect_error(advisor_r_code(s))
 }
})
