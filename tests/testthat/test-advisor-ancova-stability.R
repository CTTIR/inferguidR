.stable_ancova_template<-function() {
 e<-new.env();e$alternative<-'two.sided';e$confidence_level<-.95;e$example_group_count<-3L
 invisible(capture.output(sys.source(system.file('advisor','templates','tpl.ancova_baseline_adjusted.v1.R',package='inferguidR'),e)));e
}
test_that('weak group effects retain their analytically known ANCOVA sum of squares',{
 e<-.stable_ancova_template()
 # Balanced covariates and zero-mean residuals orthogonal to both group and baseline.
 # RSS=48, residual df=21, and the adjusted group contrast variance is 8/21.
 x<-rep(c(-3,-2,-1,1,2,3),4);noise<-rep(c(1,-2,1,1,-2,1),4);g<-rep(c('a','b'),each=12)
 for(effect in 2^c(-6,-12,-18,-24)) {
  d<-data.frame(id=seq_along(x),group=g,baseline=x,y=2+.5*x+effect*(g=='b')+noise)
  expected<-effect^2*21/8
  r<-e$analyse(d,'baseline',c('a','b'),TRUE)
  expect_gt(r$group_test$statistic,0)
  expect_equal(r$group_test$statistic/expected,1,tolerance=1e-7)
  expect_equal(r$group_test$partial_r_squared/(expected/(expected+21)),1,tolerance=1e-7)
  expect_equal(r$group_test$statistic/r$differences$statistic^2,1,tolerance=1e-13)
  expect_equal(r$group_test$p_value,r$differences$p_value,tolerance=1e-13)
  reordered<-e$analyse(d[nrow(d):1,],'baseline',c('b','a'),TRUE)
  expect_equal(reordered$group_test$statistic/expected,1,tolerance=1e-7)
  d$y<-8*d$y+16
  transformed<-e$analyse(d,'baseline',c('a','b'),TRUE)
  expect_equal(transformed$group_test$statistic/expected,1,tolerance=1e-7)
 }
})
