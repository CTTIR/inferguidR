testthat::test_that('display noncentral t coordinates agree with independent quantiles',{
 refs<-jsonlite::read_json(testthat::test_path('fixtures','report-nct-oracle.json'))$cases
 for(ref in refs){law<-list(family='t',df=ref$df,ncp=ref$ncp);q<-NULL
  testthat::expect_warning(q<-inferguidR:::`.ig_plot_law`(law,unlist(ref$probabilities),'q'),NA)
  testthat::expect_equal(q,unlist(ref$quantiles),tolerance=2e-8)
  density<-NULL;testthat::expect_warning(density<-inferguidR:::`.ig_plot_law`(law,q,'d'),NA)
  testthat::expect_equal(density,unlist(ref$densities),tolerance=2e-8)
  testthat::expect_equal(vapply(q,function(x) inferguidR:::`.ig_calc_nct_pdf`(x,ref$df,ref$ncp),0),unlist(ref$densities),tolerance=2e-8)
 }
 s<-analysis_spec('t_means_two_independent','a_priori',inputs=list(alpha=.05,power=.8,d=.5),options=list(tails='two'));pd<-power_plot_data(calculate_power(s))
 testthat::expect_identical(pd$status,'ok');testthat::expect_true(all(is.finite(pd$data$x)));testthat::expect_true(all(is.finite(pd$data$y)));testthat::expect_gt(nrow(pd$data),800)
})
