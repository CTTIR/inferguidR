test_that("variance sensitivity reaches roots beyond the generic effect scan", {
 for(tails in c("one","two")) {
  s<-analysis_spec("chi2_variance_one_sample","sensitivity",list(n_total=2,alpha=.001,power=.99),list(tails=tails))
  r<-calculate_power(s)
  expect_identical(r$status$code,"ok")
  # Independently checked with 60-digit incomplete-gamma probabilities.
  expected<-if(tails=="one") 68926.8178632028 else 77098.9106102785
  expect_equal(r$solved$ratio,expected,tolerance=1e-9)
  expect_equal(r$achieved$power,.99,tolerance=1e-10)
  expect_lt(r$diagnostics$other_branch_roots,1)
  expect_true(r$diagnostics$target_reached)
 }
})

test_that("variance sensitivity accepts the null-effect boundary consistently", {
 for(n in c(2,20,100000)) for(tails in c("one","two")) {
  r<-calculate_power(analysis_spec("chi2_variance_one_sample","sensitivity",list(n_total=n,alpha=.05,power=.05),list(tails=tails)))
  expect_identical(r$status$code,"ok")
  expect_equal(r$solved$ratio,1,tolerance=1e-9)
  expect_equal(r$achieved$power,.05,tolerance=1e-10)
 }
})

test_that("variance sensitivity preserves both branches at ordinary and large degrees of freedom", {
 for(n in c(3,50,100000)) for(tails in c("one","two")) {
  r<-calculate_power(analysis_spec("chi2_variance_one_sample","sensitivity",list(n_total=n,alpha=.001,power=.8),list(tails=tails)))
  expect_identical(r$status$code,"ok")
  expect_gt(r$solved$ratio,1);expect_lt(r$diagnostics$other_branch_roots,1)
  for(ratio in c(r$solved$ratio,r$diagnostics$other_branch_roots)) {
   f<-calculate_power(analysis_spec("chi2_variance_one_sample","post_hoc",list(n_total=n,alpha=.001,ratio=ratio),list(tails=tails)))
   expect_equal(f$achieved$power,.8,tolerance=1e-9)
  }
 }
})
