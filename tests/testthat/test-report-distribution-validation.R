.distribution_snapshot <- function(method="t_means_one_sample", inputs=list(d=.5,alpha=.05,n_total=30), options=list(tails="two")) {
 s<-analysis_spec(method,"post_hoc",inputs,options)
 create_report_snapshot(analysis_envelope(s,calculate_power(s)))
}
.distribution_rehash <- function(x) {
 x<-unclass(x);x$snapshot_hash<-NULL;x$snapshot_hash<-.ig_hash_value(x);x
}
test_that("power report import checks the retained distribution rather than just its hash", {
 x<-.distribution_snapshot()
 mutations<-list(
  function(y) {y$visualisations[[1]]$data[[1]]$y<-999;y},
  function(y) {y$visualisations[[1]]$data[[1]]$y<--1e-15;y},
  function(y) {y$visualisations[[1]]$data[[1]]$x<-999;y},
  function(y) {y$visualisations[[1]]$data[[1]]$hypothesis<-"forged";y},
  function(y) {y$visualisations[[1]]$data[[1]]$reject<-!y$visualisations[[1]]$data[[1]]$reject;y},
  function(y) {y$visualisations[[1]]$data<-y$visualisations[[1]]$data[-1];y},
  function(y) {y$visualisations[[1]]$axes$y<-"probability_mass";y},
  function(y) {y$visualisations[[1]]$caption<-"forged";y},
  function(y) {y$visualisations[[1]]$fixed_parameters$d<-.9;y},
  function(y) {y$visualisations[[1]]$critical<-list();y},
  function(y) {y$visualisations[[1]]$status<-"unavailable";y},
  function(y) {y$visualisations[[1]]$data[[1]]$extra<-1;y},
  function(y) {y$visualisations<-list();y},
  function(y) {y$visualisations<-rep(y$visualisations,2);y}
 )
 for(mutate in mutations) {
  y<-.distribution_rehash(mutate(unclass(x)))
  expect_error(.ig_report_document_from_list(y),class="inferguidR_invalid_input")
  expect_error(ig_from_json(.ig_canonical_json(y)),class="inferguidR_invalid_input")
 }
})
test_that("distribution report validation retains continuous and discrete valid displays", {
 cases<-list(
  list("f_variances_two_groups",list(ratio=1.5,n1=20,n2=30,alpha=.05),list(tails="two")),
  list("chi2_variance_one_sample",list(ratio=1.5,n_total=30,alpha=.05),list(tails="two")),
  list("t_means_one_sample",list(d=.5,alpha=.05,n_total=30),list(tails="two")),
  list("z_generic",list(mu=-2,sigma=1.5,alpha=.05),list(tails="one")),
  list("f_generic",list(df1=3,df2=40,lambda=5,alpha=.05),list()),
  list("chi2_generic",list(df=5,lambda=10,alpha=.05),list()),
  list("exact_proportion_one_sample_binomial",list(g=.2,p0=.5,n_total=20,alpha=.05),list(tails="two"))
 )
 for(args in cases) {
  x<-do.call(.distribution_snapshot,args)
  expect_identical(x$results[[1]]$status$code,"ok")
  expect_identical(x$visualisations[[1]]$status,"ok")
  expect_identical(ig_to_json(x),ig_to_json(ig_from_json(ig_to_json(x))))
 }
 x<-.distribution_snapshot()
 y<-unclass(x);y$visualisations[[1]]$data[[1]]$y<-y$visualisations[[1]]$data[[1]]$y+1e-14
 expect_s3_class(ig_from_json(.ig_canonical_json(.distribution_rehash(y))),"report_document")
})

test_that("unavailable distribution and settings-only reports remain explicit", {
 x<-.distribution_snapshot("exact_proportion_one_sample_binomial",list(g=.1,p0=.5,n_total=2100,alpha=.05),list(tails="two"))
 expect_identical(x$results[[1]]$status$code,"ok")
 expect_identical(x$visualisations[[1]]$status,"unavailable")
 expect_identical(ig_to_json(x),ig_to_json(ig_from_json(ig_to_json(x))))
 s<-analysis_spec("t_means_one_sample","post_hoc",list(d=.5,alpha=.05,n_total=30),list(tails="two"))
 x<-create_report_snapshot(analysis_envelope(s),report_kind="incomplete_settings")
 expect_length(x$visualisations,0L)
 expect_identical(ig_to_json(x),ig_to_json(ig_from_json(ig_to_json(x))))
})
