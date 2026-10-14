.operating_report <- function(counts=c(tp=40,fn=10,fp=20,tn=130),locale="en",scenario=.1) {
 s<-classification_spec(counts=counts,cutoff=.5,prevalence=scenario)
 create_report_snapshot(analysis_envelope(s,evaluate_cutoff(s)),locale=locale)
}
.operating_figures <- function(x) Filter(function(v) v$plot_type %in% c("classifier_roc_point","classifier_pr_point"),x$visualisations)
test_that("selected ROC and PR points retain aggregate and scenario meaning", {
 for(lang in c("en","de")) {
  x<-.operating_report(locale=lang);v<-.operating_figures(x)
  expect_length(v,2L)
  expect_equal(unlist(v[[1]]$data$points[[1]][c("x","y")]),c(x=2/15,y=.8),tolerance=1e-14)
  expect_equal(unlist(v[[2]]$data$points[[1]][c("x","y")]),c(x=.8,y=2/3),tolerance=1e-14)
  expect_equal(unlist(v[[2]]$data$points[[2]][c("x","y")]),c(x=.8,y=.4),tolerance=1e-14)
  expect_identical(v[[2]]$data$points[[2]]$basis,"scenario")
  expect_false(x$privacy$raw_rows_included);expect_false(x$privacy$score_curve_tables_included)
  expect_identical(ig_to_json(ig_from_json(ig_to_json(x))),ig_to_json(x))
 }
})
test_that("undefined denominators do not become zero ROC or PR coordinates", {
 for(counts in list(c(tp=0,fn=10,fp=0,tn=20),c(tp=0,fn=0,fp=3,tn=7),c(tp=3,fn=7,fp=0,tn=0))) {
  v<-.operating_figures(.operating_report(counts,scenario=NULL))
  if(sum(counts[c("tp","fn")])==0) {expect_null(v[[1]]$data$points[[1]]$y);expect_null(v[[2]]$data$points[[1]]$x)}
  if(sum(counts[c("fp","tn")])==0) expect_null(v[[1]]$data$points[[1]]$x)
  if(sum(counts[c("tp","fp")])==0) expect_null(v[[2]]$data$points[[1]]$y)
  for(fig in v) if(fig$status=="not_applicable") expect_true(nzchar(fig$reason))
 }
})
test_that("rehashing cannot forge selected point coordinates, basis or labels", {
 x<-.operating_report();idx<-which(vapply(x$visualisations,function(v) v$plot_type=="classifier_pr_point",FALSE))[1]
 mutations<-list(function(y) {y$visualisations[[idx]]$data$points[[1]]$y<-.9;y},
   function(y) {y$visualisations[[idx]]$data$points[[2]]$basis<-"observed";y},
   function(y) {y$visualisations[[idx]]$axes$y<-"AUC";y},
   function(y) {y$visualisations[[idx]]$data$points[[1]]$status<-"undefined";y},
   function(y) {y$visualisations[[idx]]$model_id<-"unknown";y},
   function(y) {y$visualisations[[idx]]$caption<-"Full curve";y},
   function(y) {y$visualisations[[idx]]$data$x_limits<-list(-1,2);y},
   function(y) {y$visualisations<-y$visualisations[-idx];y},
   function(y) {y$visualisations<-c(y$visualisations,list(y$visualisations[[idx]]));y})
 for(f in mutations) {y<-f(unclass(x));y$snapshot_hash<-NULL;y$snapshot_hash<-.ig_hash_value(y);expect_error(ig_from_json(.ig_canonical_json(y)),class="inferguidR_invalid_input")}
})
