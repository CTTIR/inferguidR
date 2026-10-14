test_that("report imports validate paired classifier comparisons", {
 s<-classification_examples("scored_synthetic")$spec
 x<-create_report_snapshot(analysis_envelope(s,evaluate_cutoff(s)))
mutations<-list(
 unknown_model=function(z) {z$model_a<-'missing';z},
 same_model=function(z) {z$model_b<-z$model_a;z},
 negative_se=function(z) {z$se<--1;z},
 invalid_p=function(z) {z$p_value<-2;z},
 reversed_ci=function(z) {z$lower<-1;z$upper<--1;z},
 wrong_difference=function(z) {z$difference<-.9;z},
 wrong_auc=function(z) {z$auc_a<-.2;z},
 unknown_status=function(z) {z$inference_status<-'invented';z},
 unknown_method=function(z) {z$method<-'invented';z},
 unknown_field=function(z) {z$raw_scores<-list(.1,.2);z},
 short_column=function(z) {z$model_a<-c('model_a','model_b');z},
 wrong_level=function(z) {z$level<-.5;z},
 wrong_z=function(z) {z$z<-99;z},
 forged_p=function(z) {z$p_value<-.4;z},
 forged_ci=function(z) {z$lower<--.9;z$upper<-.9;z},
 missing_pairs=function(z) NULL,
 duplicate_pair=function(z) lapply(z,function(v) c(v,v)))
 for (name in names(mutations)) {
  y<-unclass(x);z<-mutations[[name]](y$results[[1]]$comparisons)
  y$results[[1]]$comparisons<-z;y$comparisons[1]<-list(z)
  y$snapshot_hash<-NULL;y$snapshot_hash<-.ig_hash_value(y)
  expect_error(ig_from_json(.ig_canonical_json(y)),class="inferguidR_invalid_input",info=name)
 }
})

test_that("comparison import retains degenerate, clustered and multiple-pair results", {
 d<-data.frame(y=rep(c(1,0),each=6),a=c(8,7,7,5,4,2,6,3,3,2,1,1),b=c(5,6,4,4,3,3,5,4,2,2,2,1),cl=rep(1:6,each=2),w=rep(1:2,6))
 for(k in c("ordinary","cluster","weights","same","one_class","single_case","infeasible","three")) {
  dd<-d
  if(k=="same") dd$b<-dd$a
  if(k=="one_class") dd$y<-1
  if(k=="single_case") dd$y<-c(1,rep(0,11))
  if(k=="three") dd$c<-rev(dd$a)
  args<-list(data=dd,truth="y",positive_label=1,scores=if(k=="three") c("a","b","c") else c("a","b"),cutoff=4)
  if(k=="cluster") args$cluster<-"cl"
  if(k=="weights") args$weights<-"w"
  if(k=="infeasible") {args$cutoff<-NULL;args$rule<-cutoff_rule(selection="max_youden",constraints=list(tpr_min=1,fpr_max=0))}
  s<-do.call(classification_spec,args);r<-evaluate_cutoff(s)
  x<-create_report_snapshot(analysis_envelope(s,r),report_kind=if(r$status$code=="ok") "model_cutoff" else "incomplete_settings")
  expect_identical(ig_to_json(ig_from_json(ig_to_json(x))),ig_to_json(x),info=k)
 }
})

test_that("comparison schema rejects untyped and unexpected columns", {
 skip_if_not_installed("jsonvalidate")
 s<-classification_examples("scored_synthetic")$spec
 x<-create_report_snapshot(analysis_envelope(s,evaluate_cutoff(s)))
 expect_true(.ig_validate_schema(ig_to_json(x),"report_document")$valid)
 y<-unclass(x);y$results[[1]]$comparisons$raw_scores<-list(.1,.2)
 expect_false(.ig_validate_schema(.ig_canonical_json(y),"report_document")$valid)
 y<-unclass(x);y$results[[1]]$comparisons$se<-"small"
 expect_false(.ig_validate_schema(.ig_canonical_json(y),"report_document")$valid)
})
