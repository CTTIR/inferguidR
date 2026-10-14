test_that("rehashed classifier reports reject impossible and inconsistent aggregates", {
 s<-classification_spec(counts=c(tp=40,fn=10,fp=20,tn=130),cutoff=.5,prevalence=.1)
 x<-create_report_snapshot(analysis_envelope(s,evaluate_cutoff(s)),report_kind="model_cutoff")
mutations<-list(
 negative_cell=function(y) {y$results[[1]]$models[[1]]$operating_point$tp<--1;y},
 wrong_total=function(y) {y$results[[1]]$models[[1]]$operating_point$n_total<-1;y},
 unknown_model_field=function(y) {y$results[[1]]$models[[1]]$raw_scores<-c(.123,.456);y},
 bad_comparator=function(y) {y$results[[1]]$models[[1]]$comparator<-'invented';y},
 negative_ppv=function(y) {j<-which(unlist(y$results[[1]]$models[[1]]$metrics$metric_id)=='ppv');y$results[[1]]$models[[1]]$metrics$value[[j]]<--.1;y$metrics[[1]]$models[[1]]$metrics<-y$results[[1]]$models[[1]]$metrics;y},
 short_column=function(y) {y$results[[1]]$models[[1]]$metrics$value<-1;y$metrics[[1]]$models[[1]]$metrics<-y$results[[1]]$models[[1]]$metrics;y},
 negative_scenario=function(y) {y$results[[1]]$models[[1]]$scenario$counts$tp<--1;y$metrics[[1]]$models[[1]]$scenario<-y$results[[1]]$models[[1]]$scenario;y},
 forged_figure_model=function(y) {y$visualisations[[1]]$model_id<-'unknown';y},
 negative_figure_cell=function(y) {y$visualisations[[1]]$data$tp<--1;y},
 invalid_interval=function(y) {y$results[[1]]$models[[1]]$intervals$lower[[1]]<-.9;y$results[[1]]$models[[1]]$intervals$upper[[1]]<-.1;y},
 unknown_model_status=function(y) {y$results[[1]]$models[[1]]$status$code<-'invented';y},
 input_kind_mismatch=function(y) {y$results[[1]]$input_kind<-'scored_data';y})
 mutations$wrong_valid_ppv<-function(y) {j<-which(unlist(y$results[[1]]$models[[1]]$metrics$metric_id)=="ppv");y$results[[1]]$models[[1]]$metrics$value[[j]]<-.9;y$metrics[[1]]$models[[1]]$metrics<-y$results[[1]]$models[[1]]$metrics;y}
 mutations$unknown_metric_column<-function(y) {y$results[[1]]$models[[1]]$metrics$raw_scores<-1;y$metrics[[1]]$models[[1]]$metrics<-y$results[[1]]$models[[1]]$metrics;y}
 mutations$duplicate_metric<-function(y) {y$results[[1]]$models[[1]]$metrics$metric_id[[2]]<-"tp";y$metrics[[1]]$models[[1]]$metrics<-y$results[[1]]$models[[1]]$metrics;y}
 mutations$altered_uncertainty<-function(y) {y$uncertainty[[1]]$models[[1]]$intervals$level[[1]]<-.5;y}
 mutations$injected_figure_data<-function(y) {y$visualisations[[1]]$raw_scores<-c(.1,.2);y}
 for(mutate in mutations) {
  y<-mutate(unclass(x));y$snapshot_hash<-NULL;y$snapshot_hash<-.ig_hash_value(y)
  expect_error(.ig_report_document_from_list(y),class="inferguidR_invalid_input")
  expect_error(ig_from_json(.ig_canonical_json(y)),class="inferguidR_invalid_input")
 }
})

test_that("classifier import preserves undefined scenarios, weighted counts and clustered intervals", {
 for(counts in list(c(tp=0,fn=0,fp=0,tn=10),c(tp=10,fn=0,fp=0,tn=0),c(tp=0,fn=3,fp=7,tn=0))) for(pi in c(0,.1,1)) {
  s<-classification_spec(counts=counts,cutoff=.5,prevalence=pi)
  x<-create_report_snapshot(analysis_envelope(s,evaluate_cutoff(s)))
  expect_identical(ig_to_json(x),ig_to_json(ig_from_json(ig_to_json(x))))
 }
 d<-data.frame(y=rep(c("yes","no"),10),score=seq(.05,1,length.out=20),w=rep(1:2,10),cluster=rep(1:10,each=2))
 specs<-list(
  classification_spec(d,truth="y",positive_label="yes",scores="score",cutoff=.5,weights="w",prevalence=.1),
  classification_spec(d,truth="y",positive_label="yes",scores="score",cutoff=.5,cluster="cluster"),
  classification_spec(d,truth="y",positive_label="yes",scores="score",rule=cutoff_rule(selection="max_youden",constraints=list(tpr_min=1,fpr_max=0))))
 for(s in specs) {
  r<-evaluate_cutoff(s);x<-create_report_snapshot(analysis_envelope(s,r),report_kind=if(r$status$code=="ok") "model_cutoff" else "incomplete_settings")
  expect_identical(ig_to_json(x),ig_to_json(ig_from_json(ig_to_json(x))))
 }
})

test_that("classifier report schema rejects unknown model and metric-table fields", {
 skip_if_not_installed("jsonvalidate")
 s<-classification_examples("counts_fixed")$spec
 x<-create_report_snapshot(analysis_envelope(s,evaluate_cutoff(s)))
 expect_true(.ig_validate_schema(ig_to_json(x),"report_document")$valid)
 y<-unclass(x);y$results[[1]]$models[[1]]$raw_scores<-list(.1,.2)
 expect_false(.ig_validate_schema(.ig_canonical_json(y),"report_document")$valid)
 y<-unclass(x);y$results[[1]]$models[[1]]$metrics$unexpected<-1
 expect_false(.ig_validate_schema(.ig_canonical_json(y),"report_document")$valid)
 y<-unclass(x);y$results[[1]]$models[[1]]$metrics$value[[1]]<-"forty"
 expect_false(.ig_validate_schema(.ig_canonical_json(y),"report_document")$valid)
})
