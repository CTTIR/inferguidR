.split_fixture<-function(part,ids=TRUE,cluster=FALSE) {
 d<-classification_examples("scored_synthetic",part=part,models="model_a")$data
 if(cluster) d$cluster<-"shared_cluster"
 classification_spec(d,truth="outcome",positive_label="event",scores="model_a",id=if(ids) "id" else NULL,
  cluster=if(cluster) "cluster" else NULL,cutoff=.3,
  evaluation_split=if(part=="tuning") "held_out_tuning" else "held_out_test")
}
test_that("cutoff transfer rejects known overlap even when fingerprints differ", {
 tune<-.split_fixture("tuning");chosen<-select_cutoff(tune,"max_youden")
 all<-.split_fixture("all")
 expect_false(identical(tune$data_fingerprint,all$data_fingerprint))
 expect_error(evaluate_cutoff(all,cutoff=chosen),"share unit identifiers")
 d<-classification_examples("scored_synthetic",part="tuning",models="model_a")$data
 reverse<-classification_spec(d[nrow(d):1,],truth="outcome",positive_label="event",scores="model_a",id="id",cutoff=.3,evaluation_split="held_out_test")
 expect_error(evaluate_cutoff(reverse,cutoff=chosen),"share unit identifiers")
 ctune<-.split_fixture("tuning",cluster=TRUE);ctest<-.split_fixture("test",cluster=TRUE)
 expect_error(evaluate_cutoff(ctest,cutoff=select_cutoff(ctune,"max_youden")),"share cluster identifiers")
})
test_that("disjoint provided IDs are recorded without claiming proven independence", {
 tune<-.split_fixture("tuning");test<-.split_fixture("test");chosen<-select_cutoff(tune,"max_youden")
 r<-evaluate_cutoff(test,cutoff=chosen)
 expect_identical(r$evaluation$transfer$overlap_check$unit_identity,"no_overlap_detected")
 expect_identical(r$evaluation$transfer$overlap_check$cluster_identity,"not_available")
 expect_match(r$evaluation$statement,"declarations do not establish independence")
 expect_true("split_independence_unverified" %in% r$limitations$code)
 expect_identical(r$evaluation$label,"held_out_rule")
 # The added identity check does not reselect the operating point on test outcomes.
 direct<-evaluate_cutoff(update(test,rule=cutoff_rule(cutoff=chosen$models[[1]]$operating_point$cutoff,selection_data="held_out_tuning")))
 expect_identical(r$models,direct$models)
 expect_match(direct$evaluation$statement,"declarations do not establish independence")
 for(lang in c("en","de")) {
  report<-create_report_snapshot(analysis_envelope(classification_evaluated_spec(r),r),locale=lang,report_kind="model_cutoff")
  json<-ig_to_json(report)
  expect_identical(json,ig_to_json(ig_from_json(json)))
  expect_false(any(vapply(attr(tune,"ig_cls_data")$id,function(id) grepl(id,json,fixed=TRUE),logical(1))))
 }
})
test_that("unavailable private identities remain explicitly unchecked", {
 for(ids in c(FALSE,TRUE)) {
  tune<-.split_fixture("tuning",ids=ids);test<-.split_fixture("test",ids=FALSE)
  r<-evaluate_cutoff(test,cutoff=select_cutoff(tune,"max_youden"))
  expect_identical(r$evaluation$transfer$overlap_check$unit_identity,"not_available")
  expect_true("split_independence_unverified" %in% r$limitations$code)
 }
 chosen<-select_cutoff(.split_fixture("tuning"),"max_youden")
 imported<-ig_from_json(ig_to_json(analysis_envelope(classification_evaluated_spec(chosen),chosen)))$result;r<-evaluate_cutoff(.split_fixture("test"),cutoff=imported)
 expect_identical(r$evaluation$transfer$overlap_check$unit_identity,"not_available")
 expect_match(r$evaluation$transfer$overlap_check$scope,"does not prove independence")
 expect_true("split_independence_unverified" %in% r$limitations$code)
})
test_that("mutated attached identities cannot bypass the overlap check", {
 tune<-.split_fixture("tuning");test<-.split_fixture("test");chosen<-select_cutoff(tune,"max_youden")
 altered<-test;prep<-attr(altered,"ig_cls_data");prep$id[1]<-"changed_identity";attr(altered,"ig_cls_data")<-prep
 expect_error(evaluate_cutoff(altered,cutoff=chosen),"evaluation data no longer match")
 source_spec<-classification_evaluated_spec(chosen);prep<-attr(source_spec,"ig_cls_data");prep$id[1]<-"changed_identity";attr(source_spec,"ig_cls_data")<-prep
 attr(chosen,"ig_cls_evaluated_spec")<-source_spec
 expect_error(evaluate_cutoff(test,cutoff=chosen),"tuning data no longer match")
 attr(chosen,"ig_cls_evaluated_spec")<-test
 expect_error(evaluate_cutoff(test,cutoff=chosen),"tuning specification does not match")
})
