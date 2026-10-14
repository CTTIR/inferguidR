.metadata_reports <- function(locale="en") {
 s<-analysis_spec('t_means_one_sample','post_hoc',list(d=.5,alpha=.05,n_total=30),list(tails='two'))
 b<-boschloo_test(1,2,0,2,'greater')
 cs<-classification_examples('scored_synthetic')$spec
 cc<-classification_examples('counts_fixed')$spec
 ct<-classification_examples('binormal_equal_variance')$spec
 env<-list(power=analysis_envelope(s,calculate_power(s)),observed=analysis_envelope(b$spec,b),
   scored=analysis_envelope(cs,evaluate_cutoff(cs)),counts=analysis_envelope(cc,evaluate_cutoff(cc)),
   theory=analysis_envelope(ct,evaluate_cutoff(ct)),power_missing=analysis_envelope(s),
   observed_missing=analysis_envelope(b$spec),scored_missing=analysis_envelope(cs))
 lapply(seq_along(env),function(i) create_report_snapshot(env[[i]],locale=locale,report_kind=if(i>5) 'incomplete_settings' else 'single_analysis'))
}
.metadata_rehash<-function(x) {
 x<-unclass(x);x$snapshot_hash<-NULL;x$snapshot_hash<-.ig_hash_value(x);x
}
test_that('report metadata cannot claim unsupported content or change retained evidence',{
 mutations<-list(
  function(y){y$sweeps<-list(list(result=999));y},
  function(y){y$resources$status<-'verified';y},
  function(y){y$resources$extra<-'unsupported';y},
  function(y){y$validation[[1]]$status<-'certified';y},
  function(y){y$validation[[1]]$analysis<-2;y},
  function(y){y$validation[[1]]$evidence<-list(status='certified');y},
  function(y){y$uncertainty[[1]]$status<-'computed';y},
  function(y){y$uncertainty[[1]]$analysis<-2;y},
  function(y){y$uncertainty[[1]]$extra<-.95;y},
  function(y){y$reproducibility[[1]]$analysis<-2;y},
  function(y){y$reproducibility[[1]]$private_data_required<-!y$reproducibility[[1]]$private_data_required;y},
  function(y){y$reproducibility[[1]]$extra<-'unsupported';y})
 for(x in .metadata_reports()) {
  expect_identical(ig_to_json(x),ig_to_json(ig_from_json(ig_to_json(x))))
  for(m in mutations) {
   y<-.metadata_rehash(m(unclass(x)))
   expect_error(.ig_report_document_from_list(y),class='inferguidR_invalid_input')
   expect_error(ig_from_json(.ig_canonical_json(y)),class='inferguidR_invalid_input')
  }
 }
})
test_that('power and observed reproduction code matches the frozen result or specification',{
 reports<-.metadata_reports()
 for(i in c(1,2,6,7)) {
  y<-unclass(reports[[i]]);y$reproducibility[[1]]$r_code<-'unrelated_code()'
  expect_error(ig_from_json(.ig_canonical_json(.metadata_rehash(y))),class='inferguidR_invalid_input')
 }
 for(i in c(3,4,5,8)) {
  y<-unclass(reports[[i]]);y$reproducibility[[1]]$r_code<-list(code='not literal text')
  expect_error(ig_from_json(.ig_canonical_json(.metadata_rehash(y))),class='inferguidR_invalid_input')
 }
})
test_that('metadata validation preserves bilingual reports and exported text',{
 for(locale in c('en','de')) for(x in .metadata_reports(locale)) {
  y<-ig_from_json(ig_to_json(x))
  expect_identical(.ig_report_text(x),.ig_report_text(y))
 }
})

test_that("report schema rejects unsupported reserved sections and malformed reproduction text", {
 skip_if_not_installed("jsonvalidate")
 s<-classification_examples("counts_fixed")$spec
 x<-create_report_snapshot(analysis_envelope(s,evaluate_cutoff(s)))
 expect_true(.ig_validate_schema(ig_to_json(x),"report_document")$valid)
 mutations<-list(function(y){y$sweeps<-list(list(value=1));y},
   function(y){y$resources$status<-"verified";y},
   function(y){y$reproducibility[[1]]$r_code<-1;y},
   function(y){y$reproducibility[[1]]$unexpected<-TRUE;y})
 for(m in mutations) expect_false(.ig_validate_schema(.ig_canonical_json(m(unclass(x))),"report_document")$valid)
})
