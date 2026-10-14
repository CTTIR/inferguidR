.paired_report<-function(k=2L,identical_models=FALSE,locale='en') {
 d<-classification_examples('scored_synthetic')$data
 if(identical_models) d$model_b<-d$model_a
 if(k>2L) for(i in 3:k) d[[paste0('model_',i)]]<-pmin(1,pmax(0,d$model_a*(i/(k+1))+d$model_b*(1-i/(k+1))))
 scores<-c('model_a','model_b',if(k>2L) paste0('model_',3:k))
 s<-classification_spec(d,truth='outcome',positive_label='event',negative_label='no_event',scores=scores,id='id',score_scale='probability',cutoff=.5)
 create_report_snapshot(analysis_envelope(s,evaluate_cutoff(s)),report_kind='comparison',locale=locale)
}
.comparison_figures<-function(x) Filter(function(v) v$plot_type=='classifier_comparison',x$visualisations)
test_that('paired comparison figures retain exact orientation, intervals and shared scales', {
 x<-.paired_report(5);figs<-.comparison_figures(x)
 expect_length(figs,2L);expect_length(figs[[1]]$data$pairs,8L);expect_length(figs[[2]]$data$pairs,2L)
 expect_identical(figs[[1]]$data$x_limits,figs[[2]]$data$x_limits)
 pairs<-unlist(lapply(figs,function(v) v$data$pairs),recursive=FALSE);r<-x$results[[1]]$comparisons
 get<-function(v,j) if(is.list(v)) v[[j]] else v[j]
 for(j in seq_along(pairs)) {
  row<-pairs[[j]];expect_equal(row$pair_index,j)
  expect_identical(row$model_a,get(r$model_a,j));expect_identical(row$model_b,get(r$model_b,j))
  for(name in c('difference','lower','upper','level','inference_status')) expect_identical(row[[name]],get(r[[name]],j))
 }
 expect_match(figs[[1]]$caption,'no multiplicity adjustment',fixed=TRUE)
 expect_identical(ig_to_json(ig_from_json(ig_to_json(x))),ig_to_json(x))
})
test_that('paired figures preserve unavailable intervals without inventing precision', {
 x<-.paired_report(2,TRUE,'de');v<-.comparison_figures(x)[[1]];row<-v$data$pairs[[1]]
 expect_equal(row$difference,0);expect_null(row$lower);expect_null(row$upper)
 expect_false(row$inference_status=='ok');expect_true(nzchar(row$reason_code))
 expect_match(v$axes$x,'A minus B',fixed=TRUE)
 expect_match(v$caption,'Nicht verfuegbare Intervalle',fixed=TRUE)
})
test_that('rehashed comparison figure mutations are refused', {
 x<-.paired_report(5);idx<-which(vapply(x$visualisations,function(v) v$plot_type=='classifier_comparison',FALSE))
 reseal<-function(z) {z<-unclass(z);z$snapshot_hash<-NULL;z$snapshot_hash<-.ig_hash_value(z);.ig_canonical_json(z)}
 mutations<-list(
  function(z) {z$visualisations[[idx[1]]]$data$pairs[[1]]$difference<-.99;z},
  function(z) {z$visualisations[[idx[1]]]$data$pairs[[1]]$model_a<-'unknown';z},
  function(z) {z$visualisations[[idx[1]]]$data$x_limits<-list(-.01,.01);z},
  function(z) {z$visualisations[[idx[1]]]$data$panel_index<-99;z},
  function(z) {z$visualisations[[idx[1]]]$caption<-'Unsupported claim';z},
  function(z) {z$visualisations[[idx[1]]]$data$pairs<-rev(z$visualisations[[idx[1]]]$data$pairs);z},
  function(z) {z$visualisations<-z$visualisations[-idx[2]];z},
  function(z) {z$visualisations[[idx[2]]]<-z$visualisations[[idx[1]]];z})
 for(f in mutations) expect_error(ig_from_json(reseal(f(x))),class='inferguidR_invalid_input')
})
test_that('comparison scope still refuses a lone model or lone power analysis', {
 s<-classification_examples('counts_fixed')$spec
 expect_error(create_report_snapshot(analysis_envelope(s,evaluate_cutoff(s)),report_kind='comparison'),class='inferguidR_invalid_input')
 power<-analysis_spec('t_means_one_sample','a_priori',list(d=.5,alpha=.05,power=.8),list(tails='two'))
 expect_error(create_report_snapshot(analysis_envelope(power,calculate_power(power)),report_kind='comparison'),class='inferguidR_invalid_input')
})
