testthat::test_that('sweep rows cannot change settings, display, identity or admissibility', {
 s<-analysis_spec('t_means_one_sample','post_hoc',list(d=.5,alpha=.05,n_total=30),list(tails='two'))
 e<-analysis_envelope(s,calculate_power(s));sw<-sweep_design(s,'n_total',c(1,10,20))
 report<-function(z)create_report_snapshot(e,report_kind='incomplete_settings',sweeps=list(list(analysis=1,sweep=z)))
 for(name in c('achieved_power','achieved_beta','alpha_attained','n_total','solved_power')) {
  bad<-sw;bad[[name]][2]<-123;testthat::expect_error(report(bad),class='inferguidR_invalid_input')
 }
 bad<-sw;bad$admissible[1]<-TRUE;testthat::expect_error(report(bad),class='inferguidR_invalid_input')
 bad<-sw;bad$reason[1]<-'edited';testthat::expect_error(report(bad),class='inferguidR_invalid_input')
 bad<-sw;bad$value[1]<-Inf;testthat::expect_error(report(bad),class='inferguidR_invalid_input')
 bad<-sw;attr(bad,'base_spec_hash')<-'stale';testthat::expect_error(report(bad),class='inferguidR_invalid_input')
 bad<-sw;bad$result[[2]]<-sw$result[[3]];testthat::expect_error(report(bad),class='inferguidR_invalid_input')
 bad<-sw;bad$spec[[2]]<-sw$spec[[3]];testthat::expect_error(report(bad),class='inferguidR_invalid_input')
 fail<-sweep_design(s,'n_total',c(0,1));x<-report(fail)
 testthat::expect_identical(x$completeness,'incomplete')
 testthat::expect_identical(ig_to_json(ig_from_json(ig_to_json(x))),ig_to_json(x))
 testthat::expect_length(Filter(function(v)v$plot_type=='power_sweep',x$visualisations),0)
 sw<-sweep_design(s,'d',c(.2,.5));good<-create_report_snapshot(e,sweeps=list(list(analysis=1,sweep=sw)))
 testthat::expect_identical(good$completeness,'complete')
 for(scope in c('session','comparison')){
  x<-create_report_snapshot(list(e,e),report_kind=scope,sweeps=list(list(analysis=2,sweep=sw)))
  testthat::expect_identical(x$sweeps[[1]]$analysis,2L)
  testthat::expect_identical(ig_to_json(ig_from_json(ig_to_json(x))),ig_to_json(x))
 }
})
