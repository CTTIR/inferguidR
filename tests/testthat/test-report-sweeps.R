ns<-asNamespace('inferguidR');rehash<-function(x){x<-unclass(x);x$snapshot_hash<-NULL;x$snapshot_hash<-get('.ig_hash_value',ns)(x);get('.ig_canonical_json',ns)(x)}
testthat::test_that('sweep reports round-trip and preserve strict versioned content',{
 spec<-analysis_spec('t_means_two_independent','a_priori',inputs=list(alpha=.05,power=.8,d=.5),options=list(tails='two'));e<-analysis_envelope(spec,calculate_power(spec));sw<-sweep_design(spec,'alpha',c(.01,0,.05))
 for(lang in c('en','de')){
  x<-create_report_snapshot(e,locale=lang,report_kind='incomplete_settings',sweeps=list(list(analysis=1,sweep=sw)))
  testthat::expect_identical(x$completeness,'incomplete');testthat::expect_error(create_report_snapshot(e,sweeps=list(list(analysis=1,sweep=sw))),'incomplete_settings')
  testthat::expect_identical(x$schema_version,'inferguidR.report_document/0.2.0');j<-ig_to_json(x)
  testthat::expect_identical(ig_to_json(ig_from_json(j)),j);testthat::expect_identical(x$sweeps[[1]]$rows[[2]]$status,'invalid_input')
  figures<-Filter(function(v)v$plot_type=='power_sweep',x$visualisations);testthat::expect_length(figures,1);testthat::expect_null(figures[[1]]$data[[2]]$y)
  bad<-unclass(x);bad$sweeps[[1]]$rows[[1]]$value<-.02;testthat::expect_error(ig_from_json(rehash(bad)),class='inferguidR_invalid_input')
  bad<-unclass(x);bad$visualisations<-Filter(function(v)v$plot_type!='power_sweep',bad$visualisations);testthat::expect_error(ig_from_json(rehash(bad)),class='inferguidR_invalid_input')
  bad<-unclass(x);bad$schema_version<-'inferguidR.report_document/0.1.0';testthat::expect_error(ig_from_json(rehash(bad)),class='inferguidR_invalid_input')
  for(fmt in c('json','txt','pdf')){f<-tempfile(fileext=paste0('.',fmt));export_report(x,f,fmt,overwrite=TRUE);testthat::expect_gt(file.size(f),0)}
 }
 old<-create_report_snapshot(e);testthat::expect_identical(old$schema_version,'inferguidR.report_document/0.1.0');testthat::expect_identical(ig_to_json(ig_from_json(ig_to_json(old))),ig_to_json(old))
})

test_that("German figure captions are localised without changing stored metadata", {
 v<-list(caption="Displayed densities or probability masses describe the fixed design.")
 expect_match(.ig_report_caption(v,"de"),"Dichten")
 expect_identical(.ig_report_caption(v,"en"),v$caption)
 expect_identical(v$caption,"Displayed densities or probability masses describe the fixed design.")
})

test_that("legacy unavailable distribution stays unavailable after renderer improvements", {
 j<-paste(readLines(test_path('fixtures','report-legacy-nct-unavailable.json')),collapse='\n')
 x<-ig_from_json(j)
 expect_identical(x$visualisations[[1]]$status,'unavailable')
 expect_length(x$visualisations[[1]]$data,0)
 expect_identical(as.character(ig_to_json(x)),j)
 bad<-unclass(x);bad$visualisations[[1]]$fixed_parameters$d<-.8
 expect_error(ig_from_json(rehash(bad)),class='inferguidR_invalid_input')
 bad<-unclass(x);bad$visualisations[[1]]$data<-list(list(x=0,y=1))
 expect_error(ig_from_json(rehash(bad)),class='inferguidR_invalid_input')
})
