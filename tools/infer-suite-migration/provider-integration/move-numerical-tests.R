# Preserve numerical regression tests under their implementation owner's namespace.
base <- '/Users/rh/Documents/GitHub/CTTIR/public/inferguidR-baseline-20261031/tests/testthat'
out <- 'tools/infer-suite-migration/provider-integration/calc-tests/testthat'
dir.create(file.path(out,'fixtures'),recursive=TRUE,showWarnings=FALSE)
manifest<-list()
whole<-c('calc-fixtures','calc-functions','calc-parser','calc-trig-domain','calc-trig','classif-independent','classif-metrics','classif-theory')
for(n in whole){
 f<-paste0('test-',n,'.R');x<-readLines(file.path(base,f));
 writeLines(gsub('inferguidR::','infercalcR::',x,fixed=TRUE),file.path(out,f))
 manifest[[f]]<-data.frame(source=f,test='ALL',action='move entire numerical test file to calc')
}
# Only whitebox blocks move from mixed files. Public integration blocks stay in guid.
missing<-c('bl_all','.bl_pvalues_all','.bl_ordering','.bl_exact_order','.bl_level_tables','.bl_sup','.bl_tail_c','.bl_ambiguity_touches','.bl_set_sup','.ig_cls_wilson','.ig_cls_clopper_pearson','.ig_cls_delong_interval','.ig_cls_with_seed','.ig_cls_bootstrap_draw')
files<-c('boschloo-enumeration','boschloo-properties','boschloo-plotdata','classif-confidence-boundary','classif-uncertainty','num-dist')
for(n in files){
 f<-paste0('test-',n,'.R');txt<-readLines(file.path(base,f));p<-parse(text=txt,keep.source=TRUE);sr<-attr(p,'srcref');take<-integer();helpers<-integer();tests<-character()
 for(i in seq_along(p)){
  e<-p[[i]];is_test<-is.call(e)&&identical(e[[1]],as.name('test_that'))
  symbols<-all.names(e)
  if(is_test){
   picked<-any(missing%in%symbols)
   if(n=='num-dist')picked<-!any(c('.ig_critical_symmetric','.ig_critical_values','.ig_two_sided_prob','.ig_power_t_test')%in%symbols)
   if(picked){take<-c(take,seq.int(sr[[i]][1],sr[[i]][3]));tests<-c(tests,as.character(e[[2]]))}
  }else if(is.call(e)&&as.character(e[[1]])%in%c('<-','='))helpers<-c(helpers,seq.int(sr[[i]][1],sr[[i]][3]))
 }
 if(!length(take))next
 picked<-txt[sort(unique(c(helpers,take)))];picked<-gsub('inferguidR::','infercalcR::',picked,fixed=TRUE)
 writeLines(c('# Numerical regression blocks retained from the historical guidance suite.',picked),file.path(out,paste0('test-kernels-',n,'.R')))
 # Remove just the transferred blocks; retain original comments and test helpers.
 if(n=='boschloo-enumeration'){ for(i in seq_along(p)){e<-p[[i]];if(is.call(e)&&identical(e[[1]],as.name('<-'))&&identical(e[[2]],as.name('bl_all')))take<-union(take,seq.int(sr[[i]][1],sr[[i]][3]))}}
 writeLines(txt[-take],file.path('tests/testthat',f))
 manifest[[f]]<-data.frame(source=f,test=tests,action='move numerical whitebox test block to calc')
}
# Preserve exact helper implementations and independent input fixtures.
h<-readLines(file.path(base,'helper-reference.R'));i<-grep('^ref_reset_registry <-',h);h<-h[seq_len(i-1L)]
writeLines(h,file.path(out,'helper-legacy-reference.R'))
writeLines('th_condition_code <- function(expr) tryCatch({force(expr); NA_character_}, inferguidR_error = function(e) e$status$code)',file.path(out,'helper-legacy-condition.R'))
file.copy(file.path(base,'fixtures',c('calculator-trig-domain-independent.csv','classifier-confidence-boundary.csv')),file.path(out,'fixtures'),overwrite=TRUE)
write.csv(do.call(rbind,manifest),'tools/infer-suite-migration/provider-integration/test-moves.csv',row.names=FALSE)
