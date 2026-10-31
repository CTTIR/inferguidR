base<-'/Users/rh/Documents/GitHub/CTTIR/public/inferguidR-baseline-20261031/tests/testthat'
calc<-'tools/infer-suite-migration/provider-integration/calc-tests/testthat'
manifest<-read.csv('tools/infer-suite-migration/provider-integration/test-moves.csv',stringsAsFactors=FALSE)
count<-function(path){
 if(!file.exists(path))return(c(blocks=0L,expectations=0L))
 p<-parse(path);test<-Filter(function(e)is.call(e)&&identical(e[[1]],as.name('test_that')),as.list(p))
 c(blocks=length(test),expectations=sum(vapply(test,function(e)sum(grepl('^expect_',all.names(e,unique=FALSE))),integer(1))))
}
rows<-lapply(unique(manifest$source),function(f){
 target<-if(file.exists(file.path(calc,f)))f else sub('^test-','test-kernels-',f)
 old<-count(file.path(base,f));guid<-count(file.path('tests/testthat',f));new<-count(file.path(calc,target))
 data.frame(source=f,baseline_blocks=old[['blocks']],retained_guid_blocks=guid[['blocks']],transferred_calc_blocks=new[['blocks']],baseline_expectation_calls=old[['expectations']],retained_guid_expectation_calls=guid[['expectations']],transferred_calc_expectation_calls=new[['expectations']],expectations_preserved=old[['expectations']]==guid[['expectations']]+new[['expectations']])
})
d<-do.call(rbind,rows);write.csv(d,'tools/infer-suite-migration/provider-integration/assertion-preservation.csv',row.names=FALSE)
stopifnot(all(d$expectations_preserved));cat('PASS:',nrow(d),'historical files preserve every direct expectation call across owners.\n')
