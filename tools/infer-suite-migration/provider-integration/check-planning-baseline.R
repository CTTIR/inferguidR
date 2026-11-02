args<-commandArgs(TRUE);stopifnot(length(args)==1L)
library(inferguidR)
baseline<-readRDS(args[1]) # Trusted synthetic local baseline artifact.
now<-list()
for(v in c('less','greater','minlike','central'))now[[paste0('boschloo_power_',v)]]<-boschloo_power(5,7,.2,.7,variant=v)
now$boschloo_design<-boschloo_sample_size(.1,.8,variant='less',power=.8,n_max=30)
now$t_forward<-calculate_power(analysis_spec('t_means_two_independent','post_hoc',list(d=.5,n1=24,n2=32,alpha=.05),list(tails='two')))
now$t_inverse<-calculate_power(analysis_spec('t_means_two_independent','a_priori',list(d=.5,power=.8,alpha=.05),list(tails='two')))
normalize<-function(x){
 if(is.character(x)){x<-gsub('inferguidR::','infercalcR::',x,fixed=TRUE);x<-gsub('inferguidR Boschloo engine','infercalcR Boschloo engine',x,fixed=TRUE);return(x)}
 if(is.list(x)){a<-attributes(x);x<-lapply(x,normalize);attributes(x)<-a};x
}
for(n in names(now)){
 equal<-all.equal(normalize(now[[n]]),normalize(baseline[[n]]),tolerance=1e-13,check.attributes=TRUE)
 if(!isTRUE(equal))stop(n,': ',paste(equal,collapse='; '))
 cat('PASS complete planning baseline:',n,'\n')
}
