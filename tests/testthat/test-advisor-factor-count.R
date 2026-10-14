.pa_env<-function(){e<-new.env(parent=baseenv());invisible(capture.output(sys.source(system.file('advisor','templates','tpl.efa_common_factor_ml.v1.R',package='inferguidR'),e)));e}
test_that('SMC parallel analysis matches independent regression communalities',{
 e<-.pa_env();ref<-jsonlite::read_json(test_path('fixtures','factor-count-oracle.json'));x<-do.call(rbind,lapply(ref$x,unlist));d<-data.frame(id=seq_len(nrow(x)),x);names(d)<-c('id',paste0('v',1:6));v<-names(d)[-1]
 permutations<-array(unlist(ref$permutations),dim=c(nrow(x),6,100));r<-e$parallel_analysis(d,v,TRUE,TRUE,B=100,permutations=permutations)
 expect_identical(r$status,'ok');expect_equal(r$observed_eigenvalues,unlist(ref$observed),tolerance=1e-9)
 expect_equal(r$null_eigenvalues,do.call(rbind,lapply(ref$null,unlist)),tolerance=1e-9)
 expect_equal(r$null_thresholds,unlist(ref$thresholds),tolerance=1e-9);expect_equal(r$suggested_leading_factors,ref$leading);expect_null(r$seed)
 set.seed(729);before<-.Random.seed;kind<-RNGkind();z<-e$parallel_analysis(d,v,TRUE,TRUE,B=100)
 expect_identical(.Random.seed,before);expect_identical(RNGkind(),kind);expect_identical(z,e$parallel_analysis(d,v,TRUE,TRUE,B=100))
 expect_error(e$parallel_analysis(d,v));expect_error(e$parallel_analysis(d,v,TRUE,TRUE,B=100,quantile_probability=.99))
 dd<-d;dd$v2<-dd$v1;expect_identical(e$parallel_analysis(dd,v,TRUE,TRUE,B=100)$status,'unavailable')
 dd<-d;dd$v2<-0;expect_identical(e$parallel_analysis(dd,v,TRUE,TRUE,B=100)$status,'unavailable')
 dd<-d;dd$v1[1]<-NA;expect_error(e$parallel_analysis(dd,v,TRUE,TRUE,B=100));expect_equal(e$parallel_analysis(dd,v,TRUE,TRUE,B=100,complete_cases=TRUE)$n_excluded,1)
 bad<-permutations;bad[1,1,1]<-bad[2,1,1];expect_error(e$parallel_analysis(d,v,TRUE,TRUE,B=100,permutations=bad))
 dd<-d;dd$id[2]<-dd$id[1];expect_error(e$parallel_analysis(dd,v,TRUE,TRUE,B=100))
})
test_that('failed null permutations are not discarded',{
 e<-.pa_env();n<-24L;p<-3L;B<-100L;set.seed(2007);d<-data.frame(id=seq_len(n),v1=sample.int(n),v2=sample.int(n),v3=sample.int(n));a<-array(rep(seq_len(n),p*B),c(n,p,B))
 for(j in 1:p)a[,j,1]<-order(d[[j+1]])
 r<-e$parallel_analysis(d,c('v1','v2','v3'),TRUE,TRUE,B=B,permutations=a)
 expect_identical(r$status,'simulation_incomplete');expect_equal(r$failed_replicates,1);expect_null(r$suggested_leading_factors)
})
