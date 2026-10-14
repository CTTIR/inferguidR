.structure_env<-function(kind){e<-new.env(parent=baseenv());e$alternative<-'two.sided';e$confidence_level<-.95;invisible(capture.output(sys.source(system.file("advisor","templates",paste0("tpl.",kind,".v1.R"),package="inferguidR"),e)));e}
.structure_matrix<-function(x)do.call(rbind,lapply(x,unlist))
test_that('oblique common-factor models reproduce independent ML covariance fits',{
 ref<-jsonlite::read_json(test_path("fixtures","structure-oracle.json"))$factor;e<-.structure_env('efa_common_factor_ml');r<-e$example_result
 expect_identical(r$status,'ok');expect_equal(unname(r$uniquenesses),unlist(ref$uniqueness),tolerance=2e-5)
 expect_equal(unname(r$reproduced_correlation),.structure_matrix(ref$reproduced),tolerance=2e-5)
 expect_equal(r$objective,ref$objective,tolerance=1e-8);expect_equal(r$goodness_of_fit$statistic,ref$statistic,tolerance=1e-7)
 expect_equal(r$goodness_of_fit$p_value,ref$p_value,tolerance=1e-7);expect_equal(r$goodness_of_fit$df,ref$df)
 expect_equal(unname(diag(r$factor_correlation)),rep(1,2),tolerance=1e-12)
 expect_equal(r$reproduced_correlation,r$loadings%*%r$factor_correlation%*%t(r$loadings)+diag(r$uniquenesses))
 d<-e$example_data;vars<-paste0('v',1:6);run<-function(data=d,...)e$analyse(data,vars,2,TRUE,TRUE,...)
 expect_error(e$analyse(d,vars,2));expect_error(e$analyse(d,vars,2,TRUE,FALSE));expect_error(run(alternative='greater'))
 bad<-d;bad$v2<-bad$v1;expect_identical(run(bad)$status,'not_applicable')
 bad<-d;bad$v1[1]<-NA;expect_error(run(bad));expect_equal(run(bad,complete_cases=TRUE)$n_excluded,1)
 bad<-d;bad$v1<-bad$v1*100+1e6;expect_equal(run(bad)$uniquenesses,r$uniquenesses,tolerance=2e-5)
 expect_error(e$analyse(d,vars,5,TRUE,TRUE));expect_null(r$confidence_interval)
})
test_that('numeric cluster geometry agrees with independent Ward and separated-centre references',{
 ref<-jsonlite::read_json(test_path("fixtures","structure-oracle.json"))$ward
 for(kind in c('cluster_hierarchical_ward','cluster_kmeans')){
  e<-.structure_env(kind);r<-e$example_result;groups<-r$membership$cluster
  expect_identical(r$status,'ok');expect_equal(outer(groups,groups,'=='),.structure_matrix(ref$same_cluster))
  expect_equal(r$membership$silhouette,unlist(ref$silhouette),tolerance=1e-12)
  if(kind=='cluster_hierarchical_ward')expect_equal(r$fit$tree$height,unlist(ref$heights),tolerance=1e-12) else {
   expect_equal(r$fit$within_sum_squares,ref$kmeans_wss,tolerance=1e-12)
   expect_equal(unname(r$fit$centers),.structure_matrix(ref$kmeans_centers),tolerance=1e-12)
   for(g in 1:3)expect_equal(unname(r$fit$centers[g,]),unname(colMeans(r$fit$scaled[groups==g,,drop=FALSE])))
  }
  expect_true(all(r$stability$scores>=-1&r$stability$scores<=1));expect_equal(r$stability$failed,0)
  d<-e$example_data;run<-function(data=d,...)e$analyse(data,c('x','y'),3,TRUE,TRUE,B=20,...)
  set.seed(67);before<-.Random.seed;rr<-run();expect_identical(.Random.seed,before)
  shuffled<-run(d[rev(seq_len(nrow(d))),]);expect_identical(shuffled$membership,rr$membership);expect_equal(shuffled$stability,rr$stability)
  bad<-d;bad$x[1]<-NA;expect_error(run(bad));expect_equal(run(bad,complete_cases=TRUE)$n_excluded,1)
  bad<-d;bad$y<-1;expect_error(run(bad));expect_error(e$analyse(d,c('x','y'),3,TRUE))
  expect_equal(e$adjusted_rand(c(1,1,2,2),c(2,2,1,1)),1)
  expect_equal(e$adjusted_rand(c(1,1,2,2),c(1,2,1,2)),-.5)
 }
})
test_that('typed Gower distances and PAM agree with exhaustive independent medoid enumeration',{
 ref<-jsonlite::read_json(test_path("fixtures","structure-oracle.json"))$gower;e<-.structure_env('cluster_gower_pam');r<-e$example_result
 expect_identical(r$status,'ok');expect_equal(unname(r$dissimilarity),.structure_matrix(ref$distance),tolerance=1e-12)
 expect_equal(r$mean_dissimilarity,ref$global_mean_dissimilarity,tolerance=1e-12)
 expect_equal(outer(r$membership$cluster,r$membership$cluster,'=='),.structure_matrix(ref$same_cluster))
 expect_equal(r$mean_silhouette,ref$mean_silhouette,tolerance=1e-12)
 expect_equal(r$membership$cluster[match(r$medoid_ids,r$membership$id)],1:2)
 d<-e$example_data;vars<-names(e$example_types);weights<-setNames(rep(1,4),vars)
 run<-function(data=d,...)e$analyse(data,vars,e$example_types,weights,2,TRUE,B=20,...)
 set.seed(71);before<-.Random.seed;rr<-run();expect_identical(.Random.seed,before)
 expect_identical(run(d[rev(seq_len(nrow(d))),])$membership,rr$membership)
 bad<-d;bad$amount[1]<-NA;expect_error(run(bad));expect_identical(run(bad,missing='pairwise')$status,'ok')
 bad<-d;bad$grade<-factor(as.character(bad$grade));expect_error(run(bad))
 bad<-d;bad$flag[1]<-2;expect_error(run(bad))
 expect_error(e$analyse(d,vars,e$example_types,weights,2));expect_null(r$confidence_interval)
})
test_that('Gower pairwise weights and binary conventions match explicit independent arithmetic',{
 ref<-jsonlite::read_json(test_path("fixtures","structure-gower-edge-oracle.json"));e<-.structure_env('cluster_gower_pam')
 d<-data.frame(amount=c(0,2,4,6,NA,3),kind=factor(c('a','a','b','b','c',NA)),grade=ordered(c('low','middle','high','middle','high','low'),levels=c('low','middle','high')),
  symmetric=c(0,1,0,1,1,0),asymmetric=c(0,0,1,1,0,1))
 types<-c(amount='continuous',kind='nominal',grade='ordinal',symmetric='symmetric_binary',asymmetric='asymmetric_binary')
 r<-e$gower_fit(d,2,types,unlist(ref$weights));expect_equal(unname(r$distance),.structure_matrix(ref$distance),tolerance=1e-12)
 expect_equal(e$gower_fit(d,2,types,unlist(ref$weights)*100)$distance,r$distance,tolerance=1e-12)
})

test_that('structure advisor exports exact exploratory models and accepted alternatives',{
 routes<-list(
  efa_common_factor_ml=c(q.aim='structure',q.stage='analysis',q.structure_goal='variables',q.reduction_purpose='latent'),
  cluster_hierarchical_ward=c(q.aim='structure',q.stage='analysis',q.structure_goal='cases',q.cluster_types='all_continuous'),
  cluster_gower_pam=c(q.aim='structure',q.stage='analysis',q.structure_goal='cases',q.cluster_types='mixed'))
 for(lang in c('en','de'))for(kind in names(routes)){
  s<-advisor_replay(routes[[kind]],language=lang);r<-advisor_recommendation(s)
  expect_identical(r$preferred$analysis_variant_id,kind);expect_identical(r$r_implementation$template_status,'available')
  code<-advisor_r_code(s);expect_match(code,r$recommendation_hash,fixed=TRUE)
  e<-new.env(parent=baseenv());invisible(capture.output(eval(parse(text=code),e)));expect_identical(e$example_result$status,'ok')
  expect_error(advisor_r_code(s,alternative='greater'),'hypothesis alternative');expect_error(advisor_r_code(s,alternative='two.sided'),'hypothesis alternative')
  if(kind=='cluster_hierarchical_ward'){
   ss<-advisor_override(s,'cluster_kmeans');rr<-advisor_recommendation(ss)
   expect_identical(rr$preferred$analysis_variant_id,'cluster_kmeans');expect_identical(rr$r_implementation$template_status,'available')
   e<-new.env(parent=baseenv());invisible(capture.output(eval(parse(text=advisor_r_code(ss)),e)));expect_identical(e$example_result$status,'ok')
  }
 }
})
