e<-new.env(parent=baseenv());e$alternative<-'two.sided';e$confidence_level<-.95
invisible(capture.output(sys.source(system.file("advisor","templates","tpl.jonckheere_terpstra.v1.R",package="inferguidR"),e)))
perms<-function(x) if(length(x)==1L) matrix(x,nrow=1L) else do.call(rbind,lapply(seq_along(x),function(i) cbind(x[i],perms(x[-i]))))
test_that('integer trend scores agree with independently enumerated distributions',{
 refs<-jsonlite::fromJSON(test_path("fixtures","jonckheere-enumeration.json"),simplifyVector=FALSE)$cases
 for(ref in refs) {
  y<-unlist(ref$values);g<-rep(seq_along(ref$group_sizes),unlist(ref$group_sizes));orders<-perms(seq_along(y))
  actual<-apply(orders,1,function(i) e$jt_twice_statistic(y[i],g));tab<-table(actual)
  expect_equal(as.numeric(tab),as.numeric(unlist(ref$distribution)))
  expect_equal(names(tab),names(ref$distribution))
  expect_equal(e$jt_twice_statistic(y,g),ref$twice_observed)
  expect_equal(sum(tab),ref$labelled_permutations)
 }
})
test_that('permutation tails match a separate rank-sum implementation for the same stream',{
 d<-data.frame(id=1:6,group=rep(c('a','b','c'),c(1,2,3)),y=c(0,0,1,1,1,2));B<-999L
 # Rank-sum identity, rather than pair comparisons used by the template.
 score<-function(y) {
  groups<-split(y,factor(d$group,levels=c('a','b','c')))
  sum(vapply(combn(seq_along(groups),2,simplify=FALSE),function(i) {
   x<-groups[[i[1]]];z<-groups[[i[2]]];sum(rank(c(x,z))[-seq_along(x)])-length(z)*(length(z)+1)/2
  },0.))
 }
 old<-RNGkind();on.exit(do.call(RNGkind,as.list(old)),add=TRUE)
 RNGkind('Mersenne-Twister','Inversion','Rejection');set.seed(8128)
 obs<-score(d$y);s<-replicate(B,score(d$y[sample.int(nrow(d))]));low<-sum(s<=obs);high<-sum(s>=obs)
 for(alt in c('less','greater','two.sided')) {
  r<-e$analyse(d,c('a','b','c'),TRUE,TRUE,8128,alt,B=B)
  expected<-switch(alt,less=(low+1)/(B+1),greater=(high+1)/(B+1),two.sided=min(1,2*min(low+1,high+1)/(B+1)))
  expect_equal(r$statistic,obs);expect_equal(unname(r$simulated_exceedances),c(low,high))
  expect_equal(r$p_value,expected,tolerance=1e-15)
  expect_equal(r$simulated_replicates,B);expect_true(all(r$monte_carlo_interval>=0&r$monte_carlo_interval<=1))
  expect_lte(r$monte_carlo_interval[1],r$monte_carlo_interval[2]);expect_null(r$confidence_interval)
  expect_equal(r$minimum_p_value,if(alt=='two.sided') 2/(B+1) else 1/(B+1))
  opposite<-switch(alt,less='greater',greater='less',two.sided='two.sided')
  rev<-e$analyse(d,c('c','b','a'),TRUE,TRUE,8128,opposite,B=B)
  expect_equal(rev$p_value,r$p_value);expect_equal(rev$statistic+r$statistic,r$cross_group_pairs)
  transformed<-d;transformed$y<-10+3*d$y
  expect_equal(e$analyse(transformed,c('a','b','c'),TRUE,TRUE,8128,alt,B=B)$p_value,r$p_value)
 }
})
test_that('all tied outcomes have p one with no invented uncertainty or simulation',{
 d<-data.frame(id=1:6,group=rep(c('a','b','c'),each=2),y=1)
 for(alt in c('less','greater','two.sided')) {
  r<-e$analyse(d,c('a','b','c'),TRUE,TRUE,42,alt,B=999)
  expect_identical(r$status,'no_rank_variation');expect_equal(r$statistic,6)
  expect_equal(r$p_value,1);expect_equal(r$monte_carlo_interval,c(1,1));expect_equal(r$simulated_replicates,0)
  expect_false(r$reject);expect_equal(r$descriptive_concordance,.5)
 }
})

test_that('Monte Carlo precision bounds invert binomial tails and cover both directions',{
 d<-data.frame(id=1:6,group=rep(c('a','b','c'),each=2),y=c(1,4,2,5,3,6));B<-999L;level<-.99
 lo<-e$analyse(d,c('a','b','c'),TRUE,TRUE,43,'less',B,mc.level=1-(1-level)/2)
 hi<-e$analyse(d,c('a','b','c'),TRUE,TRUE,43,'greater',B,mc.level=1-(1-level)/2)
 both<-e$analyse(d,c('a','b','c'),TRUE,TRUE,43,'two.sided',B,mc.level=level)
 expect_equal(both$monte_carlo_interval,pmin(1,2*pmin(lo$monte_carlo_interval,hi$monte_carlo_interval)))
 for(r in list(lo,hi)) {
  k<-r$simulated_exceedances[[if(r$alternative=='less') 'lower' else 'upper']]
  expect_gt(k,0);expect_lt(k,B)
  expect_equal(stats::pbinom(k-1,B,r$monte_carlo_interval[1],lower.tail=FALSE),(1-r$monte_carlo_confidence_level)/2,tolerance=1e-12)
  expect_equal(stats::pbinom(k,B,r$monte_carlo_interval[2]),(1-r$monte_carlo_confidence_level)/2,tolerance=1e-12)
 }
})
test_that('trend validates design and restores reproducible random state',{
 d<-data.frame(id=1:6,group=rep(c('a','b','c'),each=2),y=1:6)
 run<-function(data=d,order=c('a','b','c'),prespecified=TRUE,independent=TRUE,seed=42,B=999,...) e$analyse(data,order,prespecified,independent,seed,B=B,...)
 oldkind<-RNGkind();had<-exists('.Random.seed',envir=.GlobalEnv,inherits=FALSE);if(had) oldseed<-get('.Random.seed',envir=.GlobalEnv)
 on.exit({do.call(RNGkind,as.list(oldkind));if(had) assign('.Random.seed',oldseed,envir=.GlobalEnv) else if(exists('.Random.seed',envir=.GlobalEnv,inherits=FALSE)) rm('.Random.seed',envir=.GlobalEnv)},add=TRUE)
 RNGkind("L'Ecuyer-CMRG");set.seed(17);saved<-get('.Random.seed',envir=.GlobalEnv);kind<-RNGkind()
 r<-run();expect_identical(get('.Random.seed',envir=.GlobalEnv),saved);expect_identical(RNGkind(),kind)
 expect_identical(run(),r)
 rm('.Random.seed',envir=.GlobalEnv);run();expect_false(exists('.Random.seed',envir=.GlobalEnv,inherits=FALSE));expect_identical(RNGkind(),kind)
 expect_error(run(prespecified=FALSE));expect_error(run(independent=FALSE));expect_error(run(order=c('a','a','c')))
 expect_error(run(order=c('a','b','d')));expect_error(run(seed=NA));expect_error(run(seed=.1));expect_error(run(B=998));expect_error(run(B=1e6))
 expect_error(run(conf.level=1));expect_error(run(mc.level=0));expect_error(run(alternative='increasing'))
 bad<-d;bad$id[2]<-bad$id[1];expect_error(run(bad))
 bad<-d;bad$id[2]<-NA;expect_error(run(bad))
 bad<-d;bad$y[1]<-Inf;expect_error(run(bad,complete_cases=TRUE))
 bad<-d;bad$y<-factor(bad$y);expect_error(run(bad))
 bad<-d;bad$y<-ordered(bad$y);expect_equal(run(bad)$p_value,r$p_value)
 bad<-d;bad$y[1]<-NA;expect_error(run(bad));x<-run(bad,complete_cases=TRUE)
 expect_equal(x$n_used,5);expect_equal(x$n_excluded,1);expect_equal(x$n_total,6)
 bad$y[2]<-NA;expect_error(run(bad,complete_cases=TRUE))
})

test_that('ordered rank trend exports only the selected public route and direction',{
 for(lang in c('en','de')) for(measurement in c('continuous','ordinal')) {
  a<-c(q.aim='compare',q.stage='analysis',q.outcome_type=measurement,q.unit='individuals',
   q.structure='independent',q.independence_check='no',q.n_factors='one',q.n_groups='more_than_two',
   q.covariates='no',q.estimand_continuous='rank',q.missingness='none')
  s<-advisor_override(advisor_replay(a,language=lang),'jonckheere_terpstra')
  rec<-advisor_recommendation(s)
  expect_identical(rec$preferred$analysis_variant_id,'jonckheere_terpstra')
  expect_identical(rec$r_implementation$template_status,'available')
  code<-advisor_r_code(s);env<-new.env(parent=baseenv())
  invisible(capture.output(eval(parse(text=code),env)))
  expect_identical(env$example_result$alternative,'two.sided')
  expect_identical(env$example_result$status,'ok');expect_true(is.finite(env$example_result$p_value))
  expect_match(code,rec$recommendation_hash,fixed=TRUE)
  expect_error(advisor_r_code(s,alternative='greater'),class='inferguidR_invalid_input')
  expect_error(advisor_r_code(s,alternative='less'),class='inferguidR_invalid_input')
  expect_false(rec$planning$planning_compatibility %in% c('exact_gpower_match','independently_validated_extension'))
 }
})
