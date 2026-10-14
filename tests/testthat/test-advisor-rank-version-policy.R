.rank_policy_env <- function(id) {
 e<-new.env();e$alternative<-"two.sided";e$confidence_level<-.95
 path<-system.file("advisor","templates",paste0("tpl.",id,".v1.R"),package="inferguidR")
 stopifnot(nzchar(path));invisible(capture.output(sys.source(path,e)));e
}
.rank_policy_tail <- function(lo,hi,alternative) switch(alternative,
 less=lo,greater=hi,two.sided=min(1,2*min(lo,hi)))

test_that("signed rank inference explicitly omits null differences across R versions", {
 for(id in c("signed_rank_wilcoxon_one_sample","signed_rank_wilcoxon")) {
  e<-.rank_policy_env(id)
  d<-c(-1,1,-2,3,4,0,0,0);a<-d[d!=0];ranks<-rank(abs(a));n<-length(a)
  statistic<-sum(ranks[a>0]);t<-as.numeric(table(abs(a)))
  se<-sqrt(n*(n+1)*(2*n+1)/24-sum(t^3-t)/48);delta<-statistic-n*(n+1)/4
  for(alt in c("less","greater","two.sided")) {
   correction<-switch(alt,less=-.5,greater=.5,two.sided=sign(delta)*.5)
   z<-(delta-correction)/se
   expected<-switch(alt,less=pnorm(z),greater=pnorm(z,lower.tail=FALSE),two.sided=2*pnorm(-abs(z)))
   for(mu in c(0,10)) {
    r<-e$signed_analysis(d+mu,mu,alt,.8,FALSE)
    expect_equal(r$statistic,statistic,tolerance=1e-12)
    expect_equal(r$p_value,expected,tolerance=1e-12)
    expect_equal(r$n_zero,3);expect_equal(r$n_total,8)
   }
  }
 }
})

test_that("signed rank examples retain distinct binary64 values instead of rounding ties", {
 d<-c(1,-(1+2^-30),2,-3,4,-5);n<-length(d);ranks<-rank(abs(d))
 statistic<-sum(ranks[d>0]);support<-drop(as.matrix(expand.grid(rep(list(c(0,1)),n)))%*%ranks)
 for(id in c("signed_rank_wilcoxon_one_sample","signed_rank_wilcoxon")) {
  e<-.rank_policy_env(id)
  for(alt in c("less","greater","two.sided")) for(exact in c(TRUE,FALSE)) {
   if(exact) expected<-.rank_policy_tail(mean(support<=statistic),mean(support>=statistic),alt) else {
    delta<-statistic-n*(n+1)/4;se<-sqrt(n*(n+1)*(2*n+1)/24)
    correction<-switch(alt,less=-.5,greater=.5,two.sided=sign(delta)*.5)
    z<-(delta-correction)/se
    expected<-switch(alt,less=pnorm(z),greater=pnorm(z,lower.tail=FALSE),two.sided=2*pnorm(-abs(z)))
   }
   r<-e$signed_analysis(d,0,alt,.8,exact)
   expect_equal(r$statistic,statistic,tolerance=1e-12)
   expect_equal(r$p_value,expected,tolerance=1e-12)
   expect_false(r$tied_absolute_ranks)
  }
 }
})

test_that("rank sum inference and descriptive superiority use the same unrounded ranks", {
 e<-.rank_policy_env("rank_sum_mann_whitney")
 x<-c(1,3,5);y<-c(1+2^-30,2,4);nx<-length(x);ny<-length(y);n<-nx+ny
 u<-sum(outer(x,y,">"));support<-colSums(combn(seq_len(n),nx))-nx*(nx+1)/2
 for(alt in c("less","greater","two.sided")) for(exact in c(TRUE,FALSE)) {
  if(exact) expected<-.rank_policy_tail(mean(support<=u),mean(support>=u),alt) else {
   delta<-u-nx*ny/2;se<-sqrt(nx*ny*(n+1)/12)
   correction<-switch(alt,less=-.5,greater=.5,two.sided=sign(delta)*.5)
   z<-(delta-correction)/se
   expected<-switch(alt,less=pnorm(z),greater=pnorm(z,lower.tail=FALSE),two.sided=2*pnorm(-abs(z)))
  }
  r<-e$analyse(x,y,alt,exact=exact)
  expect_equal(r$statistic,u,tolerance=1e-12)
  expect_equal(r$p_value,expected,tolerance=1e-12)
  expect_equal(r$probability_superiority,u/(nx*ny),tolerance=1e-12)
  expect_false(r$tied_ranks)
 }
})
