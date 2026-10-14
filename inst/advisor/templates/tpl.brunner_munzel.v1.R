# Independent ordinal/continuous samples: theta = P(X > Y) + 0.5 P(X = Y).
# This tests stochastic equality, not equality of medians or identical distributions.
# A t approximation is used explicitly; no automatic exact/permutation or normal fallback.
analyse <- function(x,y,alternative="two.sided",conf.level=.95) {
 stopifnot(is.numeric(x),is.numeric(y),all(is.finite(x)),all(is.finite(y)),
   is.character(alternative),length(alternative)==1L,alternative %in% c("two.sided","greater","less"),
   is.numeric(conf.level),length(conf.level)==1L,is.finite(conf.level),conf.level>0,conf.level<1)
 nx<-length(x);ny<-length(y)
 if(nx<10L||ny<10L) stop("This t-approximation teaching example requires at least 10 observations per group; smaller samples need a separately justified permutation method.")
 pooled<-rank(c(x,y),ties.method="average")
 # Cross-group placements from pooled minus within-group midranks.
 px<-(pooled[seq_len(nx)]-rank(x,ties.method="average"))/ny
 py<-1-(pooled[nx+seq_len(ny)]-rank(y,ties.method="average"))/nx
 theta<-mean(px);vx<-stats::var(px)/nx;vy<-stats::var(py)/ny;variance<-vx+vy
 if(!is.finite(variance)||variance<=0) return(list(method="Brunner-Munzel t approximation",status="unavailable",
   reason="Zero estimated placement variance: no t inference or automatic replacement is available.",
   estimate=theta,n_x=nx,n_y=ny,p_value=NULL,confidence_interval=NULL))
 se<-sqrt(variance);df<-variance^2/(vx^2/(nx-1)+vy^2/(ny-1));statistic<-(theta-.5)/se
 p<-switch(alternative,two.sided=2*stats::pt(abs(statistic),df,lower.tail=FALSE),greater=stats::pt(statistic,df,lower.tail=FALSE),less=stats::pt(statistic,df))
 critical<-stats::qt(if(alternative=="two.sided") (1-conf.level)/2 else 1-conf.level,df,lower.tail=FALSE)
 raw<-switch(alternative,two.sided=c(theta-critical*se,theta+critical*se),greater=c(theta-critical*se,Inf),less=c(-Inf,theta+critical*se))
 list(method="Brunner-Munzel t approximation",status="ok",estimate=theta,null=.5,
   estimand="P(X > Y) + 0.5 P(X = Y); greater means stochastic superiority of x",
   n_x=nx,n_y=ny,std_error=se,statistic=statistic,df=df,p_value=p,alternative=alternative,
   confidence_level=conf.level,confidence_interval=pmax(0,pmin(1,raw)),unconstrained_interval=raw,
   interval_method="Placement-variance t approximation, intersected with [0,1]; unconstrained limits retained.",
   alpha=1-conf.level,reject=p<1-conf.level,
   assumptions="Independent observations and groups; ordinal or continuous outcome; unequal shapes/spreads allowed. Small-sample accuracy is not guaranteed. No mean/median shift interpretation, automatic method substitution or data-driven direction.")
}
example_x<-(seq_len(14)*7+seq_len(14)^2)%%19/3
example_y<-(seq_len(11)*11+2*seq_len(11)^2)%%23/3
example_result<-analyse(example_x,example_y,alternative=alternative,conf.level=confidence_level)
print(example_result)
