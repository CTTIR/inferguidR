# Synthetic summaries selected by measurement meaning, not by R storage type.
# No hypothesis test or p-value is reported.
# Continuous mean intervals assume independent observations with normal errors.
# Category intervals are marginal Clopper-Pearson intervals, not simultaneous.
analyse <- function(x,measurement,conf.level=.95,event=NULL) {
 stopifnot(is.character(measurement),length(measurement)==1L,!is.na(measurement),
  measurement%in%c("continuous","binary","ordinal","nominal"),
  is.null(dim(x)),length(x)>0L,!anyNA(x),
  is.numeric(conf.level),length(conf.level)==1L,is.finite(conf.level),conf.level>0,conf.level<1)
 n<-length(x)
 if(measurement=="continuous") {
  stopifnot(is.numeric(x),all(is.finite(x)),is.null(event))
  mu<-mean(x);s<-if(n>1L) stats::sd(x) else NA_real_
  ci<-NULL;reason<-"At least two observations with resolved positive sample variance are required."
  if(n>1L && is.finite(s) && s>0 && is.finite(mu)) {
   h<-stats::qt((1+conf.level)/2,n-1)*s/sqrt(n)
   if(is.finite(h) && all(is.finite(mu+c(-h,h)))) {ci<-mu+c(-h,h);reason<-NULL}
  }
  return(list(measurement=measurement,n=n,mean=mu,sd=s,
   quantiles=stats::quantile(x,c(.25,.5,.75),type=7,names=FALSE),quantile_type=7L,
   confidence_interval=ci,interval_estimand="population mean",confidence_level=conf.level,
   interval_status=if(is.null(ci)) "not_applicable" else "ok",interval_reason=reason,p_value=NULL))
 }
 stopifnot(is.factor(x)||is.character(x)||is.logical(x),all(nzchar(as.character(x))))
 if(measurement=="binary") {
  stopifnot(is.factor(x)||is.logical(x),length(event)==1L,!is.na(event))
  declared<-if(is.logical(x)) c("FALSE","TRUE") else levels(x)
  stopifnot(length(declared)==2L,as.character(event)%in%declared)
  hits<-sum(as.character(x)==as.character(event))
  ci<-stats::binom.test(hits,n,conf.level=conf.level)$conf.int
  return(list(measurement=measurement,n=n,event=as.character(event),events=hits,proportion=hits/n,
   confidence_interval=as.numeric(ci),interval_estimand="population event probability",
   confidence_level=conf.level,interval_status="ok",p_value=NULL))
 }
 stopifnot(is.null(event))
 if(measurement=="ordinal") stopifnot(is.ordered(x))
 f<-if(is.factor(x)) x else factor(x);counts<-table(f);levels<-names(counts)
 ci<-t(vapply(as.integer(counts),function(k) as.numeric(stats::binom.test(k,n,conf.level=conf.level)$conf.int),numeric(2)))
 table<-data.frame(category=levels,count=as.integer(counts),proportion=as.numeric(counts)/n,lower=ci[,1],upper=ci[,2],stringsAsFactors=FALSE)
 quantiles<-if(measurement=="ordinal") levels[stats::quantile(as.integer(f),c(.25,.5,.75),type=1,names=FALSE)] else NULL
 list(measurement=measurement,n=n,categories=table,ordered_quantile_categories=quantiles,
  quantile_type=if(measurement=="ordinal") 1L else NULL,
  interval_estimand="each population category probability",confidence_level=conf.level,
  interval_scope="marginal, not simultaneous",p_value=NULL)
}
# With no recorded outcome type, the standalone example explicitly uses continuous data.
if(!exists("example_measurement",inherits=FALSE)) example_measurement<-"continuous"
example_x<-switch(example_measurement,
 continuous=c(1,2,4,6,7,9),
 binary=factor(c("yes","no","yes","no","yes"),levels=c("no","yes")),
 ordinal=ordered(c("low","low","mid","high"),levels=c("low","mid","high")),
 nominal=c("A","A","B","C"))
example_result<-analyse(example_x,measurement=example_measurement,conf.level=confidence_level,
 event=if(example_measurement=="binary") "yes" else NULL)
print(example_result)
