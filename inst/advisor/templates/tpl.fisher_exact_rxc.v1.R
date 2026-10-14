# Synthetic teaching example: Fisher's fixed-margin probability-ordering test.
# Independent units, mutually exclusive categories; both margins conditioned on.
# The alternative is omnibus. No odds-ratio interval is claimed for a general table.
# This teaching script limits total count to 100 and each dimension to 5.
# Larger exact problems need an explicit resource plan, not automatic approximation.
analyse <- function(counts, alternative="two.sided", conf.level=.95) {
 stopifnot(is.numeric(counts),is.matrix(counts),all(dim(counts)>=2L),all(dim(counts)<=5L),
   all(is.finite(counts)),all(counts>=0),all(counts==floor(counts)),sum(counts)<=100,
   all(rowSums(counts)>0),all(colSums(counts)>0),identical(alternative,"two.sided"),
   is.numeric(conf.level),length(conf.level)==1L,is.finite(conf.level),conf.level>0,conf.level<1)
 warnings<-character()
 fit<-tryCatch(withCallingHandlers(stats::fisher.test(counts,alternative="two.sided",
   conf.int=FALSE,simulate.p.value=FALSE,hybrid=FALSE,workspace=200000L),
   warning=function(w) {warnings<<-c(warnings,conditionMessage(w));invokeRestart("muffleWarning")}),error=function(e)e)
 if(inherits(fit,"error")) return(list(method="Fisher fixed-margin exact test",status="unavailable",
   reason=conditionMessage(fit),observed=counts,p_value=NULL,
   note="Exact calculation failed within the fixed workspace; no simulation or asymptotic substitution was used."))
 list(method=fit$method,status="ok",observed=counts,n=sum(counts),
   row_totals=rowSums(counts),column_totals=colSums(counts),
   null="Row and column independence conditional on both margins",
   ordering="Sum of probabilities of tables no more probable than the observed table",
   p_value=fit$p.value,confidence_interval=NULL,interval_status="not_computed",
   requested_confidence_level=conf.level,alpha=1-conf.level,reject=fit$p.value<1-conf.level,
   warnings=warnings,assumptions="Independent units; prespecified mutually exclusive categories; fixed-margin conditional inference.")
}
example_counts<-rbind(Group1=c(A=3,B=1,C=4),Group2=c(A=1,B=5,C=2))
example_result<-analyse(example_counts,alternative,confidence_level)
print(example_result)
