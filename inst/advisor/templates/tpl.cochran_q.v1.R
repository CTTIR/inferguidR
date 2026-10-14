# Synthetic complete-block teaching example, with one wide row per confirmed participant.
# Columns are prespecified occasions; participant IDs establish matching explicitly.
# Missing occasions require an explicit complete-case choice; excluded blocks are counted.
check_blocks <- function(data, missing, alternative, conf.level) {
 stopifnot(is.data.frame(data), "id" %in% names(data), nrow(data)>=2L,
   !anyNA(data$id),!anyDuplicated(data$id),!anyDuplicated(names(data)),
   identical(alternative,"two.sided"),is.numeric(conf.level),length(conf.level)==1L,
   is.finite(conf.level),conf.level>0,conf.level<1,
   length(missing)==1L,missing %in% c("error","complete_cases"))
 values<-data[,setdiff(names(data),"id"),drop=FALSE]
 stopifnot(ncol(values)>=3L,all(vapply(values,is.numeric,logical(1))))
 x<-as.matrix(values);stopifnot(!any(is.infinite(x)))
 keep<-stats::complete.cases(x)
 if(any(!keep) && missing=="error") stop("Missing occasions: choose and justify complete_cases explicitly.")
 stopifnot(sum(keep)>=2L)
 list(x=x[keep,,drop=FALSE],n_input=nrow(x),n_excluded=sum(!keep))
}
analyse <- function(data, alternative="two.sided", conf.level=.95, missing="error") {
 b<-check_blocks(data,missing,alternative,conf.level)
 raw<-as.matrix(data[,setdiff(names(data),"id"),drop=FALSE]);stopifnot(all(is.na(raw) | raw %in% c(0,1)))
 x<-b$x;k<-ncol(x)
 # Binary coding is explicit: 1 is the event and 0 the non-event.
 stopifnot(all(x %in% c(0,1)))
 totals<-colSums(x);rows<-rowSums(x);denominator<-sum(rows*(k-rows))
 if(denominator==0)
   return(list(status="unavailable",reason="No discordant blocks: conditional variance is zero.",
     n_input=b$n_input,n_complete=nrow(x),n_excluded=b$n_excluded,p_value=NULL))
 statistic<-(k-1)*k*sum((totals-mean(totals))^2)/denominator;df<-k-1
 p<-stats::pchisq(statistic,df,lower.tail=FALSE)
 descriptive<-list(event=1,event_counts=totals,event_proportions=colMeans(x))
 list(method="cochran_q",status="ok",n_input=b$n_input,n_complete=nrow(x),n_excluded=b$n_excluded,
   occasions=colnames(x),statistic=statistic,df=df,p_value=p,descriptive=descriptive,
   confidence_interval=NULL,interval_status="not_computed",requested_confidence_level=conf.level,
   alpha=1-conf.level,reject=p<1-conf.level,missing_strategy=missing,
   assumptions="Independent participants; matched occasions; complete blocks. Complete-case inference requires justified MCAR. Chi-squared asymptotic omnibus test; no automatic pairwise testing or method substitution.")
}
example_data<-data.frame(id=paste0("participant_",1:6),a=c(0,1,0,1,1,0), b=c(1,1,0,0,1,1), c=c(1,0,0,1,0,1))
example_result<-analyse(example_data,alternative,confidence_level)
print(example_result)
