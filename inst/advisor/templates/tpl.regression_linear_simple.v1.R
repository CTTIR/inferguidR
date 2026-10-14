# Prespecified numeric predictors and an intercept; independent units, no automatic selection.
# Coefficients describe conditional-mean associations, not causal effects without design support.
# Classical t/F inference assumes linear mean, constant residual variance and normal errors.
# Residual diagnostics do not automatically switch models or standard errors.
fit_linear <- function(data, predictors, conf.level, missing) {
 stopifnot(is.data.frame(data),!anyDuplicated(names(data)),all(c("id","y") %in% names(data)),
   !anyNA(data$id),!anyDuplicated(data$id),is.character(predictors),length(predictors)>=1L,
   !anyNA(predictors),!anyDuplicated(predictors),all(predictors %in% names(data)),
   !any(predictors %in% c("id","y")),is.numeric(conf.level),length(conf.level)==1L,
   is.finite(conf.level),conf.level>0,conf.level<1,
   length(missing)==1L,missing %in% c("error","complete_cases"))
 values<-data[,c("y",predictors),drop=FALSE]
 stopifnot(all(vapply(values,is.numeric,logical(1))),!any(is.infinite(as.matrix(values))))
 keep<-stats::complete.cases(values)
 if(any(!keep)&&missing=="error") stop("Choose and justify complete_cases explicitly before excluding incomplete rows.")
 values<-values[keep,,drop=FALSE];n<-nrow(values);p<-length(predictors)
 stopifnot(n>p+1L)
 # Local neutral names avoid treating user-supplied column names as formula code.
 names(values)<-c("y",paste0("x",seq_len(p)))
 fit<-stats::lm(stats::reformulate(names(values)[-1L],response="y"),data=values,na.action=stats::na.fail,singular.ok=FALSE)
 residual<-stats::residuals(fit);rss<-sum(residual^2);tss<-sum((values$y-mean(values$y))^2)
 if(!is.finite(rss)||!is.finite(tss)||tss<=0||rss<=64*n^2*.Machine$double.eps^2*max(tss,sum(values$y^2)))
   stop("Constant response or numerically exact fit: residual-variance inference is unavailable.")
 sm<-summary(fit);tab<-sm$coefficients;rownames(tab)<-c("(Intercept)",predictors)
 critical<-stats::qt((1-conf.level)/2,fit$df.residual,lower.tail=FALSE)
 coefficients<-data.frame(term=rownames(tab),estimate=tab[,1],std_error=tab[,2],statistic=tab[,3],
   p_value=tab[,4],lower=tab[,1]-critical*tab[,2],upper=tab[,1]+critical*tab[,2],row.names=NULL)
 f<-unname(sm$fstatistic);omnibus<-list(statistic=f[1],df_numerator=f[2],df_denominator=f[3],
   p_value=stats::pf(f[1],f[2],f[3],lower.tail=FALSE),null="All non-intercept coefficients are zero")
 list(fit=fit,values=values,n_input=nrow(data),n_excluded=sum(!keep),n=n,coefficients=coefficients,
   rss=rss,tss=tss,omnibus=omnibus,df=fit$df.residual,
   diagnostics=list(fitted=stats::fitted(fit),residuals=residual,
     leverage=stats::hatvalues(fit),cooks_distance=stats::cooks.distance(fit),
     scaled_design_condition_number=kappa(scale(as.matrix(values[,-1,drop=FALSE])),exact=TRUE)),
   assumptions="Independent units; prespecified linear mean; classical homoscedastic normal-error t/F inference. Complete-case validity needs justified missingness assumptions. No automatic robust substitution or causal interpretation.")
}
analyse <- function(data,predictor="x",alternative="two.sided",conf.level=.95,missing="error") {
 stopifnot(length(predictor)==1L,length(alternative)==1L,alternative %in% c("two.sided","greater","less"))
 b<-fit_linear(data,predictor,conf.level,missing);row<-b$coefficients[2L,];t<-row$statistic
 p<-switch(alternative,two.sided=2*stats::pt(abs(t),b$df,lower.tail=FALSE),greater=stats::pt(t,b$df,lower.tail=FALSE),less=stats::pt(t,b$df))
 critical<-stats::qt(if(alternative=="two.sided") (1-conf.level)/2 else 1-conf.level,b$df,lower.tail=FALSE)
 interval<-switch(alternative,two.sided=c(row$estimate-critical*row$std_error,row$estimate+critical*row$std_error),
   greater=c(row$estimate-critical*row$std_error,Inf),less=c(-Inf,row$estimate+critical*row$std_error))
 list(method="Simple linear regression",status="ok",n=b$n,n_input=b$n_input,n_excluded=b$n_excluded,
   predictor=predictor,estimate=row$estimate,std_error=row$std_error,statistic=t,df=b$df,p_value=p,
   alternative=alternative,confidence_interval=interval,confidence_level=conf.level,
   coefficient_table_two_sided=b$coefficients,omnibus=b$omnibus,r_squared=1-b$rss/b$tss,
   diagnostics=b$diagnostics,missing_strategy=missing,assumptions=b$assumptions,
   null="Slope is zero; greater means increasing conditional mean per unit predictor",
   alpha=1-conf.level,reject=p<1-conf.level)
}
example_data<-data.frame(id=1:12,x=1:12,y=c(3,4,4,7,6,9,8,12,10,13,12,16))
example_result<-analyse(example_data,alternative=alternative,conf.level=confidence_level)
print(example_result)
