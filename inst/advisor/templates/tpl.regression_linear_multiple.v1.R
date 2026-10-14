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
analyse <- function(data,predictors=c("x1","x2"),tested_predictors=NULL,alternative="two.sided",conf.level=.95,missing="error") {
 stopifnot(length(predictors)>=2L,identical(alternative,"two.sided"))
 b<-fit_linear(data,predictors,conf.level,missing);partial<-NULL
 if(!is.null(tested_predictors)) {
  stopifnot(is.character(tested_predictors),length(tested_predictors)>=1L,!anyNA(tested_predictors),
    !anyDuplicated(tested_predictors),all(tested_predictors %in% predictors))
  kept<-setdiff(seq_along(predictors),match(tested_predictors,predictors))
  form<-if(length(kept)) stats::reformulate(paste0("x",kept),response="y") else y~1
  reduced<-stats::lm(form,data=b$values,na.action=stats::na.fail,singular.ok=FALSE)
  delta<-sum(stats::residuals(reduced)^2)-b$rss
  # Nested least-squares sums differ only by rounding if the added effect is exactly zero.
  if(delta<0 && abs(delta)<=64*.Machine$double.eps*max(1,b$rss)) delta<-0
  stopifnot(delta>=0)
  q<-length(tested_predictors);f<-(delta/q)/(b$rss/b$df)
  partial<-list(tested_predictors=tested_predictors,statistic=f,df_numerator=q,df_denominator=b$df,
    p_value=stats::pf(f,q,b$df,lower.tail=FALSE),partial_r_squared=delta/(delta+b$rss),
    null="The prespecified tested coefficients are jointly zero, adjusted for retained predictors")
 }
 list(method="Multiple linear regression",status="ok",n=b$n,n_input=b$n_input,n_excluded=b$n_excluded,
   predictors=predictors,coefficient_table_two_sided=b$coefficients,confidence_level=conf.level,
   omnibus=b$omnibus,p_value=b$omnibus$p_value,partial_test=partial,r_squared=1-b$rss/b$tss,
   adjusted_r_squared=1-(b$rss/b$df)/(b$tss/(b$n-1)),diagnostics=b$diagnostics,
   missing_strategy=missing,assumptions=b$assumptions,alpha=1-conf.level,reject=b$omnibus$p_value<1-conf.level)
}
example_data<-data.frame(id=1:12,x1=1:12,x2=c(1,0,2,1,3,0,1,2,0,3,2,1),y=c(3,4,4,7,6,9,8,12,10,13,12,16))
example_result<-analyse(example_data,tested_predictors="x2",alternative=alternative,conf.level=confidence_level)
print(example_result)
