# Synthetic complete repeated-measures example, one row per independent subject.
# Prespecified one-way occasion effect, with routine Greenhouse-Geisser correction.
# No Mauchly pretest, covariance-model selection or automatic post-hoc comparisons.
analyse <- function(data,occasion_columns,independent_subjects=FALSE,
 alternative="two.sided",conf.level=.95,complete_cases=FALSE) {
 stopifnot(identical(alternative,"two.sided"),is.numeric(conf.level),length(conf.level)==1L,
  is.finite(conf.level),conf.level>0,conf.level<1,isTRUE(independent_subjects),
  is.data.frame(data),!anyDuplicated(names(data)),is.character(occasion_columns),
  length(occasion_columns)>=3L,length(occasion_columns)<=30L,!anyNA(occasion_columns),
  !anyDuplicated(occasion_columns),all(nzchar(occasion_columns)),!('id'%in%occasion_columns),
  all(c('id',occasion_columns)%in%names(data)),nrow(data)>=3L,nrow(data)<=10000L,
  is.atomic(data$id),is.null(dim(data$id)),!anyNA(data$id),!anyDuplicated(data$id),
  all(nzchar(as.character(data$id))),is.logical(complete_cases),length(complete_cases)==1L,!is.na(complete_cases),
  all(vapply(data[occasion_columns],function(x)is.numeric(x)&&is.null(dim(x)),logical(1))))
 y<-as.matrix(data[occasion_columns]);stopifnot(!any(is.infinite(y)))
 n_total<-nrow(y);keep<-stats::complete.cases(y)
 if(!all(keep)&&!complete_cases)stop('Incomplete subjects require an explicit complete-case choice.')
 y<-y[keep,,drop=FALSE];n<-nrow(y);k<-ncol(y);stopifnot(n>=3L)
 unavailable<-function(reason)list(status='not_applicable',reason=reason,p_value=NULL,
  confidence_interval=NULL,interval_status='not_computed',n_total=n_total,n_used=n,n_excluded=n_total-n)
 # Remove subject offsets before global scaling, retaining representable differences.
 y<-sweep(y,1L,y[,1L],'-')
 if(any(!is.finite(y)))return(unavailable('Within-subject differences exceed the numerical range.'))
 scale<-max(abs(y));if(scale==0)return(unavailable('All within-subject differences are zero.'))
 y<-y/scale;y<-y-rowMeans(y)
 occasion_mean<-colMeans(y);residual<-sweep(y,2L,occasion_mean,'-')
 sse<-sum(residual^2);ss<-n*sum(occasion_mean^2)
 if(!is.finite(sse)||sse<=0)return(unavailable('Within-subject residual variance is zero or unresolved.'))
 # Residual covariance has already been projected onto the occasion contrast space.
 V<-crossprod(residual);vscale<-max(abs(V))
 if(!is.finite(vscale)||vscale<=0)return(unavailable('The contrast covariance is unresolved.'))
 V<-V/vscale;epsilon<-sum(diag(V))^2/((k-1L)*sum(V^2));bound<-1/(k-1L)
 if(!is.finite(epsilon)||epsilon<bound-1e-12||epsilon>1+1e-12)
  return(unavailable('The Greenhouse-Geisser correction is outside its numerical bounds.'))
 epsilon<-min(1,max(bound,epsilon));df1<-k-1L;df2<-(n-1L)*df1
 F<-(ss/df1)/(sse/df2);p<-stats::pf(F,df1*epsilon,df2*epsilon,lower.tail=FALSE)
 if(!is.finite(F)||!is.finite(p))return(unavailable('The corrected F test is not numerically resolved.'))
 list(status='ok',method='One-way repeated-measures ANOVA; routine Greenhouse-Geisser correction',
  statistic=F,epsilon_gg=epsilon,df_numerator=df1*epsilon,df_denominator=df2*epsilon,p_value=p,
  uncorrected_df=c(numerator=df1,denominator=df2),partial_eta_squared=if(ss==0)0 else 1/(1+sse/ss),
  n_total=n_total,n_used=n,n_excluded=n_total-n,occasions=occasion_columns,
  alpha=1-conf.level,confidence_interval=NULL,interval_status='not_computed',requested_confidence_level=conf.level,
  interpretation='Omnibus equality of occasion means. The same F statistic uses GG-adjusted numerator and denominator degrees of freedom; the correction is applied routinely.',
  assumptions='Independent subjects with jointly normal repeated outcomes for the approximate GG F inference; one within-subject factor, no between-subject factors or covariates. Complete-case inference needs justified missingness assumptions.',
  limitations='Three to thirty occasions and at least three complete subjects. GG is an estimated finite-sample correction, not an exact covariance-free test. No post-hoc comparisons or effect-size sampling interval are computed.')
}
example_data<-data.frame(id=1:6,before=c(5,7,6,8,4,9),during=c(7,8,6,11,6,10),after=c(8,8,9,12,7,13))
example_result<-analyse(example_data,c('before','during','after'),TRUE,alternative,confidence_level)
print(example_result)
