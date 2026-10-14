# Exploratory composite scores, not latent factors or confirmatory inference.
# Choose variables and retained component count explicitly; no eigenvalue cutoff rule.
analyse <- function(data, variables=c("x1","x2","x3"), components=2L, standardise=TRUE, missing="error") {
 stopifnot(is.data.frame(data),!anyDuplicated(names(data)),"id" %in% names(data),
   !anyNA(data$id),!anyDuplicated(data$id),is.character(variables),length(variables)>=2L,
   !anyNA(variables),!anyDuplicated(variables),!"id" %in% variables,all(variables %in% names(data)),
   is.numeric(components),length(components)==1L,is.finite(components),components==floor(components),components>=1,
   is.logical(standardise),length(standardise)==1L,!is.na(standardise),
   length(missing)==1L,missing %in% c("error","complete_cases"))
 values<-data[,variables,drop=FALSE]
 stopifnot(all(vapply(values,is.numeric,logical(1))),!any(is.infinite(as.matrix(values))))
 keep<-stats::complete.cases(values)
 if(any(!keep)&&missing=="error") stop("Choose and justify complete_cases explicitly before excluding incomplete rows.")
 values<-as.matrix(values[keep,,drop=FALSE]);n<-nrow(values);p<-ncol(values)
 stopifnot(n>=3L,components<=min(n-1L,p))
 sd<-apply(values,2L,stats::sd)
 if(any(!is.finite(sd))||any(sd<=0)) stop("Constant or numerically non-finite variables cannot supply this PCA example.")
 fit<-stats::prcomp(values,center=TRUE,scale.=standardise)
 eigenvalues<-fit$sdev^2
 if(any(!is.finite(eigenvalues))||sum(eigenvalues)<=0) stop("Numerically unavailable component variances.")
 chosen<-seq_len(components);weights<-fit$rotation[,chosen,drop=FALSE];scores<-fit$x[,chosen,drop=FALSE]
 # A reproducible sign convention only; repeated eigenvalues still define a subspace, not unique axes.
 for(j in chosen) if(weights[which.max(abs(weights[,j])),j]<0) {weights[,j]<--weights[,j];scores[,j]<--scores[,j]}
 proportion<-eigenvalues/sum(eigenvalues)
 centred<-sweep(values,2L,fit$center,"-");if(standardise) centred<-sweep(centred,2L,fit$scale,"/")
 reconstruction<-scores %*% t(weights)
 residual_ss<-sum((centred-reconstruction)^2)
 list(method="Principal component analysis",status="ok",n=n,n_input=nrow(data),n_excluded=sum(!keep),
   variables=variables,components=components,standardised=standardise,center=fit$center,scale=fit$scale,
   eigenvalues=eigenvalues,variance_proportion=proportion,cumulative_variance=cumsum(proportion),
   retained_variance=sum(proportion[chosen]),weights=weights,scores=scores,row_ids=data$id[keep],
   reconstruction_residual_sum_squares=residual_ss,missing_strategy=missing,
   numerical_rank=sum(fit$sdev>max(fit$sdev)*max(n,p)*.Machine$double.eps),
   interpretation="Exploratory weighted sums; components are not latent factors. Component count is analyst-selected. No test, confidence interval or power analysis is implied. Near-equal component variances make individual axes unstable; compare their subspace. Complete-case validity needs justified missingness assumptions.")
}
example_data<-data.frame(id=1:8,x1=c(1,2,3,5,4,7,6,8),x2=c(2,1,4,3,6,5,8,7),x3=c(3,5,2,7,4,9,6,8))
example_result<-analyse(example_data,components=2L)
print(example_result)
